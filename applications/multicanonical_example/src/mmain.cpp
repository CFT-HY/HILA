////////////////////////////////////////////////////////////////////////////////
/// @author Jaakko Hällfors
////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random>

#include "hila.h"
#include "multicanonical.h"

static_assert(NDIM == 3, "NDIM must be 3 here");

Field<double> compute_local_action(Field<double> phi);
void multican_update(Field<double> &phi);
double order_parameter(Field<double> &phi);
void iterate_weights(Field<double> &phi);

// For the action we take a simple real scalar field with a
// potential with two degenerate minima.
Field<double> compute_local_action(Field<double> phi) {
    Field<double> action;
    onsites(ALL) {
        action[X] = -pow(phi[X], 2) + pow(phi[X], 4);
    }
    foralldir(d) onsites(ALL) {
        action[X] += pow(phi[X] - phi[X - d], 2);
        action[X] += pow(phi[X + d] - phi[X], 2);
    }
    return action;
}

// Construct some kind of update algorithm. Here we have a standard
// Metropolis-Hastings algorithm that creates a proposal field for the
// final multicanonical acceptance in the end.
void multican_update(Field<double> &phi, Parity PAR) {
    // Get a temporary field
    Field<double> new_phi = phi;

    // Get a field of random numbers
    Field<double> delta = 0;
    onsites(PAR) hila::gaussian_random(delta[X]);

    // Add change to temporary field
    new_phi += 0.5 * delta;

    // Get the local differences in actions
    auto old_action = compute_local_action(phi);
    auto new_action = compute_local_action(new_phi);
    onsites(PAR) {
        // compute log(exp(- delta S))
        double log_P = -(new_action[X] - old_action[X]);
        // Get a random uniform from [0,1] and reject based on that
        double log_rand = log(hila::random());
        if (log_rand > log_P) {
            new_phi[X] = phi[X];
        }
    }

    double OP_old = order_parameter(phi);
    double OP_new = order_parameter(new_phi);

    // Use the multicanonical accept_reject for a final update decision
    if (hila::muca::accept_reject(OP_old, OP_new)) {
        phi = new_phi;
    }
}

// Mean value of phi is a convenient order parameter, since its cheap
// to compute and has different values in the different minima.
double order_parameter(Field<double> &phi) {
    double OP = phi.sum();
    return OP / lattice.volume();
}

// The weight iteration consists of a simple loop
// where muca::iterate_weights is being fed with measurements
// of the order parameter. The status informs us when the iteration
// algorithm thinks its done and kills the loop. We save the weight data
// for good measure.
void iterate_weights(Field<double> &phi) {
    bool iterate_status = true;
    while (iterate_status) {
        // Perform a number of updates
        for (int i = 0; i < 25; i++) {
            multican_update(phi, ODD);
            multican_update(phi, EVEN);
        }
        // compute and feed the order parameter of the configuration to the
        // iterator
        double OP = order_parameter(phi);
        iterate_status = hila::muca::iterate_weights(OP);
    }
    hila::muca::write_weight_function(hila::muca::generate_outfile_name());
}

void multi_iterate_weights(std::vector<Field<double>> &N_phi) {
    bool iterate_status = true;
    std::vector<double> OPs(N_phi.size());
    while (iterate_status) {
        // Perform a number of updates
        for (int j = 0; j < N_phi.size(); j++) {
            for (int i = 0; i < 50; i++) {
                multican_update(N_phi[j], ODD);
                multican_update(N_phi[j], EVEN);
            }
            OPs[j] = order_parameter(N_phi[j]);
        }

        // compute and feed the order parameter of the configuration to the
        // iterator
        iterate_status = hila::muca::iterate_weights(OPs[0]);
    }
    hila::muca::write_weight_function(hila::muca::generate_outfile_name());
}

// Here we combine the above functions into a full simulation that
// iterates the weights as well as performs a run of measurements
// on the order parameter using the iterated (or loaded) weight function.
// Check the order parameter histograms to see the effect!
int main(int argc, char *argv[]) {
    hila::cmdline.add_flag("-muca", "If used, multicanonical methods are applied.");
    hila::initialize(argc, argv);
    hila::cmdline.print_help();
    // hila::finishrun();

    lattice.setup({12, 12, 12});
    hila::seed_random(0);

    hila::muca::initialise("muca_parameters");

    // Initialise a field
    Field<double> phi = 0;

    // Open a file for measurements
    bool muca = hila::cmdline.flag_present("-muca");
    std::ofstream MFile;
    if (muca) {
        hila::out0 << "Running a multicanonical simulation.\n";
        MFile.open("muca_measurements", std::ios_base::app);
    } else {
        hila::out0 << "Running a standard simulation.\n";
        MFile.open("ca_measurements", std::ios_base::app);
    }

    if (muca)
        iterate_weights(phi);

    const int N = 10000;
    // Get some measurements
    for (int i = 1; i <= N; i++) {
        double OP = order_parameter(phi);
        char buffer[1024];
        // Remember that you need the weight of each measured configuration
        // for the post-processing!
        sprintf(buffer, "%d\t%e\t%e\n", i, OP, hila::muca::weight(OP));
        if (hila::myrank() == 0)
            MFile << std::string(buffer);

        if (i % (N / 100) == 0)
            hila::out0 << 100 * i / N << "% done\n";
        // Perform a set of updates
        for (int i = 0; i < 50; i++) {
            multican_update(phi, ODD);
            multican_update(phi, EVEN);
        }
    }

    MFile.close();
    hila::finishrun();
}
