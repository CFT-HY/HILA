/**
 * Test deep halo gathers (HALO_DEPTH > 1).
 *
 * Field value is an exact function of the global coordinates.  After
 * f.gather_depth(d, k) the element at site_layout neighbour_offset(i, d, m), m = 1..k,
 * must equal the value at coordinate(i) + m*d (modded, sign flipped across an
 * antiperiodic boundary).  Also tests the gather bookkeeping: nn halo already valid,
 * repeated gathers, shallow -> deep, invalidation by field changes.
 *
 * Lattice size is read from file "parameters" (key "lattice size").
 */

#include "hila.h"

#if HALO_DEPTH < 2
#error "test_halo_depth: compile with HALO_DEPTH >= 2"
#endif

static int64_t errors = 0;
static int reported = 0;

static void report(const std::string &s) {
    errors++;
    if (reported++ < 10)
        hila::out << "rank " << hila::myrank() << ": " << s << '\n';
}

// value of the field at global coordinate c
static double value_at(const CoordinateVector &c, const CoordinateVector &mult, double add) {
    return (double)c.dot(mult) + add;
}

// check halo layers up to depth k of direction d, field f = value_at(X) + add
static void check(const Field<double> &f, Direction d, int k, const CoordinateVector &mult,
                  double add, bool antiperiodic, const std::string &what) {
    const auto &layout = lattice->mynode.layout;
    const CoordinateVector lsz = lattice.size();
    int dim = abs(d);
    for (unsigned i = 0; i < lattice->mynode.volume; i++) {
        CoordinateVector c = layout.coordinates(i);
        for (int m = 1; m <= k; m++) {
            unsigned idx = layout.neighbour_offset(i, d, m);
            CoordinateVector cg = c + m * d;
            double sign = 1;
            if (antiperiodic && (cg[dim] < 0 || cg[dim] >= lsz[dim]))
                sign = -1;
            cg = cg.mod(lsz);
            double expect = sign * value_at(cg, mult, add);
            double v = f.get_value_at(idx);
            if (v != expect)
                report(what + ": dir " + std::to_string((int)d) + " k " + std::to_string(k) +
                       " m " + std::to_string(m) + " site " + std::to_string(i) + " index " +
                       std::to_string(idx) + " value " + std::to_string(v) + " expected " +
                       std::to_string(expect));
        }
    }
}

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

    hila::input par("parameters");
    CoordinateVector lsize = par.get("lattice size");
    par.close();

    lattice.setup(lsize);

    CoordinateVector mult;
    {
        int64_t m = 1;
        foralldir (d) {
            mult[d] = m;
            m *= 37;
        }
    }

    Field<double> f;
    double add = 0.5;
    f[ALL] = X.coordinates().dot(mult) + add;

    // 1: all directions and depths, fresh field each time
    for (int k = 2; k <= HALO_DEPTH; k++) {
        for (Direction d = (Direction)0; d < NDIRS; ++d) {
            f[ALL] = X.coordinates().dot(mult) + add; // invalidates halos
            f.gather_depth(d, k);
            check(f, d, k, mult, add, false, "fresh");
        }
    }

    // 2: nn halo already valid -> layer 0 not sent
    add = 1.5;
    f[ALL] = X.coordinates().dot(mult) + add;
    for (Direction d = (Direction)0; d < NDIRS; ++d) {
        f.gather(d, ALL);
        f.gather_depth(d, HALO_DEPTH);
        check(f, d, HALO_DEPTH, mult, add, false, "nn first");
    }

    // 3: repeated gather is avoided, shallower too.  Deep gather also makes nn valid
    {
        int64_t done = hila::n_gather_done;
        for (Direction d = (Direction)0; d < NDIRS; ++d) {
            f.gather_depth(d, HALO_DEPTH);
            f.gather_depth(d, 2);
            f.gather(d, ALL);
        }
        if (hila::n_gather_done != done)
            report("repeated gathers were not avoided");
    }

    // 4: shallow started, then deeper; and start without wait followed by a change
    add = 2.5;
    f[ALL] = X.coordinates().dot(mult) + add;
    for (Direction d = (Direction)0; d < NDIRS; ++d) {
        f.start_gather_depth(d, 2);
        f.gather_depth(d, HALO_DEPTH);
        check(f, d, HALO_DEPTH, mult, add, false, "shallow->deep");
    }
    for (Direction d = (Direction)0; d < NDIRS; ++d)
        f.start_gather_depth(d, HALO_DEPTH);
    add = 3.5;
    f[EVEN] = X.coordinates().dot(mult) + add; // drops ongoing comms, invalidates
    f[ODD] = X.coordinates().dot(mult) + add;
    for (Direction d = (Direction)0; d < NDIRS; ++d) {
        if (f.is_gathered_depth(d, 2) && !lattice->mynode.layout.wrap[abs(d)])
            report("halo still marked valid after change");
        f.gather_depth(d, HALO_DEPTH);
        check(f, d, HALO_DEPTH, mult, add, false, "after change");
    }

    // 5: several directions on the fly at the same time
    add = 4.5;
    f[ALL] = X.coordinates().dot(mult) + add;
    for (Direction d = (Direction)0; d < NDIRS; ++d)
        f.start_gather_depth(d, HALO_DEPTH);
    for (Direction d = (Direction)0; d < NDIRS; ++d) {
        f.wait_gather_depth(d, HALO_DEPTH);
        check(f, d, HALO_DEPTH, mult, add, false, "simultaneous");
    }

#ifdef SPECIAL_BOUNDARY_CONDITIONS
    // 6: antiperiodic b.c. in divided directions
    {
        Field<double> g;
        foralldir (d)
            if (!lattice->mynode.layout.wrap[d])
                g.set_boundary_condition(d, hila::bc::ANTIPERIODIC);
        add = 5.5;
        g[ALL] = X.coordinates().dot(mult) + add;
        for (Direction d = (Direction)0; d < NDIRS; ++d) {
            bool ap = !lattice->mynode.layout.wrap[abs(d)];
            g.gather_depth(d, HALO_DEPTH);
            check(g, d, HALO_DEPTH, mult, add, ap, "antiperiodic");
        }
    }
#endif

    int64_t total = errors;
    hila::reduce_node_sum(total);

    hila::out0 << "halo_depth test, HALO_DEPTH " << HALO_DEPTH << ", lattice " << lsize << ", "
               << hila::number_of_nodes() << " nodes, node divisions "
               << lattice->nodes.n_divisions << ": "
               << (total == 0 ? "PASSED" : "FAILED, " + std::to_string(total) + " errors")
               << '\n';

    hila::finishrun();
}
