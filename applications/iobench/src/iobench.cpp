/**
 * @file iobench.cpp
 * @brief Benchmark binary Field I/O: parallel (MPI-IO) vs. serial (rank 0) path.
 *
 * Writes and reads an SU(3) double precision gauge configuration with
 * GaugeField::config_write() / config_read(), first with hila::parallel_io = true,
 * then (optionally) with the serial rank 0 path.  Each read is checked to reproduce
 * the written configuration exactly.
 *
 * Parameters (file "parameters"):
 *     lattice size     32, 32, 32, 32    # NDIM values
 *     repeats          3                 # timed write+read cycles per mode
 *     serial           yes               # also time the serial path (yes/no)
 *     file name        iobench.cfg       # put this on the file system to be tested
 *
 * Run e.g. "mpirun -np 8 ./build/iobench".
 *
 * On Lustre the striping of the file is the main factor for parallel write performance.
 * Files get the striping of the directory when created ("lfs getstripe <dir>").  E.g. on
 * LUMI-g, 64^4 lattice on 8 nodes: the default progressive layout (1 -> 4 -> 8 stripes)
 * gave 2.4 GB/s write, "lfs setstripe -c 16 -S 4M <dir>" 6.2 GB/s.  With Cray MPICH,
 * striping can also be set with MPI-IO hints, e.g.
 *     export MPICH_MPIIO_HINTS="*:striping_factor=16:striping_unit=4194304"
 *     export MPICH_MPIIO_HINTS_DISPLAY=1     # show hints in use
 *     export MPICH_MPIIO_STATS=1             # I/O statistics
 * Note: the read rate measured here can be inflated by caching, because the file was
 * just written by the same job.
 */

#include "hila.h"

struct io_timing {
    double min_write = 1e100, sum_write = 0, min_read = 1e100, sum_read = 0;
    int n = 0;
};

io_timing time_io(GaugeField<SU<3, double>> &U, const std::string &fname, int repeats,
                  bool parallel) {

    io_timing t;
    GaugeField<SU<3, double>> V;

    hila::parallel_io = parallel;

    for (int rep = 0; rep < repeats; rep++) {

        hila::synchronize();
        double t0 = hila::gettime();
        U.config_write(fname);
        hila::synchronize();
        double tw = hila::gettime() - t0;

        foralldir (d)
            V[d] = 0;

        hila::synchronize();
        t0 = hila::gettime();
        V.config_read(fname);
        hila::synchronize();
        double tr = hila::gettime() - t0;

        // check that the config was read back exactly
        double diff = 0;
        foralldir (d) {
            onsites (ALL) {
                diff += (V[d][X] - U[d][X]).squarenorm();
            }
        }
        if (diff != 0) {
            hila::out0 << "ERROR: configuration read back differs from written, diff " << diff
                       << '\n';
            hila::terminate(1);
        }

        hila::out0 << (parallel ? "  parallel" : "  serial  ") << "  rep " << rep << ":  write "
                   << tw << " s,  read " << tr << " s\n";

        t.min_write = std::min(t.min_write, tw);
        t.min_read = std::min(t.min_read, tr);
        t.sum_write += tw;
        t.sum_read += tr;
        t.n++;
    }
    return t;
}

void report(const char *name, const io_timing &t, double gb) {
    hila::out0 << name << "  write: best " << t.min_write << " s (" << gb / t.min_write
               << " GB/s), mean " << t.sum_write / t.n << " s   read: best " << t.min_read
               << " s (" << gb / t.min_read << " GB/s), mean " << t.sum_read / t.n << " s\n";
}

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

    hila::input par("parameters");
    CoordinateVector lsize = par.get("lattice size");
    int repeats = par.get("repeats");
    bool do_serial = (par.get_item("serial", {"yes", "no"}) == 0);
    std::string fname = par.get("file name");
    par.close();

    lattice.setup(lsize);
    hila::seed_random(1);

    GaugeField<SU<3, double>> U;
    foralldir (d)
        U[d].gaussian_random();

    double gb = (double)NDIM * lattice.volume() * sizeof(SU<3, double>) / 1e9;
    hila::out0 << "I/O benchmark: SU(3) double gauge configuration, " << gb << " GB, "
               << hila::number_of_nodes() << " ranks, file " << fname << '\n';

    io_timing tpar = time_io(U, fname, repeats, true);
    io_timing tser;
    if (do_serial)
        tser = time_io(U, fname, repeats, false);

    hila::print_dashed_line();
    report("PARALLEL", tpar, gb);
    if (do_serial) {
        report("SERIAL  ", tser, gb);
        hila::out0 << "Speedup: write " << tser.min_write / tpar.min_write << ", read "
                   << tser.min_read / tpar.min_read << '\n';
    }
    hila::print_dashed_line();

    if_rank0 () {
        std::remove(fname.c_str());
    }

    hila::finishrun();
}
