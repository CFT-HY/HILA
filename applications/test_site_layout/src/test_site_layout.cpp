/**
 * Test plumbing/site_layout.h: on-the-fly coordinates, site indices and neighbour
 * indices must agree with the lattice arrays (coordinates, neighb, halo offsets)
 * on every site and direction.
 *
 * Lattice size is read from file "parameters" (key "lattice size").
 */

#include "hila.h"
#include "plumbing/site_layout.h"

#if defined(GPU_OVERLAP_COMM) || defined(SUBNODE_LAYOUT) || !defined(EVEN_SITES_FIRST)
#error "test_site_layout: only EVEN_SITES_FIRST layout without GPU_OVERLAP_COMM is supported"
#endif

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

    hila::input par("parameters");
    CoordinateVector lsize = par.get("lattice size");
    par.close();

    lattice.setup(lsize);

    const auto &node = lattice->mynode;

    hila::site_layout_struct layout;
    layout.setup(node.min, node.size, lattice->nodes.n_divisions);

    int64_t errors = 0;
    int reported = 0;
    auto report = [&](const std::string &s) {
        errors++;
        if (reported++ < 10)
            hila::out << "rank " << hila::myrank() << ": " << s << '\n';
    };

    if (layout.volume != node.volume)
        report("volume " + std::to_string(layout.volume) + " vs " +
               std::to_string(node.volume));
    if (layout.evensites != node.evensites)
        report("evensites " + std::to_string(layout.evensites) + " vs " +
               std::to_string(node.evensites));
    if (layout.field_alloc_size != node.field_alloc_size)
        report("field_alloc_size " + std::to_string(layout.field_alloc_size) + " vs " +
               std::to_string(node.field_alloc_size));

    for (int d = 0; d < NDIRS; d++) {
        const auto &fn = lattice->nn_comminfo[d].from_node;
        if (fn.sites > 0) {
            if (layout.halo_offset[d] != fn.buffer || layout.halo_evensites[d] != fn.evensites)
                report("halo dir " + std::to_string(d) + " offset " +
                       std::to_string(layout.halo_offset[d]) + " vs " +
                       std::to_string(fn.buffer) + ", evensites " +
                       std::to_string(layout.halo_evensites[d]) + " vs " +
                       std::to_string(fn.evensites));
        }
    }

    for (unsigned i = 0; i < node.volume; i++) {
        CoordinateVector c = layout.coordinates(i);
        if (c != lattice->coordinates(i)) {
            report("site " + std::to_string(i) + " coordinates differ");
            continue;
        }
        if (layout.site_index(c) != i)
            report("site " + std::to_string(i) + " site_index gives " +
                   std::to_string(layout.site_index(c)));
        if (layout.site_parity(i) != c.parity())
            report("site " + std::to_string(i) + " parity differs");

        auto s = layout.site(i);
        for (int d = 0; d < NDIRS; d++) {
            unsigned nb = layout.neighbour(s, d);
            if (nb != lattice->neighb[d][i])
                report("site " + std::to_string(i) + " dir " + std::to_string(d) +
                       " neighbour " + std::to_string(nb) + " vs " +
                       std::to_string(lattice->neighb[d][i]));
            if (layout.is_off_node(s, d) != (lattice->neighb[d][i] >= node.volume))
                report("site " + std::to_string(i) + " dir " + std::to_string(d) +
                       " is_off_node differs");
        }
    }

    int64_t total = errors;
    hila::reduce_node_sum(total);

    hila::out0 << "site_layout test, lattice " << lsize << ", " << hila::number_of_nodes()
               << " nodes, node divisions " << lattice->nodes.n_divisions << ": "
               << (total == 0 ? "PASSED" : "FAILED, " + std::to_string(total) + " errors")
               << '\n';

    hila::finishrun();
}
