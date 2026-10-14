/**
 * Test plumbing/site_layout.h: on-the-fly coordinates, site indices and neighbour
 * indices must agree with the lattice arrays (coordinates, neighb, halo offsets)
 * on every site and direction.  Also tests neighbour_offset() and deep halo layers
 * (HALO_DEPTH > 1) for depths 1..4, independently of the compiled HALO_DEPTH.
 *
 * With LAYOUT_ON_THE_FLY (always with EVEN_SITES_FIRST=0) there are no lattice arrays:
 * coordinates are checked against direct lexicographic decoding, on-node neighbours against
 * lattice site_index(), and halo indices for being one-to-one and in their ranges.
 *
 * Lattice size is read from file "parameters" (key "lattice size").
 */

#include "hila.h"
#include "plumbing/site_layout.h"

#if defined(GPU_OVERLAP_COMM) || defined(SUBNODE_LAYOUT)
#error "test_site_layout: GPU_OVERLAP_COMM or SUBNODE_LAYOUT layouts are not supported"
#endif

#ifndef LAYOUT_ON_THE_FLY
#define HAVE_LATTICE_ARRAYS
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

    // fast_divisor must agree with integer division for all 32-bit numerators
    {
        std::vector<unsigned> divisors;
        for (unsigned d = 1; d <= 2048; d++)
            divisors.push_back(d);
        for (unsigned d : {4095u, 65535u, 65536u, 1000003u, 0x7fffffffu, 0x80000000u,
                           0xfffffffeu, 0xffffffffu})
            divisors.push_back(d);
        uint64_t rng = 12345;
        for (unsigned d : divisors) {
            hila::fast_divisor fd;
            fd.set(d);
            auto check = [&](unsigned n) {
                if (fd.div(n) != n / d)
                    report("fast_divisor " + std::to_string(n) + " / " + std::to_string(d) +
                           " gives " + std::to_string(fd.div(n)));
            };
            for (unsigned n : {0u, 1u, d - 1, d, d + 1, 2 * d - 1, 0xffffffffu, 0xfffffffeu,
                               0xffffffffu - 0xffffffffu % d, 0xffffffffu - 0xffffffffu % d - 1})
                check(n);
            for (int i = 0; i < 2000; i++) {
                rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
                check((unsigned)(rng >> 32));
            }
        }
    }

    // incremental advance() must reproduce site(idx) for all idx
    {
        auto s = layout.site(0);
        for (unsigned i = 1; i < node.volume; i++) {
            layout.advance(s, i);
            auto r = layout.site(i);
            bool same = (s.lex == r.lex && s.parity == r.parity);
            foralldir (e)
                if (s.lc[e] != r.lc[e])
                    same = false;
            if (!same) {
                report("advance() to site " + std::to_string(i) + " differs from site()");
                s = r;
            }
        }
    }

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
#ifdef HAVE_LATTICE_ARRAYS
        CoordinateVector cref = lattice->coordinates(i);
#else
        // lexicographic order of sites, or lex order of even sites then odd sites
        CoordinateVector cref;
        {
#ifdef EVEN_SITES_FIRST
            // i-th site of its parity in lex order: search (slow but simple)
            static std::vector<CoordinateVector> esf_order;
            if (esf_order.empty()) {
                esf_order.resize(node.volume);
                unsigned ne = 0, no = node.evensites;
                for (unsigned lex = 0; lex < node.volume; lex++) {
                    CoordinateVector cc;
                    unsigned r = lex;
                    foralldir (e) {
                        cc[e] = r % node.size[e] + node.min[e];
                        r /= node.size[e];
                    }
                    if (cc.parity() == EVEN)
                        esf_order[ne++] = cc;
                    else
                        esf_order[no++] = cc;
                }
            }
            cref = esf_order[i];
#else
            unsigned r = i;
            foralldir (e) {
                cref[e] = r % node.size[e] + node.min[e];
                r /= node.size[e];
            }
#endif
        }
#endif
        if (c != cref) {
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
#ifdef HAVE_LATTICE_ARRAYS
            if (nb != lattice->neighb[d][i])
                report("site " + std::to_string(i) + " dir " + std::to_string(d) +
                       " neighbour " + std::to_string(nb) + " vs " +
                       std::to_string(lattice->neighb[d][i]));
            if (layout.is_off_node(s, d) != (lattice->neighb[d][i] >= node.volume))
                report("site " + std::to_string(i) + " dir " + std::to_string(d) +
                       " is_off_node differs");
#else
            CoordinateVector cn = (c + (Direction)d).mod(lattice.size());
            int dim = layout.dim_of(d);
            bool off = !layout.wrap[dim] && !lattice->is_on_mynode(cn);
            // with 2 node divisions the neighbour can be on node, but via halo
            if (!layout.wrap[dim] && (is_up_dir((Direction)d) ? s.lc[dim] == node.size[dim] - 1
                                                              : s.lc[dim] == 0))
                off = true;
            if (layout.is_off_node(s, d) != off)
                report("site " + std::to_string(i) + " dir " + std::to_string(d) +
                       " is_off_node differs");
            if (!off && nb != lattice->site_index(cn))
                report("site " + std::to_string(i) + " dir " + std::to_string(d) +
                       " neighbour " + std::to_string(nb) + " vs " +
                       std::to_string(lattice->site_index(cn)));
            if (off && (nb < layout.halo_offset[d] ||
                        nb >= layout.halo_offset[d] + layout.face_size[dim]))
                report("site " + std::to_string(i) + " dir " + std::to_string(d) +
                       " halo index " + std::to_string(nb) + " outside its range");
#endif
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // neighbour_offset() and deep halo layers, for depths 1..4.
    //  - on-node targets must agree with lattice site_index (incl. wrap on node)
    //  - layer 0 halo must agree with following neighb[] k-1 steps and then 1 step
    //  - (direction, global coordinate) <-> halo index must be one-to-one, indices
    //    inside [halo_offset/deep_offset ranges] and below deep_halo_end
    ///////////////////////////////////////////////////////////////////////////

    const CoordinateVector lsz = lattice.size();
    for (int depth = 1; depth <= 4; depth++) {
        bool ok = true;
        foralldir (d)
            if (lattice->nodes.n_divisions[d] > 1 && node.size[d] < depth)
                ok = false;
        if (!ok)
            continue;

        hila::site_layout_struct lay;
        lay.setup(node.min, node.size, lattice->nodes.n_divisions, depth);

        if (depth == HALO_DEPTH && lay.field_alloc_size != node.field_alloc_size)
            report("depth " + std::to_string(depth) + " field_alloc_size " +
                   std::to_string(lay.field_alloc_size) + " vs " +
                   std::to_string(node.field_alloc_size));
#ifdef EVEN_SITES_FIRST
        // nn halos do not depend on the depth
        if (lay.halo_end != layout.halo_end)
            report("depth " + std::to_string(depth) + " halo_end changed");
#else
        // all layers together after each direction
        if (lay.halo_end != lay.deep_halo_end)
            report("depth " + std::to_string(depth) + " halo_end != deep_halo_end");
#endif

        // halo index -> (direction, global coordinate), and reverse
        std::vector<int> owner_dir(lay.deep_halo_end, -1);
        std::vector<CoordinateVector> owner_c(lay.deep_halo_end);
        size_t n_filled = 0;

        for (unsigned i = 0; i < node.volume; i++) {
            CoordinateVector c = lay.coordinates(i);
            auto s = lay.site(i);
            for (int d = 0; d < NDIRS; d++) {
                Direction dir = (Direction)d;
                int dim = abs(dir);
                for (int k = 1; k <= depth; k++) {
                    unsigned nb = lay.neighbour_offset(s, d, k);
#ifdef EVEN_SITES_FIRST
                    // nn halo indices do not depend on the depth
                    const auto &nn_layout = layout;
#else
                    const auto &nn_layout = lay;
#endif
                    if (k == 1 && nb != nn_layout.neighbour(s, d))
                        report("depth " + std::to_string(depth) + " site " + std::to_string(i) +
                               " dir " + std::to_string(d) + " k=1 differs from neighbour()");

                    CoordinateVector cg = c + k * dir;
                    cg = cg.mod(lsz);

                    CoordinateVector lc = cg - node.min;
                    bool on_node = true;
                    foralldir (e)
                        if (lc[e] < 0 || lc[e] >= node.size[e])
                            on_node = false;
                    // with 2 divisions, wrapped target may be back on node, but it is
                    // still fetched through the halo
                    if (lay.wrap[dim]) {
                        if (!on_node || nb != lattice->site_index(cg))
                            report("depth " + std::to_string(depth) + " site " +
                                   std::to_string(i) + " dir " + std::to_string(d) +
                                   " k " + std::to_string(k) + " wrapped index " +
                                   std::to_string(nb));
                        continue;
                    }

                    int x = c[dim] - node.min[dim] + (is_up_dir(dir) ? k : -k);
                    if (x >= 0 && x < node.size[dim]) {
                        if (nb != lattice->site_index(cg))
                            report("depth " + std::to_string(depth) + " site " +
                                   std::to_string(i) + " dir " + std::to_string(d) +
                                   " k " + std::to_string(k) + " on-node index " +
                                   std::to_string(nb) + " vs " +
                                   std::to_string(lattice->site_index(cg)));
                        continue;
                    }

                    int l = (x < 0) ? -1 - x : x - node.size[dim];
                    if (l == 0) {
                        unsigned j = i;
                        for (int m = 0; m < k; m++)
#ifdef HAVE_LATTICE_ARRAYS
                            j = lattice->neighb[d][j];
#else
                            j = lay.neighbour(j, d);
#endif
                        if (nb != j)
                            report("depth " + std::to_string(depth) + " site " +
                                   std::to_string(i) + " dir " + std::to_string(d) +
                                   " k " + std::to_string(k) + " nn halo " +
                                   std::to_string(nb) + " vs " + std::to_string(j));
                    } else {
                        if (nb < lay.deep_offset[d] ||
                            nb >= lay.deep_offset[d] + (depth - 1) * lay.face_size[dim])
                            report("depth " + std::to_string(depth) + " site " +
                                   std::to_string(i) + " dir " + std::to_string(d) +
                                   " k " + std::to_string(k) + " deep index " +
                                   std::to_string(nb) + " outside its range");
                    }

                    if (nb < node.volume || nb >= lay.deep_halo_end) {
                        report("depth " + std::to_string(depth) + " halo index " +
                               std::to_string(nb) + " out of range");
                        continue;
                    }
                    if (owner_dir[nb] < 0) {
                        owner_dir[nb] = d;
                        owner_c[nb] = cg;
                        n_filled++;
                    } else if (owner_dir[nb] != d || owner_c[nb] != cg) {
                        report("depth " + std::to_string(depth) + " halo index " +
                               std::to_string(nb) + " used for two different sites");
                    }
                }
            }
        }

        // every halo element must be used exactly by one (dir, coordinate):
        // count of distinct targets == number of halo elements
        if (n_filled != lay.deep_halo_end - node.volume)
            report("depth " + std::to_string(depth) + " halo elements filled " +
                   std::to_string(n_filled) + " vs " +
                   std::to_string(lay.deep_halo_end - node.volume));
    }

    int64_t total = errors;
    hila::reduce_node_sum(total);

    hila::out0 << "site_layout test, lattice " << lsize << ", " << hila::number_of_nodes()
               << " nodes, node divisions " << lattice->nodes.n_divisions << ": "
               << (total == 0 ? "PASSED" : "FAILED, " + std::to_string(total) + " errors")
               << '\n';

    hila::finishrun();
}
