#include "hila.h"
#include "catch.hpp"
#include "catch_main.hpp"
#include "read_write_snapshots.hpp"
#include <vector>
#include <map>
#include <cmath>

/**
 * @brief If static data is used in multiple tests, it is defined here
 * 
 */
struct SnapshotStruct {
    int ls = lattice_size;
    //Adjacent_nodes_list depicts what rank is adjacent to given rank ordered in directions {e_x,e_y,e_z,-e_x,-e_y,-e_z}
    std::vector<std::vector<int>> adjacent_nodes_list_8 = {
        {1, 2, 4, 1, 2, 4}, {0, 3, 5, 0, 3, 5}, //z - direction front x{0}{1}
        {3, 0, 6, 3, 0, 6}, {2, 1, 7, 2, 1, 7}, //                    x{2}{3}
                                                //                      y  y 
        {5, 6, 0, 5, 6, 0}, {4, 7, 1, 4, 7, 1}, //z - direction back  x{4}{5}
        {7, 4, 2, 7, 4, 2}, {6, 5, 3, 6, 5, 3}  //                    x{6}{7}
                                                //                      y  y 
    };
    std::vector<std::vector<int>> adjacent_nodes_list_4 = {
        {0, 1, 2, 0, 1, 2}, {1, 0, 3, 1, 0, 3},
        {2, 3, 0, 2, 3, 0}, {3, 2, 1, 3, 2, 1}
    }; 
    std::vector<std::vector<int>> adjacent_nodes_list_2 = {{0, 0, 1, 0, 0, 1}, {1, 1, 0, 1, 1, 0}};

    std::vector<std::vector<int>> adjacent_nodes_list_1 = {{0, 0, 0, 0, 0, 0}};

    std::map<int, std::vector<std::vector<int>>> adjacent_nodes_dict = {
        { 1, adjacent_nodes_list_1 },
        { 2, adjacent_nodes_list_2 },
        { 4, adjacent_nodes_list_4 },
        { 8, adjacent_nodes_list_8 }
    };

    std::vector<CoordinateVector> node_first_coordinate_1 = {{0, 0, 0}};
    std::vector<CoordinateVector> node_first_coordinate_2 = {{0, 0, 0}, {0, 0, ls/2}};
    std::vector<CoordinateVector> node_first_coordinate_4 = {
        {   0,   0,   0},{   0,ls/2,   0},
        {   0,   0,ls/2},{   0,ls/2,ls/2}
    }; 
    std::vector<CoordinateVector> node_first_coordinate_8 = {
        {   0,   0,   0}, {ls/2,   0,   0}, //z - direction front x{0}{1}
        {   0,ls/2,   0}, {ls/2,ls/2,   0}, //                    x{2}{3}
                                            //                      y  y 
        {   0,   0,ls/2}, {ls/2,   0,ls/2}, //z - direction back  x{4}{5}
        {   0,ls/2,ls/2}, {ls/2,ls/2,ls/2}  //                    x{6}{7}
                                            //                      y  y 
    };

    std::map<int, std::vector<CoordinateVector>> node_first_coordinate = {
        { 1, node_first_coordinate_1 },
        { 2, node_first_coordinate_2 },
        { 4, node_first_coordinate_4 },
        { 8, node_first_coordinate_8 }
    };
};

class TestLattice : public SnapshotStruct {
  public:
    lattice_struct lat = lattice;
    int num_nodes;
    int my_rank;
    int total_lattice_size = lattice_size * lattice_size * lattice_size;
    std::string snapshot_dir = "./snapshots/test_lattice";

    TestLattice() {
        MPI_Comm_size(lattice->mpi_comm_lat, &num_nodes);
        MPI_Comm_rank(lattice->mpi_comm_lat, &my_rank);

        if (num_nodes != 2 && num_nodes != 4 && num_nodes != 8) {
            hila::out0 << "Incorrect number of ranks exiting \n";
            exit(0);
        }
    }

    template<typename T>
    void generate_snapshot(std::vector<T> data, std::string filename) {
        Snapshot temp_snapshot(data, filename);
        temp_snapshot.write_to_file();
    }

    std::vector<int> adjacent_nodes(int direction=1) {
        std::vector<int> nodes_list;
        foralldir (d)
            nodes_list.push_back(lattice.mynode.nn[d]);
        foralldir (d)
            nodes_list.push_back(lattice.mynode.nn[-d]);
        return nodes_list;
    }

    std::vector<int> nn_comminfo_rank() {
        std::vector<int> nodes_list;
        foralldir(d)
            nodes_list.push_back(lattice.nn_comminfo[ d].to_node.rank);
        foralldir(d)
            nodes_list.push_back(lattice.nn_comminfo[-d].to_node.rank);
        return nodes_list;
    }

    std::vector<int> nn_comminfo_n_sites(Parity par) {
        std::vector<int> n_sites_list;
        foralldir(d)
            n_sites_list.push_back(lattice.nn_comminfo[ d].from_node.n_sites(par));
        foralldir(d)
            n_sites_list.push_back(lattice.nn_comminfo[-d].from_node.n_sites(par));
        return n_sites_list;
    }

    std::vector<int> nn_comminfo_offset(Parity par) {
        std::vector<int> offset_list;
        foralldir(d)
            offset_list.push_back(lattice.nn_comminfo[ d].from_node.offset(par));
        foralldir(d)
            offset_list.push_back(lattice.nn_comminfo[-d].from_node.offset(par));
        return offset_list;
    }

};

TEST_CASE_METHOD(TestLattice, "Lattice routines", "[MPI][.]") {
    SECTION("Test size routiens") {
        REQUIRE(lattice.volume() == lattice_size * lattice_size * lattice_size);
        REQUIRE(lattice.size(e_x) == lattice_size);
        REQUIRE(lattice.size(e_y) == lattice_size);
        REQUIRE(lattice.size(e_z) == lattice_size);
        REQUIRE(lattice.size(0) == lattice_size);
        REQUIRE(lattice.size(1) == lattice_size);
        REQUIRE(lattice.size(2) == lattice_size);
        REQUIRE(lattice.size() == CoordinateVector({lattice_size,lattice_size,lattice_size}));
    }
    SECTION("Test coordinate and looping methods") {
        std::map<int,CoordinateVector> local_coordinate_dict = {
            {1,CoordinateVector{127,127,127}},
            {2,CoordinateVector{127,127, 63}},
            {4,CoordinateVector{127, 63, 63}},
            {8,CoordinateVector{ 63, 63, 63}}
        };
        CoordinateVector first_coordinate = node_first_coordinate[num_nodes][my_rank];
        REQUIRE(lattice.coordinates(0)==first_coordinate);
        REQUIRE(lattice.coordinate(0,e_x)==first_coordinate[e_x]);
        REQUIRE(lattice.coordinate(0,e_y)==first_coordinate[e_y]);
        REQUIRE(lattice.coordinate(0,e_z)==first_coordinate[e_z]);
        REQUIRE(lattice.site_parity(0) == Parity(EVEN));
        REQUIRE(lattice.site_parity(lattice.mynode.volume/2) == Parity(ODD));
        REQUIRE(lattice.loop_begin(ODD) == std::round(total_lattice_size / num_nodes / 2));
        REQUIRE(lattice.loop_begin(EVEN) == 0);
        REQUIRE(lattice.loop_end(ODD) == std::round(total_lattice_size / num_nodes));
        REQUIRE(lattice.loop_end(EVEN) == std::round(total_lattice_size / num_nodes / 2));
        REQUIRE(lattice.global_coordinates(0) == CoordinateVector({0,0,0}));
        REQUIRE(lattice.global_coordinates(total_lattice_size-1) == CoordinateVector({ls-1,ls-1,ls-1}));
        REQUIRE(lattice.local_coordinates(0) == CoordinateVector(0));
        REQUIRE(lattice.local_coordinates(total_lattice_size/num_nodes-1) == local_coordinate_dict[num_nodes]);
    }
}

TEST_CASE_METHOD(TestLattice, "Node information", "[MPI][.]") {
    SECTION("Node specific information: struct node_struct") {
        INFO("Test that MPI node recognizes correct rank")

        REQUIRE(lattice.mynode.rank == my_rank);
        REQUIRE(lattice.node_rank() == my_rank);
        REQUIRE(lattice.n_nodes() == num_nodes);

        GIVEN("Total lattice size " << total_lattice_size) {
            THEN("Each node should have a site total of total_lattice_size/number_of_nodes "
                "and even/odd site total of total_lattice_size/number_of_nodes/2")
            REQUIRE(lattice.mynode.volume == lattice.mynode.volume());

            REQUIRE(lattice.mynode.volume == std::round(total_lattice_size / num_nodes));
            REQUIRE(lattice.mynode.evensites == std::round(total_lattice_size / num_nodes / 2));
            REQUIRE(lattice.mynode.oddsites == std::round(total_lattice_size / num_nodes / 2));
        }
    }
    SECTION("Shared node information: struct allnodes") {

        REQUIRE(lattice->nodes.number == num_nodes);

        std::map<int, CoordinateVector> n_division_for_nodes = {
            { 1, CoordinateVector({1,1,1}) },
            { 2, CoordinateVector({1,1,2}) },
            { 4, CoordinateVector({1,2,2}) },
            { 8, CoordinateVector({2,2,2}) }
        };
        std::map<int, CoordinateVector> max_size_for_node = {
            { 1, CoordinateVector({ls  ,ls  ,ls  }) },
            { 2, CoordinateVector({ls  ,ls  ,ls/2}) },
            { 4, CoordinateVector({ls  ,ls/2,ls/2}) },
            { 8, CoordinateVector({ls/2,ls/2,ls/2}) }
        };

        REQUIRE(lattice->nodes.n_divisions == n_division_for_nodes[num_nodes]);
        REQUIRE(lattice->nodes.max_size == max_size_for_node[num_nodes]);

    }
    SECTION("Communication node information: struct comm_node_struct") {
        std::string filename = "/test_nn_comminfo_indicies_z_dir_" + std::to_string(num_nodes) + 
                            "_nodes_" + std::to_string(lattice_size) + "_size.txt";
        Snapshot<int> dummy_data_halo_indicies(snapshot_dir + filename);
        std::vector<int> sitelist_snapshot = dummy_data_halo_indicies.read_from_file();

        SECTION("Testing nearest node communication ranks") {  

            REQUIRE(nn_comminfo_rank() == adjacent_nodes_dict[num_nodes][my_rank]); 
        }
        SECTION("Testing indicies of halo coordinates. Choosing z direction, since there will always be a z-dir split") {
            std::vector<int> coord_indicies;
            int lower=0;
            int upper = lattice.nn_comminfo[e_z].to_node.sites;
            for (int i = lower; i < upper; i++) {
                coord_indicies.push_back(lattice.nn_comminfo[e_z].to_node.site_index(i,ALL));
            }

            REQUIRE(coord_indicies == sitelist_snapshot);
        }
        SECTION("Testing sitelist in comm_node_struct") {
            int size;
            const unsigned *sitelist = lattice.nn_comminfo[e_z].to_node.get_sitelist(ALL,size);
            std::vector<int> sitelist_vector(sitelist, sitelist + size);
            REQUIRE(sitelist_vector == sitelist_snapshot);
        }
        SECTION("Testing n_sites and offset functions") {
            std::map<int, std::vector<int>> n_sites_dict = {
                {2, {    0,    0,16384,    0,    0,16384}},
                {4, {    0, 8192, 8192,    0, 8192, 8192}},
                {8, { 4096, 4096, 4096, 4096, 4096, 4096}}
            };
            std::map<int, std::vector<int>> n_sites_dict_ODD_EVEN = {
                {2, {    0,    0, 8192,    0,    0, 8192}},
                {4, {    0, 4096, 4096,    0, 4096, 4096}},
                {8, { 2048, 2048, 2048, 2048, 2048, 2048}}
            };
            std::map<int, std::vector<int>> offset_dict = {
                {2, { 1048576, 1048576, 1048576, 1081344, 1081344, 1064960}},
                {4, { 524288, 524288, 532480, 557056, 548864, 540672 }},
                {8, { 262144, 266240, 270336, 282624, 278528, 274432}}
            };
            std::map<int, std::vector<int>> offset_dict_ODD = {
                {2, { 1048576 , 1048576 , 1056768 , 1081344, 1081344 , 1073152  }},
                {4, { 524288, 528384, 536576, 557056, 552960, 544768 }},
                {8, { 264192 , 268288 , 272384 , 284672 , 280576 , 276480  }}
            };
            REQUIRE(nn_comminfo_n_sites(ALL) == n_sites_dict[num_nodes]);
            REQUIRE(nn_comminfo_n_sites(EVEN) == n_sites_dict_ODD_EVEN[num_nodes]);
            REQUIRE(nn_comminfo_n_sites(ODD) == n_sites_dict_ODD_EVEN[num_nodes]);

            REQUIRE(nn_comminfo_offset(ALL) == offset_dict[num_nodes]);
            REQUIRE(nn_comminfo_offset(EVEN) == offset_dict[num_nodes]);
            REQUIRE(nn_comminfo_offset(ODD) == offset_dict_ODD[num_nodes]);
        }
    }
}

TEST_CASE_METHOD(TestLattice, "Lattice split", "[MPI][.]") {

    if (num_nodes == 2) {
        GIVEN("A 3D lattice split into 2 nodes") {
            THEN("The lattice will be split in the z direction and adjacent node in "
                 "y-x-direction will be itself") {                
                REQUIRE(adjacent_nodes() == adjacent_nodes_dict[num_nodes][my_rank]);
                REQUIRE(lattice.mynode.min == node_first_coordinate[num_nodes][my_rank]);           
                REQUIRE(lattice.mynode.size == CoordinateVector({ls,ls,ls/2}));
            }
        }
    }
    if (num_nodes == 4) {
        GIVEN("A 3D lattice split into 4 nodes") {
            THEN("The lattice will be split in the y-z direction and adjacent node in "
                 "x-direction will be itself") {
                REQUIRE(adjacent_nodes() == adjacent_nodes_dict[num_nodes][my_rank]);
                REQUIRE(lattice.mynode.min == node_first_coordinate[num_nodes][my_rank]);
                REQUIRE(lattice.mynode.size == CoordinateVector({ls,ls/2,ls/2}));
            }
        }
    }
    if (num_nodes == 8) {
        GIVEN("A 3D lattice split into 8 nodes") {
            THEN("The lattice will be split in all direction evenly") {

                REQUIRE(adjacent_nodes() == adjacent_nodes_dict[num_nodes][my_rank]);
                REQUIRE(lattice.mynode.min == node_first_coordinate[num_nodes][my_rank]);
                REQUIRE(lattice.mynode.size == CoordinateVector({ls/2,ls/2,ls/2}));
            }
        }
    }
}