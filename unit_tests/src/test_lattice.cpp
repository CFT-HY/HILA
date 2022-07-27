#include "hila.h"
#include "catch.hpp"
#include "catch_main.hpp"
#include "read_write_snapshots.hpp"
#include <vector>
#include <map>
#include <cmath>

class TestLattice {
  public:
    lattice_struct lat = *lattice;
    int ls = lattice_size;
    int num_nodes;
    int my_rank;
    int total_lattice_size = lattice_size * lattice_size * lattice_size;
    std::string snapshot_dir = "./snapshots/test_lattice";

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

    std::map<int, std::vector<std::vector<int>>> adjacent_nodes_list = {
        { 1, adjacent_nodes_list_1 },
        { 2, adjacent_nodes_list_2 },
        { 4, adjacent_nodes_list_4 },
        { 8, adjacent_nodes_list_8 }
    };

    TestLattice() {
        MPI_Comm_size(lattice->mpi_comm_lat, &num_nodes);
        MPI_Comm_rank(lattice->mpi_comm_lat, &my_rank);
        if (num_nodes != 2 && num_nodes != 4 && num_nodes != 8) {
            output0 << "Incorrect number of ranks exiting \n";
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
            nodes_list.push_back(lattice->mynode.nn[d]);
        foralldir (d)
            nodes_list.push_back(lattice->mynode.nn[-d]);
        return nodes_list;
    }

    std::vector<int> nn_comminfo_rank() {
        std::vector<int> nodes_list;
        foralldir(d)
            nodes_list.push_back(lattice->nn_comminfo[ d].to_node.rank);
        foralldir(d)
            nodes_list.push_back(lattice->nn_comminfo[-d].to_node.rank);
        return nodes_list;
    }

    void test() {
        std::vector<int> coord_indicies;
        int lower=0;
        int upper = lattice->nn_comminfo[3].to_node.sites;
        // int lower=  2048+1024;
        // int upper= 2048+1024+32;
        for (int i = lower; i < upper; i++) {
            coord_indicies.push_back(lattice->nn_comminfo[2].to_node.site_index(i,ALL));
        }
        std::string filename = "/test_nn_comminfo_indicies_z_dir_" + std::to_string(num_nodes) + 
                               "_nodes_" + std::to_string(lattice_size) + "_size.txt";

        generate_snapshot(coord_indicies,snapshot_dir + filename);
    }

    std::vector<int> generate_halo_indicies(int dir) {
        std::vector<int> coord_indicies;
        int lower=0;
        int upper = lattice->nn_comminfo[dir].to_node.sites;
        // int lower=  2048+1024;
        // int upper= 2048+1024+32;
        for (int i = lower; i < upper; i++) {
            coord_indicies.push_back(lattice->nn_comminfo[dir].to_node.site_index(i,ALL));
        }
        return coord_indicies;
    }

    // void print() {
    //     if(my_rank == 1) {
    //        // std::cout << lattice->coordinates(lattice->nn_comminfo[0].to_node.site_index(1,ALL)) << '\n';
    //         std::cout << lattice->nn_comminfo[0].to_node.rank << '\n';

    
    //     }     
    //     if(my_rank==0) {
    //         //for (size_t i = 0; i < lattice->nn_comminfo[e_z].from_node.sites; i++)
    //         int lower=0;
    //         int upper = lattice->nn_comminfo[3].to_node.sites;
    //         // int lower=  2048+1024;
    //         // int upper= 2048+1024+32;
    //         for (int i = lower; i < upper; i++)
    //         {
    //             output0 << i << " " << lattice->nn_comminfo[3].to_node.site_index(i,ALL) << '\n'
    //                                 << lattice->coordinates(lattice->nn_comminfo[3].to_node.site_index(i,ALL)+0) << "\n"
    //                                 << lattice->coordinates(lattice->nn_comminfo[3].to_node.site_index(i,ALL)+1) << "\n" 
    //                                 << lattice->coordinates(lattice->nn_comminfo[3].to_node.site_index(i,ALL)+2) << '\n'
    //                                 << lattice->coordinates(lattice->nn_comminfo[3].to_node.site_index(i,ALL)+3) << '\n'
    //                                 << lattice->coordinates(lattice->nn_comminfo[3].to_node.site_index(i,ALL)+4) << '\n'
    //                                 << lattice->coordinates(lattice->nn_comminfo[3].to_node.site_index(i,ALL)+5) << '\n'
    //                                 << lattice->coordinates(lattice->nn_comminfo[3].to_node.site_index(i,ALL)+6) << '\n'
    //                                 << lattice->coordinates(lattice->nn_comminfo[3].to_node.site_index(i,ALL)+7) << '\n';
    //         }
    //     }
    //}

};

TEST_CASE("Size", "[MPI][.]") {
    REQUIRE(lattice->volume() == lattice_size * lattice_size * lattice_size);
}

TEST_CASE_METHOD(TestLattice, "Node information", "[MPI][.]") {
    SECTION("Node specific information") {
        INFO("Test that MPI node recognizes correct rank")
        REQUIRE(lattice->mynode.rank == my_rank);
        GIVEN("Total lattice size " << total_lattice_size) {
            THEN("Each node should have a site total of total_lattice_size/number_of_nodes "
                "and even/odd site total of total_lattice_size/number_of_nodes/2")
            REQUIRE(lattice->mynode.sites == lattice->mynode.volume());
            REQUIRE(lattice->mynode.sites == std::round(total_lattice_size / num_nodes));
            REQUIRE(lattice->mynode.evensites == std::round(total_lattice_size / num_nodes / 2));
            REQUIRE(lattice->mynode.oddsites == std::round(total_lattice_size / num_nodes / 2));
        }
    }
    SECTION("Shared node information") {
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
    SECTION("Communication node information") {
        INFO("Testing nearest node communication ranks") {  
            REQUIRE(nn_comminfo_rank() == adjacent_nodes_list[num_nodes][my_rank]); 
        }
        INFO("Testing indicies of halo coordinates. Choosing z direction, since there will always be a z-dir split") {
            std::vector<int> coord_indicies;
            int lower=0;
            int upper = lattice->nn_comminfo[e_z].to_node.sites;
            for (int i = lower; i < upper; i++) {
                coord_indicies.push_back(lattice->nn_comminfo[e_z].to_node.site_index(i,ALL));
            }
            std::string filename = "/test_nn_comminfo_indicies_z_dir_" + std::to_string(num_nodes) + 
                                "_nodes_" + std::to_string(lattice_size) + "_size.txt";
            Snapshot<int> dummy_data_halo_indicies(snapshot_dir + filename);
            REQUIRE(coord_indicies == dummy_data_halo_indicies.read_from_file());
        }
    }
}

TEST_CASE_METHOD(TestLattice, "Lattice split", "[MPI][.]") {

    if (num_nodes == 2) {
        GIVEN("A 3D lattice split into 2 nodes") {
            THEN("The lattice will be split in the z direction and adjacent node in "
                 "y-x-direction will be itself") {
                std::vector<CoordinateVector> node_first_coordinate = {{0, 0, 0}, {0, 0, ls/2}};
                
                REQUIRE(adjacent_nodes() == adjacent_nodes_list[num_nodes][my_rank]);
                REQUIRE(lattice->mynode.min == node_first_coordinate[my_rank]);           
                REQUIRE(lattice->mynode.size == CoordinateVector({ls,ls,ls/2}));
            }
        }
    }
    if (num_nodes == 4) {
        GIVEN("A 3D lattice split into 4 nodes") {
            THEN("The lattice will be split in the y-z direction and adjacent node in "
                 "x-direction will be itself") {

                std::vector<CoordinateVector> node_first_coordinate = {
                    {   0,   0,   0},{   0,ls/2,   0},
                    {   0,   0,ls/2},{   0,ls/2,ls/2}}; 

                REQUIRE(adjacent_nodes() == adjacent_nodes_list[num_nodes][my_rank]);
                REQUIRE(lattice->mynode.min == node_first_coordinate[my_rank]);
                REQUIRE(lattice->mynode.size == CoordinateVector({ls,ls/2,ls/2}));
            }
        }
    }
    if (num_nodes == 8) {
        GIVEN("A 3D lattice split into 8 nodes") {
            THEN("The lattice will be split in all direction evenly") {

                std::vector<CoordinateVector> node_first_coordinate = {
                    {   0,   0,   0}, {ls/2,   0,   0}, //z - direction front x{0}{1}
                    {   0,ls/2,   0}, {ls/2,ls/2,   0}, //                    x{2}{3}
                                                        //                      y  y 
                    {   0,   0,ls/2}, {ls/2,   0,ls/2}, //z - direction back  x{4}{5}
                    {   0,ls/2,ls/2}, {ls/2,ls/2,ls/2}  //                    x{6}{7}
                                                        //                      y  y 
                };

                REQUIRE(adjacent_nodes() == adjacent_nodes_list[num_nodes][my_rank]);
                REQUIRE(lattice->mynode.min == node_first_coordinate[my_rank]);
                REQUIRE(lattice->mynode.size == CoordinateVector({ls/2,ls/2,ls/2}));
            }
        }
    }
}