#include "hila.h"
#include "catch.hpp"
#include "catch_main.hpp"
#include <vector>
#include <map>

class TestLattice {
  public:
    lattice_struct lat = *lattice;
    int num_nodes;
    int my_rank;
    int total_lattice_size = lattice_size * lattice_size * lattice_size;
    TestLattice() {
        MPI_Comm_size(lattice->mpi_comm_lat, &num_nodes);
        MPI_Comm_rank(lattice->mpi_comm_lat, &my_rank);

        if (num_nodes != 2 && num_nodes != 4 && num_nodes != 16) {
            output0 << "Incorrect number of ranks exiting \n";
            exit(0);
        }
    }
    std::vector<int> adjacent_nodes() {
        std::vector<int> nodes_list;
        foralldir (d)
            nodes_list.push_back(lattice->mynode.nn[d]);
        return nodes_list;
    }

    void get_nearest_node_each_direction() {
        if(my_rank==0) {
            std::cout << lattice->nn_comminfo[0].to_node.rank << '\n';
            std::cout << lattice->nn_comminfo[1].to_node.rank << '\n';
            std::cout << lattice->nn_comminfo[2].to_node.rank << '\n';
            std::cout << lattice->nn_comminfo[3].to_node.rank << '\n';
            std::cout << lattice->nn_comminfo[4].to_node.rank << '\n';
            std::cout << lattice->nn_comminfo[5].to_node.rank << '\n';
            std::cout << lattice->nn_comminfo[6].to_node.rank << '\n';

        };
    }

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
            REQUIRE(lattice->mynode.sites == total_lattice_size / num_nodes);
            REQUIRE(lattice->mynode.evensites == total_lattice_size / num_nodes / 2);
            REQUIRE(lattice->mynode.oddsites == total_lattice_size / num_nodes / 2);
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
            { 1, CoordinateVector({128,128,128}) },
            { 2, CoordinateVector({128,128, 64}) },
            { 4, CoordinateVector({128, 64, 64}) },
            { 8, CoordinateVector({ 64, 64, 64}) }
        };
        REQUIRE(lattice->nodes.n_divisions == n_division_for_nodes[num_nodes]);
        REQUIRE(lattice->nodes.max_size == max_size_for_node[num_nodes]);

    }
    SECTION("Communication node information") {

        // if(num_nodes==4){
        //     std::vector<std::vector<int>> adjacent_nodes_list = {
        //             {0, 1, 2}, {1, 0, 3},
        //             {2, 3, 0}, {3, 2, 1}}; 
        //     std::vector<int> nn_in_each_dir = get_nearest_node_each_direction();
        // }
        get_nearest_node_each_direction();
    }
}

TEST_CASE_METHOD(TestLattice, "Lattice split", "[MPI][.]") {

    if (num_nodes == 2) {
        GIVEN("A 3D lattice split into 2 nodes") {
            THEN("The lattice will be split in the z direction and adjacent node in "
                 "y-x-direction will be itself") {
                std::vector<std::vector<int>> adjacent_nodes_list = {{0, 0, 1}, {1, 1, 0}};
                std::vector<CoordinateVector> node_first_coordinate = {{0, 0, 0}, {0, 0, 64}};
                
                REQUIRE(adjacent_nodes() == adjacent_nodes_list[my_rank]);
                REQUIRE(lattice->mynode.min == node_first_coordinate[my_rank]);           
                REQUIRE(lattice->mynode.size == CoordinateVector({128,128,64}));
            }
        }
    }
    if (num_nodes == 4) {
        GIVEN("A 3D lattice split into 4 nodes") {
            THEN("The lattice will be split in the y-z direction and adjacent node in "
                 "x-direction will be itself") {
                std::vector<std::vector<int>> adjacent_nodes_list = {
                    {0, 1, 2}, {1, 0, 3},
                    {2, 3, 0}, {3, 2, 1}}; 

                std::vector<CoordinateVector> node_first_coordinate = {
                    {  0,  0,  0},{  0, 64,  0},
                    {  0,  0, 64},{  0, 64, 64}}; 

                REQUIRE(adjacent_nodes() == adjacent_nodes_list[my_rank]);
                REQUIRE(lattice->mynode.min == node_first_coordinate[my_rank]);
                REQUIRE(lattice->mynode.size == CoordinateVector({128,64,64}));
            }
        }
    }
    if (num_nodes == 8) {
        GIVEN("A 3D lattice split into 8 nodes") {
            THEN("The lattice will be split in all direction evenly") {
                std::vector<std::vector<int>> adjacent_nodes_list = {
                    {1, 2, 4}, {0, 3, 5}, //z - direction front x{0}{1}
                    {3, 0, 6}, {2, 1, 7}, //                    x{2}{3}
                                          //                      y  y 
                    {5, 6, 0}, {4, 7, 1}, //z - direction back  x{4}{5}
                    {7, 4, 2}, {6, 5, 3}  //                    x{6}{7}
                                          //                      y  y 
                };

                std::vector<CoordinateVector> node_first_coordinate = {
                    {  0,  0,  0}, { 64,  0,  0}, //z - direction front x{0}{1}
                    {  0, 64,  0}, { 64, 64,  0}, //                    x{2}{3}
                                                //                      y  y 
                    {  0,  0, 64}, { 64,  0, 64}, //z - direction back  x{4}{5}
                    {  0, 64, 64}, { 64, 64, 64}  //                    x{6}{7}
                                                //                      y  y 
                };

                REQUIRE(adjacent_nodes() == adjacent_nodes_list[my_rank]);
                REQUIRE(lattice->mynode.min == node_first_coordinate[my_rank]);
                REQUIRE(lattice->mynode.size == CoordinateVector({64,64,64}));
            }
        }
    }
}