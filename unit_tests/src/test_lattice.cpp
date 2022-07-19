#include "hila.h"
#include "catch.hpp"
#include "catch_main.hpp"
#include <vector>

class TestLattice {
  public:
    lattice_struct lat = *lattice;
    int num_nodes;
    int my_rank;
    int total_lattice_size = lattice_size * lattice_size * lattice_size;
    TestLattice() {
        MPI_Comm_size(lattice->mpi_comm_lat, &num_nodes);
        MPI_Comm_rank(lattice->mpi_comm_lat, &my_rank);
    }
    std::vector<int> adjacent_nodes() {
        std::vector<int> nodes_list;
        foralldir (d)
            nodes_list.push_back(lattice->mynode.nn[d]);
        return nodes_list;
    }
};

TEST_CASE("Size", "[MPI][.]") {
    REQUIRE(lattice->volume() == lattice_size * lattice_size * lattice_size);
}

TEST_CASE_METHOD(TestLattice, "Node information", "[MPI][.]") {
    INFO("Test that MPI node recognizes correct rank")
    REQUIRE(lattice->mynode.rank == my_rank);
    GIVEN("Total lattice size" << total_lattice_size) {
        THEN("Each node should have a site total of total_lattice_size/number_of_nodes "
             "and even/odd site total of total_lattice_size/number_of_nodes/2")
        REQUIRE(lattice->mynode.sites == total_lattice_size / num_nodes);
        REQUIRE(lattice->mynode.evensites == total_lattice_size / num_nodes / 2);
        REQUIRE(lattice->mynode.oddsites == total_lattice_size / num_nodes / 2);
    }
    std::cout << lattice->mynode.min - e_x << lattice->mynode.size << " " << my_rank
              << '\n';
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