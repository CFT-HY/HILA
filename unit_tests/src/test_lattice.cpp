#include "hila.h"
#include "catch.hpp"
#include "catch_main.hpp"
#include <vector>
#include <map>
#include "read_write_snapshots.hpp"
class TestLattice {
  public:
    lattice_struct lat = *lattice;
    int num_nodes;
    int my_rank;
    int total_lattice_size = lattice_size * lattice_size * lattice_size;

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
        {2, 3, 0, 2, 3, 0}, {3, 2, 1, 3, 2, 1}}; 

    std::vector<std::vector<int>> adjacent_nodes_list_2 = {{0, 0, 1, 0, 0, 1}, {1, 1, 0, 1, 1, 0}};

    std::map<int,std::vector<std::vector<int>>> adjacent_nodes_list = {
        {2,adjacent_nodes_list_2},
        {4,adjacent_nodes_list_4},
        {8,adjacent_nodes_list_8},
    };

    TestLattice() {
        MPI_Comm_size(lattice->mpi_comm_lat, &num_nodes);
        MPI_Comm_rank(lattice->mpi_comm_lat, &my_rank);

        if (num_nodes != 2 && num_nodes != 4 && num_nodes != 8) {
            output0 << "Incorrect number of ranks exiting \n";
            exit(0);
        }
    }
    std::vector<int> adjacent_nodes() {
        std::vector<int> nodes_list;
        foralldir(d)
            nodes_list.push_back(lattice->mynode.nn[ d]);
        foralldir(d)
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

    void generate_snapshot() {
        
    }
    void test() {
        std::vector<int> test_data = {1,2,3,4,5,6,7};
        std::vector<int> return_data;
        Snapshots test_data_snapshot(test_data,"./src/snapshots/test_data.txt");
        test_data_snapshot.write_to_file();
        test_data_snapshot.clear_data();
        return_data = test_data_snapshot.read_from_file();
        for(int elem : return_data) {
            std::cout << elem << '\n';
        }

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
        REQUIRE(nn_comminfo_rank() == adjacent_nodes_list[num_nodes][my_rank]); 
        if (my_rank == 0) TestLattice::test();
    }
}

TEST_CASE_METHOD(TestLattice, "Lattice split", "[MPI][.]") {

    if (num_nodes == 2) {
        GIVEN("A 3D lattice split into 2 nodes") {
            THEN("The lattice will be split in the z direction and adjacent node in "
                 "y-x-direction will be itself") {
                std::vector<CoordinateVector> node_first_coordinate = {{0, 0, 0}, {0, 0, 64}};
                
                REQUIRE(adjacent_nodes() == adjacent_nodes_list[num_nodes][my_rank]);
                REQUIRE(lattice->mynode.min == node_first_coordinate[my_rank]);           
                REQUIRE(lattice->mynode.size == CoordinateVector({128,128,64}));
            }
        }
    }
    if (num_nodes == 4) {
        GIVEN("A 3D lattice split into 4 nodes") {
            THEN("The lattice will be split in the y-z direction and adjacent node in "
                 "x-direction will be itself") {

                std::vector<CoordinateVector> node_first_coordinate = {
                    {  0,  0,  0},{  0, 64,  0},
                    {  0,  0, 64},{  0, 64, 64}}; 

                REQUIRE(adjacent_nodes() == adjacent_nodes_list[num_nodes][my_rank]);
                REQUIRE(lattice->mynode.min == node_first_coordinate[my_rank]);
                REQUIRE(lattice->mynode.size == CoordinateVector({128,64,64}));
            }
        }
    }
    if (num_nodes == 8) {
        GIVEN("A 3D lattice split into 8 nodes") {
            THEN("The lattice will be split in all direction evenly") {

                std::vector<CoordinateVector> node_first_coordinate = {
                    {  0,  0,  0}, { 64,  0,  0}, //z - direction front x{0}{1}
                    {  0, 64,  0}, { 64, 64,  0}, //                    x{2}{3}
                                                //                      y  y 
                    {  0,  0, 64}, { 64,  0, 64}, //z - direction back  x{4}{5}
                    {  0, 64, 64}, { 64, 64, 64}  //                    x{6}{7}
                                                //                      y  y 
                };

                REQUIRE(adjacent_nodes() == adjacent_nodes_list[num_nodes][my_rank]);
                REQUIRE(lattice->mynode.min == node_first_coordinate[my_rank]);
                REQUIRE(lattice->mynode.size == CoordinateVector({64,64,64}));
            }
        }
    }
}