#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "hila.h"
#include "catch_main.hpp"
//#include <mpi.h>

int main(int argc, char *argv[]) {
    char p0[] = "catch_main";
    char *fake_argv[] = {p0, NULL};

    hila::initialize(1, fake_argv);
    //lattice_size=128 defined in catch_main.hpp 
    lattice->setup({lattice_size,lattice_size,lattice_size});

    hila::seed_random(0);

    int result = Catch::Session().run(argc, argv);

    hila::finishrun();

    return result;
}


