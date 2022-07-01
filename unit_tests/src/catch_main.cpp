#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "hila.h"



int main( int argc, char* argv[] ) {
    // global setup...
    char *hila_argument[] = {NULL};
    hila::initialize(0, hila_argument);
    lattice->setup({128,128,128}); 
    hila::seed_random(0);

    int result = Catch::Session().run( argc, argv);
    hila::finishrun();
    // global clean-up...

    return result;

}