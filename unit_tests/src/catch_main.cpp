#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "hila.h"

int main( int argc, char* argv[] ) {
    // global setup...
    char p0[] = "catch_main";
    char *fake_argv[] = { p0, NULL };
    hila::initialize(1, fake_argv);
    lattice->setup({128,128,128}); 
    hila::seed_random(0);

    int result = Catch::Session().run(argc,argv);
    hila::finishrun();
    // global clean-up...

    return result;

}