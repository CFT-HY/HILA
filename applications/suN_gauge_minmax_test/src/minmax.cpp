#include "hila.h"

#include <random>


static_assert(NDIM == 4, "NDIM must be 4 here");

using ftype=float;
using mygroup=SU<NCOLOR,ftype>;

int main(int argc, char *argv[]) {

    // initialize system
    hila::initialize(argc,argv);

    int lsize=12;

    lattice.setup({lsize, lsize, lsize, lsize});

    // Random numbers are used here - use time to seed
    hila::seed_random(0);

    // lattice field
    mygroup S=0;
    S.random();
    GaugeField<mygroup> U;
    Field<ftype> maxabsU;
    foralldir(d) {
        onsites(ALL) U[d][X]=S;

        U[d][{(int)d,1,1,1}]=1.0;
    }

    ftype val1,val2,val3,val4,val5,val6;
    CoordinateVector loc1, loc2, loc3, loc4, loc5, loc6;
    foralldir(d) {
        onsites(ALL) maxabsU[X]=U[d][X].max_abs();

        val1=maxabsU.min(ODD, loc1);
        val2=maxabsU.min(EVEN,loc2);
        val3=maxabsU.max(ODD, loc3);
        val4=maxabsU.max(EVEN,loc4);
        val5=maxabsU.min(ALL, loc5);
        val6=maxabsU.max(ALL, loc6);

        hila::out0<<"dir "<<(int)d<<":\n";
        hila::out0<<"Min value of ODD sites "<<val1<<" at location: "<<loc1<<'\n';
        hila::out0<<"Min value of EVEN sites  "<<val2<<" at location: "<<loc2<<'\n';
        hila::out0<<"Max value of ODD sites "<<val3<<" at location: "<<loc3<<'\n';
        hila::out0<<"Max value of EVEN sites  "<<val4<<" at location: "<<loc4<<'\n';
        hila::out0<<"Min value of ALL sites  "<<val5<<" at location: "<<loc5<<'\n';
        hila::out0<<"Max value of ALL sites  "<<val6<<" at location: "<<loc6<<'\n';
        hila::out0<<"\n\n";

    }
    


    //hila::out0 << "Reduction test " << reduced << "\n";
    hila::finishrun();
    return 0;
}
