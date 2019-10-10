#include "SUN.h"

// Direct output to stdout
std::ostream &hila::output = std::cout;
std::ostream &output = std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;


// Define some parameters for the simulation
double beta = 8;
int n_measurements=100;
int n_updates_per_measurement=10;
long seed = 123456;
int NX=8, NY=8, NZ=8, NT=8;
int VOLUME = NX*NY*NZ*NT;




field<matrix<N,N,cmplx<double>>>
calc_staples( field<matrix<N,N,cmplx<double>>> U[NDIM], direction dir)
{
  /* Calculate the sum of staples connected to links in direction
   * dir 
   */
  field<matrix<N,N,cmplx<double>>> down_staple, staple_sum;
  staple_sum[ALL] = 0;
  foralldir(d2){
    direction dir2 = (direction)d2;
    //Calculate the down side staple.
    //This will be communicated up.
    down_staple[ALL] = U[dir2][X].conjugate()
                     * U[dir][X]
                     * U[dir2][X+dir];
    // Forward staple
    staple_sum[ALL] += U[dir2][X+dir]
                     * U[dir][X+dir2].conjugate()
                     * U[dir2][X].conjugate();
    // Add the two staples together
    staple_sum[ALL] += down_staple[X - dir2];
  }
  return staple_sum;
}

 
template<typename T>
void update(
  T &U, const T &staple,
  double beta
){
  monte( U, staple, beta );
}

void update(
  matrix<2,2,cmplx<double>> &U,
  const matrix<2,2,cmplx<double>> &staple,
  double beta
){
  matrix<2,2,cmplx<double>> temp = -beta*staple;
  KennedyPendleton( U, temp );
}


int main()
{
  // Basic setup
  lattice->setup( NX, NY, NZ, NT );
  // Define a field
  field<matrix<N,N,cmplx<double>>> U[NDIM];

  seed_mersenne( seed );

  /* "Warm up" the rng generator */
  for( int i=0; i<543210; i++ ) mersenne();
  
  // Set to 1
  foralldir(d) {
    U[d][ALL] = 1;
  }

  // Run update-measure loop
  for( int i=0; i<n_measurements; i++ ){

    // Run a number of updates
    for(int j=0; j<n_updates_per_measurement; j++){
      foralldir(d) {
        // update direction dir
        direction dir = (direction)d;
        // First we need the staple sum
        field<matrix<N,N,cmplx<double>>> staple = calc_staples(U, dir);

        // Now update, first even then odd
        parity p = EVEN;
        for( int par=0; par < 2; par++ ){
          onsites(p){
            update( U[dir][X], staple[X], beta );
          }
          p = opp_parity(p);
        }
      }
    }

    // Measure plaquette
    double Plaq=0;
    foralldir(d1) foralldir(d2) if(d1 != d2){
      direction dir1 = (direction)d1, dir2 = (direction)d2;
      onsites(ALL){
        matrix<N,N,cmplx<double>> temp;
        temp =  U[dir1][X] * U[dir2][X+dir1];
        temp *= U[dir1][X+dir2].conjugate();
        temp *= U[dir2][X].conjugate();
        Plaq += 1-temp.trace().re/N;
      }
    }
    printf("Plaquette %f\n", Plaq/(VOLUME*NDIM*(NDIM-1)));
  }
  
  return 0;
}
