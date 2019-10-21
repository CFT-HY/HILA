#include "bench.h"

#define N 3
#define n_runs_multiplier 1


void dirac_naive(
    field<matrix<N,N, cmplx<double>> > &matrix1,
    field<matrix<1,N, cmplx<double>> > &vector1,
    field<matrix<1,N, cmplx<double>> > &vector2);



int main(){
    int n_runs;
    double msecs;
    clock_t init, end;
    double timing;
    double sum;


    bench_setup();

    field<matrix<N,N, cmplx<double>> > matrix1;
    field<matrix<N,N, cmplx<double>> > matrix2;
    field<matrix<N,N, cmplx<double>> > matrix3;
    field<matrix<1,N, cmplx<double>> > vector1;
    field<matrix<1,N, cmplx<double>> > vector2;

    matrix1[ALL] = 1; 
    matrix2[ALL] = 1;

    // Time MATRIX * MATRIX
    n_runs = n_runs_multiplier*10;

    init = clock();

    for( int i=0; i<n_runs; i++){
        matrix3[ALL] = matrix1[X]*matrix2[X];
    }

    end = clock();

    timing = ((double)(end - init)) *1000.0 / (CLOCKS_PER_SEC) / n_runs;
    printf("Matrix * Matrix: %g ms \n", timing);

    matrix1[ALL] = 1; 
    onsites(ALL){
        for(int i=0; i<N; i++){
            vector1[X].c[0][1]=1;
        }
    }


    // Time VECTOR * MATRIX
    n_runs = n_runs_multiplier*10;

    init = clock();

    for( int i=0; i<n_runs; i++){
        vector1[ALL] = vector1[X]*matrix1[X];
    }

    end = clock();

    timing = ((double)(end - init)) *1000 / (CLOCKS_PER_SEC) / n_runs;
    output0 << "Vector * Matrix: " << timing << " ms \n";


    // Time VECTOR * MATRIX
    n_runs = n_runs_multiplier*1;
    
    sum=0;
    init = clock();

    for( int i=0; i<n_runs; i++){
        onsites(ALL){
            sum += vector1[X].norm_sq();
        }
    }

    end = clock();

    timing = ((double)(end - init)) *1000 / (CLOCKS_PER_SEC) / n_runs;
    output0 << "Vector square sum: " << timing << " ms \n";




    // Time naive Dirac operator 
    n_runs = n_runs_multiplier;

    init = clock();

    for( int i=0; i<n_runs; i++){
        dirac_naive(matrix1, vector1, vector2);
    }

    end = clock();

    timing = ((double)(end - init)) *1000.0 / (CLOCKS_PER_SEC) / n_runs;
    printf("Dirac: %g ms \n", timing);


    // Conjugate gradient step 
    n_runs = n_runs_multiplier;

    {
        field<matrix<1,N, cmplx<double>> > r, rnew, p, Dp;

        init = clock();

        for( int i=0; i<n_runs; i++){
            r[ALL] = vector1[X];
            p[ALL] = vector1[X];
            onsites(ALL){
                 for(int i=0; i<N; i++){
                    vector2[X].c[0][i] = 0;
                 }
            }
            
            dirac_naive(matrix1, p, Dp);

            double pDDp = 0;
            double rr = 0;
            onsites(ALL){
                double rr_temp = r[X].norm_sq();
                double pDDp_temp = Dp[X].norm_sq();
                pDDp += pDDp_temp;
                rr += rr_temp;
            }

            double alpha = rr / pDDp;

            vector2[ALL] = r[X] + alpha*p[X];
            r[ALL] = r[X] - alpha*Dp[X];

            double rrnew = 0;
            onsites(ALL){
                double rr_temp = r[X].norm_sq();
                rrnew += rr_temp;
            }

            double beta = rrnew/rr;
            
            p[ALL] = r[X] + beta*p[X];

        }
        
        end = clock();
    }

    timing = ((double)(end - init)) *1000.0 / (CLOCKS_PER_SEC) / n_runs;
    printf("CG: %g ms \n", timing);


    return 0;
}




void dirac_naive(
    field<matrix<N,N, cmplx<double>> > &matrix1,
    field<matrix<1,N, cmplx<double>> > &vector1,
    field<matrix<1,N, cmplx<double>> > &vector2){
    static field<double> eta[NDIM]; // The staggered phase
    static bool initialized = false;

    double mass = 0.1;

    if(!initialized){
        foralldir(d){
            onsites(ALL){
                location l = coordinates(X);
                int sumcoord = 0;
                for(int d2=0;d2<d;d2++){
                    sumcoord += l[d];
                }
                if( sumcoord %2 ){
                    eta[d][X] = 1;
                } else {
                    eta[d][X] =-1;
                }
            }
        }
        initialized = true;
    }
    

    vector2[ALL] = mass * vector1[X];

    foralldir(d){
        // Positive directions: get the vector and multiply by matrix stored here
        direction dir = (direction)d;
        vector2[ALL] += 0.5*eta[d][X]*vector1[X+dir]*matrix1[X];

        // Negative directions: get both form neighbour
        dir = opp_dir( (direction)d );
        direction dir2 = opp_dir( (direction)d );
        vector2[ALL] -= 0.5*eta[d][X]*vector1[X+dir]*matrix1[X+dir2].conjugate();
    }
}


