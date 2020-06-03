/***************************************************************
 *  mersenne.h   
 *  for the inline version of the mersenne generator
 */

#define MERSENNE_N 624

extern int mersenne_i;
extern double mersenne_array[MERSENNE_N];

#define mersenne() ( mersenne_i > 0 ? mersenne_array[--mersenne_i] : mersenne_generate() )

void seed_mersenne(long);
double mersenne_generate();
