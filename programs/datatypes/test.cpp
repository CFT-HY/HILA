#include "sun.h"
#include "../plumbing/defs.h"

int main(){
	seed_mersenne(123);
	SU2<float> a, b, c;
	SU_vector<3, float> A, B;
	a.random();
	b.random();
	c.random();
	c = (a*conj(b)); 
	SU3<float> C(A, B);
	return 0;
}
