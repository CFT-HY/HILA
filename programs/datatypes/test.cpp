#include "sun.h"
#include "../plumbing/defs.h"

int main(){
	seed_mersenne(123);
	SU2<float> a, b, c;
	a.random();
	b.random();
	c.random();
	c = (a*conjugate<float>(b));
	return 0;
}
