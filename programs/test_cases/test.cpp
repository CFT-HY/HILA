#include "../plumbing/field.h"

#pragma transformer loop_function
double h();

int main(){
    
  field<double> a;
  field<double> b[2];
    
  foralldir(d){
    onsites(ALL){
      a[X] = b[1][X+d];
      a[X] += b[0][X+d];
    }
  }


  double sum=0;
  onsites(ALL){
    a[X] = b[0][X+YDOWN];
    sum += a[X];
  }

  return 0;
}


