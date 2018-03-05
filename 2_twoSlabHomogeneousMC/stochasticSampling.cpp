#include <iostream>
#include <math.h>

int main(){
  int n = 1e6;
  float sumMu = 0;
  for ( int i = 0; i < n; ++i ){
    float xi = (float)rand()/RAND_MAX;
    float innerP = pow(16*xi*xi - 16*xi + 5, 0.5 );
    float first = pow( innerP - 4*xi + 2, 0.33333 );
    float mu = first - (1.0/first);
    sumMu += mu;

  }
  std::cout << sumMu / ( n * 1.0 ) << std::endl;
  return 0;
}

