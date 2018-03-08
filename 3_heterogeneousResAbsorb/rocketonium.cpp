#include <iostream>
#include <vector>
#include <math.h>



int main(){
  std::vector<double> e = { 999.7, 999.9, 1000, 1000.1, 1000.3, 1999.7, 1999.9, 2000, 2000.1, 2000.3 };
  
  // Infinite Dilution
  double n1 =  5000.0  * ( e[1] - e[0] ) - 4998500.0  * log( e[1] / e[0] );
  double n2 =  20000.0 * ( e[2] - e[1] ) - 19996000.0 * log( e[2] / e[1] );
  double n3 = -20000.0 * ( e[3] - e[2] ) + 20004000.0 * log( e[3] / e[2] );
  double n4 = -5000.0  * ( e[4] - e[3] ) + 5001500.0  * log( e[4] / e[3] );

  double n5 =  5000.0  * ( e[6] - e[5] ) - 9998500.0  * log( e[6] / e[5] );
  double n6 =  20000.0 * ( e[7] - e[6] ) - 39996000.0 * log( e[7] / e[6] );
  double n7 = -20000.0 * ( e[8] - e[7] ) + 40004000.0 * log( e[8] / e[7] );
  double n8 = -5000.0  * ( e[9] - e[8] ) + 10001500.0 * log( e[9] / e[8] );

  double den = log( 5000.0 / 1.0 );

  double sigma_g = (n1+n2+n3+n4+n5+n6+n7+n8) / den;
  // Narrow Resonance  
  double a = 5000.0, b = -4998500.0, c = 50.0, d = 5000.0;
  n1 = ( (c+d)*log(a*e[1]+b+c+d) + b*log(e[1]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[0]+b+c+d) + b*log(e[0]) ) * ( (c+d) / (b+c+d) );
  a = 20000.0, b = -19996000.0;
  n2 = ( (c+d)*log(a*e[2]+b+c+d) + b*log(e[2]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[1]+b+c+d) + b*log(e[1]) ) * ( (c+d) / (b+c+d) );
  a = -20000.0, b = 20004000.0;
  n3 = ( (c+d)*log(a*e[3]+b+c+d) + b*log(e[3]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[2]+b+c+d) + b*log(e[2]) ) * ( (c+d) / (b+c+d) );
  a = -5000.0, b = 5001500.0;
  n4 = ( (c+d)*log(a*e[4]+b+c+d) + b*log(e[4]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[3]+b+c+d) + b*log(e[3]) ) * ( (c+d) / (b+c+d) );
  a = 5000.0, b = -9998500.0;
  n5 = ( (c+d)*log(a*e[6]+b+c+d) + b*log(e[6]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[5]+b+c+d) + b*log(e[5]) ) * ( (c+d) / (b+c+d) );
  a = 20000.0, b = -39996000.0;
  n6 = ( (c+d)*log(a*e[7]+b+c+d) + b*log(e[7]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[6]+b+c+d) + b*log(e[6]) ) * ( (c+d) / (b+c+d) );
  a = -20000.0, b = 40004000.0;
  n7 = ( (c+d)*log(a*e[8]+b+c+d) + b*log(e[8]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[7]+b+c+d) + b*log(e[7]) ) * ( (c+d) / (b+c+d) );
  a = -5000.0, b = 10001500.0;
  n8 = ( (c+d)*log(a*e[9]+b+c+d) + b*log(e[9]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[8]+b+c+d) + b*log(e[8]) ) * ( (c+d) / (b+c+d) );

  double numerator = n1+n2+n3+n4+n5+n6+n7+n8;

  a = 5000.0, b = -4998500.0, c = 50.0, d = 5000.0;
  n1 = ( (c+d)*( log(e[1]) - log(a*e[1]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[0]) - log(a*e[0]+b+c+d) ) ) / ( b + c + d );
  a = 20000.0, b = -19996000.0;
  n2 = ( (c+d)*( log(e[2]) - log(a*e[2]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[1]) - log(a*e[1]+b+c+d) ) ) / ( b + c + d );
  a = -20000.0, b = 20004000.0;
  n3 = ( (c+d)*( log(e[3]) - log(a*e[3]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[2]) - log(a*e[2]+b+c+d) ) ) / ( b + c + d );
  a = -5000.0, b = 5001500.0;
  n4 = ( (c+d)*( log(e[4]) - log(a*e[4]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[3]) - log(a*e[3]+b+c+d) ) ) / ( b + c + d );
  a = 5000.0, b = -9998500.0;
  n5 = ( (c+d)*( log(e[6]) - log(a*e[6]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[5]) - log(a*e[5]+b+c+d) ) ) / ( b + c + d );
  a = 20000.0, b = -39996000.0;
  n6 = ( (c+d)*( log(e[7]) - log(a*e[7]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[6]) - log(a*e[6]+b+c+d) ) ) / ( b + c + d );
  a = -20000.0, b = 40004000.0;
  n7 = ( (c+d)*( log(e[8]) - log(a*e[8]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[7]) - log(a*e[7]+b+c+d) ) ) / ( b + c + d );
  a = -5000.0, b = 10001500.0;
  n8 = ( (c+d)*( log(e[9]) - log(a*e[9]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[8]) - log(a*e[8]+b+c+d) ) ) / ( b + c + d );


  double denominator = n1+n2+n3+n4+n5+n6+n7+n8+log(999.7/1.0)+log(1999.7/1000.3)+log(5000.0/2000.3);

  std::cout << "NR, ratio 10 gives us a group xs of " << numerator/denominator << std::endl;

  // Wide Resonance  
  a = 5000.0, b = -4998500.0, c = 0.0, d = 5000.0;
  n1 = ( (c+d)*log(a*e[1]+b+c+d) + b*log(e[1]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[0]+b+c+d) + b*log(e[0]) ) * ( (c+d) / (b+c+d) );
  a = 20000.0, b = -19996000.0;
  n2 = ( (c+d)*log(a*e[2]+b+c+d) + b*log(e[2]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[1]+b+c+d) + b*log(e[1]) ) * ( (c+d) / (b+c+d) );
  a = -20000.0, b = 20004000.0;
  n3 = ( (c+d)*log(a*e[3]+b+c+d) + b*log(e[3]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[2]+b+c+d) + b*log(e[2]) ) * ( (c+d) / (b+c+d) );
  a = -5000.0, b = 5001500.0;
  n4 = ( (c+d)*log(a*e[4]+b+c+d) + b*log(e[4]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[3]+b+c+d) + b*log(e[3]) ) * ( (c+d) / (b+c+d) );
  a = 5000.0, b = -9998500.0;
  n5 = ( (c+d)*log(a*e[6]+b+c+d) + b*log(e[6]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[5]+b+c+d) + b*log(e[5]) ) * ( (c+d) / (b+c+d) );
  a = 20000.0, b = -39996000.0;
  n6 = ( (c+d)*log(a*e[7]+b+c+d) + b*log(e[7]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[6]+b+c+d) + b*log(e[6]) ) * ( (c+d) / (b+c+d) );
  a = -20000.0, b = 40004000.0;
  n7 = ( (c+d)*log(a*e[8]+b+c+d) + b*log(e[8]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[7]+b+c+d) + b*log(e[7]) ) * ( (c+d) / (b+c+d) );
  a = -5000.0, b = 10001500.0;
  n8 = ( (c+d)*log(a*e[9]+b+c+d) + b*log(e[9]) ) * ( (c+d) / (b+c+d) ) -
       ( (c+d)*log(a*e[8]+b+c+d) + b*log(e[8]) ) * ( (c+d) / (b+c+d) );

  numerator = n1+n2+n3+n4+n5+n6+n7+n8;


  a = 5000.0, b = -4998500.0, c = 0.0, d = 5000.0;
  n1 = ( (c+d)*( log(e[1]) - log(a*e[1]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[0]) - log(a*e[0]+b+c+d) ) ) / ( b + c + d );
  a = 20000.0, b = -19996000.0;
  n2 = ( (c+d)*( log(e[2]) - log(a*e[2]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[1]) - log(a*e[1]+b+c+d) ) ) / ( b + c + d );
  a = -20000.0, b = 20004000.0;
  n3 = ( (c+d)*( log(e[3]) - log(a*e[3]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[2]) - log(a*e[2]+b+c+d) ) ) / ( b + c + d );
  a = -5000.0, b = 5001500.0;
  n4 = ( (c+d)*( log(e[4]) - log(a*e[4]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[3]) - log(a*e[3]+b+c+d) ) ) / ( b + c + d );
  a = 5000.0, b = -9998500.0;
  n5 = ( (c+d)*( log(e[6]) - log(a*e[6]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[5]) - log(a*e[5]+b+c+d) ) ) / ( b + c + d );
  a = 20000.0, b = -39996000.0;
  n6 = ( (c+d)*( log(e[7]) - log(a*e[7]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[6]) - log(a*e[6]+b+c+d) ) ) / ( b + c + d );
  a = -20000.0, b = 40004000.0;
  n7 = ( (c+d)*( log(e[8]) - log(a*e[8]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[7]) - log(a*e[7]+b+c+d) ) ) / ( b + c + d );
  a = -5000.0, b = 10001500.0;
  n8 = ( (c+d)*( log(e[9]) - log(a*e[9]+b+c+d) ) ) / ( b + c + d ) -
       ( (c+d)*( log(e[8]) - log(a*e[8]+b+c+d) ) ) / ( b + c + d );


  denominator = n1+n2+n3+n4+n5+n6+n7+n8+log(999.7/1.0)+log(1999.7/1000.3)+log(5000.0/2000.3);

  std::cout << "WR, ratio 10 gives us a group xs of " << numerator/denominator << std::endl;





  return 0;

}

