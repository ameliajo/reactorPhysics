#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <iomanip>  

/* Two slab problem, assuming isotropic scatting in center of mass, with a 
 * uniformly distributed isotropic unit source in first slab. Vacuum boundary
 * conditions, and 1 energy group for the neutrons.
 */


struct Neutron {
  float x, mu, phi, w;
  Neutron( float x, float mu, float phi, float w ): 
    x(x), mu(mu), phi(phi), w(w) {}
};


int main() {
  float xsMacroTotal1 = 1,   xsMacroAbsor1 = 0.5, 
        xsMacroTotal2 = 1.5, xsMacroAbsor2 = 1.2;
  float A = 12;


  std::ofstream myfile;
  myfile.open ("data.txt");
  int numNeutrons = 10;

  for ( int i = 0; i < numNeutrons; ++i ){ 
    // Initialize 
    float x   = 2*(float)rand()/RAND_MAX;      // Uniform  0 -> 1
    float mu  = 2*(float)rand()/RAND_MAX-1;    // Uniform -1 -> 1
    float phi = 2*M_PI*(float)rand()/RAND_MAX; // Uniform  0 -> 2pi
    float w   = 2*(float)rand()/RAND_MAX-1;
     //myfile  << std::setw(15) << std::left << x << 
     //           std::setw(15) << std::left << mu << 
     //           std::setw(15) << std::left << phi << 
     //           std::setw(15) << std::left << phi << std::endl;
 
    Neutron n(x,mu,phi,w);
    while ( true ){
      //std::cout << "sample travel distance" << std::endl;
      float xsi = (float)rand()/RAND_MAX;
      float s = -log(xsi)/xsMacroTotal1;
      float mu_lab = ( 1 + A*n.mu ) / pow( A*A + 2*A*n.mu + 1, 0.5 );
      n.w = mu_lab*n.w - pow( 1-mu_lab*mu_lab, 0.5)*pow( 1-n.w*n.w,0.5)*cos(n.phi);
      if ( n.x < 2 and n.x + n.w*x > 2 ){
        std::cout << "going from 1 --> 2" << "     " << w << std::endl;
        n.x = 2+1e-5;
      }
      else if ( n.x > 2 and n.x + n.w*x < 2 ){
        std::cout << "going from 2 --> 1" << "     " << w << std::endl;
        n.x = 2-1e-5;
      }
      else if ( n.x + n.w*x < 0 ){
        std::cout << "leak below\n" << std::endl;
        break;
      }
      else if ( n.x + n.w*x > 4 ){
        std::cout << "leak above\n" << std::endl;
        break;
      }
      else {
        float total = n.x < 2 ? xsMacroTotal1 : xsMacroTotal2;
        float absor = n.x < 2 ? xsMacroAbsor1 : xsMacroAbsor2;
        xsi = (float)rand()/RAND_MAX;
        if ( xsi*total < absor ){
          if ( n.x < 2 ){ std::cout << "Absorption in 1\n" << std::endl; }
          else          { std::cout << "Absorption in 2\n" << std::endl; }
          break;
        }
        else {
          if ( n.x < 2 ){ std::cout << "Scattering in 1" << std::endl; }
          else          { std::cout << "Scattering in 2" << std::endl; }

          n.mu  = 2*(float)rand()/RAND_MAX-1;    // Uniform -1 -> 1
          n.phi = 2*M_PI*(float)rand()/RAND_MAX; // Uniform  0 -> 2pi
        }

        n.x += n.w*s;

      }
      
      //std::cout << "x position:  " << n.x << std::endl;
      //myfile << wprime << std::endl;
    } 
  }
  myfile.close();


  return 0;
}
