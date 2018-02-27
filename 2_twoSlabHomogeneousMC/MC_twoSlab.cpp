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
  float x, phi, mu;
  Neutron( float x, float phi, float mu ): 
    x(x), phi(phi), mu(mu) {}
};


int main() {
  float xsMacroTotal1 = 1,   xsMacroAbsor1 = 0.5, 
        xsMacroTotal2 = 1.5, xsMacroAbsor2 = 1.2;

  std::ofstream myfile;
  myfile.open ("data.txt");
  int numNeutrons = 1;

  for ( int i = 0; i < numNeutrons; ++i ){ 
    // Initialize 
    float x   = 2*(float)rand()/RAND_MAX;      // Uniform  0 -> 1
    float mu  = 2*(float)rand()/RAND_MAX-1;    // Uniform -1 -> 1
    float phi = 2*M_PI*(float)rand()/RAND_MAX; // Uniform  0 -> 2pi
     //myfile  << std::setw(15) << std::left << x << 
     //           std::setw(15) << std::left << mu << 
     //           std::setw(15) << std::left << phi << 
     //           std::setw(15) << std::left << phi << std::endl;
 
    Neutron n(x,mu,phi);
    while ( true ){
      //std::cout << "sample travel distance" << std::endl;
      float xsi = (float)rand()/RAND_MAX;
      float s = -log(xsi)/xsMacroTotal1;

      break;
    } 
  }
  myfile.close();

  /*
    while ( n.getEnergy() > 1 ){

      // Get the macroscopic cross section of U-238 at this energy
      sigmaU238_absorption = 5;
      sigmaU238_elastic = 2;
      sigmaU238 = sigmaU238_absorption + sigmaU238_elastic;
      sigmaH1   = 20;
      NH1 = 10;
      NU238 = 1;
      total = sigmaU238 * NU238 + sigmaH1 * NH1;
      float r1 = (float)rand()/RAND_MAX;
      if ( r1*total < sigmaU238 * NU238 ){ 
        float r2 = (float)rand()/RAND_MAX;
        if ( r2*sigmaU238 < sigmaU238_absorption ){ 
          numAbsorbed += 1;
          break;
        }
        else {
          collideNeutronsU238( n ); 
        }
      }
      else {
        collideNeutronsH1( n );  
      }

      dataFile2<< n.getEnergy() << std::endl;

    } // while energy is above 10 eV
  } // for neutron loop

  std::cout << "Num. neutrons absorbed in U-238: " << numAbsorbed << std::endl;
  std::cout << "%    neutrons absorbed in U-238: " << 100.0*numAbsorbed/numNeutrons << std::endl;
  dataFile2.close();


  */
  return 0;
}
