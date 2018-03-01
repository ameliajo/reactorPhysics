#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <iomanip>  
#include<string>


/* Two slab problem, assuming isotropic scatting in center of mass, with a 
 * uniformly distributed isotropic unit source in first slab. Vacuum boundary
 * conditions, and 1 energy group for the neutrons.
 */


struct Neutron {
  float x, mu, phi, w;
  std::vector<float> distances;  // First index is for distance traveled
  std::vector<int>   collisions; // in first material, second index is for 
                                 // distances traveled in second material
  Neutron( float x, float mu, float phi, float w, 
    std::vector<float> distances = {0,0}, std::vector<int> collisions = {0,0} ): 
    x(x), mu(mu), phi(phi), w(w), distances(distances), collisions(collisions) {}
};


int main() {
  float xsMacroTotal1 = 1,   xsMacroAbsor1 = 0.5, 
        xsMacroTotal2 = 1.5, xsMacroAbsor2 = 1.2;
  float A = 12;
  
  //srand (4);


  std::ofstream myfile;
  myfile.open ("data.txt");
  int numNeutrons = 1e4;

  std::vector<float> pathLength1 ( numNeutrons, 0 ), 
                     pathLength2 ( numNeutrons, 0 );
  std::vector<int> collisions1 ( numNeutrons, 0 ), collisions2 ( numNeutrons, 0 );


  for ( int i = 0; i < numNeutrons; ++i ){ 
    // Initialize 
    float x   = 2*(float)rand()/RAND_MAX;      // Uniform  0 -> 1
    float mu  = 2*(float)rand()/RAND_MAX-1;    // Uniform -1 -> 1
    float phi = 2*M_PI*(float)rand()/RAND_MAX; // Uniform  0 -> 2pi
    float w   = 2*(float)rand()/RAND_MAX-1;
 
    Neutron n(x,mu,phi,w);
    while ( true ){
      float xsi = (float)rand()/RAND_MAX;
      float s = -log(xsi)/xsMacroTotal1;

      if ( n.x < 2 and n.x + n.w*s > 2 ){
//        std::cout << "going from 1 --> 2" << "     " << n.x+n.w*s << std::endl;
        n.distances[0] += abs( 2 - n.x );
        n.x = 2+1e-5;
      }
      else if ( n.x > 2 and n.x + n.w*s < 2 ){
//        std::cout << "going from 2 --> 1" << "     " << n.x+n.w*s << std::endl;
        n.distances[1] += abs( n.x - 2 );
        n.x = 2-1e-5;
      }
      else if ( n.x + n.w*s < 0 or n.x + n.w*s > 4 ){
        if ( n.x + n.w*s < 0 ){
          n.distances[0] += abs( n.x );
//          std::cout << "leak left\n" << std::endl;
        }
        else {
          n.distances[1] += abs( 4 - n.x );
//          std::cout << "leak right\n" << std::endl;
        }
        break;
      }
      else {
        float total = n.x < 2 ? xsMacroTotal1 : xsMacroTotal2;
        float absor = n.x < 2 ? xsMacroAbsor1 : xsMacroAbsor2;
        if ( n.x < 2  ){ n.distances[0] += abs( n.w*s ); }
        if ( n.x < 2  ){ n.collisions[0] += 1; }
        if ( n.x >= 2 ){ n.distances[1] += abs( n.w*s ); }
        if ( n.x >= 2 ){ n.collisions[1] += 1; }
        xsi = (float)rand()/RAND_MAX;
        if ( xsi*total < absor ){
//          if ( n.x < 2 ){ std::cout << "Absorption in 1\n" << std::endl; }
//          else          { std::cout << "Absorption in 2\n" << std::endl; }
          break;
        }
        else {
//          if ( n.x < 2 ){ std::cout << "Scattering in 1" << std::endl; }
//          else          { std::cout << "Scattering in 2" << std::endl; }

          n.mu  = 2*(float)rand()/RAND_MAX-1;    // Uniform -1 -> 1

          float mu_lab = ( 1 + A*n.mu ) / pow( A*A + 2*A*n.mu + 1, 0.5 );
          n.w = mu_lab*n.w - pow( 1-mu_lab*mu_lab, 0.5)*pow( 1-n.w*n.w,0.5)*cos(n.phi);
          n.phi = 2*M_PI*(float)rand()/RAND_MAX; // Uniform  0 -> 2pi
        }

        n.x += n.w*s;

      }
      
      //std::cout << "x position:  " << n.x << std::endl;
      //myfile << wprime << std::endl;
    } 
    pathLength1[i] = n.distances[0];
    pathLength2[i] = n.distances[1];
    collisions1[i] = n.collisions[0];
    collisions2[i] = n.collisions[1];
  }
  myfile.close();

  float totalDist1, totalDist2, totalCollide1, totalCollide2;
  for ( size_t i = 0; i < numNeutrons; ++i ){
    std::cout << pathLength1[i] << "         " << pathLength2[i] << std::endl;
    totalDist1    += pathLength1[i]; totalDist2    += pathLength2[i];
    totalCollide1 += collisions1[i]; totalCollide2 += collisions2[i];
  }

  std::cout << "flux in 1 is " << std::setw(15) << 
    totalDist1 / (numNeutrons*2.0)  << std::setw(15) <<
    totalCollide1 / (2.0*xsMacroTotal1*numNeutrons) << std::endl;
  std::cout << "flux in 2 is " << std::setw(15) << 
    totalDist2 / (numNeutrons*4.0)  << std::setw(15)<< 
    totalCollide2 / (4.0*xsMacroTotal2*numNeutrons) << std::endl;

  return 0;
}
