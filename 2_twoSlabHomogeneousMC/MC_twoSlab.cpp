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

int getCellNum ( float x ){
  float spacingforRes = 0.1;
  return int( x / spacingforRes );
}

int getMaterial ( float x ){
  return x < 2 ? 1 : 2;
}

struct Neutron {
  float x, mu, phi, w;
  std::vector<float> distances;  // First index is for distance traveled
  std::vector<int>   collisions; // in first material, second index is for 
                                 // distances traveled in second material
                                 
  int material;
  
  Neutron( float x, float mu, float phi, float w, int numSlits,
    int material = 1 ): 
    x(x), mu(mu), phi(phi), w(w), 
    distances (std::vector<float> (numSlits,0)), 
    collisions(std::vector<int>   (numSlits,0)),
    material(material) {}
};


int main() {
  float xsMacroTotal1 = 1,   xsMacroAbsor1 = 0.5, 
        xsMacroTotal2 = 1.5, xsMacroAbsor2 = 1.2;
  float A = 1;
  
  float spacingforRes = 0.1;
  int numSlits = round(6 / spacingforRes);

  int numNeutrons = 1e7;

  std::vector<float> totalDistances  ( numSlits, 0 ), totalDistancesSq  ( numSlits, 0 );
  std::vector<int>   totalCollisions ( numSlits, 0 ), totalCollisionsSq ( numSlits, 0 );

  for ( int i = 0; i < numNeutrons; ++i ){ 
    // Initialize 
    float x   = 2*(float)rand()/RAND_MAX;      // Uniform  0 -> 2
    float mu  = 2*(float)rand()/RAND_MAX-1;    // Uniform -1 -> 1
    float phi = M_PI*(float)rand()/RAND_MAX; // Uniform  0 -> 2pi
    float w   = cos(phi);
 
    
    //std::cout << "--------------------" << std::endl;
    Neutron n(x,mu,phi,w,numSlits);

    while ( true ){
      float xsi = (float)rand()/RAND_MAX;
      float s = -log(xsi)/xsMacroTotal1;

      if ( n.material == 1 and getMaterial(n.x + n.w*s) == 2 ){                  // Going 1 --> 2
     //   std::cout << "1--->2" << std::endl;
        n.material = 2;
        int cellI = getCellNum( n.x );
        int cellF = getCellNum( 2.0001 );
        int cellM = cellI;
     //   std::cout << "goes from " << n.x << "   to    " << n.x+n.w*s<< std::endl;
     //   std::cout << "goes from " << cellI << "   to    " << cellF << std::endl;
        while (cellM != cellF ){
          if ( cellM == cellI ){
            n.distances[cellM] += abs(( getCellNum(n.x) * spacingforRes - n.x )/n.w);
          } 
          else {
            n.distances[cellM] += abs(spacingforRes / n.w); 
          }
          cellM = cellI < cellF ? cellM + 1 : cellM - 1 ;
        }
        /*
        std::cout << n.distances[14] << std::endl;
        std::cout << n.distances[15] << std::endl;
        std::cout << n.distances[16] << std::endl;
        std::cout << n.distances[17] << std::endl;
        std::cout << n.distances[18] << std::endl;
        std::cout << n.distances[19] << std::endl;
        std::cout << n.distances[20] << std::endl;
        std::cout << n.distances[21] << std::endl;
        */

        n.x = 2.00001;


      }
      else if ( n.material == 2 and getMaterial(n.x + n.w*s) == 1 ){       // Going 2 --> 1
      //  std::cout << "2--->1" << std::endl;
        n.material = 1;
        int cellI = getCellNum( n.x );
        int cellF = getCellNum( 1.9999 );
        int cellM = cellI;
     //   std::cout << "goes from " << n.x << "   to    " << n.x+n.w*s<< std::endl;
     //   std::cout << "goes from " << cellI << "   to    " << cellF << std::endl;

        while (cellM != cellF ){
          if ( cellM == cellI ){
            n.distances[cellM] += abs(( getCellNum(n.x) * spacingforRes - n.x )/n.w);
          } 
          else {
            n.distances[cellM] += abs(spacingforRes / n.w); 
          }
          cellM = cellI < cellF ? cellM + 1 : cellM - 1 ;
        }
        
        n.x = 1.99999;

      }
      else if ( n.x + n.w*s < 0 or n.x + n.w*s > 6 ){      
        if ( n.x + n.w*s < 0 ){
          //std::cout << "Leak left" << std::endl;
          int cellI = getCellNum(n.x);
         // std::cout << n.x << "   " << cellI << std::endl;
         // std::cout << n.distances[0] << std::endl;
         // std::cout << n.distances[1] << std::endl;
        //  std::cout << n.distances[2] << std::endl;
          n.distances[getCellNum(n.x)] += abs( (n.x - ( getCellNum(n.x) * spacingforRes ))/n.w ); 
          for ( int i = 0; i < getCellNum(n.x); ++i ){
            n.distances[i] += abs( spacingforRes / n.w );
          }
      //    std::cout << n.distances[0] << std::endl;
      //    std::cout << n.distances[1] << std::endl;
      //    std::cout << n.distances[2] << std::endl;

         // if ( cellI != 0 ){
         // return 0;
         // }
        }
        else {
      //    std::cout << "Leak right" << std::endl;
          
//          n.distances[numSlits-1] += abs( (4.0 - n.x)/n.w );                // Leak right
        int cellI = getCellNum( n.x );
        int cellF = numSlits-1;
        int cellM = cellI;

        while (cellM != cellF ){
          if ( cellM == cellI ){
              n.distances[cellM] += abs(( (getCellNum( n.x )+1.0)* spacingforRes - n.x )/n.w);
          } 
          else {
            n.distances[cellM] += abs(spacingforRes / n.w); 
          }
          cellM = cellI < cellF ? cellM + 1 : cellM - 1 ;
        }
        
        if ( cellI != cellF ){ 
            n.distances[cellF] += abs(spacingforRes/n.w);
        }
        else {
          n.distances[cellI] +=  s;
        }

 
        }
        break;
      }
      else {
 //       std::cout << "collision    " << n.x << "     " << n.x+n.w*s<< std::endl;
        float total = n.material == 1 ? xsMacroTotal1 : xsMacroTotal2;
        float absor = n.material == 1 ? xsMacroAbsor1 : xsMacroAbsor2;

  //      std::cout << "goes from " << n.x << "   to    " << n.x+n.w*s << std::endl;
        int cellI = getCellNum( n.x );
        int cellF = getCellNum( n.x + n.w*s );
        int cellM = cellI;

        while (cellM != cellF ){
          if ( cellM == cellI ){
            if ( cellI < cellF ){ 
   //           std::cout << getCellNum(n.x)*spacingforRes << "   " <<  getCellNum( n.x)* spacingforRes - n.x << std::endl;
              n.distances[cellM] += abs(( (getCellNum( n.x )+1.0)* spacingforRes - n.x )/n.w);
            }
            else {
              n.distances[cellM] += abs(( n.x - getCellNum(n.x) * spacingforRes )/n.w);
            }
          } 
          else {
            n.distances[cellM] += abs(spacingforRes / n.w); 
          }
          cellM = cellI < cellF ? cellM + 1 : cellM - 1 ;
        }
        
        if ( cellI != cellF ){ 
          if ( cellI < cellF ){
            n.distances[cellF] += abs(( n.x+n.w*s - cellF * spacingforRes )/n.w);
          } 
          else {
            n.distances[cellF] += abs(( (cellF+1) * spacingforRes - ( n.x+n.w*s) )/n.w);
          } 
        }
        else {
          n.distances[cellI] +=  s;
        }

        n.collisions[cellF] += 1;
    //    std::cout << n.distances[19] << std::endl;
    //    std::cout << n.distances[20] << std::endl;
    //    std::cout << n.distances[21] << std::endl;
    //    std::cout << n.distances[22] << std::endl;

        /*
        if ( cellI == 20 ){ 
          std::ofstream out;
          out.open("collisions.txt");
          for ( size_t i = 0; i < totalCollisions.size(); ++i ){
            out << std::setw(10) << i*spacingforRes << std::setw(15) <<  
            // 1.0*totalCollisions[i] << std::setw(15) <<   
            1.0*totalCollisions[i]/numNeutrons << std::setw(15) <<   
            // 1.0*totalDistances[i]<< std::endl; 
            1.0*totalDistances[i]/(numNeutrons*6.0) << std::endl; 
          }
          out.close();
          return 0;
        }
        */



        xsi = (float)rand()/RAND_MAX;
        n.x += n.w*s;
        if ( xsi*total < absor ){
     //     std::cout << "Absorb!" << std::endl;
          break;
        }
        else {
      //    std::cout << "Scatter!" << std::endl;

          n.mu  = 2*(float)rand()/RAND_MAX-1;    // Uniform -1 -> 1

          float mu_lab = ( 1 + A*n.mu ) / pow( A*A + 2*A*n.mu + 1, 0.5 );
          n.w = mu_lab*n.w - pow( 1-mu_lab*mu_lab, 0.5)*pow( 1-n.w*n.w,0.5)*cos(n.phi);
          n.phi = M_PI*(float)rand()/RAND_MAX; // Uniform  0 -> pi
        }


       // std::cout << n.x << std::endl;
      }
      
    } 
    for ( size_t i = 0; i < n.distances.size(); ++i ){
      totalDistances[i]    += n.distances[i];
      totalDistancesSq[i]  += n.distances[i] * n.distances[i];
      totalCollisions[i]   += n.collisions[i];
      totalCollisionsSq[i] += n.collisions[i] * n.collisions[i];
    }

  }
  std::ofstream out;
  out.open("collisions.txt");
  for ( size_t i = 0; i < totalCollisions.size(); ++i ){
    float xbarCollision = 1.0*totalCollisions[i]/(numNeutrons*1.0) ;
    float xbarDistance = 1.0*totalDistances[i]/(numNeutrons) ;
    out << std::setw(10) << i*spacingforRes << std::setw(15) <<  
     // 1.0*totalCollisions[i] << std::setw(15) <<   
      1.0*totalCollisions[i]/(numNeutrons*6.0) << std::setw(15) <<   
     // 1.0*totalDistances[i]<< std::endl; 
      1.0*totalDistances[i]/(numNeutrons*6.0) << std::setw(15) << 
      sqrt( (1.0/(numNeutrons-1)) * ( 1.0*totalCollisionsSq[i] / (1.0*numNeutrons) -  xbarCollision*xbarCollision ) ) << std::setw(15) << 

      sqrt( ( 1.0 / ( numNeutrons - 1 ) ) * ( 1.0*totalDistancesSq[i] / (1.0*numNeutrons) -  xbarDistance*xbarDistance ) ) << 
      std::endl; 
  }
  out.close();
 



  return 0;


  //myfile.close();
}
