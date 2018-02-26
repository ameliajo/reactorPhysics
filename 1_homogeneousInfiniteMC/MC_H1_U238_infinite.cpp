#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>



class Neutron {
  private:
    float energy;
  public:
    // Constructor
    Neutron( float energy ) : energy(energy) { }
    float getEnergy() { return energy; }
    void  setEnergy( float newEnergy ){ energy = newEnergy; }
  };

void collideNeutronsH1( Neutron& n ){
  float r = (float)rand()/RAND_MAX;
  float newEnergy = r*n.getEnergy();
  n.setEnergy( newEnergy );
}

void collideNeutronsU238( Neutron& n ){
  float r = (float)rand()/RAND_MAX;
  float alpha = (237.0*237.0/(239.0*239.0));
  float newEnergy = n.getEnergy()*( alpha - 1.0 ) * r + n.getEnergy();
  n.setEnergy( newEnergy );
}





int main() {
  float sigmaU238, sigmaH1, NU238, NH1, total, sigmaU238_absorption, sigmaU238_elastic;

  std::ifstream uResData;
  uResData.open("uraniumRes1000.txt");
  std::ifstream is("uraniumRes1000.txt");
  std::istream_iterator<float> start(is), end;
  std::vector<float> uData(start, end);
  std::vector<float>::const_iterator beginning = uData.begin();
  std::vector<float> energyList (beginning,      beginning+1000),
                     elasticList(beginning+1000, beginning+2000),
                     gammaList  (beginning+2000, beginning+3000);


  int numNeutrons = 1e4;

  std::vector<int> energyTallies( 1000, 0 );
  int j = 0;


  std::ofstream dataFile2;
  dataFile2.open("data_H_U.txt");


  int numAbsorbed = 0;
  for ( int i = 0; i < numNeutrons; ++i ){ 
    Neutron n(1000);

    while ( n.getEnergy() > 1 ){

      // Get the macroscopic cross section of U-238 at this energy
      sigmaU238_absorption = gammaList[int(n.getEnergy())-1];
      sigmaU238_elastic = elasticList[int(n.getEnergy())-1];
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


  return 0;
}
