#include <iostream>
#include <vector>
#include <fstream>


class Neutron {
  private:
    float energy;
  public:
    // Parameterized Constructor
    Neutron( float energy ) : energy(energy) { }
    
    float getEnergy() { return energy; }
    void setEnergy( double newEnergy ){ energy = newEnergy; }
  };

void collideNeutrons( Neutron& n ){
  float r = (float)rand()/RAND_MAX;
  float newEnergy = r*n.getEnergy();
  n.setEnergy( newEnergy );
}




int main() {
  // This is a quick script for neutrons colliding in infinite H1 material 
  // to show neutron slowing down via elastic scattering. Neutrons are cut off
  // once they reach below 1 eV.
  std::ofstream dataFile;
  dataFile.open("data.txt");

  int numNeutrons = 8e3;

  std::vector<int> energyTallies( 1000, 0 );
  int j = 0;

  for ( int i = 0; i < numNeutrons; ++i ){ 
    Neutron n(1000);
    while ( n.getEnergy() > 1 ){


      float r = (float)rand()/RAND_MAX;

      if ( r < 1.0 ){  // Should we even collide? This is dependent on the 
                       // macroscopic xs which I don't have but I did kind of 
                       // want to add it anyway
        collideNeutrons( n );  

      } // if we should collide

      dataFile << n.getEnergy() << std::endl;

    } // while energy is above 10 eV
  } // for neutron loop

  dataFile.close();

  return 0;
}
