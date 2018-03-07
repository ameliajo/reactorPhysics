#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

// This uses the SLBW model that I made in #1, and computes the numerical 
// resonance integrals and group cross sections for 
//
//     T (K)     |     6-10 eV     |      10-25 eV       |      25-50 eV
// -----------------------------------------------------------------------
//      300      |   RI    XS      |     RI     XS       |    RI    XS
//     1000      |   RI    XS      |     RI     XS       |    RI    XS
// -----------------------------------------------------------------------




int main() {
  

  int vecSize = 5000;
  std::vector<double> energy(vecSize), capture300(vecSize), capture1000(vecSize);
  std::ifstream in("capture.txt");
  int i = 0;
  while(!in.eof()){
    in >> energy[i];
    in >> capture300[i];
    in >> capture1000[i];
    ++i;
  }

  // i = 0 corresponds 
  std::vector<double> integral_sigma_flux { 0, 0, 0 };
  std::vector<double> resonance_integral{ 0, 0, 0 };
  for ( int i = 1; i < vecSize; ++i ){
    // now we're going to integrate from eLeft --> eRight
    double eLeft  = energy[i-1], 
           eRight = energy[i],
           xs     = ( capture300[i-1] + capture300[i] ) / 2.0,
           toAdd  = xs * log( eRight / eLeft );
    if      ( eLeft <= 10 ){ integral_sigma_flux[0] += toAdd; }
    else if ( eLeft <= 25 ){ integral_sigma_flux[1] += toAdd; }
    else                   { integral_sigma_flux[2] += toAdd; }

    toAdd = log( eRight / eLeft );
    if      ( eLeft <= 10 ){ resonance_integral[0] += toAdd; }
    else if ( eLeft <= 25 ){ resonance_integral[1] += toAdd; }
    else                   { resonance_integral[2] += toAdd; }

  } 
  std::cout << "---------  " << resonance_integral[0] << "     " << resonance_integral[0] / integral_sigma_flux[0]<< std::endl;
  std::cout << "---------  " << resonance_integral[1] << "     " << resonance_integral[0] / integral_sigma_flux[1]<< std::endl;
  std::cout << "---------  " << resonance_integral[2] << "     " << resonance_integral[0] / integral_sigma_flux[2]<< std::endl;
  return 0;

}
