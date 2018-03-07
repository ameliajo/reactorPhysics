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
  
  /*
  std::ifstream uResData;
  uResData.open("capture.txt");
  std::ifstream is("capture.txt");
  std::istream_iterator<float> start(is), end;
  std::vector<float> uData(start, end);
  std::vector<float>::const_iterator beginning = uData.begin();
  std::vector<float> energyList (beginning,      beginning+1000),
                     elasticList(beginning+1000, beginning+2000),
                     gammaList  (beginning+2000, beginning+3000);
                     */

  std::vector<double> energy(10000), capture300(10000), capture1000(10000);
  std::ifstream in("capture.txt");
  int i = 0;
  while(!in.eof()){
    in >> energy[i];
    in >> capture300[i];
    in >> capture1000[i];
    ++i;
  }

  // i = 0 corresponds 

  double integral_sigma_flux_6_to_10  = 0, integral_sigma_flux_10_to_25 = 0,
         integral_sigma_flux_25_to_50 = 0;
  double resonance_integral_6_to_10 = 0, resonance_integral_10_to_25 = 0,
         resonance_integral_25_to_50 = 0;
  for ( int i = 1; i < 10000; ++i ){
    // now we're going to integrate from eLeft --> eRight
    double eLeft  = energy[i-1];
    double eRight = energy[i];
    double xs = ( capture300[i-1] + capture300[i] ) / 2.0;
    double toAdd = xs * log( eRight / eLeft );
    if      ( eLeft <= 10 ){ integral_sigma_flux_6_to_10  += toAdd; }
    else if ( eLeft <= 25 ){ integral_sigma_flux_10_to_25 += toAdd; }
    else                   { integral_sigma_flux_25_to_50 += toAdd; }

    toAdd = log( eRight / eLeft );
    if      ( eLeft <= 10 ){ resonance_integral_6_to_10  += toAdd; }
    else if ( eLeft <= 25 ){ resonance_integral_10_to_25 += toAdd; }
    else                   { resonance_integral_25_to_50 += toAdd; }


  } 
  std::cout << resonance_integral_6_to_10 / integral_sigma_flux_6_to_10 << std::endl;
  std::cout << integral_sigma_flux_10_to_25 << std::endl;
  std::cout << integral_sigma_flux_25_to_50 << std::endl;
  return 0;

}
