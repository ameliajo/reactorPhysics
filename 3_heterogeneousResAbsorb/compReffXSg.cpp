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



auto narrowResonanceFlux( double sigmaTotU238, double sigmaDilution ){
  double sigmaPotU238 = 11.2934;
  return ( sigmaPotU238 + sigmaDilution ) / ( sigmaTotU238 + sigmaDilution );
}

auto wideResonanceFlux( double sigmaAbsU238, double sigmaDilution ){
  return ( sigmaDilution ) / ( sigmaAbsU238 + sigmaDilution );
}

auto intermediateResonanceFlux( double sigmaAbsU238, double sigmaDilution ){
  double lambda = 0.5;
  double sigmaPotU238 = 11.2934;
  return ( lambda * sigmaPotU238 + sigmaDilution ) / 
         ( sigmaAbsU238 +  lambda * sigmaPotU238 + sigmaDilution );
}




  
auto calcResIntegralXS(){
  bool infiniteDilution = false;
  bool narrowResonance = true;
  bool wideResonance = true;

  int vecSize = 5000;
  std::vector<double> energy(vecSize), capture300(vecSize), capture1000(vecSize), total300(vecSize), total1000(vecSize);
  std::ifstream in("capture.txt");
  int i = 0;
  while(!in.eof()){
    in >> energy[i];
    in >> capture300[i];
    in >> capture1000[i];
    in >> total300[i];
    in >> total1000[i];
    ++i;
  }

  std::vector<double> integral_sigma_flux_300  { 0, 0, 0 }, 
                      integral_sigma_flux_1000 { 0, 0, 0 };
  std::vector<double> resonance_integral_300   { 0, 0, 0 },
                      resonance_integral_1000  { 0, 0, 0 };
  double eLeft, eRight, sigmaAbsU238_300, sigmaAbsU238_1000, intPieceXS300,
         intPieceXS1000,resIntPiece300,resIntPiece1000, sigmaTotU238_300,
         sigmaTotU238_1000, sigmaDilution;
  sigmaDilution = 2000;
  for ( int i = 1; i < vecSize; ++i ){
    if ( infiniteDilution ){
      // THIS IS FOR INFINITE DILUTION
      // now we're going to integrate from eLeft --> eRight
      eLeft  = energy[i-1], 
      eRight = energy[i],
      sigmaAbsU238_300  = ( capture300[i-1]  + capture300[i]  ) / 2.0,
      sigmaAbsU238_1000 = ( capture1000[i-1] + capture1000[i] ) / 2.0,
      intPieceXS300  = sigmaAbsU238_300  * log( eRight / eLeft ),
      intPieceXS1000 = sigmaAbsU238_1000 * log( eRight / eLeft ),
      resIntPiece300  = log( eRight / eLeft ),
      resIntPiece1000 = log( eRight / eLeft );

    } 
    else if ( narrowResonance ){
      // THIS IS FOR NARROW RESONANCE APPROXIMATION
      // now we're going to integrate from eLeft --> eRight
      eLeft  = energy[i-1], 
      eRight = energy[i],
      sigmaTotU238_300  = ( total300[i-1]    + total300[i]   ) / 2.0,
      sigmaTotU238_1000 = ( total1000[i-1]   + total1000[i]   ) / 2.0,
      sigmaAbsU238_300  = ( capture300[i-1]  + capture300[i] ) / 2.0,
      sigmaAbsU238_1000 = ( capture1000[i-1] + capture1000[i] ) / 2.0;

      double flux300  = narrowResonanceFlux( sigmaTotU238_300,  sigmaDilution ); 
      double flux1000 = narrowResonanceFlux( sigmaTotU238_1000, sigmaDilution ); 

      intPieceXS300   = flux300  * sigmaAbsU238_300  * log( eRight / eLeft ),
      intPieceXS1000  = flux1000 * sigmaAbsU238_1000 * log( eRight / eLeft ),
      resIntPiece300  = flux300  * log( eRight / eLeft ),
      resIntPiece1000 = flux1000 * log( eRight / eLeft );

    }
    else if ( wideResonance ){
      // THIS IS FOR WIDE RESONANCE APPROXIMATION
      // now we're going to integrate from eLeft --> eRight
      eLeft  = energy[i-1], 
      eRight = energy[i],
      sigmaTotU238_300  = ( total300[i-1]    + total300[i]   ) / 2.0,
      sigmaTotU238_1000 = ( total1000[i-1]   + total1000[i]   ) / 2.0,
      sigmaAbsU238_300  = ( capture300[i-1]  + capture300[i] ) / 2.0,
      sigmaAbsU238_1000 = ( capture1000[i-1] + capture1000[i] ) / 2.0;

      double flux300  = wideResonanceFlux( sigmaAbsU238_300,  sigmaDilution ); 
      double flux1000 = wideResonanceFlux( sigmaAbsU238_1000, sigmaDilution ); 

      intPieceXS300   = flux300  * sigmaAbsU238_300  * log( eRight / eLeft ),
      intPieceXS1000  = flux1000 * sigmaAbsU238_1000 * log( eRight / eLeft ),
      resIntPiece300  = flux300  * log( eRight / eLeft ),
      resIntPiece1000 = flux1000 * log( eRight / eLeft );

    }
    else {
      // THIS IS FOR INTERMEDIATE RESONANCE APPROXIMATION
      // now we're going to integrate from eLeft --> eRight
      eLeft  = energy[i-1], 
      eRight = energy[i],
      sigmaTotU238_300  = ( total300[i-1]    + total300[i]   ) / 2.0,
      sigmaTotU238_1000 = ( total1000[i-1]   + total1000[i]   ) / 2.0,
      sigmaAbsU238_300  = ( capture300[i-1]  + capture300[i] ) / 2.0,
      sigmaAbsU238_1000 = ( capture1000[i-1] + capture1000[i] ) / 2.0;

      double flux300  = intermediateResonanceFlux( sigmaAbsU238_300,  sigmaDilution ); 
      double flux1000 = intermediateResonanceFlux( sigmaAbsU238_1000, sigmaDilution ); 

      intPieceXS300   = flux300  * sigmaAbsU238_300  * log( eRight / eLeft ),
      intPieceXS1000  = flux1000 * sigmaAbsU238_1000 * log( eRight / eLeft ),
      resIntPiece300  = flux300  * log( eRight / eLeft ),
      resIntPiece1000 = flux1000 * log( eRight / eLeft );

    }
      int index;
      if      ( eLeft <= 10 ){ index = 0; } 
      else if ( eLeft <= 25 ){ index = 1; }
      else                   { index = 2; }
      integral_sigma_flux_300[index]  += intPieceXS300; 
      integral_sigma_flux_1000[index] += intPieceXS1000; 
      resonance_integral_300[index]   += resIntPiece300; 
      resonance_integral_1000[index]  += resIntPiece1000; 






  } 
  std::cout << "---------  " << resonance_integral_300[0] << "     " << resonance_integral_300[0] / integral_sigma_flux_300[0]<< std::endl;
  std::cout << "---------  " << resonance_integral_300[1] << "     " << resonance_integral_300[0] / integral_sigma_flux_300[1]<< std::endl;
  std::cout << "---------  " << resonance_integral_300[2] << "     " << resonance_integral_300[0] / integral_sigma_flux_300[2]<< std::endl;

}





int main() {
  calcResIntegralXS();

  return 0;

}
