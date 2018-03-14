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
  bool infiniteDilution = true;
  bool narrowResonance = false;
  bool wideResonance = false;

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

  std::vector<double> numerator_300  { 0, 0, 0 }, 
                      numerator_1000 { 0, 0, 0 };
  std::vector<double> denominator_300   { 0, 0, 0 },
                      denominator_1000  { 0, 0, 0 };
  double eLeft, eRight, sigmaAbsU238_300, sigmaAbsU238_1000, numeratorPiece300,
         numeratorPiece1000,denominatorPiece300,denominatorPiece1000, sigmaTotU238_300,
         sigmaTotU238_1000, sigmaDilution;
  sigmaDilution = 20;
  for ( int i = 1; i < vecSize; ++i ){
    if ( infiniteDilution ){
      // THIS IS FOR INFINITE DILUTION
      // now we're going to integrate from eLeft --> eRight
      eLeft  = energy[i-1], 
      eRight = energy[i],
      sigmaAbsU238_300  = ( capture300[i-1]  + capture300[i]  ) / 2.0,
      sigmaAbsU238_1000 = ( capture1000[i-1] + capture1000[i] ) / 2.0,
      numeratorPiece300  = sigmaAbsU238_300  * log( eRight / eLeft ),
      numeratorPiece1000 = sigmaAbsU238_1000 * log( eRight / eLeft ),
      denominatorPiece300  = log( eRight / eLeft ),
      denominatorPiece1000 = log( eRight / eLeft );

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

      numeratorPiece300   = flux300  * sigmaAbsU238_300  * log( eRight / eLeft ),
      numeratorPiece1000  = flux1000 * sigmaAbsU238_1000 * log( eRight / eLeft ),
      denominatorPiece300  = flux300  * log( eRight / eLeft ),
      denominatorPiece1000 = flux1000 * log( eRight / eLeft );

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

      numeratorPiece300   = flux300  * sigmaAbsU238_300  * log( eRight / eLeft ),
      numeratorPiece1000  = flux1000 * sigmaAbsU238_1000 * log( eRight / eLeft ),
      denominatorPiece300  = flux300  * log( eRight / eLeft ),
      denominatorPiece1000 = flux1000 * log( eRight / eLeft );

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

      numeratorPiece300   = flux300  * sigmaAbsU238_300  * log( eRight / eLeft ),
      numeratorPiece1000  = flux1000 * sigmaAbsU238_1000 * log( eRight / eLeft ),
      denominatorPiece300  = flux300  * log( eRight / eLeft ),
      denominatorPiece1000 = flux1000 * log( eRight / eLeft );

    }
      int index;
      if      ( eLeft <= 10 ){ index = 0; } 
      else if ( eLeft <= 25 ){ index = 1; }
      else                   { index = 2; }
      numerator_300[index]  += numeratorPiece300; 
      numerator_1000[index] += numeratorPiece1000; 
      denominator_300[index]   += denominatorPiece300; 
      denominator_1000[index]  += denominatorPiece1000; 






  } 
  std::cout << "--------- 6-10  at 300 " << numerator_300[0] << "     " << 1.0/ ( denominator_300[0] / numerator_300[0] ) << std::endl;
  std::cout << "--------- 10-25 at 300 " << numerator_300[1] << "     " << 1.0/ (denominator_300[1] / numerator_300[1] ) << std::endl;
  std::cout << "--------- 25-50 at 300 " << numerator_300[2] << "     " << 1.0/ (denominator_300[2] / numerator_300[2] ) << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "--------- 6-10  at 1000 " << numerator_1000[0] << "     " << 1.0/ ( denominator_1000[0] / numerator_1000[0] ) << std::endl;
  std::cout << "--------- 10-25 at 1000 " << numerator_1000[1] << "     " << 1.0/ (denominator_1000[1] / numerator_1000[1] ) << std::endl;
  std::cout << "--------- 25-50 at 1000 " << numerator_1000[2] << "     " << 1.0/ (denominator_1000[2] / numerator_1000[2] ) << std::endl;
  std::cout << std::endl;




}





int main() {
  calcResIntegralXS();
  return 0;
}
