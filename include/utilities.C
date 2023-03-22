// ---------------------------------------------------------
// Content: Ultility functions uses for laserball MC
// Authors:  Martti Nirkko and Sammy Valder, University of Sussex (2022)
// ---------------------------------------------------------

// C++ stuff
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

// ROOT stuff
#include <TCanvas.h>
#include <TColor.h>
#include <TF2.h>
#include <TEllipse.h>
#include <TError.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TVector2.h>
#include <TVector3.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TStyle.h>

// Helpful functions

// Include frequently used functions from 'std' (entire namespace is bad practice)
using std::cerr;
using std::cout;
using std::endl;
using std::flush;

// -----------------------------------------------------------------------------
// Run time parameters
const int MORE_OUTPUT = 1;                  // additional plots for testing

// Global constants
const double pi = TMath::Pi();
const TVector3 e1(1,0,0);
const TVector3 e2(0,1,0);
const TVector3 e3(0,0,1);

// Physical constants
const double LAMBDA   = 5.0e-4;             // mm -> TELLIE wavelength
const double N_WATER  = 1.33772;            // at 500 nm -> http://www.philiplaven.com/p20.html
const double C_VACUUM = 299.792458;         // mm/ns (detector units)
const double C_WATER  = C_VACUUM/N_WATER;   // mm/ns

// -----------------------------------------------------------------------------
// Initialise functions
std::string printVector(const TVector3&);
void printProgress(int, int);
void GetRotationAngles(const TVector3&, double&, double&);
void DrawCircle(const TVector3&, double, TVector3**, int);
double getFWHM(TH1*);

// -----------------------------------------------------------------------------
/// Display vector as a std::string
std::string printVector(const TVector3& v) {
  std::string out;
  if (v.Mag() < 10) out = Form("(%.3f, %.3f, %.3f)", v.X(),  v.Y(), v.Z());
  else              out = Form("(%.1f, %.1f, %.1f)", v.X(),  v.Y(), v.Z());
  return out.c_str();
}

// -----------------------------------------------------------------------------
/// Display progress bar within a loop (can't have any other output in loop!)
void printProgress(int it, int n) {
  int div = (n - (n % 100))/100;              // 1% division
  if (it % div != 0 && it != n-1) return;
  float prog = (float)it/n;
  int barWidth = 70;
  cout << "[";
  int pos = barWidth * prog;
  for (int i=0; i<barWidth; ++i) {
    if (i < pos)                cout << "=";  // processed
    else if (pos+1 == barWidth) cout << "=";  // reached 100%
    else if (i == pos)          cout << ">";  // processing
    else                        cout << " ";  // not yet processed
  }
  int perc = (int)round(100.*prog);
  if (it < n-1) cout << "] " << perc << "%\r" << flush;
  else          cout << "] " << perc << "%\r" << endl;
}

// -----------------------------------------------------------------------------
/// Get rotation angles to rotate from given direction to xyz-frame
void GetRotationAngles(const TVector3& vec, double &rot_Z, double &rot_X) {

  // Get rotation angle from z-axis to central direction
  rot_Z = acos(e3*(vec.Unit()));         // rotate counter-clockwise, towards x-axis by [0,pi]
  
  // Get rotation angle from x-axis to central direction (projected to transverse plane)
  int sign = 0;
  if (vec[1] > 0) sign = +1;
  else            sign = -1;
  TVector3 vect(vec[0],vec[1],0);
  rot_X = sign*acos(e1*(vect.Unit()));   // rotate around z-axis by [-pi,pi]
}


// -----------------------------------------------------------------------------
/// Get FWHM of histogram
double getFWHM(TH1* hist) {
  int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2.);
  int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2.);
  double fwhm = hist->GetXaxis()->GetBinUpEdge(bin2) - hist->GetXaxis()->GetBinLowEdge(bin1);
  return fwhm;
}

