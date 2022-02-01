#include <omp.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMath.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TH1.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <TF1.h>
#include <sstream>
#include <TUUID.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>

#include "DatabasePDG2.h"
#include "gen.h"

// ############################################################
//  execution modes:
//  1) calculation of polarization
//     #define PLOTS commented out in src/gen.h
//  2) plots of invariant combinations of beta derivatives
//     #define PLOTS UNcommented out in src/gen.h
// ############################################################

using namespace std;
int getNlines(char *filename);

int ranseed;

extern "C" {
void getranseedcpp_(int *seed) { *seed = ranseed; }
}

// ########## MAIN block ##################

int main(int argc, char **argv) {
 // command-line parameters

	 if(argc!=3 && argc !=5){
	 cout<< "INVALID SINTAX!"<<endl;
	 cout<<"use './calc <surface_file> <output_file>' to compute Lambda polarization at decoupling."<<endl;
	 cout<<"use './calc <surface_file> <output_file> <pdg_mother> <pdg_son>' to compute Lambda polarization inherited from the decay mother->Lambda+son."<<endl;
	 exit(1);
	 }

 char surface_file[200], output_file[200];
 strcpy(surface_file, argv[1]);
 strcpy(output_file, argv[2]);
 //========= particle database init
 DatabasePDG2 *database = new DatabasePDG2("Tb/ptl3.data", "Tb/dky3.mar.data");
 database->LoadData();
 //	database->SetMassRange(0.01, 10.0); //-------without PHOTONS
 //	database->SetWidthRange(0., 10.0);
 database->SortParticlesByMass();
 database->CorrectBranching();
 database->DumpData();
 cout << " pion index = " << database->GetPionIndex() << endl;
 gen::database = database;

 #ifdef PLOTS
 TApplication theApp("App", &argc, argv);
 gStyle->SetOptStat(kFALSE);
 #endif
 // ========== generator init
 gen::initCalc();
 gen::load(surface_file, getNlines(surface_file));
 #ifndef PLOTS
 gen::calcEP1();
 switch(argc){
	 case 3:
		gen::doCalculations();
		gen::outputPolarization(output_file);
	break;
	case 5:
		gen::Polfromdecay(atoi(argv[3]),atoi(argv[4]), output_file); //allowed decays 3212,22 or 3224,211 
	break;
	default:
	cout<<"ERROR: something went wrong! Check the arguments of main."<<endl;
}
 #else
 gen::calcInvariantQuantities();
 #endif

 // ========== trees & files
 time_t start, end;
 time(&start);

 time(&end);
 float diff2 = difftime(end, start);
 cout << "Execution time = " << diff2 << " [sec]" << endl;
 #ifdef PLOTS
 theApp.Run();
 #endif
 return 0;
}

// auxiliary function to get the number of lines
int getNlines(char *filename) {
 ifstream fin(filename);
 if (!fin) {
  cout << "getNlines: cannot open file: " << filename << endl;
  exit(1);
 }
 string line;
 int nlines = 0;
 while (fin.good()) {
  getline(fin, line);
  nlines++;
 };
 fin.close();
 return nlines - 1;
}
