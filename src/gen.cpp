#include <omp.h>
#include <TGraph.h>
#include <TMath.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TF1.h>

#include "DatabasePDG2.h"
#include "gen.h"
#include "particle.h"
#include "const.h"

using namespace std ;

#define _USE_MATH_DEFINES
const double C_Feq = (pow(0.5/M_PI/hbarC,3)) ;

namespace gen{


int Nelem ;
double *ntherm, dvMax, dsigmaMax ;
TRandom3 *rnd ;
DatabasePDG2 *database ;
int NPART ;

struct element {
 double tau, x, y, eta ;
 double u[4] ;
 double dsigma[4] ;
 double T, mub, muq, mus ;
 double pi[10] ;
 double Pi ;
} ;

element *surf ;
int *npart ;               // number of generated particles in each event

const double pmax = 10.0 ;
const double rapmax = 6.0 ;

const double c1 = pow(1./2./hbarC/TMath::Pi(),3.0) ;
double *cumulantDensity ; // particle densities (thermal). Seems to be redundant, but needed for fast generation
double totalDensity ; // sum of all thermal densities


// ######## load the elements
void load(char *filename, int N)
{
 double dV, vEff=0.0, vEffOld=0.0, dvEff, dvEffOld ;
 int nfail=0, ncut=0 ;
 TLorentzVector dsigma ;
 Nelem = N ;
 surf = new element [Nelem] ;

 cout<<"reading "<<N<<" lines from  "<<filename<<"\n" ;
 ifstream fin(filename) ;
 if(!fin){ cout << "cannot read file " << filename << endl ; exit(1) ; }
 dvMax=0. ;
 dsigmaMax=0. ;
 // ---- reading loop
 string line ;
 istringstream instream ;
 cout<<"1?: failbit="<<instream.fail()<<endl ;
 for(int n=0; n<Nelem; n++){
   getline(fin, line) ;
   instream.str(line) ;
   instream.seekg(0) ;
   instream.clear() ; // does not work with gcc 4.1 otherwise
    instream>>surf[n].tau>>surf[n].x>>surf[n].y>>surf[n].eta
      >>surf[n].dsigma[0]>>surf[n].dsigma[1]>>surf[n].dsigma[2]>>surf[n].dsigma[3]
      >>surf[n].u[0]>>surf[n].u[1]>>surf[n].u[2]>>surf[n].u[3]
      >>surf[n].T>>surf[n].mub>>surf[n].muq>>surf[n].mus>>dvEff ;
      if(surf[n].muq>0.12){ surf[n].muq=0.12 ; // omit charge ch.pot. for test
	ncut++ ;
      }
      if(surf[n].muq<-0.12){ surf[n].muq=-0.12 ; // omit charge ch.pot. for test
	ncut++ ;
      }

   if(instream.fail()){ cout<<"reading failed at line "<<n<<"; exiting\n" ; exit(1) ; }
   // calculate in the old way
   dvEffOld = surf[n].dsigma[0]*surf[n].u[0]+surf[n].dsigma[1]*surf[n].u[1]+
   surf[n].dsigma[2]*surf[n].u[2]+surf[n].dsigma[3]*surf[n].u[3] ;
   vEffOld += dvEffOld ;
   if(dvEffOld<0.0){
     //cout<<"!!! dvOld!=dV " << dvEffOld <<"  " << dV << "  " << surf[n].tau <<endl ;
     nfail++ ;
   }
   //if(nfail==100) exit(1) ;
   // ---- boost
   //dsigma.SetXYZT(-surf[n].dsigma[1],-surf[n].dsigma[2],-surf[n].dsigma[3],surf[n].dsigma[0]) ;
   //dsigma.Boost(-surf[n].u[1]/surf[n].u[0],-surf[n].u[2]/surf[n].u[0],-surf[n].u[3]/surf[n].u[0]) ;
   //// ######################################################################
   //// ###  now and later surf.dsigma are boosted to the fluid rest frame  ##
   //// ######################################################################
   //surf[n].dsigma[0] = dsigma.T() ;
   //surf[n].dsigma[1] = -dsigma.X() ;
   //surf[n].dsigma[2] = -dsigma.Y() ;
   //surf[n].dsigma[3] = -dsigma.Z() ;
   //dvEff = surf[n].dsigma[0] ;
   //vEff += dvEff ;
   //if(dvMax<dvEff) dvMax = dvEff ;
   //// maximal value of the weight max(W) = max(dsigma_0+|\vec dsigma_i|)   for equilibrium DFs
   //if(dsigma.T()+dsigma.Rho()>dsigmaMax) dsigmaMax = dsigma.T()+dsigma.Rho() ;
  }
 //cout<<"..done.\n" ;
 //cout<<"Veff = "<<vEff<<"  dvMax = "<<dvMax<<endl ;
 //cout<<"Veff(old) = "<<vEffOld<<endl ;
 cout<<"failed elements = "<<nfail<<endl ;
 cout<<"mu_cut elements = "<<ncut<<endl ;
// ---- prepare some stuff to calculate thermal densities
 NPART=database->GetNParticles() ;
 cout<<"NPART="<<NPART<<endl ;
 cout<<"dsigmaMax="<<dsigmaMax<<endl ;
}

double ffthermal(double *x, double *par)
{
  double &T = par[0] ;
  double &mu = par[1] ;
  double &mass = par[2] ;
  double &stat = par[3] ;
  return x[0]*x[0]/( exp((sqrt(x[0]*x[0]+mass*mass)-mu)/T) - stat ) ;
}

void doCalculations(char *out_file)
{
 const double gmumu [4] = {1., -1., -1., -1.} ;
 const double mass = 0.1396; // pion
 vector<double> pt, phiP;
 for(double pti=0.1; pti<2.5; pti+=0.1){
  pt.push_back(pti);
 }
 for(double phi=0.0; phi<2.0*M_PI; phi+=M_PI/100.0){
  phiP.push_back(phi);
 }
 double dndp [pt.size()][phiP.size()];
 for(int ip=0; ip<pt.size(); ip++)
 for(int iphi=0; iphi<phiP.size(); iphi++)
  dndp[ip][iphi] = 0.0;
 
 ofstream fout (out_file);
 if(!fout) { cout<<"I/O error with " << out_file << endl; exit(1); }
 
 for(int iel=0; iel<Nelem; iel++){ // loop over all elements
  for(int ip=0; ip<pt.size(); ip++)
  for(int iphi=0; iphi<phiP.size(); iphi++){
  double mT = sqrt(mass*mass + pt[ip]*pt[ip]);
  double p [4] = {mT, pt[ip]*cos(phiP[iphi]), pt[ip]*sin(phiP[iphi]), 0};
  double pds = 0., pu = 0.;
  for(int mu=0; mu<4; mu++){
   pds += p[mu]*surf[iel].dsigma[mu];
   pu += p[mu]*surf[iel].u[mu]*gmumu[mu];
  }
  dndp[ip][iphi] += c1 * pds / ( exp(pu/surf[iel].T) - 1.0);
 }
 } // loop over all elements
 // Fourier coefficients
  const int nhar = 6;
  double vn [nhar][pt.size()];
  for(int ip=0; ip<pt.size(); ip++)
  for(int ihar=0; ihar<nhar; ihar++){
   vn[ihar][ip] = 0.0;
  }

  // Fourier integrals
  for(int ip=0; ip<pt.size(); ip++)
  for(int iphi=0; iphi<phiP.size(); iphi++)
  for(int ihar=0; ihar<nhar; ihar++){
   vn[ihar][ip] += dndp[ip][iphi] * cos(ihar * phiP[iphi]) / (double)(phiP.size());
  }
  // normalisations
  for(int ip=0; ip<pt.size(); ip++)
  for(int ihar=1; ihar<nhar; ihar++){
   vn[ihar][ip] = vn[ihar][ip] / vn[0][ip];
  }
 // writing results to file
 for(int ip=0; ip<pt.size(); ip++){
  fout << setw(14) << pt[ip];
  for(int ihar=0; ihar<nhar; ihar++)
   fout << setw(14) << vn[ihar][ip];
  fout << endl;
 }
 fout.close();
 cout << "calculation done.\n" << endl ;
}

} // end namespace gen
