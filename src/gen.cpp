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

using namespace std;

#define _USE_MATH_DEFINES
const double C_Feq = (pow(0.5 / M_PI / hbarC, 3));

// Levi-Civita symbols
// i,j,k,l = 0...3
int levi(int i, int j, int k, int l)
{
 if((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l)) return 0;
 else return ( (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)/12 );
}

namespace gen {

int Nelem;
double *ntherm, dvMax, dsigmaMax;
TRandom3 *rnd;
DatabasePDG2 *database;
ParticlePDG2 *particle; // chosen particle sort for the polarization calc
int NPART;

struct element {
 double tau, x, y, eta;
 double u[4];
 double dsigma[4];
 double T, mub, muq, mus;
 double dbeta [4][4];
};

element *surf;
vector<double> px, py, rap;
vector<vector<vector<double> > > Pi_num; // numerator of Eq. 34
vector<vector<double> > Pi_den; // denominator of Eq. 34
int nhydros;

const double c1 = pow(1. / 2. / hbarC / TMath::Pi(), 3.0);

// ######## load the elements
void load(char *filename, int N) {
 double dV, vEff = 0.0, vEffOld = 0.0, dvEff, dvEffOld;
 int nfail = 0, ncut = 0;
 TLorentzVector dsigma;
 Nelem = N;
 surf = new element[Nelem];

 cout << "reading " << N << " lines from  " << filename << "\n";
 ifstream fin(filename);
 if (!fin) {
  cout << "cannot read file " << filename << endl;
  exit(1);
 }
 dvMax = 0.;
 dsigmaMax = 0.;
 // ---- reading loop
 string line;
 istringstream instream;
 for (int n = 0; n < Nelem; n++) {
  getline(fin, line);
  instream.str(line);
  instream.seekg(0);
  instream.clear();  // does not work with gcc 4.1 otherwise
  instream >> surf[n].tau >> surf[n].x >> surf[n].y >> surf[n].eta >>
      surf[n].dsigma[0] >> surf[n].dsigma[1] >> surf[n].dsigma[2] >>
      surf[n].dsigma[3] >> surf[n].u[0] >> surf[n].u[1] >> surf[n].u[2] >>
      surf[n].u[3] >> surf[n].T >> surf[n].mub >> surf[n].muq >> surf[n].mus;
  for(int i=0; i<4; i++)
  for(int j=0; j<4; j++)
   instream >> surf[n].dbeta[i][j];
  //if (surf[n].muq > 0.12) {
   //surf[n].muq = 0.12;  // omit charge ch.pot. for test
   //ncut++;
  //}
  //if (surf[n].muq < -0.12) {
   //surf[n].muq = -0.12;  // omit charge ch.pot. for test
   //ncut++;
  //}

  if (instream.fail()) {
   cout << "reading failed at line " << n << "; exiting\n";
   exit(1);
  }
  // calculate in the old way
  dvEffOld =
      surf[n].dsigma[0] * surf[n].u[0] + surf[n].dsigma[1] * surf[n].u[1] +
      surf[n].dsigma[2] * surf[n].u[2] + surf[n].dsigma[3] * surf[n].u[3];
  vEffOld += dvEffOld;
  if (dvEffOld < 0.0) {
   // cout<<"!!! dvOld!=dV " << dvEffOld <<"  " << dV << "  " << surf[n].tau
   // <<endl ;
   nfail++;
  }
  // if(nfail==100) exit(1) ;
  // ---- boost
  // dsigma.SetXYZT(-surf[n].dsigma[1],-surf[n].dsigma[2],-surf[n].dsigma[3],surf[n].dsigma[0])
  // ;
  // dsigma.Boost(-surf[n].u[1]/surf[n].u[0],-surf[n].u[2]/surf[n].u[0],-surf[n].u[3]/surf[n].u[0])
  // ;
  //// ######################################################################
  //// ###  now and later surf.dsigma are boosted to the fluid rest frame  ##
  //// ######################################################################
  // surf[n].dsigma[0] = dsigma.T() ;
  // surf[n].dsigma[1] = -dsigma.X() ;
  // surf[n].dsigma[2] = -dsigma.Y() ;
  // surf[n].dsigma[3] = -dsigma.Z() ;
  // dvEff = surf[n].dsigma[0] ;
  // vEff += dvEff ;
  // if(dvMax<dvEff) dvMax = dvEff ;
  //// maximal value of the weight max(W) = max(dsigma_0+|\vec dsigma_i|)   for
  ///equilibrium DFs
  // if(dsigma.T()+dsigma.Rho()>dsigmaMax) dsigmaMax = dsigma.T()+dsigma.Rho() ;
 }
 // cout<<"..done.\n" ;
 // cout<<"Veff = "<<vEff<<"  dvMax = "<<dvMax<<endl ;
 // cout<<"Veff(old) = "<<vEffOld<<endl ;
 // cout<<"failed elements = "<<nfail<<endl ;
 // cout<<"mu_cut elements = "<<ncut<<endl ;
 //// ---- prepare some stuff to calculate thermal densities
 // NPART=database->GetNParticles() ;
 // cout<<"NPART="<<NPART<<endl ;
 // cout<<"dsigmaMax="<<dsigmaMax<<endl ;
}

void initCalc() {
 for (double _px = -4.0; _px < 4.05; _px += 0.2) {
  px.push_back(_px);
 }
 for (double _py = -4.0; _py < 4.05; _py += 0.2) {
  py.push_back(_py);
 }
 for (double _rap = -1.0; _rap < 1.05; _rap += 0.4) {
  rap.push_back(_rap);
 }
 Pi_num.resize(px.size());
 Pi_den.resize(px.size());
 for (int i = 0; i < Pi_num.size(); i++) {
  Pi_num[i].resize(py.size());
  Pi_den[i].resize(py.size());
  for (int j = 0; j < Pi_num[i].size(); j++) {
   Pi_den[i][j] = 0.0;
   Pi_num[i][j].resize(4);
   for(int k=0; k<4; k++)
    Pi_num[i][j][k] = 0.0;
  }
 }
 nhydros = 0;
}

double ffthermal(double *x, double *par) {
 double &T = par[0];
 double &mu = par[1];
 double &mass = par[2];
 double &stat = par[3];
 return x[0] * x[0] / (exp((sqrt(x[0] * x[0] + mass * mass) - mu) / T) - stat);
}

void doCalculations() {
 cout <<"call doCalculations()\n";
 const double gmumu[4] = {1., -1., -1., -1.};
 particle = database->GetPDGParticle(3122);
 const double mass = particle->GetMass();  // pion
 const double baryonCharge = particle->GetBaryonNumber();
 const double electricCharge = particle->GetElectricCharge();
 const double strangeness = particle->GetStrangeness();
 cout << "calculations for: " << particle->GetName() << ", charges = "
  << baryonCharge << "  " << electricCharge << "  " << strangeness << endl;
 int nFermiFail = 0; // number of elements where nf>1.0
 int nBadElem = 0;
 for (int iel = 0; iel < Nelem; iel++) {  // loop over all elements
  if(fabs(surf[iel].dbeta[0][0])>1000.0) nBadElem++;
  //if(fabs(surf[iel].dbeta[0][0])>1000.0) continue;
  for (int ipx = 0; ipx < px.size(); ipx++)
   for (int ipy = 0; ipy < py.size(); ipy++) {
    double mT = sqrt(mass * mass + px[ipx] * px[ipx] + py[ipy] * py[ipy]);
    double p[4] = {mT, px[ipx], py[ipy], 0};
    double p_[4] = {mT, -px[ipx], -py[ipy], 0};
    double pds = 0., pu = 0.;
    for (int mu = 0; mu < 4; mu++) {
     pds += p[mu] * surf[iel].dsigma[mu];
     pu += p[mu] * surf[iel].u[mu] * gmumu[mu];
    }
    const double mutot = surf[iel].mub * baryonCharge
      + surf[iel].muq * electricCharge + surf[iel].mus * strangeness;
    const double nf = c1 / (exp( (pu - mutot) / surf[iel].T) + 1.0);
    if(nf > 1.0) nFermiFail++;
    Pi_den[ipx][ipy] += pds * nf ;
    for(int mu=0; mu<4; mu++)
     for(int nu=0; nu<4; nu++)
      for(int rh=0; rh<4; rh++)
       for(int sg=0; sg<4; sg++)
        Pi_num[ipx][ipy][mu] += pds * nf * (1. - nf) * levi(mu, nu, rh, sg)
                                * p_[sg] * surf[iel].dbeta[nu][rh];
   }
 }  // loop over all elements
 delete[] surf;
 cout << "doCalculations: total, bad = " << setw(12) << Nelem << setw(12) << nBadElem << endl;
 cout << "number of elements*pT configurations where nf>1.0: " << nFermiFail
  << endl;
}

void calcBasicQuantities() {
 cout << "call calcBasicQuantities()\n";
 const double gmumu[4] = {1., -1., -1., -1.};
 particle = database->GetPDGParticle(3122);
 const double mass = particle->GetMass();  // pion
 const double baryonCharge = particle->GetBaryonNumber();
 const double electricCharge = particle->GetElectricCharge();
 const double strangeness = particle->GetStrangeness();
 double den=0., num_v2 = 0., num_pt=0.;
 vector<double> den_v1 (rap.size(), 0.);
 vector<double> num_v1 (rap.size(), 0.);
 for (int iel = 0; iel < Nelem; iel++) {  // loop over all elements
  for (int ipx = 0; ipx < px.size(); ipx++)
   for (int ipy = 0; ipy < py.size(); ipy++) {
    double mT = sqrt(mass * mass + px[ipx] * px[ipx] + py[ipy] * py[ipy]);
    double p[4] = {mT, px[ipx], py[ipy], 0};
    double pds = 0., pu = 0.;
    for (int mu = 0; mu < 4; mu++) {
     pds += p[mu] * surf[iel].dsigma[mu];
     pu += p[mu] * surf[iel].u[mu] * gmumu[mu];
    }
    const double mutot = surf[iel].mub * baryonCharge
      + surf[iel].muq * electricCharge + surf[iel].mus * strangeness;
    const double nf = c1 / (exp( (pu - mutot) / surf[iel].T) - 1.0);
    //if(nf > 1.0) nFermiFail++;
    den += pds * nf;
    num_v2 += (px[ipx] * px[ipx] - py[ipy] * py[ipy])/
     (px[ipx] * px[ipx] + py[ipy] * py[ipy]) * pds * nf;
    num_pt += sqrt(px[ipx] * px[ipx] + py[ipy] * py[ipy]) * pds * nf;
    // *** v_1 calc
    for (int irap = 0; irap < rap.size(); irap++) {
     double p[4] = {mT*cosh(rap[irap]), px[ipx], py[ipy], mT*sinh(rap[irap])};
     double pds = 0., pu = 0.;
     for (int mu = 0; mu < 4; mu++) {
      pds += p[mu] * surf[iel].dsigma[mu];
      pu += p[mu] * surf[iel].u[mu] * gmumu[mu];
     }
     const double nf = c1 / (exp( (pu - mutot) / surf[iel].T) - 1.0);
     den_v1[irap] += pds * nf;
     num_v1[irap] += px[ipx] / sqrt(px[ipx] * px[ipx] + py[ipy] * py[ipy])
      * pds * nf;
    } // rapidity loop for v_1
   }
 }  // loop over all elements
 cout << "*** calcBasicQuantities:\n";
 cout << "v_2:   " << num_v2/den << "  num, den = " << num_v2 << "  " << den << endl;
 cout << "<p_T>: " << num_pt/den << "  num, den = " << num_pt << "  " << den << endl;
 cout << "rapidity-dependent v_1:\n";
 for (int irap = 0; irap < rap.size(); irap++)
  cout << setw(14) << num_v1[irap]/den_v1[irap];
 cout << endl;
}

void sumPtPolarization() {
 double Pi [4] = {0.}, norm=0.0;
 for (int ipx = 0; ipx < px.size(); ipx++)
  for (int ipy = 0; ipy < py.size(); ipy++) {
   for(int mu=0; mu<4; mu++)
    Pi[mu] += Pi_num[ipx][ipy][mu] * hbarC / (8.0 * particle->GetMass());
   norm += Pi_den[ipx][ipy];
 }
 cout << "*** polarization components:\n";
 for(int mu=0; mu<4; mu++)
  cout << setw(14) << Pi[mu]/norm;
 cout << endl;
}

void outputPolarization(char *out_file) {
 ofstream fout(out_file);
 if (!fout) {
  cout << "I/O error with " << out_file << endl;
  exit(1);
 }
 for (int ipx = 0; ipx < px.size(); ipx++)
  for (int ipy = 0; ipy < py.size(); ipy++) {
   fout << setw(14) << px[ipx] << setw(14) << py[ipy]
     << setw(14) << Pi_den[ipx][ipy];
   for(int mu=0; mu<4; mu++)
    fout << setw(14) << Pi_num[ipx][ipy][mu] * hbarC / (8.0 * particle->GetMass());
   fout << endl;
 }
 fout.close();
}

}  // end namespace gen
