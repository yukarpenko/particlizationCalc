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

int leviC(int i, int j, int k, int l)
// contravariant (upper indexes) Levi-Civita tensor in Minkowski coordinates
{
 const double gmumu[4] = {1., -1., -1., -1.};
 if((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l)) return 0;
 else return -( (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)/12 )*
   gmumu[i]*gmumu[j]*gmumu[k]*gmumu[l];
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
 double du [4][4];
 double dT [4];
};

element *surf;
vector<double> pT, phi;
vector<vector<vector<double> > > Pi_num_gradT, Pi_num_accel, Pi_num_vorti; // numerator of Eq. 34
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
   instream >> surf[n].du[i][j];
  for(int i=0; i<4; i++)
   instream >> surf[n].dT[i];
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
 for (double _pT = 0.0; _pT < 3.01; _pT += 0.2) {
  pT.push_back(_pT);
 }
 for (double _phi = 0.0; _phi < 2.0*M_PI-1e-5; _phi += M_PI/20) {
  phi.push_back(_phi);
 }
 Pi_num_gradT.resize(pT.size());
 Pi_num_accel.resize(pT.size());
 Pi_num_vorti.resize(pT.size());
 Pi_den.resize(pT.size());
 for (int i = 0; i < Pi_den.size(); i++) {
  Pi_num_gradT[i].resize(phi.size());
  Pi_num_accel[i].resize(phi.size());
  Pi_num_vorti[i].resize(phi.size());
  Pi_den[i].resize(phi.size());
  for (int j = 0; j < Pi_den[i].size(); j++) {
   Pi_den[i][j] = 0.0;
   Pi_num_gradT[i][j].resize(4);
   Pi_num_accel[i][j].resize(4);
   Pi_num_vorti[i][j].resize(4);
   for(int k=0; k<4; k++) {
    Pi_num_gradT[i][j][k] = 0.0;
    Pi_num_accel[i][j][k] = 0.0;
    Pi_num_vorti[i][j][k] = 0.0;
   }
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
 double Qx1=0., Qy1=0., Qx2=0., Qy2=0.;
 for (int iel = 0; iel < Nelem; iel++) {  // loop over all elements
  if(fabs(surf[iel].du[0][0])>1000.0) nBadElem++;
  //if(fabs(surf[iel].dbeta[0][0])>1000.0) continue;
  if(surf[iel].T<1e-8) continue;
  for (int ipt = 0; ipt < pT.size(); ipt++)
   for (int iphi = 0; iphi < phi.size(); iphi++) {
    double mT = sqrt(mass * mass + pT[ipt] * pT[ipt]);
    double p[4] = {mT, pT[ipt]*cos(phi[iphi]), pT[ipt]*sin(phi[iphi]), 0};
    double p_[4] = {mT, -pT[ipt]*cos(phi[iphi]), -pT[ipt]*sin(phi[iphi]), 0};
    double u_[4];
    for(int i=0; i<4; i++) u_[i] = surf[iel].u[i] * gmumu[i];
    double pds = 0., pu = 0.;
    for (int mu = 0; mu < 4; mu++) {
     pds += p[mu] * surf[iel].dsigma[mu];
     pu += p[mu] * surf[iel].u[mu] * gmumu[mu];
    }
    const double mutot = surf[iel].mub * baryonCharge
      + surf[iel].muq * electricCharge + surf[iel].mus * strangeness;
    const double nf = c1 / (exp( (pu - mutot) / surf[iel].T) + 1.0);
    if(nf > 1.0) nFermiFail++;
    Pi_den[ipt][iphi] += pds * nf ;
    double omega[4] = {0.}, A_[4] = {0.};
    for(int mu=0; mu<4; mu++)
     for(int rh=0; rh<4; rh++)
      for(int sg=0; sg<4; sg++)
       for(int ta=0; ta<4; ta++) {
        omega[mu] += 0.5 * leviC(mu, rh, sg, ta) * surf[iel].du[rh][sg] * u_[ta];
       }
    for(int mu=0; mu<4; mu++)
     for(int rh=0; rh<4; rh++) {
      A_[mu] += surf[iel].u[rh] * surf[iel].du[rh][mu];
     }
    double sum1[4] = {0.}, sum2[4] = {0.};
    for(int mu=0; mu<4; mu++)
     for(int rh=0; rh<4; rh++)
      for(int sg=0; sg<4; sg++)
       for(int ta=0; ta<4; ta++){
        sum1[mu] += leviC(mu, rh, sg, ta) * surf[iel].dT[rh] * u_[sg] * p_[ta];
        sum2[mu] += leviC(mu, rh, sg, ta) * A_[sg] * u_[rh] * p_[ta] / surf[iel].T;
        if(sum2[mu]!=sum2[mu]) {
         cout << "sum2 is naan\n";
         exit(77);
        }
       }
    double wp = 0.;
    for(int mu=0; mu<4; mu++) {
     wp += p_[mu] * omega[mu];
    }
    const double pref = pds * nf * (1. - nf) * hbarC / (8.0 * particle->GetMass());
    for(int mu=0; mu<4; mu++) {
     Pi_num_gradT[ipt][iphi][mu] += pref * sum1[mu];
     Pi_num_accel[ipt][iphi][mu] += pref * sum2[mu];
     Pi_num_vorti[ipt][iphi][mu] += pref * 2./surf[iel].T *
      (omega[mu] * pu - wp * surf[iel].u[mu]);
     }
    Qx1 += p[1] * pds * nf;
    Qy1 += p[2] * pds * nf;
    Qx2 += (p[1]*p[1] - p[2]*p[2])/(pT[ipt]+1e-10) * pds * nf;
    Qy2 += (p[1]*p[2])/(pT[ipt]+1e-10) * pds * nf;
   }
 }  // loop over all elements
 delete[] surf;
 cout << "doCalculations: total, bad = " << setw(12) << Nelem << setw(12) << nBadElem << endl;
 cout << "number of elements*pT configurations where nf>1.0: " << nFermiFail
  << endl;
 cout << "event_plane_vectors: " << Qx1 << "  " << Qy1 << "  "
   << Qx2 << "  " << Qy2 << endl;
}

void calcEP1() {
 const double gmumu[4] = {1., -1., -1., -1.};
 particle = database->GetPDGParticle(2112);
 const double mass = particle->GetMass();  // pion
 int nBadElem = 0;
 double Qx1=0., Qy1=0., Qx2=0., Qy2=0.;
 const double coshY = cosh(1.0);
 const double sinhY = sinh(1.0);
 const double pT = 1.0;
 int nFFail = 0;
 for (int iel = 0; iel < Nelem; iel++) {  // loop over all elements
  if(fabs(surf[iel].du[0][0])>1000.0) nBadElem++;
  //if(fabs(surf[iel].dbeta[0][0])>1000.0) continue;
  for (int iphi = 0; iphi < phi.size(); iphi++) {
   double mT = sqrt(mass * mass + pT * pT);
   const double sin_phi = sin(phi[iphi]);
   const double cos_phi = cos(phi[iphi]);
   double p1[4] = {mT*coshY, pT*cos_phi, pT*sin_phi, mT*sinhY};
   double p2[4] = {mT*coshY, pT*cos_phi, pT*sin_phi, -mT*sinhY};
   double pds1 = 0., pds2 = 0., pu1 = 0., pu2 = 0.;
   for (int mu = 0; mu < 4; mu++) {
    pds1 += p1[mu] * surf[iel].dsigma[mu];
    pu1 += p1[mu] * surf[iel].u[mu] * gmumu[mu];
    pds2 += p2[mu] * surf[iel].dsigma[mu];
    pu2 += p2[mu] * surf[iel].u[mu] * gmumu[mu];
   }
   const double f1 = c1 * exp( - pu1 / surf[iel].T );
   const double f2 = c1 * exp( - pu2 / surf[iel].T );
   if(f1 > 1.0) nFFail++;
   Qx1 += p1[1] * pds1 * f1;
   Qy1 += p1[2] * pds1 * f1;
   Qx2 += p2[1] * pds2 * f2;
   Qy2 += p2[2] * pds2 * f2;
  }
 }  // loop over all elements
 cout << "EP1_vectors: " << Qx1 << "  " << Qy1 << "  "
   << Qx2 << "  " << Qy2 << endl;
 cout << "EP_angles: " << atan2(Qy1, Qx1) << "  " << atan2(Qy2, Qx2) << endl;
}

void outputPolarization(char *out_file) {
 ofstream fout(out_file);
 if (!fout) {
  cout << "I/O error with " << out_file << endl;
  exit(1);
 }
 for (int ipt = 0; ipt < pT.size(); ipt++)
  for (int iphi = 0; iphi < phi.size(); iphi++) {
   fout << setw(14) << pT[ipt] << setw(14) << phi[iphi]
     << setw(14) << Pi_den[ipt][iphi];
   for(int mu=0; mu<4; mu++)
    fout << setw(14) << (Pi_num_gradT[ipt][iphi][mu] +
      Pi_num_accel[ipt][iphi][mu] + Pi_num_vorti[ipt][iphi][mu]);
   for(int mu=0; mu<4; mu++)
    fout << setw(14) << Pi_num_gradT[ipt][iphi][mu];
   for(int mu=0; mu<4; mu++)
    fout << setw(14) << Pi_num_accel[ipt][iphi][mu];
   for(int mu=0; mu<4; mu++)
    fout << setw(14) << Pi_num_vorti[ipt][iphi][mu];
   fout << endl;
 }
 fout.close();
 char dim_file [200];
 strcpy(dim_file, out_file);
 strcat(dim_file, ".dim");
 ofstream fdim(dim_file);
 fdim << pT.size() << "  " << phi.size() << endl;
 fdim.close();
}

}  // end namespace gen
