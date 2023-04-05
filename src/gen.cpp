#include <omp.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
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
#include <TH1D.h>

///////
#include <chrono>
using namespace std::chrono;
///////

/////CLASS FOR PROGRESS BAR/////
#include <iomanip>
#include <ostream>
#include <string>

class progress_bar
{
    static const auto overhead = sizeof " [100%]";

    std::ostream& os;
    const std::size_t bar_width;
    std::string message;
    const std::string full_bar;

 public:
    progress_bar(std::ostream& os, std::size_t line_width,
                 std::string message_, const char symbol = '#')
        : os{os},
          bar_width{line_width - overhead},
          message{std::move(message_)},
          full_bar{std::string(bar_width, symbol) + std::string(bar_width, ' ')}
    {
        if (message.size()+1 >= bar_width || message.find('\n') != message.npos) {
            os << message << '\n';
            message.clear();
        } else {
            message += ' ';
        }
        write(0.0);
    }

    // not copyable
    progress_bar(const progress_bar&) = delete;
    progress_bar& operator=(const progress_bar&) = delete;

    ~progress_bar()
    {
        os << '\n';
    }

    void write(double fraction,double time=0,bool total=false, std::string time_unit = "seconds");
};

void progress_bar::write(double fraction,double time,bool total,std::string time_unit)
{
    // clamp fraction to valid range [0,1]
    if (fraction < 0)
        fraction = 0;
    else if (fraction > 1)
        fraction = 1;

    auto width = bar_width - message.size();
    auto offset = bar_width - static_cast<unsigned>(width * fraction);
       
       if(total){
               os << '\r' << message;
               os.write(full_bar.data() + offset, width);
               os << " [" << std::setw(3) << static_cast<int>(100) << "%] total execution time: "<<time<< " "<<time_unit <<std::flush;
               }
       else{
    os << '\r' << message;
    os.write(full_bar.data() + offset, width);
    os << " [" << std::setw(3) << static_cast<int>(100*fraction) << "%] single step in "<<time<< " "<<time_unit <<std::flush;
}
}
////////////////////////////////


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
vector<double> pT, phi;
vector<vector<vector<double> > > Pi_num; // numerator of eq. 1 of 2103.14621
vector<vector<vector<double> > > Pi_num_exact; // (minus) numerator of eq. 3 of 2103.14621
vector<vector<double> > Pi_den; // denominator of eq. 1 of 2103.14621
int nhydros;
TCanvas *plotSymm, *plotAsymm, *plotMod;
TH1D *histMod, *histSymm, *histAsymm;

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
 for (double _pT = 0.0; _pT < 6.01; _pT += 0.2) {
  pT.push_back(_pT);
 }
 for (double _phi = 0.0; _phi < 2.0*M_PI-1e-5; _phi += M_PI/20) {
  phi.push_back(_phi);
 }
 Pi_num.resize(pT.size());
 Pi_num_exact.resize(pT.size());
 Pi_den.resize(pT.size());
 for (int i = 0; i < Pi_num.size(); i++) {
  Pi_num[i].resize(phi.size());
  Pi_num_exact[i].resize(phi.size());
  Pi_den[i].resize(phi.size());
  for (int j = 0; j < Pi_num[i].size(); j++) {
   Pi_den[i][j] = 0.0;
   Pi_num[i][j].resize(4);
   Pi_num_exact[i][j].resize(4);
   for(int k=0; k<4; k++) {
    Pi_num[i][j][k] = 0.0;
    Pi_num_exact[i][j][k] = 0.0;
   }
  }
 }
 nhydros = 0;
 #ifdef PLOTS
 plotSymm = new TCanvas("plotSymm","symmetric derivatives");
 plotAsymm = new TCanvas("plotAsymm","Asymmetric derivatives");
 plotMod = new TCanvas("plotMod","Derivatives");
 histSymm = new TH1D("histSymm", "histSymm", 100, 0., 5.0);
 histAsymm = new TH1D("histAsymm", "histAsymm", 100, -1.2, 0.2);
 histMod = new TH1D("histMod", "histMod", 100, 0., 2.0);
 #endif
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
 const double tvect[4] = {1.,0., 0., 0.};
 particle = database->GetPDGParticle(3122);//3122=Lambda 3212=sigma0 3224=sigmastar
 const double mass = particle->GetMass();  // get particle's mass
 const double baryonCharge = particle->GetBaryonNumber();
 const double electricCharge = particle->GetElectricCharge();
 const double strangeness = particle->GetStrangeness();
 cout << "calculations for: " << particle->GetName() << ", charges = "
  << baryonCharge << "  " << electricCharge << "  " << strangeness << endl;
 int nFermiFail = 0; // number of elements where nf>1.0
 int nBadElem = 0;
 double Qx1=0., Qy1=0., Qx2=0., Qy2=0.;
 
 progress_bar progress(std::clog, 70u, "Working");
 auto tot_start = high_resolution_clock::now();

 for (int iel = 0; iel < Nelem; iel++) {  // loop over all elements of the freeze-out hypersurface
	 auto start = high_resolution_clock::now();
  if(fabs(surf[iel].dbeta[0][0])>1000.0) nBadElem++;
  //if(fabs(surf[iel].dbeta[0][0])>1000.0) continue;
  for (int ipt = 0; ipt < pT.size(); ipt++)
   for (int iphi = 0; iphi < phi.size(); iphi++) {
    double mT = sqrt(mass * mass + pT[ipt] * pT[ipt]);
    double p[4] = {mT, pT[ipt]*cos(phi[iphi]), pT[ipt]*sin(phi[iphi]), 0};
    double p_[4] = {mT, -pT[ipt]*cos(phi[iphi]), -pT[ipt]*sin(phi[iphi]), 0};
    double s[4] = {0.,0.,0.,0.};
    double pds = 0., pu = 0., s_sq=0.;
    
    for(int mu=0; mu<4; mu++)
     for(int nu=0; nu<4; nu++)
      for(int rh=0; rh<4; rh++)
       for(int sg=0; sg<3; sg++){ //sg=3 is zero because p_[3]=0. This speeds up the program
		   s[mu] += levi(mu, nu, rh, sg)
                                * p_[sg] * surf[iel].dbeta[nu][rh]/(2*particle->GetMass());
		   }
    for (int mu = 0; mu < 4; mu++) {
     pds += p[mu] * surf[iel].dsigma[mu];
     pu += p[mu] * surf[iel].u[mu] * gmumu[mu];
     s_sq += s[mu]*s[mu]*gmumu[mu];
    }

	const double mutot = surf[iel].mub * baryonCharge
      + surf[iel].muq * electricCharge + surf[iel].mus * strangeness;
    const double nf = c1 / (exp( (pu - mutot) / surf[iel].T) + 1.0);
    
	if(nf > 1.0) nFermiFail++;
    Pi_den[ipt][iphi] += pds * nf ;
    for(int mu=0; mu<4; mu++){
		Pi_num[ipt][iphi][mu] += pds * nf *  s[mu] *(1. - nf/c1)/4.;
		Pi_num_exact[ipt][iphi][mu] += 0.5*pds * nf *  (s[mu]/sqrt(-s_sq)) * sinh(sqrt(-s_sq)*0.5)/(cosh(sqrt(-s_sq)*0.5)+exp(-(pu - mutot)/surf[iel].T));
		}

    Qx1 += p[1] * pds * nf;
    Qy1 += p[2] * pds * nf;
    Qx2 += (p[1]*p[1] - p[2]*p[2])/(pT[ipt]+1e-10) * pds * nf;
    Qy2 += (p[1]*p[2])/(pT[ipt]+1e-10) * pds * nf;
   }
   auto stop = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>(stop - start);
   progress.write(double(iel)/(Nelem),duration.count(),false,"microseconds");
 }  // loop over all elements
 auto tot_stop = high_resolution_clock::now();
 auto tot_duration = duration_cast<minutes>(tot_stop - tot_start);
 progress.write(1,tot_duration.count(),true,"minutes");
  delete[] surf;
 cout << endl<<"doCalculations: total, bad = " << setw(12) << Nelem << setw(12) << nBadElem << endl;
 cout << "number of elements*pT configurations where nf>1.0: " << nFermiFail
  << endl;
 cout << "event_plane_vectors: " << Qx1 << "  " << Qy1 << "  "
   << Qx2 << "  " << Qy2 << endl;
}

void calcInvariantQuantities() {
 const double gmumu[4] = {1., -1., -1., -1.};
 const double tvect[4] = {1.,0., 0., 0.};
 int nBadElem = 0;
 for (int iel = 0; iel < Nelem; iel++) {  // loop over all elements
  if(fabs(surf[iel].dbeta[0][0])>1000.0) nBadElem++;
  double symm_deriv = 0.0, asymm_deriv = 0.0, mod_deriv = 0.0;
  for(int mu=0; mu<4; mu++)
   for(int nu=0; nu<4; nu++) {
    symm_deriv += 0.25*pow(hbarC*(surf[iel].dbeta[mu][nu]+surf[iel].dbeta[nu][mu]),2)
		    *gmumu[mu]*gmumu[nu];
    asymm_deriv += 0.25*pow(hbarC*(surf[iel].dbeta[mu][nu]-surf[iel].dbeta[nu][mu]),2)
		    *gmumu[mu]*gmumu[nu];
    mod_deriv += pow(hbarC*(surf[iel].dbeta[mu][nu]),2)
		    *gmumu[mu]*gmumu[nu];
    //for(int rh=0; rh<4; rh++)
     //for(int sg=0; sg<4; sg++) {
      //// computing the 'standard' polarization expression
      //Pi_num[ipt][iphi][mu] += pds * nf * (1. - nf) * levi(mu, nu, rh, sg)
                              //* p_[sg] * surf[iel].dbeta[nu][rh];
      //// computing the extra 'xi' term for the polarization
      //for(int ta=0; ta<4; ta++)
       //Pi_num_xi[ipt][iphi][mu] += pds * nf * (1. - nf) * levi(mu, nu, rh, sg)
                   //* p_[sg] * p[ta] / p[0] * tvect[nu]
                   //* ( surf[iel].dbeta[rh][ta] + surf[iel].dbeta[ta][rh]);
     //} // mu-nu-rho-sigma loop
    } // mu-nu loop
    if(fabs(surf[iel].eta)<0.5) {
     //sqroot1 = sqrt(sqroot1);
     if(symm_deriv!=symm_deriv)
      cout << "symm_deriv=nan\n";
     histSymm->Fill(symm_deriv);
     histAsymm->Fill(asymm_deriv);
     histMod->Fill(mod_deriv);
   }
 }  // loop over all elements
 histSymm->Scale(1./histSymm->Integral(),"width");
 histAsymm->Scale(1./histAsymm->Integral(),"width");
 histMod->Scale(1./histMod->Integral(),"width");
 plotSymm->cd();
 histSymm->Draw("H");
 plotAsymm->cd();
 histAsymm->Draw("H");
 plotMod->cd();
 histMod->Draw("H");
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
  if(fabs(surf[iel].dbeta[0][0])>1000.0) nBadElem++;
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
    fout << setw(14) << Pi_num[ipt][iphi][mu] * hbarC ;// / (8.0 * particle->GetMass()); THE MASS IS NOW TAKEN INTO ACCOUNT IN doCalculations
   for(int mu=0; mu<4; mu++)
    fout << setw(14) << Pi_num_exact[ipt][iphi][mu] * hbarC; // / (8.0 * particle->GetMass()); 
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


void Polfromdecay(int pdg_mother, int pdg_second_son,char *out_file){
	//Computes the contribution to the lambda polarization from the decay of pdg_mother. The result is written as a function of p in out_file
	ofstream fout(out_file);
	 if (!fout) {
		 cout << "I/O error with " << out_file << endl;
		 exit(1);
	 }
	const double gmumu[4] = {1., -1., -1., -1.};
	const double tvect[4] = {1.,0., 0., 0.};
	particle = database->GetPDGParticle(3122); // son of which to compute the polarization
	const double mass = particle->GetMass(); 
	
	//fetch information concerning mother and second son
	ParticlePDG2 *mother = database->GetPDGParticle(pdg_mother); //By default only Sigma* and Sigma0 are included
	const double mother_mass = mother->GetMass();  
	const double mother_baryonCharge = mother->GetBaryonNumber();
	const double mother_electricCharge = mother->GetElectricCharge();
	const double mother_strangeness = mother->GetStrangeness();
	
	ParticlePDG2 *second_son = database->GetPDGParticle(pdg_second_son);
	const double second_son_mass = second_son->GetMass();  
	//////////////////////////////////////////////////////////
	if(pdg_mother==3212 && pdg_second_son==22 ||pdg_mother==3224 && pdg_second_son==211)//ALLOWED DECAYS
	cout << "calculations for the decay: " << mother->GetName() << " to " << particle->GetName() << " and " << second_son->GetName() << endl;
	
	else{
		cout<<"ERROR: I don't know about the decay " <<mother->GetName() << " to "<< particle->GetName() << " and " << second_son->GetName()<< "!"<<endl;
		exit(1);
		}
	
	int nFermiFail = 0; // number of elements where nf>1.0
	int nBadElem = 0;
	progress_bar progress(std::clog, 70u, "Working");
    double counter = 0;
    auto tot_start = high_resolution_clock::now();
	//Loop over Lambda momenta
	//What follows computes formula 31 of 1905.03123 as a function of pt and theta for partilces at midrapidity. The values of pt and theta are the same as in doCalculations
	for (int ipt = 0; ipt < pT.size(); ipt++){
		for (int iphi = 0; iphi < phi.size(); iphi++) {
			double result_den=0;
			double result_num[3]={0,0,0};
			auto start = high_resolution_clock::now();
			
			//Loop for integral in solid angle: theta
			double dangle = M_PI/10; 
			#pragma omp parallel for num_threads(5) reduction(+:result_den,result_num)
			for(int ith_int=1; ith_int<10; ith_int++){ 
				double ith = dangle*ith_int;
				double mT = sqrt(mass * mass + pT[ipt] * pT[ipt]); // Lambda energy in Labframe
				double p[4] = {mT, pT[ipt]*cos(phi[iphi]), pT[ipt]*sin(phi[iphi]), 0}; //Lambda momentum in Labframe
				//Variables usefull for the decays
				double E_Lambda_rest = (mother_mass*mother_mass + mass*mass - second_son_mass*second_son_mass)/(2*mother_mass); 
				double p_rest_abs = sqrt(E_Lambda_rest*E_Lambda_rest - mass*mass);
				double p_rest[4];
				double jacobian; //jacobian from eq. 30 of 1905.03123
				double Energy_mother;
				double P_mother[4];
				double P_mother_[4];
				double integralphi_den=0;
				double integralphi_num[3]={0,0,0};
				
				//Loop for integral in solid angle: phi
				for(double iph=0; iph<2.0*M_PI+1e-5; iph+=dangle){
					double rest_frame_Pi[3] = {0,0,0}; //polarization in the rest frame of the mother
					double polarization_mother_tot[4] = {0,0,0,0}; //polarization from both vorticity and thermal shear
					double S_vect[3]={0,0,0};//the particular vector depending on the decay. For details 1905.03123
					double spectrum = 0;
					double polarization_mother[4]={0,0,0,0};
					double polarization_mother_Xi[4]={0,0,0,0};
					
					//the momentum of Lambda in mother rest frame. p_rest[0] is irrelevant for the caluclations
					p_rest[0] = 0;
					p_rest[1] = p_rest_abs*sin(ith)*cos(iph);
					p_rest[2] = p_rest_abs*sin(iph)*sin(ith);
					p_rest[3] = p_rest_abs*cos(ith); 
					
					jacobian = pow(mother_mass,3)*pow(mT + E_Lambda_rest,2)*
									(pow(mT + E_Lambda_rest,2)-(mass*mass+mT*E_Lambda_rest+p[1]*p_rest[1]+p[2]*p_rest[2]+p[3]*p_rest[3]))/
									(E_Lambda_rest*pow(mass*mass+mT*E_Lambda_rest+p[1]*p_rest[1]+p[2]*p_rest[2]+p[3]*p_rest[3],3));
					
					//momentum of the mother
					//upper indices
					P_mother[1] = (p[1]-p_rest[1])*2*mother_mass*(mT + E_Lambda_rest)/
								(pow(mT + E_Lambda_rest,2)-pow(p[1]-p_rest[1],2)-pow(p[2]-p_rest[2],2)-pow(p[3]-p_rest[3],2));
					P_mother[2] = (p[2]-p_rest[2])*2*mother_mass*(mT + E_Lambda_rest)/
								(pow(mT + E_Lambda_rest,2)-pow(p[1]-p_rest[1],2)-pow(p[2]-p_rest[2],2)-pow(p[3]-p_rest[3],2));
					P_mother[3] = (p[3]-p_rest[3])*2*mother_mass*(mT + E_Lambda_rest)/
								(pow(mT + E_Lambda_rest,2)-pow(p[1]-p_rest[1],2)-pow(p[2]-p_rest[2],2)-pow(p[3]-p_rest[3],2));
					Energy_mother = sqrt(mother_mass*mother_mass+pow(P_mother[1],2)+pow(P_mother[2],2)+pow(P_mother[3],2));
					P_mother[0] = Energy_mother; 
					//lower indices
					P_mother_[1] = -P_mother[1];
					P_mother_[2] = -P_mother[2];
					P_mother_[3] = -P_mother[3];
					P_mother_[0] = Energy_mother;
					
					for(int mu=0;mu<4;mu++){
						polarization_mother_tot[mu] = 0;
						}
					for(int mu=0;mu<3;mu++){
						rest_frame_Pi[mu] = 0;
						S_vect[mu] = 0;
						}
					//INTEGRAL OVER HYPERSURFACE	
					for (int iel = 0; iel < Nelem; iel++) {  		
						double mutot = surf[iel].mub * mother_baryonCharge
								+ surf[iel].muq * mother_electricCharge + surf[iel].mus * mother_strangeness;
						double pds = 0;
						double pu = 0;	

						for (int mu = 0; mu < 4; mu++) {
							 pds += P_mother[mu] * surf[iel].dsigma[mu];
							 pu += P_mother[mu] * surf[iel].u[mu] * gmumu[mu];
						}

						double nf =  c1 / (exp( (pu - mutot) / surf[iel].T) + 1.0); 
						if(nf > 1.0) nFermiFail++;
						spectrum += pds * nf ;
						for(int mu=0; mu<4; mu++)
							for(int nu=0; nu<4; nu++)
								for(int rh=0; rh<4; rh++)
									for(int sg=0; sg<4; sg++) {
										// computing the 'standard' polarization expression
										polarization_mother[mu] += pds * nf * (1. - nf) * levi(mu, nu, rh, sg)
																* P_mother_[sg] * surf[iel].dbeta[nu][rh];
										// computing the extra 'xi' term for the polarization
									if(nu==0)
									for(int ta=0; ta<4; ta++)
											polarization_mother_Xi[mu] += pds * nf * (1. - nf) * levi(mu, nu, rh, sg)
											* P_mother_[sg] * P_mother[ta] / P_mother[0] * tvect[nu]
											* ( surf[iel].dbeta[rh][ta] + surf[iel].dbeta[ta][rh]);
									}		
						}//end hypersurface integral
									
						for(int mu=0;mu<4;mu++){
							polarization_mother_tot[mu]+=hbarC*(polarization_mother[mu]-polarization_mother_Xi[mu])/(8*mother_mass);
						}
										
						//boost to the rest frame of the mother: the inverse standard boost is given by {{e/m, -(px/m), -(py/m), -(pz/m)}, {-(px/m), 1 + px^2/(e m + m^2), (px py)/(e m + m^2), (px pz)/(e m + m^2)}, {-(py/m), (px py)/(e m + m^2), 1 + py^2/(e m + m^2), (py pz)/(e m + m^2)}, {-(pz/m), (px pz)/(e m + m^2), (py pz)/(e m + m^2), 1 + pz^2/(e m + m^2)}}
						rest_frame_Pi[0] = -((polarization_mother_tot[0]*P_mother[1])/mother_mass) + 
								polarization_mother_tot[1]*(1 + pow(P_mother[1],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
								(polarization_mother_tot[2]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
								(polarization_mother_tot[3]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
								
						rest_frame_Pi[1] = -((polarization_mother_tot[0]*P_mother[2])/mother_mass) + 
								polarization_mother_tot[2]*(1 + pow(P_mother[2],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
								(polarization_mother_tot[1]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
								(polarization_mother_tot[3]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
					
						rest_frame_Pi[2] = -((polarization_mother_tot[0]*P_mother[3])/mother_mass) + 
								polarization_mother_tot[3]*(1 + pow(P_mother[3],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
								(polarization_mother_tot[1]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
								(polarization_mother_tot[2]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
						
						//computation of the vector S depending on the decay
						switch(pdg_mother){
							case 3212: //Sigma0 to Lambda and photon
							for(int mu=0;mu<3;mu++){ 
								S_vect[mu] = -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_Pi[0]*p_rest[1]+rest_frame_Pi[1]*p_rest[2]+rest_frame_Pi[2]*p_rest[3]); 
								}
							break;
							case 3224: //Sigma* to Lambda and pion
							for(int mu=0;mu<3;mu++){
								S_vect[mu] =2*rest_frame_Pi[mu] -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_Pi[0]*p_rest[1]+rest_frame_Pi[1]*p_rest[2]+rest_frame_Pi[2]*p_rest[3]); 
								}
							break;
							default:
							cout <<"I don't know about this decay!"<< endl;
							exit(1);
							}
						
						//compute integrals in phi and theta  
						integralphi_den += dangle*spectrum*jacobian/P_mother[0];
						for(int mu=0;mu<3;mu++){
							integralphi_num[mu] += dangle*jacobian*S_vect[mu]/P_mother[0];
							}
					} //end loop in iph
					//~ integraltheta_den += dangle*sin(ith)*integralphi_den;
					//~ for(int mu=0;mu<3;mu++){
						//~ integraltheta_num[mu] += dangle*sin(ith)*integralphi_num[mu];
						//~ }
					//~ result_den += integraltheta_den;
					//~ for(int mu=0;mu<3;mu++){
								//~ result_num[mu] += integraltheta_num[mu];
								//~ }
					result_den += dangle*sin(ith)*integralphi_den;
					for(int mu=0;mu<3;mu++){
							result_num[mu] += dangle*sin(ith)*integralphi_num[mu];
							}
				} //end loop in ith

				
				auto stop = high_resolution_clock::now();
				auto duration = duration_cast<seconds>(stop - start);
				/////////////////////////////////////////////////////////////////////////////////////////////////
				//write results in output file
				fout << setw(14) << pT[ipt] << setw(14) << phi[iphi]
				<< setw(14) << result_den << setw(14)  << result_num[0]
				<< setw(14) << result_num[1] << setw(14) << result_num[2] << endl;
				 if(ipt==0){ //for pT=0 the results are always the same, I write them all at once without recomputing the polarization. This saves 39 iterations.
					 for(int i=1; i<phi.size();i++){
							fout << setw(14) << pT[ipt] << setw(14) << phi[i]
								<< setw(14) << result_den << setw(14)  << result_num[0]
								<< setw(14) << result_num[1] << setw(14) << result_num[2] << endl;
						 }
						 iphi = phi.size()+1;
						 counter = phi.size()-1;
					 }
				progress.write((counter)/(pT.size()*phi.size()),duration.count());
				counter++;
	   }
   }
		auto tot_stop = high_resolution_clock::now();
		auto tot_duration = duration_cast<hours>(tot_stop - tot_start);
		progress.write(1,tot_duration.count(),true,"hours");
	    fout.close();
	}
	
void Pol_tablewriter(int pdg_mother,int pdg_second_son,char* out_file){
	//THIS FUNCTION WAS NOT CAREFULLY TESTED AND CAN BE CONSIDERED IN A BETA VERSION, USE IT AT YOUR OWN RISK!
	//write the polarization of the particle with id pdg_mother in the output file. This is done in view of the decay mother->lambda+second_son
	ofstream fout(out_file);
	 if (!fout) {
		 cout << "I/O error with " << out_file << endl;
		 exit(1);
	 }
	const double gmumu[4] = {1., -1., -1., -1.};
	const double tvect[4] = {1.,0., 0., 0.};
	particle = database->GetPDGParticle(3122);
	const double mass = particle->GetMass(); 
	
	//fetch information concerning mother and second son
	ParticlePDG2 *mother = database->GetPDGParticle(pdg_mother);
	const double mother_mass = mother->GetMass();  
	const double mother_baryonCharge = mother->GetBaryonNumber();
	const double mother_electricCharge = mother->GetElectricCharge();
	const double mother_strangeness = mother->GetStrangeness();
	
	ParticlePDG2 *second_son = database->GetPDGParticle(pdg_second_son);
	const double second_son_mass = second_son->GetMass();
	 cout << "calculations for the mother S vector in the decay: " << mother->GetName()<< " to "<<particle->GetName()<<" and "<<second_son->GetName() << endl;
	 int nFermiFail = 0; // number of elements where nf>1.0
	 int nBadElem = 0;
	 for (int ipt = 0; ipt < pT.size(); ipt++){
		for (int iphi = 0; iphi < phi.size(); iphi++) {
			auto start = high_resolution_clock::now();
			double mT = sqrt(mass * mass + pT[ipt] * pT[ipt]); // Lambda energy
			double p[4] = {mT, pT[ipt]*cos(phi[iphi]), pT[ipt]*sin(phi[iphi]), 0}; //Lambda momentum in Labframe
			double p_[4] = {mT, -pT[ipt]*cos(phi[iphi]), -pT[ipt]*sin(phi[iphi]), 0};
			//Variables usefull for the decays
			double E_Lambda_rest = (mother_mass*mother_mass + mass*mass - second_son_mass*second_son_mass)/(2*mother_mass);
			double p_rest_abs = sqrt(E_Lambda_rest*E_Lambda_rest - mass*mass);
			double Pimother[4]={0,0,0,0};
			double Pimother_Xi[4]={0,0,0,0};
			double Pimother_den=0;
			double p_rest[4];
			double P_motherx;
			double P_mothery;
			double P_motherz;
			double Energy_mother;
			double P_mother[4];
			double P_mother_[4];
			double spectrum = 0;


			double dangle = M_PI/10; 
			for(double ith=0; ith<M_PI-1e-5; ith+=dangle){
				for(double iph=0; iph<2.0*M_PI-1e-5; iph+=dangle){
					p_rest[0] = 0;
					p_rest[1] = p_rest_abs*sin(ith)*cos(iph);
					p_rest[2] = p_rest_abs*sin(iph)*sin(ith);
					p_rest[3] = p_rest_abs*cos(ith); //the momentum of Lambda in mother rest frame. p_rest[0] is unimportant
					P_mother[1] = (p[1]-p_rest[1])*2*mother_mass*(mT + E_Lambda_rest)/(pow(mT + E_Lambda_rest,2)-pow(p[1]-p_rest[1],2)-pow(p[2]-p_rest[2],2)-pow(p[3]-p_rest[3],2));
					P_mother[2] = (p[2]-p_rest[2])*2*mother_mass*(mT + E_Lambda_rest)/(pow(mT + E_Lambda_rest,2)-pow(p[1]-p_rest[1],2)-pow(p[2]-p_rest[2],2)-pow(p[3]-p_rest[3],2));
					P_mother[3] = (p[3]-p_rest[3])*2*mother_mass*(mT + E_Lambda_rest)/(pow(mT + E_Lambda_rest,2)-pow(p[1]-p_rest[1],2)-pow(p[2]-p_rest[2],2)-pow(p[3]-p_rest[3],2));
					Energy_mother = sqrt(mother_mass*mother_mass+pow(P_motherx,2)+pow(P_mothery,2)+pow(P_motherz,2));
					P_mother[0] = Energy_mother; //momentum of the mother
					P_mother_[1] = -P_mother[1];
					P_mother_[2] = -P_mother[2];
					P_mother_[3] = -P_mother[3];
					P_mother_[0] = Energy_mother;
					for(int mu=0;mu<4;mu++){
						Pimother[mu] = 0;
						Pimother_Xi[mu] = 0;
						}
					spectrum = 0;
					
					//compute the polarization for fixed momentum
					#pragma omp parallel for num_threads(5) reduction(+:spectrum,Pimother,Pimother_Xi)
					for (int iel = 0; iel < Nelem; iel++) {  
						if(fabs(surf[iel].dbeta[0][0])>1000.0) nBadElem++;
						double pds = 0;
						double pu = 0;
						double polarization_mother[4];
						double polarization_mother_Xi[4];
						for (int mu = 0; mu < 4; mu++) {
							 pds += P_mother[mu] * surf[iel].dsigma[mu];
							 pu += P_mother[mu] * surf[iel].u[mu] * gmumu[mu];
							 polarization_mother[mu]=0;
							 polarization_mother_Xi[mu]=0;
							}
						double mutot = surf[iel].mub * mother_baryonCharge
						+ surf[iel].muq * mother_electricCharge + surf[iel].mus * mother_strangeness;
						double nf =  c1 / (exp( (pu - mutot) / surf[iel].T) + 1.0); 
						if(nf > 1.0) nFermiFail++;
						spectrum += pds * nf ; 
						for(int mu=0; mu<4; mu++)
							for(int nu=0; nu<4; nu++)
								for(int rh=0; rh<4; rh++)
									for(int sg=0; sg<4; sg++) {
										// computing the 'standard' polarization expression
										polarization_mother[mu] += pds * nf * (1. - nf) * levi(mu, nu, rh, sg)
																* P_mother_[sg] * surf[iel].dbeta[nu][rh];
										// computing the extra 'xi' term for the polarization
										if(nu==0)
									for(int ta=0; ta<4; ta++)
											polarization_mother_Xi[mu] += pds * nf * (1. - nf) * levi(mu, nu, rh, sg)
											* P_mother_[sg] * P_mother[ta] / P_mother[0] * tvect[nu]
											* ( surf[iel].dbeta[rh][ta] + surf[iel].dbeta[ta][rh]);
									}		
									
									
									for(int mu=0;mu<4;mu++){
										Pimother[mu]+=polarization_mother[mu];
										Pimother_Xi[mu]+=-polarization_mother_Xi[mu];
										}				
						}//end loop over the hypersurface
						// boost to restframe
						double rest_frame_Pi[3], rest_frame_PiXi[3],S_vect[3],S_vectXi[3];
						
						rest_frame_Pi[0] = -((Pimother[0]*P_mother[1])/mother_mass) + 
								Pimother[1]*(1 + pow(P_mother[1],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
								(Pimother[2]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
								(Pimother[3]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
								
						rest_frame_Pi[1] = -((Pimother[0]*P_mother[2])/mother_mass) + 
								Pimother[2]*(1 + pow(P_mother[2],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
								(Pimother[1]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
								(Pimother[3]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
					
						rest_frame_Pi[2] = -((Pimother[0]*P_mother[3])/mother_mass) + 
								Pimother[3]*(1 + pow(P_mother[3],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
								(Pimother[1]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
								(Pimother[2]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
										
						rest_frame_PiXi[0] = -((Pimother_Xi[0]*P_mother[1])/mother_mass) + 
								Pimother_Xi[1]*(1 + pow(P_mother[1],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
								(Pimother_Xi[2]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
								(Pimother_Xi[3]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
								
						rest_frame_PiXi[1] = -((Pimother_Xi[0]*P_mother[2])/mother_mass) + 
								Pimother_Xi[2]*(1 + pow(P_mother[2],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
								(Pimother_Xi[1]*P_mother[1]*P_mother[2])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
								(Pimother_Xi[3]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
					
						rest_frame_PiXi[2] = -((Pimother_Xi[0]*P_mother[3])/mother_mass) + 
								Pimother_Xi[3]*(1 + pow(P_mother[3],2)/(P_mother[0]*mother_mass + pow(mother_mass,2))) + 
								(Pimother_Xi[1]*P_mother[1]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2)) + 
								(Pimother_Xi[2]*P_mother[2]*P_mother[3])/(P_mother[0]*mother_mass + pow(mother_mass,2));
						for(int i=0;i<3;i++){// NOT SURE ABOUT THE "-" FOR THE XI TERM
							rest_frame_Pi[i]*=hbarC/(8*mother_mass*spectrum);
							rest_frame_PiXi[i]*=-hbarC/(8*mother_mass*spectrum);
							}
						
						//computation of the vector S depending on the decay
						switch(pdg_mother){
							case 3212: //Sigma0 to Lambda and photon
							for(int mu=0;mu<3;mu++){ 
								S_vect[mu] = -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_Pi[0]*p_rest[1]+rest_frame_Pi[1]*p_rest[2]+rest_frame_Pi[2]*p_rest[3]); 
								S_vectXi[mu] = -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_PiXi[0]*p_rest[1]+rest_frame_PiXi[1]*p_rest[2]+rest_frame_PiXi[2]*p_rest[3]); 
								}
							break;
							case 3224: //Sigma* to Lambda and pion
							for(int mu=0;mu<3;mu++){
								S_vect[mu] =2*rest_frame_Pi[mu] -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_Pi[0]*p_rest[1]+rest_frame_Pi[1]*p_rest[2]+rest_frame_Pi[2]*p_rest[3]); 
								S_vectXi[mu] =2*rest_frame_PiXi[mu] -(p_rest[mu+1]/pow(p_rest_abs,2))*(rest_frame_PiXi[0]*p_rest[1]+rest_frame_PiXi[1]*p_rest[2]+rest_frame_PiXi[2]*p_rest[3]); 
								}
							break;
							default:
							cout <<"I don't know about this decay!"<< endl;
							exit(1);
							}
						//write table: pt phi theta phi Pi_rest[mu] Pi_Xi_rest[mu]
						fout<<pT[ipt]<<setw(14)<<phi[iphi]<<setw(14) << ith << setw(14) << iph
				<< setw(14)   << S_vect[0]<< setw(14)  << S_vect[1]<<setw(14)  << S_vect[2]
				<<setw(14)  << S_vectXi[0]<<setw(14)  << S_vectXi[1]<<setw(14)  << S_vectXi[2]<< endl;
					}//end loop iph
				}//end loop ith
				auto stop = high_resolution_clock::now();
				auto duration = duration_cast<seconds>(stop - start);
				cout<< ipt+1<<" of " << pT.size()<<" and "<<iphi+1<<" of "<< phi.size()<< " completed in "<< duration.count()<<" seconds"<<endl;
	
	}
}
}

}  // end namespace gen
