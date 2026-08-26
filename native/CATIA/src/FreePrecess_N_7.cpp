//
// June 3, 2009:
//    Modified by G. Bouvignies:
//                  - to take into account the temperature dependance of some variables;
//                    TODO: needs to be checked by Flemming for the non-trivial variables
//                  - to use function for calculating kex, pb and delatO.
//

#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>
namespace ublas = boost::numeric::ublas;
using boost::to_lower_copy;

void Catia::FreePrecess_N_7(ublas::matrix<double>& G,Dataset& dset, int Atom){
  /*
    NameOfBasis: N_7

    The basis is:
    { E,
    [ Iph(x) ,Iph(y) ,Iph(z)
      ] x { Site A, Site B }
    }

    An inphase basis, where relaxation rates are calculated from the
    S2tc (spectral density function) and other spectral parameters

  */
  std::vector<std::string> basis=dset.vpar("basis");
  //
  G.resize(7,7);
  G.clear();
  // Now fetch the parameters from the different Dataset/Atom parameters
  const double kex = Kex(dset);
  const double pb = Pb(dset);
  //
  std::string Nucl=dset._nucleus;
  double gammaS=0.;
  double gammaI=0.;
  double r_is=0.;
  double hbar=6.626075e-34;
  if(Nucl=="N"){
    gammaS=_gammaN;
    r_is=1.02E-10;
  } else {
    std::cerr<<" The basis set: N_7 only allows nitrogen as the primary nucleus"<<std::endl;
    std::cerr<<" Function: .FreePrecess_N_7()\n;"<<std::endl;
    Abort(1);
  };
  //Store it for later use
  dset._gamma=gammaS;
  //
  std::string cn=to_lower_copy(dset.vpar("couplednucleus")[1]);
  if(cn.find('h')<cn.npos){
    gammaI=_gammaH;
  } else {
    std::cerr<<" The basis set: N_7 only allows hidrogen as the coupled nucleus"<<std::endl;
    std::cerr<<" Function: .FreePrecess_N_7()\n;"<<std::endl;
    Abort(1);
  };
  //
  char line[MAX_STRING_LENGTH];
  sprintf(line, "_%.0f", dset._sfrq);
  std::string fieldMarker(line);
  ClearBuf(line,sizeof(line));
  //
  std::string temperatureMarker("");
  if(_multipleTemperatures) {
    sprintf(line, "_%.1f", dset._temperature);
    temperatureMarker = line;
    ClearBuf(line,sizeof(line));
  }
  std::string marker(fieldMarker+temperatureMarker);
  //
  //Check that the parameters are there.
  std::vector<std::string> RequiredParam;
  RequiredParam.push_back("R1iph"+marker); // 1/s
  RequiredParam.push_back("CSA");               // ppm
  RequiredParam.push_back("S2tc");              // nanoseconds - This is indeed (5/8){4J(0)+3J(wN)}=(5/8){4*2/5*S2tc + .. }
  RequiredParam.push_back("Omega"+temperatureMarker);             // ppm             // ppm
  RequiredParam.push_back("DeltaS2tc");         // nanoseconds
  for(unsigned int i=0;i<RequiredParam.size();i++){
    if(!(BResolveParam(LocalParam[Atom],RequiredParam[i]))){
      std::cerr<<" The parameters:"<<RequiredParam[i]<<" is required by the basisset";
      std::cerr<<" "<<basis[1]<<"\n but is not provided for atom "<<AtomNumber2AtomName(Atom)<<"\n";
      std::cerr<<" please provide in the LocalParameter set\n";
      std::cerr<<" Function .FreePrecess_N_NH_13()\n"<<std::endl;
      Abort(1);
    };
  };
  //
  double R1iph_a=ResolveParam(LocalParam[Atom],"R1iph"+marker);   // 15N R1 (1/s)    \ These are population
  double CSA=ResolveParam(LocalParam[Atom],"CSA");    // CSA (ppm)
  double S2tc=ResolveParam(LocalParam[Atom],"S2tc");      // Order parameter x correlation time (no unit)
  double DeltaS2tc=ResolveParam(LocalParam[Atom],"DeltaS2tc"); //Delta(orderparameter) (no unit)
  double Omega=((gammaS)/_gammaH)*dset._sfrq*2*DFH_PI*  //offset -> 'OMEGA' is in ppm
    (ResolveParam(LocalParam[Atom],"Omega"+temperatureMarker)-dset._xcar);
  double DeltaO=((gammaS)/_gammaH)*dset._sfrq*2*DFH_PI* // Change in chemical shift, B-A (ppm)
    DeltaOmega(dset,Atom);
  //
  // We follow the numeclarture from the 'exchange-free' paper.
  double B0=dset._sfrq*1E6*2.*DFH_PI/_gammaH;
  double cc=B0*_gammaN*CSA*1e-6/sqrt(3.);
  double dd=(1.e-7)*hbar*_gammaH*_gammaN*pow(r_is,-3.)/(DFH_PI*2.);
  //
  //
  double R1iph[2]={ R1iph_a,R1iph_a};
  //
  double JwS=R1iph[0]/(cc*cc+dd*dd*0.75);
  //
  double R0iph[2] ={ (dd*dd/8.+cc*cc/6.)*( (8./5.)*S2tc*1e-9 + 3.*JwS),
  		     (dd*dd/8.+cc*cc/6.)*( (8./5.)*(S2tc+DeltaS2tc)*1e-9 + 3.*JwS)
  };
  // Relaxations
  G(1,1) = R0iph[0];
  G(2,2) = R0iph[0];
  G(3,3) = R1iph[0];
  G(3,0) = -(1-pb)*(R1iph[0])*gammaS/(_gammaH*2.);
  G(4,4) = R0iph[1];
  G(5,5) = R0iph[1];
  G(6,6) = R1iph[1];
  G(6,0) = -(1-pb)*(R1iph[1])*gammaS/(_gammaH*2.);
  //
  // Chemical shift
  G(1,2) =  (Omega);
  G(2,1) = -(Omega);
  G(4,5) =  (Omega+DeltaO);
  G(5,4) = -(Omega+DeltaO);

  //
  // kex additions
  for(unsigned int i=0;i<3;i++){
    G(i+1,i+1) = G(i+1,i+1)+kex*pb;
    G(i+4,i+1) = G(i+4,i+1)-kex*pb;
    //
    G(i+4,i+4) = G(i+4,i+4)+kex*(1-pb);
    G(i+1,i+4) = G(i+1,i+4)-kex*(1-pb);
  };
  /*
  std::cerr<<" R2iph "<<R0iph[0]<<"\t"<<R0iph[1]<<std::endl;
  G.print("%9.3f");
  exit(10);
  */
};
