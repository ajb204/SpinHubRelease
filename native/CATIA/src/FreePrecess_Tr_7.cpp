//
// June 3, 2009:
//    Modified by G. Bouvignies:
//                  - to take into account the temperature dependance of some variables;
//                  - to use function for calculating kex, pb and delatO.
//

#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>
namespace ublas = boost::numeric::ublas;
using boost::to_lower_copy;

void Catia::FreePrecess_Tr_7(ublas::matrix<double>& G,Dataset& dset, int Atom){
  /*
    NameOfBasis: Iph_7

    The basis is:
    { E,
      [ Tr(x) ,Tr(y) ,Tr(z)] x { Site A, Site B }
    }

    Thus - a 7 element basis
  */
  std::vector<std::string> basis=dset.vpar("basis");
  //
  G.resize(7,7);
  G.clear();
  // Now fetch the parameters from the different Dataset/Atom parameters
  //
  const double kex = Kex(dset);
  const double pb = Pb(dset);
  //
  std::string Nucl=dset._nucleus;
  double gammaS=0.;
  double gammaI=0.;
  double delta_csa=0;
  double r_is=0.;
  double hbar=6.626075e-34;
  if(Nucl=="N"){
    gammaS=_gammaN;
    delta_csa=-172E-6;
    r_is=1.02E-10;
  } else if (Nucl=="C"){
    gammaS=_gammaC;
    delta_csa=-25E-6;
    r_is=1.098e-10; // From Kowalewski, JMR, 2002, p171-177
  } else if (Nucl=="H"){
    gammaS=_gammaH;
    delta_csa=-10E-6;
    r_is=1.02E-10;
  } else if (Nucl=="D"){
    gammaS=_gammaD;
    delta_csa=0.;
    r_is=1.0E-10;
  } else {
    std::cerr<<" Could not resolve the nucleus of type "<<Nucl<<"\n";
    std::cerr<<" in Dataset "<<dset._id<<"\n";
    std::cerr<<" Functions: FreePrecess_Tr_7()\n";
    Abort(1);
  };
  //Store it for later use
  dset._gamma=gammaS;
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
  RequiredParam.push_back("R1iph"+marker);
  RequiredParam.push_back("R1aph"+marker);
  RequiredParam.push_back("R0Tr"+marker);
  RequiredParam.push_back("Omega"+temperatureMarker);
  //RequiredParam.push_back("JIS");
  RequiredParam.push_back("DeltaJ");
  for(unsigned int i=0;i<RequiredParam.size();i++){
    if(!(BResolveParam(LocalParam[Atom],RequiredParam[i]))){
      std::cerr<<" The parameters:"<<RequiredParam[i]<<" is required by the basisset";
      std::cerr<<" "<<basis[1]<<"\n but is not provided for atom "<<AtomNumber2AtomName(Atom)<<"\n";
      std::cerr<<" please provide in the LocalParameter set\n";
      std::cerr<<" Function FreePrecess_Tr_7()\n"<<std::endl;
      Abort(1);
    };
  };
  if( !(
	 BResolveParam(LocalParam[Atom],"JIS") ||
	 BResolveParam(LocalParam[Atom],"JIS"+marker)
	)){
    std::cerr<<" The parameter JIS or JIS"<<marker<<" is required by the basisset";
    std::cerr<<" "<<basis[1]<<"\n but is not provided for atom "<<AtomNumber2AtomName(Atom)<<"\n";
    std::cerr<<" please provide one of these parameters in the LocalParameter set\n";
    std::cerr<<" Function FreePrecess_Tr_7()\n"<<std::endl;
    Abort(1);
  };
  double R1iph=ResolveParam(LocalParam[Atom],"R1iph"+marker);
  double R1aph=ResolveParam(LocalParam[Atom],"R1aph"+marker);
  double R0Tr=ResolveParam(LocalParam[Atom],"R0Tr"+marker);
  double JIS=0.;
  if(BResolveParam(LocalParam[Atom],"JIS"+marker)){
    JIS=ResolveParam(LocalParam[Atom],"JIS"+marker);
  } else {
    JIS=ResolveParam(LocalParam[Atom],"JIS");
  };
  //  double JIS=ResolveParam(LocalParam[Atom],"JIS");
  double DeltaJ=ResolveParam(LocalParam[Atom],"DeltaJ");
  if(dset.Bvpar("deltajscaling")){
    double djs=atof(dset.vpar("deltajscaling")[0].c_str());
    DeltaJ=DeltaJ*djs;
    sprintf(line, "Scaling used for DeltaJ, i.e., DeltaJ(%s)=%g*DeltaJ", dset._id.c_str(),djs);
    std::string note(line);
    ClearBuf(line,sizeof(line));
    LocalNotes[Atom]["DeltaJ@"+dset._id]=note;
  };
  double Omega=(gammaS/_gammaH)*dset._sfrq*2*DFH_PI*
    (ResolveParam(LocalParam[Atom],"Omega"+temperatureMarker)-dset._xcar);
  double DeltaO=(gammaS/_gammaH)*dset._sfrq*2*DFH_PI*
    DeltaOmega(dset,Atom);
  //
  gammaI=_gammaH;
  double EtaZ=0;
  bool estimate_etaZ=false;
  for(unsigned int i=0;i<basis.size();i++){
    if(to_lower_copy(basis[i])=="estimate_etaz"){
      estimate_etaZ=true;
      break;
    };
  };
  if(estimate_etaZ){
    //We estimate from R1iph and assume
    //CSA of -172 for Nitrogen, 25 for Carbon, and -10ppm for H
    double B0=dset._sfrq*1E6*2.*DFH_PI/_gammaH;
    double cc=B0*gammaS*delta_csa/pow(3.,0.5);
    double dd=(1.e-7)*hbar*gammaI*gammaS*pow(r_is,-3.)/(DFH_PI*2.);
    double JwS=R1iph/(cc*cc+dd*dd*0.75);
    EtaZ=-pow(3.,0.5)*cc*dd*JwS;
  };
  //
  // Approximate value for equlibrium value
  // site A
  //
  // Auto Relaxations
  G(1,1) = R0Tr;
  G(2,2) = R0Tr;
  G(3,3) = (R1iph+R1aph)/2.-EtaZ;
  G(3,0) = -(1-pb)*(R1iph-EtaZ )*gammaS/(gammaI*2.);
  //
  // Omega
  G(1,2) = -JIS*DFH_PI+(Omega);
  G(2,1) =  JIS*DFH_PI-(Omega);

  //
  // site B
  //
  G(4,4) = R0Tr;
  G(5,5) = R0Tr;
  G(6,6) = (R1iph+R1aph)/2.-EtaZ;
  G(6,0) = -pb*(R1iph-EtaZ )*gammaS/(gammaI*2.);
  //
  G(4,5) = -(JIS+DeltaJ)*DFH_PI+(Omega+DeltaO);
  G(5,4) =  (JIS+DeltaJ)*DFH_PI-(Omega+DeltaO);
  //
  // kex additions
  for(unsigned int i=0;i<3;i++){
    G(i+1,i+1) = G(i+1,i+1)+kex*pb;
    G(i+4,i+1) = G(i+4,i+1)-kex*pb;
    //
    G(i+4,i+4) = G(i+4,i+4)+kex*(1-pb);
    G(i+1,i+4) = G(i+1,i+4)-kex*(1-pb);
  };
  //G.print("%10.3f");
  //exit(10);
};
