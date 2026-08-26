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

void Catia::FreePrecess_N_NH_13(ublas::matrix<double>& G,Dataset& dset, int Atom){
  /*
    NameOfBasis: N_NH_13

    The basis is:
    { E,
    [ Tr(x) ,Tr(y) ,Tr(z),
     ATr(x),ATr(y),ATr(z)
    ] x { Site A, Site B }
    }

    where Tr(x) is the Trosy x component = Ix-2IxSz
    ATr(x) is the AntiTrosy x component = Ix+2IxSz

    Thus - a 13 element basis
  */
  std::vector<std::string> basis=dset.vpar("basis");
  //
  G.resize(13,13);
  G.clear();
  // Now fetch the parameters from the different Dataset/Atom parameters
  const double kex = Kex(dset);
  const double pb = Pb(dset);
  //
  std::string Nucl=dset._nucleus;
  double gammaS=0.;
  double gammaI=0.;
  double r_is=0.;
  double hbar=6.626075e-34/(2.*DFH_PI);
  if(Nucl=="N"){
    gammaS=_gammaN;
    r_is=1.02E-10;
  } else {
    std::cerr<<" The basis set: N_NH_13 only allows nitrogen as the primary nucleus"<<std::endl;
    std::cerr<<" Function: .FreePrecess_N_NH_13()\n;"<<std::endl;
    Abort(1);
  };
  //Store it for later use
  dset._gamma=gammaS;
  //
  std::string cn=to_lower_copy(dset.vpar("couplednucleus")[1]);
  if(cn.find('h')<cn.npos){
    gammaI=_gammaH;
  } else {
    std::cerr<<" The basis set: N_NH_13 only allows hidrogen as the coupled nucleus"<<std::endl;
    std::cerr<<" Function: .FreePrecess_N_NH_13()\n;"<<std::endl;
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
  RequiredParam.push_back("R1aph"+marker); // 1/s
  RequiredParam.push_back("DeltaSF");           // 1/s  This is the change in Spin Flip rate
  RequiredParam.push_back("CSA");               // ppm
  RequiredParam.push_back("S2tc");              // nanoseconds - This is indeed (5/8){4J(0)}=(5/8){4*2/5*S2tc }
  RequiredParam.push_back("Omega"+temperatureMarker);             // ppm             // ppm
  RequiredParam.push_back("DeltaS2tc");         // nanoseconds
  RequiredParam.push_back("DeltaJ");            // Hertz
  for(unsigned int i=0;i<RequiredParam.size();i++){
    if(!(BResolveParam(LocalParam[Atom],RequiredParam[i]))){
      std::cerr<<" The parameters:"<<RequiredParam[i]<<" is required by the basisset";
      std::cerr<<" "<<basis[1]<<"\n but is not provided for atom "<<AtomNumber2AtomName(Atom)<<"\n";
      std::cerr<<" please provide in the LocalParameter set\n";
      std::cerr<<" Function .FreePrecess_N_NH_13()\n"<<std::endl;
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
    std::cerr<<" Function .FreePrecess_N_NH_13()\n"<<std::endl;
    Abort(1);
  };
  //
  double R1iph_a=ResolveParam(LocalParam[Atom],"R1iph"+marker);   // 15N R1 (1/s)    \ These are population
  double R1aph_a=ResolveParam(LocalParam[Atom],"R1aph"+marker);   // R1(2NzHz) (1/s) / weighted averages
  double CSA=ResolveParam(LocalParam[Atom],"CSA");    // CSA (ppm)
  double S2tc=ResolveParam(LocalParam[Atom],"S2tc");      // Order parameter x correlation time (no unit)
  double DeltaS2tc=ResolveParam(LocalParam[Atom],"DeltaS2tc"); //Delta(orderparameter) (no unit)
  double Omega=((gammaS)/_gammaH)*dset._sfrq*2*DFH_PI*  //offset -> 'OMEGA' is in ppm
    (ResolveParam(LocalParam[Atom],"Omega"+temperatureMarker)-dset._xcar);
  double DeltaO=((gammaS)/_gammaH)*dset._sfrq*2*DFH_PI* // Change in chemical shift, B-A (ppm)
    DeltaOmega(dset,Atom);
  double DeltaSF=ResolveParam(LocalParam[Atom],"DeltaSF"); // change in SpinFlip rate (B-A).
  double JIS=0.;                                      // Scalar coupling
  if(BResolveParam(LocalParam[Atom],"JIS"+marker)){
    JIS=ResolveParam(LocalParam[Atom],"JIS"+marker);
  } else {
    JIS=ResolveParam(LocalParam[Atom],"JIS");
  };
  //
  double DeltaJ=ResolveParam(LocalParam[Atom],"DeltaJ"); //Change in coupling (scalar or dipolar)
  if(dset.Bvpar("deltajscaling")){
    double djs=atof(dset.vpar("deltajscaling")[1].c_str());
    DeltaJ=DeltaJ*djs;
    sprintf(line, "Scaling used for DeltaJ, i.e., DeltaJ(%s)=%g*DeltaJ", dset._id.c_str(),djs);
    std::string note(line);
    ClearBuf(line,sizeof(line));
    LocalNotes[Atom]["DeltaJ@"+dset._id]=note;
  };
  //
  // We follow the numeclarture from the 'exchange-free' paper.
  double B0=dset._sfrq*1E6*2.*DFH_PI/_gammaH;
  double dd=(1.e-7)*hbar*_gammaH*_gammaN*pow(r_is,-3.);
  double cc=B0*_gammaN*CSA*1e-6/sqrt(3.);
  //
  double phi_CSA_DD= 22.*DFH_PI/180.; // angle between CSA and N-H vector is assumed 22 degrees
  //
  //
  //
  double R1iph[2]={ R1iph_a,R1iph_a};
  double R1aph[2]={ R1aph_a - pb*DeltaSF,
		    R1aph_a + (1-pb)*DeltaSF
  };
  double JwS=R1iph[0]/(cc*cc+dd*dd*0.75);
  double EtaZ[2] ={ -pow(3.,0.5)*cc*dd*0.5*(3*pow(cos(phi_CSA_DD),2.)-1.)*JwS,
		    -pow(3.,0.5)*cc*dd*0.5*(3*pow(cos(phi_CSA_DD),2.)-1.)*JwS
  }; // we assume that this does not change between gound and excited state
  //
  double EtaXY[2]={ -sqrt(3.)/6*cc*dd*  0.5*(3*pow(cos(phi_CSA_DD),2.)-1.) *( (8./5.)*S2tc*1e-9 + 3.*JwS),
		    -sqrt(3.)/6*cc*dd*  0.5*(3*pow(cos(phi_CSA_DD),2.)-1.) *( (8./5.)*(S2tc+DeltaS2tc)*1e-9 + 3.*JwS)
  };
  //
  // Matrix elements
  double R0T[2] ={ (dd*dd/8.+cc*cc/6.)*( (8./5.)*S2tc*1e-9+ 3.*JwS )           -EtaXY[0]+ 0.5*(R1aph[0]-R1iph[0]),
		   (dd*dd/8.+cc*cc/6.)*( (8./5.)*(S2tc+DeltaS2tc)*1e-9+3.*JwS) -EtaXY[1]+ 0.5*(R1aph[1]-R1iph[1])
  };
  double R0AT[2]={ (dd*dd/8.+cc*cc/6.)*( (8./5.)*S2tc*1e-9+ 3.*JwS )           +EtaXY[0]+ 0.5*(R1aph[0]-R1iph[1]),
		   (dd*dd/8.+cc*cc/6.)*( (8./5.)*(S2tc+DeltaS2tc)*1e-9+3.*JwS) +EtaXY[1]+ 0.5*(R1aph[1]-R1iph[1])
  };
  double R1T[2]={ (R1iph[0]+R1aph[0])/2.-EtaZ[0],
		  (R1iph[1]+R1aph[1])/2.-EtaZ[1]
  };
  double R1AT[2]={(R1iph[0]+R1aph[0])/2.+EtaZ[0],
		  (R1iph[1]+R1aph[1])/2.+EtaZ[1]
  };
  // Auto Relaxations
  //
  G(1,1) = R0T[0];
  G(2,2) = R0T[0];
  G(3,3) = R1T[0];
  G(3,0) = -(1-pb)*(R1iph[0]-EtaZ[0])*gammaS/(_gammaH*2.);
  G(4,4) = R0AT[0];
  G(5,5) = R0AT[0];
  G(6,6) = R1AT[0];
  G(6,0) = -(1-pb)*(R1iph[0]+EtaZ[0] )*gammaS/(_gammaH*2.);
  //
  G(4,1) = (R1iph[0]-R1aph[0])/2.;
  G(5,2) = (R1iph[0]-R1aph[0])/2.;
  G(6,3) = (R1iph[0]-R1aph[0])/2.;
  //
  G(1,4) = (R1iph[0]-R1aph[0])/2.;
  G(2,5) = (R1iph[0]-R1aph[0])/2.;
  G(3,6) = (R1iph[0]-R1aph[0])/2.;
  //
  // J couplings
  G(1,2) =  (JIS)*DFH_PI+Omega;
  G(2,1) = -(JIS)*DFH_PI-Omega;
  G(4,5) = -(JIS)*DFH_PI+Omega;
  G(5,4) =  (JIS)*DFH_PI-Omega;
  //
  // site B
  //
  // Auto Relaxations
  G(7,7) = R0T[1];
  G(8,8) = R0T[1];
  G(9,9) = R1T[1];
  G(9,0) = -pb*gammaS*(R1iph[1]-EtaZ[1])/(2.*_gammaH);
  G(10,10) = R0AT[1];
  G(11,11) = R0AT[1];
  G(12,12) = R1AT[1];
  G(12,0) = -pb*gammaS*(R1iph[1]+EtaZ[1])/(2.*_gammaH);
  //
  G(7,10) = (R1iph[1]-R1aph[1])/2.;
  G(8,11) = (R1iph[1]-R1aph[1])/2.;
  G(9,12) = (R1iph[1]-R1aph[1])/2.;
  //
  G(10,7) = (R1iph[1]-R1aph[1])/2.;
  G(11,8) = (R1iph[1]-R1aph[1])/2.;
  G(12,9) = (R1iph[1]-R1aph[1])/2.;
  //
  // J couplings
  G( 7,8) =  (JIS+DeltaJ)*DFH_PI+(Omega+DeltaO);
  G( 8,7) = -(JIS+DeltaJ)*DFH_PI-(Omega+DeltaO);
  G(10,11) = -(JIS+DeltaJ)*DFH_PI+(Omega+DeltaO);
  G(11,10) =  (JIS+DeltaJ)*DFH_PI-(Omega+DeltaO);
  //
  // kex additions
  for(unsigned int i=0;i<6;i++){
    G(i+1,i+1) = G(i+1,i+1)+kex*pb;
    G(i+7,i+1) = G(i+7,i+1)-kex*pb;
    //
    G(i+7,i+7) = G(i+7,i+7)+kex*(1-pb);
    G(i+1,i+7) = G(i+1,i+7)-kex*(1-pb);
  };
  /*
  std::cerr<<" R1iph "<<R1iph[0]<<"\t"<<R1iph[1]<<"\n";
  std::cerr<<" R1aph "<<R1aph[0]<<"\t"<<R1aph[1]<<"\n";
  std::cerr<<" EtaZ  "<<EtaZ[0] <<"\t"<<EtaZ[1] <<"\n";
  std::cerr<<" EtaXY "<<EtaXY[0]<<"\t"<<EtaXY[1]<<"\n";
  std::cerr<<" R0T   "<<R0T[0]  <<"\t"<<R0T[1]  <<"\n";
  std::cerr<<" R0AT  "<<R0AT[0] <<"\t"<<R0AT[1] <<"\n";
  std::cerr<<" (dd/8) S2tc "<<(dd*dd/8.)*(8./5.)*S2tc*1e-9<<std::endl;
  std::cerr<<" (cc/6) S2tc "<<(cc*cc/6.)*(8./5.)*S2tc*1e-9<<std::endl;
  std::cerr<<" dd          "<<dd<<std::endl;
  std::cerr<<" cc          "<<cc<<std::endl;
  std::cerr<<std::endl;
  G.print("%10.3f");
  exit(10);
  */

};
