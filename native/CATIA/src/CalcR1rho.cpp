#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>
#include <StringMethods.h>

#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

//#include <MatrixExponential.hpp>
#include <MatrixExponentialLapack.hpp>

using namespace boost::numeric::ublas;

namespace ublas = boost::numeric::ublas;

using ublas::prod;

typedef double real_t;
typedef std::complex<real_t> cmplx_t;
typedef ublas::vector<cmplx_t> cvec_t;
typedef ublas::matrix<real_t, ublas::row_major> rrmat_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> rcmat_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> ccmat_t;

/*
 This sequence is the seqfil="r1rh"

 Essentially it is build up of
Take Iz, apply a rotation, apply the spin lock on x, reverse the rotation, read off Iz.


 */

//#return the angle in radians
double Catia::GetTheta(double w1,double offset) 
{
  double th;
  if(w1==0)
    th=0.0;
  else if(offset==0 and w1!=0)
    th=90.0/180.0*DFH_PI;
  else
    th=atan(w1/fabs(offset));
  
  if(offset>0.0)
    th=th*-1.0;
  return th;
}

void Catia::FreePrecess_R1rho(ublas::matrix<double>& G,Dataset& dset, int Atom){
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
    delta_csa=-25E-6; //Calpha
    r_is=1.098e-10; // From Kowalewski, JMR, 2002, p171-177
  } else if (Nucl=="H"){
    gammaS=_gammaH;
    delta_csa=-10E-6;
    r_is=1.02E-10;
  } else if (Nucl=="D"){
    gammaS=_gammaD;
    delta_csa=0.;
    r_is=1.0E-10;
  } else if (Nucl=="F"){
    gammaS=_gammaF;
    delta_csa=0.;
    r_is=1.0E-10;
  } else {
    std::cerr<<" Could not resolve the nucleus of type "<<Nucl<<"\n";
    std::cerr<<" in Dataset "<<dset._id<<"\n";
    std::cerr<<" Functions: FreePrecess_Iph_7()\n";
    Abort(1);
  };
  //Store it for later use
  dset._gamma=gammaS;
  //
  char line[MAX_STRING_LENGTH];
  //char line[100];
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
  RequiredParam.push_back("R0iph"+marker);
  RequiredParam.push_back("Omega"+temperatureMarker);
  for(unsigned int i=0;i<RequiredParam.size();i++){
    if(!(BResolveParam(LocalParam[Atom],RequiredParam[i]))){
      std::cerr<<" The parameters:"<<RequiredParam[i]<<" is required by the basisset";
      std::cerr<<" "<<basis[1]<<"\n but is not provided for atom "<<AtomNumber2AtomName(Atom)<<"\n";
      std::cerr<<" please provide in the LocalParameter set\n";
      std::cerr<<" Function FreePrecess_Iph_7()\n"<<std::endl;
      Abort(1);
    };
  };
  double R1iph=fabs(ResolveParam(LocalParam[Atom],"R1iph"+marker));
  double R0iph=fabs(ResolveParam(LocalParam[Atom],"R0iph"+marker));
  double Omega=(gammaS/_gammaH)*dset._sfrq*2*DFH_PI*
    (ResolveParam(LocalParam[Atom],"Omega"+temperatureMarker)-dset._xcar);
  double DeltaO=(gammaS/_gammaH)*dset._sfrq*2*DFH_PI*
    DeltaOmega(dset,Atom);
  //
  gammaI=_gammaH;
  //
  // Approximate value for equlibrium value
  // site A
  //
  // Auto Relaxations
  G(1,1) = R0iph;
  G(2,2) = R0iph;
  G(3,3) = R1iph;
  G(3,0) = -(1-pb)*(R1iph)*gammaS/(gammaI*2.);
  //
  // Omega
  //G(1,2) =  (cs_offset); //turn off offset
  //G(2,1) = -(cs_offset);
  //
  // site B
  //
  G(4,4) = R0iph;
  G(5,5) = R0iph;
  G(6,6) = R1iph;
  G(6,0) = -pb*(R1iph)*gammaS/(gammaI*2.);
  //
  //G(4,5) =  (cs_offset+DeltaO); //turn off offset
  //G(5,4) = -(cs_offset+DeltaO);
  G(4,5) =  (DeltaO);
  G(5,4) = -(DeltaO);
  
  // kex additions
  for(unsigned int i=0;i<3;i++){
    G(i+1,i+1) = G(i+1,i+1)+kex*pb;
    G(i+4,i+1) = G(i+4,i+1)-kex*pb;
    //
    G(i+4,i+4) = G(i+4,i+4)+kex*(1-pb);
    G(i+1,i+4) = G(i+1,i+4)-kex*(1-pb);
  };
}


//For main spin lock.
//filter out eigen frequencies with imaginary components
void Catia::PulseX_R1rhoSpinLock(rrmat_t* P,ublas::matrix<double>& G,double w1,double time){


  P->resize(7,7);
  P->clear();
(*P)=G;
   // calculate the pulsing field.
  //double w1 = DFH_PI / (2*pw);
  //apply field to x
  //rrmat_t P = ublas::zero_matrix<double>(G.size1(), G.size2());
  for (unsigned int k = 0; k < 2; k++) {
    //inphase magnetisation only considered
    (*P)(2 + 3 * k, 3 + 3 * k) += w1;
    (*P)(3 + 3 * k, 2 + 3 * k) += -w1;
  }
  //(*P) += G;
  //P *= -time;
  //P = ublas::expm_pad(P);
  
  //now take eigenvalues, filter out those with non zero imaginary components,
  //and make propagator for specified time.

  //UNTESTED!!! //NEED TO CHECK THIS!!!
  expm_eig(*P,G,time);


  /*
  const unsigned int n = P->size1();
  cvec_t eigValVect(n);
  ccmat_t lEigVectMat(1, 1);
  ccmat_t rEigVectMat(n, n);
  ccmat_t freeEvolMatColMaj((*P));
	    
  lapack::geev('N', 'V',freeEvolMatColMaj, eigValVect, lEigVectMat, rEigVectMat, lapack::optimal_workspace());
		
  ccmat_t rEigVectMatInv(rEigVectMat);
  ublas::vector<int> ipiv(n); // pivot vector
  lapack::getrf(rEigVectMatInv, ipiv); // no lu_factor() alias for getrf() available
  lapack::getri(rEigVectMatInv, ipiv); // no lu_invert() alias for getrf() available
  
  ublas::zero_matrix<cmplx_t, ublas::column_major> zeroMat(n);
  ccmat_t tempDiagMat(n, n);


  // we should allocate space on G 
  tempDiagMat = zeroMat;
  for (unsigned int j = 0; j < n; ++j) {
    if(fabs(imag(eigValVect(j)))<1E-6)//include only if eigenvector is real
      tempDiagMat(j, j) = exp(-time * eigValVect(j));
  }
  (*P) = real(prod(rEigVectMat, ccmat_t(prod(tempDiagMat, rEigVectMatInv))));
  */
}


void Catia::PulseY_R1rho(ublas::matrix<double>&P,ublas::matrix<double>& G,double pw,double time){
  P.resize(7,7);
  P.clear();
   
  // calculate the pulsing field.
  double w1=DFH_PI / (2*pw); //get B1
  //apply field to x
  //rrmat_t P = ublas::zero_matrix<double>(G.size1(), G.size2());
  for (unsigned int k = 0; k < 2; k++) {
    //inphase
    P(1 + 3 * k, 3 + 3 * k) = -w1;
    P(3 + 3 * k, 1 + 3 * k) = w1;
  }
  P += G;
  //P *= -time;
  //P = ublas::expm_pad(P); //could probably be made real
  expm_eig(P,P,time);
}



//add the offset to the free precession matrix
void Catia::FreePrecess_R1rhoOffset(ublas::matrix<double>& P,ublas::matrix<double>& G,double cs_offset){
  P.resize(7,7);
  P.clear();

  P=G;

  // Omega
  P(1,2) += -(cs_offset);
  P(2,1) += +(cs_offset);
  //
  // site B
  P(4,5) += -(cs_offset);
  P(5,4) += +(cs_offset);
  

}
  

 void ShowCoh(ublas::vector<double> sigma)
 {
   std::cout << sigma[0] << std::endl;
   std::cout << sigma[1] << std::endl;
   std::cout << sigma[2] << std::endl;
   std::cout << sigma[3] << std::endl;
   std::cout << sigma[4] << std::endl;
   std::cout << sigma[5] << std::endl;
   std::cout << sigma[6] << std::endl;
   std::cout << " " << std::endl;
   

 }

//#make relaxation matrix for nicolai's formula
double Catia::GetExchangeInducedShift(Dataset& dset, int Atom){
  double gammaS=_gammaN;
  double DeltaO=(gammaS/_gammaH)*dset._sfrq*2*DFH_PI*DeltaOmega(dset,Atom);
  const double kex = Kex(dset);
  const double pb = Pb(dset);



  double DeltaR2=0.0;



  /*
  ublas::matrix<std::complex<double> >L; //declare new matricies for propagators
  L.resize(2,2);
  L.clear();
  double k_ab=pb*kex;
  double k_ba=(1-pb)*kex;

  //#add chemical shift evolution
  L(0,0)=0.0;
  L(1,1)=std::complex<double> (0,-DeltaO);//  #dOmega in rads-1
  
  //#add intrinsic relaxation
  L(0,0) -= 0.0;
  L(1,1) -= DeltaR2;
  
  //#add exchange
  L(0, 0) -= k_ab;
  L(1, 0) += k_ab;
  L(1, 1) -= k_ba;
  L(0, 1) += k_ba;

  const unsigned int n = L.size1();
  
  cvec_t eigValVect(n);
  ccmat_t lEigVectMat(1, 1);
  ccmat_t rEigVectMat(n, n);
  ccmat_t freeEvolMatColMaj(L);

  lapack::geev('N', 'V', freeEvolMatColMaj, eigValVect, lEigVectMat, rEigVectMat, lapack::optimal_workspace());
  //return -1*eig_values[0].imag; 
  return -1.*imag(eigValVect[0]);//#return shift in rad s-1
  */
  double pa=(1-pb);
  double keg=kex*(1-pb);
  double kge=kex*pb;
  //deltaR2=R2e-R2g;
  //#########################################################################
  //  #get the real and imaginary components of the exchange induced shift
    double g1=2*DeltaO*(DeltaR2+keg-kge);//                   #same as carver richards zeta
    double g2=pow(DeltaR2+keg-kge,2)+4*keg*kge-pow(DeltaO,2);//   #same as carver richards psi
    //double g3=cos(0.5*atan2(g1,g2))*pow(g1*g1+g2*g2,1/4.0);//   #trig faster than square roots
    double g4=sin(0.5*atan2(g1,g2))*pow(g1*g1+g2*g2,1./4.0);//   #trig faster than square roots
    //#########################################################################
    //#time independent factors
    //std::cout << -1*imag(eigValVect[0]) << " " << 0.5*(DeltaO-g4) << std::endl;
    return 0.5*(DeltaO-g4);

  
	
}

  //void Catia::CalcR1rho_baldwin(double gB1,double offset,double dw,double kex,double Pb,double R1,double R2,double R2E)

//Implementation of Baldwin R1rho formula 2012
void Catia::CalcR1rho_baldwin(Dataset& dset, int GlobalAtom) {

  // GlobalAtom is the AtomNumber of the Global set of Atoms;
  int LocalAtom = -1;
  for (unsigned int i = 0; i < dset._localToGlobalAtomIndex.size(); i++) {
    if (GlobalAtom == dset._localToGlobalAtomIndex[i]) {
      LocalAtom = i;
      break;
    }
  }
  if (LocalAtom == -1) {
    std::cerr << " Atom :" << AtomNumber2AtomName(GlobalAtom);
    std::cerr << " is not found in the Dataset named :" << dset._id;
    std::cerr << std::endl;
    Abort(1);
  }
  
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
  double w1 = atof(dset.vpar("w1")[1].c_str())*2*DFH_PI; //get pwx_cp from dataset file

  double gammaS=_gammaN;
  double DeltaO=(gammaS/_gammaH)*dset._sfrq*2*DFH_PI*DeltaOmega(dset,GlobalAtom); //get DeltaO in rad s-1
  double R1=fabs(ResolveParam(LocalParam[GlobalAtom],"R1iph"+marker));
  double R2=fabs(ResolveParam(LocalParam[GlobalAtom],"R0iph"+marker));
  double R2E=R2;

  const double kex = Kex(dset);
  const double pb = Pb(dset);

  double Pa=(1-pb);
  double R=R2-R1;
  double dr=R2E-R2;
  double f1p=+   Pa*pb*DeltaO*DeltaO ; // #sin
  double f2= dr*Pa+2*kex+w1*w1/kex;    // #sin 

  double shift=GetExchangeInducedShift(dset,GlobalAtom); //get exchange induced shift

  /*
  double time=10E-3; //fake time for r1rho decay rate constant
  ublas::vector<double> Sigma(7); // #get equilibrium magnetisiation
  CalcCoherence("Sz", Sigma, dset); //the dset contains the basis information
  double InitMag = MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), Detect("Sz", 0,dset, Sigma), Detect("Sz", 1, dset, Sigma));
  ublas::matrix<double> Gf; //space for free precession matrix
  ublas::matrix<double> P1,P2,Go; //declare new matricies for propagators
  rrmat_t P_SL;  //matrix for the spin lock.
  FreePrecess_R1rho(Gf, dset, GlobalAtom); //get onresonance free precession matrix (offset zero)
  */

  for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) //for each offset
    {
      //std::cout << " " << std::endl;
      double cs_offset;
      if(DeltaO>0) //adjust exchange induced shifted with sign of deltaO
	cs_offset=dset.ncyc[LocalAtom][ncycC]*2*DFH_PI - shift;
      else
	cs_offset=dset.ncyc[LocalAtom][ncycC]*2*DFH_PI + shift;	    
      

      double angle=GetTheta(w1,cs_offset)/DFH_PI*180.;  //get the flip angle in degrees

      //# making offset go in opposite direction for spectrometer
      double deltaA=-cs_offset;
      double deltaB=-cs_offset+DeltaO;
      double deltaBar=-cs_offset+pb*DeltaO;
      
      double theta=atan(w1/deltaBar);
      double  thetaflip=atan(w1/deltaA);
      
      // #  for approximation:
      // #  deltaBar=deltaA
      // #  theta=thetaflip
      double OmegaA=sqrt(w1*w1+deltaA*deltaA);
      double OmegaB=sqrt(w1*w1+deltaB*deltaB);
      double OmegaE=sqrt(w1*w1+deltaBar*deltaBar);
      
      //# Palmer original terms (not dependent on relaxation difference)
      double f2p=+   kex*kex+w1*w1 + pow(deltaB*deltaA/deltaBar,2); //#cos
      double dp =+   kex*kex + pow(OmegaA*OmegaB/OmegaE,2) ;      //     #Palmer denominator
      //# terms affected by relaxation difference
      double f1=  (OmegaA*OmegaA + kex*kex + Pa*kex*dr )*pb;//  #sin # BIG TERM (2)    
      double f3= 3*pb*kex+( 2*Pa*kex+w1*w1/kex + dr + dr*pow(pb*kex/OmegaA,2) )*1/(pow(sin(thetaflip),2));//  #fits neither sin nor cos
      
      //# Coefficients
      double R2ex=(f1p*kex+dr*f1) / (dp+dr*f3*pow(sin(theta),2));
      double CR1=(f2p+(f1p+dr*(f3-f2))*pow(tan(theta),2.)) / (dp+dr*f3*pow(sin(theta),2.));
      double CR2=(dp/pow(sin(theta),2)-f2p/pow(tan(theta),2)-f1p+dr*f2) / (dp+dr*f3*pow(sin(theta),2));
  
      //# composite expression
      double R1rhoNew=CR1*R1*pow(cos(theta),2.)+(CR2*R2+R2ex)*pow(sin(theta),2.);
      dset.R2_calc[LocalAtom][ncycC] = R1rhoNew;
      
      /*
      PulseY_R1rho(P1,Gf,1E-6,angle/90.*1E-6); //onresonance flip down propagator to take to angle
      PulseY_R1rho(P2,Gf,1E-6,-angle/90.*1E-6); //on resonance flip up propagator to take back to z
      FreePrecess_R1rhoOffset(Go,Gf,cs_offset); //adjust offset (Go)
      PulseX_R1rhoSpinLock(&P_SL,Go,w1,time); //get spin lock propagator
      Sigma.clear(); //Initialise magnetisation
      CalcCoherence("Sz", Sigma, dset); //the dset contains the basis information
      propagate(P1   ,Sigma); //apply flip down (+yangle)
      propagate(P_SL ,Sigma); //apply the spin lock 
      propagate(P2   ,Sigma); //apply flip up (-yangle)
      //ShowCoh(Sigma);
      double Intensity[2];
      Intensity[0] = Detect("Sz", 0, dset, Sigma); //get ground state signal intensity
      Intensity[1] = Detect("Sz", 1, dset, Sigma); //get ground state signal intensity
      //dset.R2_calc[LocalAtom][ncycC] = -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),Intensity[0], Intensity[1]) / InitMag) / time; 
      //dset.R2_calc[LocalAtom][ncycC] = -log(Intensity[0]/(1-pb))/time;
      //std::cout << log(Intensity[0]/(1-pb))/(time) << std::endl;
      //if(cs_offset==0)
      //  Abort(1);
      */
      //std::cout << R1rhoNew << " " << -log(Intensity[0]/(1-pb))/time << " " << -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),Intensity[0], Intensity[1]) / InitMag) / time << std::endl;
    }
  
}



void Catia::CalcR1rho(Dataset& dset, int GlobalAtom) {
	// GlobalAtom is the AtomNumber of the Global set of Atoms;
	int LocalAtom = -1;
	for (unsigned int i = 0; i < dset._localToGlobalAtomIndex.size(); i++) {
		if (GlobalAtom == dset._localToGlobalAtomIndex[i]) {
			LocalAtom = i;
			break;
		}
	}
	if (LocalAtom == -1) {
		std::cerr << " Atom :" << AtomNumber2AtomName(GlobalAtom);
		std::cerr << " is not found in the Dataset named :" << dset._id;
		std::cerr << std::endl;
		Abort(1);
	}
	///////////////////////////////////////////////////
	double w1 = atof(dset.vpar("w1")[1].c_str())*2*DFH_PI; //get pwx_cp from dataset file
	//double time_T2 = atof(dset.vpar("time_t2")[1].c_str()); //get pwx_cp from dataset file
	double gammaS=_gammaN;
	double DeltaO=(gammaS/_gammaH)*dset._sfrq*2*DFH_PI*DeltaOmega(dset,GlobalAtom); //get DeltaO in rad s-1
	const double pb = Pb(dset);

	double time=10E-3; //fake time for r1rho decay rate constant
	double shift=GetExchangeInducedShift(dset,GlobalAtom); //get exchange induced shift

	ublas::vector<double> Sigma(7); // #get equilibrium magnetisiation
	CalcCoherence("Sz", Sigma, dset); //the dset contains the basis information
	double InitMag = MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom), Detect("Sz", 0,dset, Sigma), Detect("Sz", 1, dset, Sigma));

	
	ublas::matrix<double> Gf; //space for free precession matrix
	ublas::matrix<double> P1,P2,Go; //declare new matricies for propagators
	rrmat_t P_SL;  //matrix for the spin lock.

	FreePrecess_R1rho(Gf, dset, GlobalAtom); //get onresonance free precession matrix (offset zero)
	for (unsigned int ncycC = 0; ncycC < dset.ncyc[LocalAtom].size(); ncycC++) //for each offset
	  {
	    //std::cout << " " << std::endl;
	    double cs_offset;
	    if(DeltaO>0) //adjust exchange induced shifted with sign of deltaO
	      cs_offset=dset.ncyc[LocalAtom][ncycC]*2*DFH_PI - shift;
	    else
	      cs_offset=dset.ncyc[LocalAtom][ncycC]*2*DFH_PI + shift;	    


	    double angle=GetTheta(w1,cs_offset)/DFH_PI*180.;  //get the flip angle in degrees
	    //std::cout << angle << " " << cs_offset << std::endl;

	    PulseY_R1rho(P1,Gf,1E-6,angle/90.*1E-6); //onresonance flip down propagator to take to angle
	    PulseY_R1rho(P2,Gf,1E-6,-angle/90.*1E-6); //on resonance flip up propagator to take back to z
	    FreePrecess_R1rhoOffset(Go,Gf,cs_offset); //adjust offset (Go)
	    PulseX_R1rhoSpinLock(&P_SL,Go,w1,time); //get spin lock propagator

	    Sigma.clear(); //Initialise magnetisation
	    CalcCoherence("Sz", Sigma, dset); //the dset contains the basis information

	    propagate(P1   ,Sigma); //apply flip down (+yangle)
	    propagate(P_SL ,Sigma); //apply the spin lock 
	    propagate(P2   ,Sigma); //apply flip up (-yangle)
	    //ShowCoh(Sigma);

	    double Intensity[2];
	    Intensity[0] = Detect("Sz", 0, dset, Sigma); //get ground state signal intensity
	    Intensity[1] = Detect("Sz", 1, dset, Sigma); //get ground state signal intensity
	    dset.R2_calc[LocalAtom][ncycC] = -log(MajorPeakIntensity(Pb(dset), Kex(dset), (dset._gamma / _gammaH) * dset._sfrq * 2 * DFH_PI * DeltaOmega(dset, GlobalAtom),Intensity[0], Intensity[1]) / InitMag) / time; 
	    //dset.R2_calc[LocalAtom][ncycC] = -log(Intensity[0]/(1-pb))/time;
	    //std::cout << log(Intensity[0]/(1-pb))/(time) << std::endl;
	    //if(cs_offset==0)
	    //  Abort(1);
	    
	  }
	
	//delete[] Gfree;
	
	//  exit(10);
}


					    
