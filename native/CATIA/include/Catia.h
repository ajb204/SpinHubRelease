/*
 CPMG and Trosy Intelligent Analysis (CATIA)

 D. Flemming Hansen
 flemming@pound.med.utoronto.ca,
 Aug. 2007 - Aug. 2008

 27 May 2009: Modified by Guillaume Bouvignies
 - The selection of the atoms now follows the perl regular expression conventions.
 ex: "+." -> "*"; "[0-9]+N" -> "*N" where "*" is any number
 - The function ReadDataSet can now be called several times, in particular when new atoms has been added to Atoms.

 */

#ifndef CATIA_H
#define CATIA_H

#include <standard.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <nr3.h>

#ifdef DFH_NRUTIL
#define NRANSI

//  extern "C" {
//    #include <nr.h>
//    #include <nrutil.h>
//  };
#endif

#ifdef DFH_MPI
#include <mpi.h>
#endif

#ifndef DFH_PI
#define DFH_PI 3.1415926535
#endif

#ifndef CATIA_VERSION
#define CATIA_VERSION 0.60
#endif



namespace ublas = boost::numeric::ublas;
typedef double real_t;
typedef ublas::matrix<real_t, ublas::row_major> rrmat_t;

// Forward declarations
class Dataset;

class Catia {

public:
	// Constructors
	Catia();

	// Detstructor
	~Catia();

	//Main Calculations/minimization
	//
	double Enorm(Dataset&, int); // Ready for the NR routines.

	double Minimize(std::vector<int>, std::vector<int>, std::string, bool&);
	double Minimize(std::vector<int>, std::vector<std::string>, std::string, bool&);
	double Minimize(std::string, std::vector<int>, std::string, bool&);
	double Minimize(std::vector<int>, std::string, std::string, bool&);
	double Minimize(std::string, std::string, std::string, bool&);
	ublas::matrix<double> LastFitCovar;
	void CalcR2(Dataset&, int);
	//
	// Book-keeping functions and internal IO
	int AtomName2AtomNumber(std::string);
	std::string AtomNumber2AtomName(int) const;
	//
	// Parameter handling
	void FreeGlobalParam(std::string, int); //(name,free)
	void FreeLocalParam(std::vector<int>, std::string, int); //(Atoms,ParamName,free);
	void FreeLocalParam(unsigned int, std::string, int); //(AtomNumber,ParamName,free);
	void FreeLocalParam(std::string, std::string, int); //(AtomName,ParamName,free);
	void FreeLocalParam(std::vector<std::string>, std::string, int); //(AtomNames,ParamName,free);
	//
	// Set/Get value of a specific parameter
	void SetGlobalParam(std::string, double); //(ParamName,NewValue)
	double GetGlobalParam(std::string); //(ParamName)
	void SetLocalParam(std::string, std::string, double); //(AtomName,ParamName,NewValue);
	double GetLocalParam(std::string, std::string); //(AtomName,ParamName)
	void SetInternalFixThres(double);
	void SetRateType(std::string);
	void SetMultipleTemperatures(bool multipleTemperatures) {
		_multipleTemperatures = multipleTemperatures;
	}
	void SetDeltaOmegaType(const std::string);
	//
	// IO functions
	int ReadDataset(std::string); // Return the number assigned
	void ReadParam_Global(std::string);
	void ReadParam_Local(std::string);
	void ReadParam(std::string, std::string, int, int); //(paramname,filename,Col[Name],Col[Val]  (Columns start from zero)
	void SortData();
	//
	std::string PrintParam(std::string); // sofar for "all" / "*"
	std::string PrintParam(bool, int, std::string); // global,GlobalAtomNumber,Name,
	std::string PrintArrheniusParam(void) const;
	//
	void PrintData(int); //Atomnumber
	void PrintData(std::string, std::ostream&); // AtomName,outstream
	ublas::matrix<double> CorrelationM(std::vector<std::string>);
	//
	// Maps
	std::map<std::string, int> LastFitParamName;
	std::map<std::string, double> GlobalParam; //includes kex, pb etc.
	std::map<std::string, int> GlobalParamF; //0: free, 1: fixed
	std::map<std::string, double> GlobalParamE; //uncertainty (-1: not fitted, 0: fixed )
	//
	std::vector<std::map<std::string, double> > LocalParam; //Include DeltaOmega,R2inf,..
	std::vector<std::map<std::string, int> > LocalParamF; //0: free, 1: fixed
	std::vector<std::map<std::string, double> > LocalParamE; //uncertainty
	std::vector<std::map<std::string, std::string> > LocalNotes; // Can store information that will be printed out
	// in the final output file.

	std::map<int, std::vector<int> > X2dset_atom; //This is a kind-of wrapper map for helping the numerical recepies functions out.
	// They only take _one_ x coordinates, which we here map into an array,
	// Thus, x={0,1,2,3,4....}, for each datapoint, which then maps
	// to the different parameters (Dataset,atomnumber,ncyc)
	//
	// Variables
	std::vector<Dataset> Datasets;
	std::map<std::string, int> Atoms;
	std::vector<int> Atoms2Fit;
	int CalcDeriv; // 0: No derivatives;
	//                1: Fast, (f(a1,..,ai*(1.+1E-4),..,an)-f(a1,...,ai,...,an))/1E-4*ai
	//                2: Slow, use the NR subroutine (Not implemented yet)
	//
	// we define the Levernberg-Marquart object class as derived from catia.
#include <fitmrq.h>

	//

private:

	void CalcR2drv(VecDoub_I, VecDoub&, MatDoub&, const std::vector<double*>, const std::vector<int*>, const std::vector<int>, const std::vector<std::string>);

	//
	// NMR sequences
	void CalcR2CW_CPMG(Dataset&, int);
	// No decoupling!!! .. only cpmg part
	// Sy -> CPMG(Y) - 180(x,-x) - CPMG(Y) - 90x -> detect Sz
	void CalcR2Fast_CPMG(Dataset&, int); // Averaged Liouvillian
	void CalcR2Trosy_CPMG(Dataset&, int); // PNAS paper (Also works for alpha C)
	void CalcR2AntiTrosy_CPMG(Dataset&, int); // PNAS paper (Also works for alpha C)
	void CalcR2TrosySel_CPMG(Dataset&, int); // As above, but selects Trosy before CPMG (initial condition)
	void CalcR2AntiTrosySel_CPMG(Dataset&, int); // As above, but selects AntiTrosy before CPMG
	void CalcR2Trosy_CPMG_vo(Dataset&, int); // With the XY sequence (primarily used for amide protons).
	void CalcR2PE_CPMG(Dataset&, int); // (CPMG)y - tau - 180 - tau - (CPMG)y
	void CalcR2AP_CPMG(Dataset&, int); // (CPMG)y - 180x - (CPMG)y
	void CalcR2CW_3st_CPMG(Dataset&, int); // (CPMG)y - 180x - (CPMG)y
	void CalcR2CW_CPMG_Tr(Dataset&, int);
	// Set init conditions  - optionally from input.
	// CPMG(y) - 180x,-x CPMG(y)
	// Select trosy component.
	void CalcR2CW_CPMG_ATr(Dataset&, int); // Same as CW_CPMG_Tr - just select anti-trosy in the end!
	void CalcR2COsc_CPMG(Dataset&, int);
	// the carbonyl cpmg, with the sidechain (CO;ASP;ASN) inversion filter.
	// CPMG(x) - 90(y) - delay(taucc) - 180(y) - delay(taucc) - 90(y) - CPMG(x) -
	void CalcR2CW_CPMG_CH3(Dataset&, int); // The methyl group cpmg, with selection of the four lines.
	void CalcR2PE_CPMG_CH3(Dataset&, int);
	// The methyl group cpmg, with selection of the four lines.
	// .. with palmer element in the middle.
	void CalcR1_Sz(Dataset&, int); // The T1 experiment: Create +/- Nz and apply proton pulses every 5 ms.
	void CalcR1_SzIz(Dataset&, int); // The NzHz experiment: 2NzHz -> T/4 - 180(H) - T/4 - 180(N) - T/4 - 180(H) - T/4 -> detect
	void CalcR1rho(Dataset&, int); // The NzHz experiment: 2NzHz -> T/4 - 180(H) - T/4 - 180(N) - T/4 - 180(H) - T/4 -> detect
	void CalcR1rho_baldwin(Dataset&, int); // The NzHz experiment: 2NzHz -> T/4 - 180(H) - T/4 - 180(N) - T/4 - 180(H) - T/4 -> detect

	double GetTheta(double w1,double offset); //#return the angle in radians
	void FreePrecess_R1rho(ublas::matrix<double>& G,Dataset& dset, int Atom);
	void PulseX_R1rhoSpinLock(rrmat_t* P,ublas::matrix<double>& G,double w1,double time);
	void PulseY_R1rho(ublas::matrix<double>&P,ublas::matrix<double>& G,double pw,double time);
	void FreePrecess_R1rhoOffset(ublas::matrix<double>& P,ublas::matrix<double>& G,double cs_offset);
	double GetExchangeInducedShift(Dataset& dset, int Atom);
	//
	//
	// Matrix diagonalization
	void CalcMatrix(ublas::matrix<double>&, Dataset&, int, std::string, double);//Diagonalizes the evoution matrices, including pulses etc.
	void CalcMatrix(ublas::matrix<double>*, Dataset&, int, std::string, std::vector<double>);
	void CalcMatrix(ublas::matrix<std::complex<double> >*, Dataset&, int, std::string, std::vector<double>);
	//
	//
	//Calculates the FreePrecessing Liovillian (brance out to different bases)
	void FreePrecess(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess(ublas::matrix<std::complex<double> >&, Dataset&, int);
	void FreePrecess_Fast_2(ublas::matrix<std::complex<double> >&, Dataset&, int);
	void FreePrecess_Fast_3(ublas::matrix<std::complex<double> >&, Dataset&, int);
	void FreePrecess_TrATr_13(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_IphAph_13(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_IphAph_13_deltaR2(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_IphAph_13_deltaR2_simple(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_Iph_7(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_Aph_7(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_Tr_7(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_C_CH3_25(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_N_NH_13(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_N_7(ublas::matrix<double>&, Dataset&, int);
	void FreePrecess_3st_Iph_10(ublas::matrix<double>&, Dataset&, int);
	//
	//Detecting the signal, etc.
	void CalcCoherence(std::string, ublas::vector<double>&, Dataset); //coherence,Vector,dset
	void CalcCoherence_3st(std::string, ublas::vector<double>& Eq, Dataset); //coherence,Vector,dset
	double MajorPeakIntensity(double, double, double, double, double);
	double MajorPeakIntensity3States(double, double, double, double, double, double, double, double, double, double);
	double Detect(std::string, int, Dataset, ublas::vector<double>&);
	//
	//Calculate kex
	double Kex(const Dataset&) const;
	double Kex_ab_3st(const Dataset&) const;
	double Kex_ac_3st(const Dataset&) const;
	double Kex_bc_3st(const Dataset&) const;
	double KexStandard(const Dataset&) const;
	double KexStandard_3st_ab(const Dataset&) const;
	double KexStandard_3st_ac(const Dataset&) const;
	double KexStandard_3st_bc(const Dataset&) const;
	double KexArrhenius(const Dataset&) const;
	double KexArrhenius_3st_ab(const Dataset&) const;
	double KexArrhenius_3st_ac(const Dataset&) const;
	double KexArrhenius_3st_bc(const Dataset&) const;
	double KexKabKba() const;
	//
	//Calculate pb
	double Pb(const Dataset&) const;
	double Pb_3st(const Dataset&) const;
	double Pc_3st(const Dataset&) const;
	double PbStandard(const Dataset&) const;
	double PcStandard(const Dataset&) const;
	double PbArrhenius(const Dataset&) const;
	double PbArrhenius_3st(const Dataset&) const;
	double PcArrhenius_3st(const Dataset&) const;
	double PbKabKba() const;
	//
	//Calculate deltaOmega
	double DeltaOmega(const Dataset&, const int) const;
	double DeltaOmega_ab(const Dataset&, const int) const;
	double DeltaOmega_ac(const Dataset&, const int) const;
	double DeltaOmegaStandard(const Dataset&, const int) const;
	double DeltaOmegaStandard_ab(const Dataset&, const int) const;
	double DeltaOmegaStandard_ac(const Dataset&, const int) const;
	double DeltaOmegaLinear(const Dataset&, const int) const;
	double DeltaOmegaLinear_ab(const Dataset&, const int) const;
	double DeltaOmegaLinear_ac(const Dataset&, const int) const;
	double DeltaOmegaHarmonic(const Dataset&, const int) const;
	//
	void AddAtom(std::string);
	//
	std::string PrintArrheniusParam2st(void) const;
	std::string PrintArrheniusParam3st(void) const;
	//
	double ResolveParam(const std::map<std::string, double>&, const std::string&) const;
	bool BResolveParam(const std::map<std::string, double>&, const std::string&) const; //Same as ResolveParam, but returns a bool to indicate if parameter exists
	//
	void HasMultipleTemperatures(void);
	//
	const double _gammaH;
	const double _gammaN;
	const double _gammaC;
	const double _gammaD;
	const double _gammaF;
	//
	double _fixParamLimit;
	bool _stopFitting; //Catch a CTRL-C and stop fitting (Not implemented yet)
	//
	// Temperature Dependence Activated?
	bool _multipleTemperatures;
	std::string _rateType;
	//

};

#endif
