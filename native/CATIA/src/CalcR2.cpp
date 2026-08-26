#include <Catia.h>
#include <Dataset.h>
#include <Abort.h>

void Catia::CalcR2(Dataset& dset, int Atom) {
  /*
	 Calculates the R2 for all the experimental points
	 in the Dataset #Dataset and atom Atom
	 */
	std::string seqfil = dset.vpar("seqfil")[1];
	//Now call the sequence
	if (seqfil == "fast_cpmg") {
		CalcR2Fast_CPMG(dset, Atom);
	} else if (seqfil == "cw_cpmg") {
		CalcR2CW_CPMG(dset, Atom);
	} else if (seqfil == "trosy_cpmg") {
		CalcR2Trosy_CPMG(dset, Atom);
	} else if (seqfil == "trosy_cpmg_select") {
		CalcR2TrosySel_CPMG(dset, Atom);
	} else if (seqfil == "antitrosy_cpmg") {
		CalcR2AntiTrosy_CPMG(dset, Atom);
	} else if (seqfil == "antitrosy_cpmg_select") {
		CalcR2AntiTrosySel_CPMG(dset, Atom);
	} else if (seqfil == "trosy_cpmg_vo") {
		CalcR2Trosy_CPMG_vo(dset, Atom);
	} else if (seqfil == "pe_cpmg") {
		CalcR2PE_CPMG(dset, Atom);
	} else if (seqfil == "ap_cpmg") {
		CalcR2AP_CPMG(dset, Atom);
	} else if (seqfil == "cw_cpmg_tr") {
		CalcR2CW_CPMG_Tr(dset, Atom);
	} else if (seqfil == "cw_cpmg_atr") {
		CalcR2CW_CPMG_ATr(dset, Atom);
	} else if (seqfil == "cosc_cpmg") {
		CalcR2COsc_CPMG(dset, Atom);
	} else if (seqfil == "cw_cpmg_ch3") {
		CalcR2CW_CPMG_CH3(dset, Atom);
	} else if (seqfil == "pe_cpmg_ch3") {
		CalcR2PE_CPMG_CH3(dset, Atom);
	} else if (seqfil == "cw_3st_cpmg") {
		CalcR2CW_3st_CPMG(dset, Atom);
	} else if (seqfil == "r1rho") {
		CalcR1rho(dset, Atom);
	} else if (seqfil == "r1rho_baldwin") {
		CalcR1rho_baldwin(dset, Atom);
	} else {
		std::cerr << " No sequence compiled for seqfil=" << seqfil << "\n";
		std::cerr << " Function: CalcR2();\n";

		Abort(1);
	}
}
