#include <standard.h>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>

typedef std::pair<std::string, std::string> pairSS;
typedef std::pair<std::string, int> pairSI;
typedef std::pair<std::string, std::vector<std::string> > pairSVs;

using std::string;
using std::vector;
using std::map;
using std::ifstream;
using std::cerr;
using std::endl;

using boost::to_lower_copy;
using boost::trim_copy;
using boost::trim_if;
using boost::trim;
using boost::erase_all;
using boost::contains;
using boost::is_any_of;
using boost::replace_all;

map<string, vector<string> > ParseDataSetFile(const string&);
void ReadInitialIntensities(const vector<string>&);
void RemoveComments(string&);
string ParseName(string);
vector<string> ParseValues(string);
void CheckInitialIntensityFormat(const vector<string>&);

int Catia::ReadDataset(string inputFileName) {

	vector<string> initialIntensityFormat;

	Dataset out;

	out._inputFileName = inputFileName;
	map<string, vector<string> > parameters(ParseDataSetFile(out._inputFileName));

	//Now read the parameters
	//
	out.ReadId(parameters);
	out.ReadNucleus(parameters);
	out.ReadIntensityFileFormat(parameters);
	out.ReadMinError(parameters);
	out.ReadDataDirectory(parameters);
	out.ReadData(parameters);

	BOOST_FOREACH(pairSI atomNameToLocalAtomNumber, out._atomNameToLocalAtomNumber) {
		string atomName(atomNameToLocalAtomNumber.first);
		//See if we have already data for this residue?
		if (AtomName2AtomNumber(atomName) == -1) {
			AddAtom(atomName);
		} else {
			//Check that we have defined this atom with the same nucleus type in the other Datasets.
			//type.
			BOOST_FOREACH (Dataset dataset, Datasets) {
				BOOST_FOREACH (int atom, dataset._localToGlobalAtomIndex) {
					if (atom == AtomName2AtomNumber(atomName)) {
						if (!(dataset._nucleus == out._nucleus)) {
							cerr << " Mismatch nucleus type for Atom with name: " << atomName << endl;
							cerr << " While reading Dataset \'" << out._id << "\', which has defined nucleus type of " << out._nucleus << endl;
							cerr << " Previously, in Dataset \'" << dataset._id << "\' this atom was defined with type: " << dataset._nucleus << endl;
							cerr << " Rading Dataset from file:" << out._inputFileName << " aborted\n" << endl;
							Abort(1);
						}
					}
				}
			}
		}
		out._localToGlobalAtomIndex[atomNameToLocalAtomNumber.second] = AtomName2AtomNumber(atomName);
	}

	BOOST_FOREACH(pairSVs element, parameters) {
		string name(element.first);
		vector<string> values(element.second);

		if (name == "initialintensity") {

			initialIntensityFormat = values;
			CheckInitialIntensityFormat(initialIntensityFormat);
			out.ReadInitialIntensities(initialIntensityFormat);

		} else if (name == "sfrq") {

			out._sfrq = atof(values[0].c_str());

		} else if (name == "xcar") {

			out._xcar = atof(values[0].c_str());

		} else if (name == "temperature") {

			out._temperature = atof(values[0].c_str());

		} else if (name == "deltaomegatype") {

			out.SetDeltaOmegaTempDep(values[0]);

		} else if (!contains("id;format;data", name)) { //Scoop the parameters into the pulseseq_param{ //Scoop the things into the pulseseq_param
			out.PutInProcpar(name,values);
		}
	}
	//
	//
	//Now do some checking before submitting the Dataset
	if (!(out._localToGlobalAtomIndex.size() == out.ncyc.size() && out._localToGlobalAtomIndex.size() == out.R2_exp.size() && out._localToGlobalAtomIndex.size() == out.R2_esd.size())) {
		cerr << " There is a mismatch in the size of the arrays in the Dataset\n";
		cerr << out._id << "\n";
		cerr << " .Atom.size()  =" << out._localToGlobalAtomIndex.size() << "\n";
		cerr << " .ncyc.size()  =" << out.ncyc.size() << "\n";
		cerr << " .R2_exp.size()=" << out.R2_exp.size() << "\n";
		cerr << " .R2_esd.size()=" << out.R2_esd.size() << "\n";
		Abort(1);
	}

	//
	// Read the initial intensities
	//

	Datasets.push_back(out);

	HasMultipleTemperatures();

	return Datasets.size();

}


void Catia::HasMultipleTemperatures(void) {

	double noTemperature = 1000.0;
	double oldTemperature = noTemperature;

	_multipleTemperatures = false;

	BOOST_FOREACH(Dataset dset, Datasets) {
		if (oldTemperature == noTemperature) {
			oldTemperature = dset._temperature;
		} else if (oldTemperature != dset._temperature) {
			_multipleTemperatures = true;
		}
	}
	std::cout << _multipleTemperatures << std::endl;
}


map<string, vector<string> > ParseDataSetFile(const string& fileName) {

	map<string, vector<string> > result;

	ifstream infile(fileName.c_str());

	if (infile.is_open()) {
		string line;

		while (getline(infile, line)) {

            trim(line);

            if (IsComment(line) or line.empty()) {
				continue;
			}

			RemoveComments(line);

			if (!contains(line, string("="))) {
				continue;
			}

			if (contains(line, string("("))) { // We have an arried parameter
				while (!contains(line, string(")"))) {
					string subLine;
					getline(infile, subLine);
					if (IsComment(subLine)) {
						continue;
					}
					RemoveComments(subLine);
					line += subLine;
				}
			}

			result[ParseName(line)] = ParseValues(line);

		}

	} else {
	  //cerr << " Could not open the inputfile " << infile << "\n";
		cerr << " Could not open the inputfile "  << "\n";
		cerr << " Function .ReadDataset()\n";
		cerr << endl;
		Abort(1);
	}

	return result;
}


void RemoveComments(string& line) {
	vector<string> tokens;
	boost::split(tokens, line, is_any_of("#"));
	line = tokens[0];
	trim(line);
}

string ParseName(string fullLine) {
	vector<string> tokens;
	boost::split(tokens, fullLine, is_any_of("="));
	return to_lower_copy(trim_copy(tokens[0]));
}

vector<string> ParseValues(string fullLine) {
	vector<string> tokens;

	boost::split(tokens, fullLine, is_any_of("="));

	vector<string> values;

	if (contains(tokens[1], "[")) {
		trim_if(tokens[1], is_any_of(" ()[];\t"));
		replace_all(tokens[1], "];[", "|");
		boost::split(values, tokens[1], is_any_of("|"));
	} else {
		trim_if(tokens[1], is_any_of(" ()[];\t"));
		boost::split(values, tokens[1], is_any_of(";"));
	}

	return values;
}

void CheckInitialIntensityFormat(const vector<string>& values) {

	if (values.size() != 4) {
		cerr << " Wrong format of InitialIntensity line\n";
		cerr << " USAGE: InitialIntensity=(FileName;col[Name];col[init_Trosy];col[init_AntiTrosy])  \n";
		cerr << " file " << __FILE__ << " at line " << __LINE__ << "\n";
		cerr << " Function .ReadDataset();" << endl;
		cerr << " PROGRAM ABORTED" << endl;
		Abort(1);
	}

}
