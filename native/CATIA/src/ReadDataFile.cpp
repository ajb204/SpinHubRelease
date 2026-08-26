#include <standard.h>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>

using std::string;
using std::vector;
using std::valarray;
using std::map;
using std::ifstream;
using std::cerr;
using std::endl;
using std::istringstream;
using std::ws;

using boost::to_lower_copy;
using boost::trim_copy;
using boost::trim;
using boost::trim_if;
using boost::erase_all;
using boost::contains;
using boost::is_any_of;
using boost::replace_all;


void Dataset::ReadDataFile(const string& inputFileName, const string& resName) {

	vector<vector<double> > theData;

	vector<double> EmptyArray;

	for (unsigned int i = 0; i < _intensityFileFormat.size(); i++) {
		theData.push_back(EmptyArray);
	}

	ifstream ifs(inputFileName.c_str());

	if (!ifs) {
		cerr << " Could not open the inputfile " << inputFileName << "\n";
		cerr << " Function .ReadDataFile()\n";
		Abort(1);
	}

	char line[MAX_STRING_LENGTH];

	while (!ifs.eof()) {

		// First we read a line a strip according to #
		ifs.getline(line, sizeof(line));
		Tab2Space(line, sizeof(line));

		if (line[0] == '#') {
			continue;
		}

		string l(line);

		vector<string> its = split(l, "#");

		l = its[0];

		istringstream iss(l);

		if (its[0].length() < 1) {
			continue;
		}

		vector<double> RowVector;

		int RowCounter = 0;

		while (!iss.eof()) {
			RowVector.push_back(0.);
			iss >> RowVector[RowCounter++] >> ws;
		}


		for (unsigned int i = 0; i < _intensityFileFormat.size(); i++) {
			theData[i].push_back(RowVector[_intensityFileFormat[i]]);
		}

	}

	//Move From Vector to ValArray

	valarray<double> EmptyValArray;
	vector<valarray<double> > OutData;

	for (unsigned int i = 0; i < theData.size(); i++) {

		EmptyValArray.resize(theData[i].size());

		for (unsigned int j = 0; j < theData[i].size(); j++) {
			EmptyValArray[j] = theData[i][j];
		}

		OutData.push_back(EmptyValArray);

	}

	_localToGlobalAtomIndex.push_back(0);
	_atomNameToLocalAtomNumber[resName] = _localToGlobalAtomIndex.size() - 1;
	_localAtomNumberToAtomName[_localToGlobalAtomIndex.size() - 1] = resName;
	ncyc.push_back(OutData[0]);
	R2_exp.push_back(OutData[1]);
	R2_esd.push_back(OutData[2]);
	R2_calc.push_back(OutData[1]); // store the observed

	CorrectIntensityError(resName);

}

void Dataset::CorrectIntensityError(const string& resName) {
	unsigned int localAtomNumber = _atomNameToLocalAtomNumber[resName];
	for (unsigned int i = 0; i < R2_esd[localAtomNumber].size(); i++) {

		double minimumErrorPercentage = _minError[0] * 0.01 * R2_exp[localAtomNumber][i];
		if(R2_esd[localAtomNumber][i] < minimumErrorPercentage) {
			R2_esd[localAtomNumber][i] = minimumErrorPercentage;
		}

		double minimumErrorThresold = _minError[1];
		if(R2_esd[localAtomNumber][i] < minimumErrorThresold) {
			R2_esd[localAtomNumber][i] = minimumErrorThresold;
		}
	}
}

void Dataset::ReadId(const map<string, vector<string> >& parameters) {

	map<string, vector<string> >::const_iterator parameterIt = parameters.find("id");

	if (parameterIt != parameters.end()) {

		vector<string> values = parameterIt->second;
		string id(values[0]);
		replace_all(id, "\t", "_");
		replace_all(id, " ", "_");
		_id = id;

	} else {

		cerr << " The Dataset inputfile " << _inputFileName << " must include an ID statement\n";
		cerr << " usage: id = dataset name, e.g.\n";
		cerr << " id=N15 CW CPMG @ 800\n";
		cerr << " Function ReadDataset();\n";
		cerr << "\n" << endl;
		Abort(1);

	}

}

void Dataset::ReadNucleus(const map<string, vector<string> >& parameters) {

	map<string, vector<string> >::const_iterator parameterIt = parameters.find("nucleus");

	if (parameterIt != parameters.end()) {

		vector<string> values = parameterIt->second;
		string nucleus(to_lower_copy(values[0]));

		if (contains(nucleus, string("d"))) {
			_nucleus = "D";
		} else if (contains(nucleus, string("n"))) {
			_nucleus = "N";
		} else if (contains(nucleus, string("c"))) {
			_nucleus = "C";
		} else if (contains(nucleus, string("h"))) {
			if (nucleus == "H2") {
				_nucleus = "D";
			} else {
				_nucleus = "H";
			}
		} else if (contains(nucleus, string("f"))) {
			_nucleus = "F";
		} else {
			cerr << " Nucleus=" << _nucleus << endl;
			cerr << " Nucleus type " << values[0] << " is unknown" << endl;
			cerr << " Function ReadDataset();" << endl;
			Abort(1);
		}

	} else {

		cerr << " The Dataset inputfile " << _inputFileName << " must include a nucleus statement\n";
		cerr << " usage: nucleus={N,C,H,D,F}, e.g.\n";
		cerr << " nucleus=N\n";
		cerr << " Function ReadDataset();\n";
		cerr << "\n" << endl;
		Abort(1);

	}
}

void Dataset::ReadMinError(const map<string, vector<string> >& parameters) {

	map<string, vector<string> >::const_iterator parameterIt = parameters.find("minerror");

	if (parameterIt != parameters.end()) {
		//Check what unit we have:
		BOOST_FOREACH(string value, parameterIt->second)	{
			if (contains(value, "%")) {
				erase_all(value, "%");
				_minError[0] = atof(value.c_str());
			} else if (contains(value, "/s")) {
				erase_all(value, "/s");
				_minError[1] = atof(value.c_str());
			} else {
				cerr << " Only valid units for minerror are % and /s\n";
				cerr << " Function .ReadDataset()\n";
				cerr << endl;
				Abort(1);
			}
		}
	}

}

void Dataset::ReadIntensityFileFormat(const map<string, vector<string> >& parameters) {

	map<string, vector<string> >::const_iterator parameterIt = parameters.find("format");

	bool parameterFound = (parameterIt != parameters.end());
	bool parameterOfCorrectSize = ((parameterIt->second).size() == 3);

	if (parameterFound && parameterOfCorrectSize) {

		vector<string> values = parameterIt->second;

		_intensityFileFormat.clear();
		BOOST_FOREACH (string value, values) {
			_intensityFileFormat.push_back(atoi(value.c_str()));
		}

	} else if (!parameterOfCorrectSize) {

		cerr << " Format parameter is not allowed" << endl;
		cerr << " usage: format = (col[ncyc];col[R2];col[esd(R2)])" << endl;
		cerr << " Function ReadDataset();" << endl;
		Abort(1);

	} else {

		cerr << " The Dataset inputfile " << _inputFileName << " must include a format statement\n";
		cerr << " usage: format=(col[ncyc];col[R2];col[esd(R2)]), e.g.\n";
		cerr << " format=(0;1;2)\n";
		cerr << " Function ReadDataset();\n";
		cerr << "\n" << endl;
		Abort(1);

	}
}

void Dataset::ReadDataDirectory(const map<string, vector<string> >& parameters) {

	map<string, vector<string> >::const_iterator parameterIt = parameters.find("datadirectory");

	if (parameterIt != parameters.end()) {
		_dataDirectory = (parameterIt->second)[0];
	} else {
		_dataDirectory = "./";
	}
}

void Dataset::ReadData(const map<string, vector<string> >& parameters) {

	map<string, vector<string> >::const_iterator parameterIt = parameters.find("data");

	if (parameterIt != parameters.end()) {

		vector<string> values = parameterIt->second;
		vector<string> allPreviousAtomNames;

		BOOST_FOREACH(string value, values)	{

			vector<string> tokens;
			boost::split(tokens, value, is_any_of(";"));
			string atomName = trim_copy(tokens[0]);
			string fileName = trim_copy(tokens[1]);

			//
			//Check that this atom is not already there
			BOOST_FOREACH(string previousAtomName, allPreviousAtomNames) {
				if(atomName == previousAtomName) {
					cerr << " It seems that you have multiple definitions for the atomname:" << atomName << endl;
					cerr << " observed while reading the Dataset " << _id << endl;
					cerr << " Function .ReadDataset();" << endl;
					cerr << endl;
					Abort(1);
				}
			}
			allPreviousAtomNames.push_back(atomName);

			ReadDataFile(_dataDirectory + fileName, atomName);

		}

	} else {
		cerr << " The Dataset inputfile " << _inputFileName << " must include a data statement\n";
		cerr << " usage: data=([Name,FileName],....,[Name,FileName]), e.g.\n";
		cerr << " data=([14N,14N-HN.out.cpmg],[17N,17N-HN.out.cpmg])\n";
		cerr << " Function ReadDataset();\n";
		cerr << endl << endl;
		Abort(1);
	}
}

void Dataset::ReadInitialIntensities(const vector<string>& initialIntensityFormat) {

	vector<valarray<double> > initialIntensities; //Trosy,AntiTrosy intensity.
	initialIntensities.resize(_localToGlobalAtomIndex.size());
	BOOST_FOREACH (int atom, _localToGlobalAtomIndex)	{
		initialIntensities[atom].resize(2);
		initialIntensities[atom][0] = 0.;
		initialIntensities[atom][1] = 0.;
	}
	string initialIntensityFileName = initialIntensityFormat[0];
	const unsigned int nameCol = atoi(initialIntensityFormat[1].c_str());
	const unsigned int trosyCol = atoi(initialIntensityFormat[2].c_str());
	const unsigned int aTrosyCol = atoi(initialIntensityFormat[3].c_str());
	//
	ifstream initialIntensityFile(initialIntensityFileName.c_str());
	//

	if (initialIntensityFile.is_open()) {
		string aLine;
		while (getline(initialIntensityFile, aLine)) {
			replace_all(aLine, "\t", " ");
			vector<string> tokens;
			boost::split(tokens, aLine, is_any_of(" "));
			// do we have a comment line ?
			if (aLine.length() < 1 || tokens.size() < 2 || IsComment(aLine)) {
				continue;
			}
			if (!(nameCol < tokens.size() || trosyCol < tokens.size() || aTrosyCol < tokens.size())) {
				cerr << " You are trying to access column number: " << nameCol << "," << trosyCol << " and " << aTrosyCol << "\n";
				cerr << " in the file " << initialIntensityFileName << " but only " << tokens.size() << " columns are available\n";
				cerr << " Function .ReadDataset @ InitialIntensity();\n";
				cerr << endl;
				Abort(1);
			}
			string theAtomName = trim_copy(tokens[nameCol]);

			if (_atomNameToLocalAtomNumber.find(theAtomName) != _atomNameToLocalAtomNumber.end() ) {
				int localAtomIndex = _atomNameToLocalAtomNumber[theAtomName];
				valarray<double> theInitialIntensities;
				theInitialIntensities.resize(2);
				theInitialIntensities[0] = atof(tokens[trosyCol].c_str());
				theInitialIntensities[1] = atof(tokens[aTrosyCol].c_str());
				initialIntensities[localAtomIndex] = theInitialIntensities;
			}
		}
	} else {
		cerr << " Could not open the inputfile " << initialIntensityFileName << "\n";
		cerr << " Function .ReadDataset() @ InitialIntensity\n";
		cerr << endl;
		Abort(1);
	}

	//Check that we found intensities for all nuclei;
	BOOST_FOREACH (int localAtomIndex, _localToGlobalAtomIndex) {
		if (sqrt(pow(initialIntensities[localAtomIndex][0], 2.) + pow(initialIntensities[localAtomIndex][1], 2.)) < 1E-10) {
			cerr << " WARNING: No initial intensity found (or intensities < 1e-10) for"<<endl;
			cerr << " Atom: " << _localAtomNumberToAtomName[localAtomIndex] << " in the Dataset " << _id << endl;
			cerr << " Initial intensity set to (1.,-1), i.e., AntiPhase\n" << endl;
			initialIntensities[localAtomIndex][0] = 1.;
			initialIntensities[localAtomIndex][1] = -1.;
		}
	}
	_haveInitIntensity = true;
	initialIntensities = initialIntensities;

}



void Dataset::PutInProcpar(const string& name, const vector<string>& values) {

	vector<string> ida;
	ida.push_back(name);
	BOOST_FOREACH(string value, values)
	ida	.push_back(to_lower_copy(value));

	_procpar.push_back(ida);

}

