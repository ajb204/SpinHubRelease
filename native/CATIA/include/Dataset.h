/*
 *  Dataset.h
 *
 *
 *  Created by Guillaume Bouvignies on 04/06/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DATASET_H
#define DATASET_H

#include <standard.h>

class Dataset {

public:

	Dataset();
	~Dataset();

	//
	// Public methods
	//
	std::vector<std::string> vpar(const std::string) const; //view a parameter from the procpar list
	bool Bvpar(const std::string&); //True if parameter exists
	void PutInProcpar(const std::string&, const std::vector<std::string>&);
	void ReadId(const std::map<std::string, std::vector<std::string> >&);
	void ReadNucleus(const std::map<std::string, std::vector<std::string> >&);
	void ReadIntensityFileFormat(const std::map<std::string, std::vector<std::string> >&);
	void ReadMinError(const std::map<std::string, std::vector<std::string> >&);
	void ReadDataDirectory(const std::map<std::string, std::vector<std::string> >&);
	void ReadDataFile(const std::string&, const std::string&);
	void ReadInitialIntensities(const std::vector<std::string>&);
	void ReadData(const std::map<std::string, std::vector<std::string> >&);

	void SetDeltaOmegaTempDep(const std::string deltaOmegaTempDep) {
		_deltaOmegaTempDep = deltaOmegaTempDep;
	}
	const std::string& DeltaOmegaTempDep() const {
		return _deltaOmegaTempDep;
	}
	bool DeltaOmegaTempDepIs(const std::string deltaOmegaTempDep) const {
		return _deltaOmegaTempDep == deltaOmegaTempDep;
	}

	//
	// Data members
	//
	std::string _inputFileName;
	std::string _id; //
	double _sfrq; // Spectrometer 1H frequency
	double _temperature; // Temperature
	double _xcar; // carrier offset
	std::vector<std::vector<std::string> > _procpar; //
	std::string _nucleus; // N,C,H, etc ..
	double _gamma; // Gyro. mag. ratio.
	double _minError[2]; // minimum error (pct,/s)
	std::vector<int> _localToGlobalAtomIndex; // Atoms included in this dataset
	std::map<std::string, int> _atomNameToLocalAtomNumber; // Simular to Atoms, but maps atomname to number
	std::map<int, std::string> _localAtomNumberToAtomName; // just the inverse of _atomNameToLocalAtomNumber
	bool _initialized; // Whether first and initial R2 calculation has been done
	// This includes all the checking etc, later on in the real
	//   ... minimization this checking is skipped.
	bool _haveInitIntensity; // Do we have the initial intensities - otherwise start with Eq.
	std::vector<int> _intensityFileFormat;
	std::string _dataDirectory;

	//
	// Main data!
	//
	std::vector<std::valarray<double> > ncyc;
	std::vector<std::valarray<double> > R2_exp;
	std::vector<std::valarray<double> > R2_esd;
	std::vector<std::valarray<double> > R2_calc;

	//
	// Initial intensities.
	//
	std::vector<std::valarray<double> > initialIntensities;

private:

	void CorrectIntensityError(const std::string&);
	std::string _deltaOmegaTempDep;
	std::string _sfrqStr;
	std::string _temperatureStr;

};

#endif
