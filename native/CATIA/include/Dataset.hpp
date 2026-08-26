/*
 * Dataset.hpp
 *
 *  Created on: May 18, 2010
 *      Author: guillaume
 */

#ifndef DATASET_HPP_
#define DATASET_HPP_


class Dataset {

public:
	Dataset();
	Dataset(const std::string inputFileName);
	~Dataset();

	void ReadFile(const std::string inputFileName);

	// Setter and Getter for _deltaOmegaTempDep
	void SetDeltaOmegaTempDep(const std::string deltaOmegaTempDep) {
		_deltaOmegaTempDep = deltaOmegaTempDep;
	}
	const std::string& DeltaOmegaTempDep() const {
		return _deltaOmegaTempDep;
	}
	bool DeltaOmegaTempDepIs(const std::string deltaOmegaTempDep) const {
		return (_deltaOmegaTempDep == deltaOmegaTempDep);
	}


private:
	// Name of the input dataset file
	std::string _inputFileName;

	// ID of the dataset
	std::string _id;

	// Temperature dependence of delta_omega values
	std::string _deltaOmegaTempDep;

	// Minimum error (% of the R2,eff , /s)
	double _minError[2];

	// Map of local atom indexes to the global atom indexes
	std::vector<int> _localToGlobalAtomIndex;

	// Map of atom names to global atom indexes
	std::map<std::string, int> _atomNameToLocalAtomNumber;

	// Map of global atom indexes to atom names
	std::map<int, std::string> _localAtomNumberToAtomName;

	// Whether first and initial R2 calculation has been done
	// This includes all the checking etc, later on in the real
	//   ... minimization this checking is skipped.
	bool _initialized;

	// Do we have the initial intensities - otherwise start with Eq.
	bool _haveInitIntensity;

	// Format of intensity data file
	std::vector<int> _intensityFileFormat;

	// Directory of the intensity data files
	std::string _dataDirectory;


};


#endif /* DATASET_HPP_ */
