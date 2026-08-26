/*
 * DatasetInventory.hpp
 *
 *  Created on: May 18, 2010
 *      Author: guillaume
 */

#ifndef DATASETINVENTORY_HPP_
#define DATASETINVENTORY_HPP_

// Forward declarations
class Dataset;

class DatasetInventory {

public:
public:

	DatasetInventory();  // constructor - no return type
	~DatasetInventory(); // destructor - no return type

	void AddDataset(const std::string inputFileName);

private:
	void CheckConsistency();

private:
	std::vector<Dataset> Datasets;

};


#endif /* DATASETINVENTORY_HPP_ */
