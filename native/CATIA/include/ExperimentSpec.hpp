/*
 * ExperimentSpec.hpp
 *
 *  Created on: May 18, 2010
 *      Author: guillaume
 */

#ifndef EXPERIMENTSPEC_HPP_
#define EXPERIMENTSPEC_HPP_

class ExperimentSpec {

public:
	ExperimentSpec();
	~ExperimentSpec();

private:
	// Spectrometer 1H frequency
	double _sfrq;

	// Studied nucleus - N, C, H, etc.
	std::string _nucleus;

	// Gyromagnetic ratio of the studied nucleus
	double _gamma;

	// Temperature of the sample
	double _temperature;

	// Carrier offset
	double _xcar;

	// Additional Parameters
	std::vector<std::vector<std::string> > _procpar;

};

#endif /* EXPERIMENTSPEC_HPP_ */
