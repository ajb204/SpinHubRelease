#include <standard.h>
#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>
#include <Abort.h>

int returnval = 0;

namespace ublas = boost::numeric::ublas;


std::string Convert(std::string is) {
	char line[MAX_STRING_LENGTH];
	char oline[MAX_STRING_LENGTH];
	ClearBuf(line, sizeof(line));
	sprintf(line, "%s", is.c_str());
	int cc = 0;
	int tabs = 0;
	for (unsigned int i = 0; i < (is.size() + tabs); i++) {
		// enter
		if (line[i + tabs] == '\\' && line[1 + i + tabs] == 'n') {
			oline[cc] = '\n';
			cc++;
			tabs++;
		} else if (line[i + tabs] == '\\' && line[1 + i + tabs] == 'b') {
			oline[cc] = '\b';
			cc++;
			tabs++;
		} else if (line[i + tabs] == '\\' && line[1 + i + tabs] == 't') {
			oline[cc] = '\t';
			cc++;
			tabs++;
		} else if (line[i + tabs] == '\\' && line[1 + i + tabs] == 'r') {
			oline[cc] = '\r';
			cc++;
			tabs++;
		} else {
			oline[cc] = line[i + tabs];
			cc++;
		}
	}
	oline[cc] = 0;
	return std::string(oline);
}

std::string lowercase(std::string l) {
	std::transform(l.begin(), l.end(), // source
			l.begin(), // destination
			tolower); // operation
	return l;
}

void ErrorLine(std::string line) {
	std::cerr << " Error while reading the commandline:\n --->";
	std::cerr << line << "<-----\n";
	std::cerr << " Syntax error\n" << std::endl;
	returnval = 1;
}

void WriteHelp(std::string func) {
	if (func.substr(0, 3) == "all") {
		std::cout << " Usage: help([Function])\n" << std::endl;
		std::cout << " Functions available:\n";
		std::cout
				<< "  - ChiSq              : Calculate ChiSq for Datasets/atoms\n";
		std::cout
				<< "  - CorrMatrix         : Prints the correlation matrix of global parameters\n";
		std::cout
				<< "  - exit(N)            : Abort the program, with return value N \n";
		std::cout << "  - echo               : Write text to screen \n";
		std::cout
				<< "  - FreeLocalParam     : Fix or set free local parameters\n";
		std::cout
				<< "  - FreeGlobalParam    : Fix or set free global parameters\n";
		std::cout << "  - Minimize           : Perform minimization\n";
		std::cout
				<< "  - PrintData          : Print experimental data v.s. fitted data\n";
		std::cout << "  - PrintParam         : Print fitted parameters\n";
		std::cout << "  - ReadDataset        : Read another Dataset\n";
		std::cout
				<< "  - ReadParam          : Read parameters from an input file\n";
		std::cout << "  - ReadParam_Global   : Read global parameters\n";
		std::cout << "  - ReadParam_Local    : Read local paramters\n";
		std::cerr
				<< "  - SetLocalParam      : Set the value of a local parameter\n";
		std::cout
				<< "  - SetInternalFixThres: Set the threshold for automatic fixing\n";
		std::cout
				<< "  - quit(N)            : Abort the program with return value N \n";
		std::cout << std::endl;
		return;
	}
	//
	std::cout << "\n\n Usage: ";
	if (func == "ChiSq") {
		std::cout << " ChiSq([DataSet];[AtomName])\n";
		std::cout << "   DataSet : ID of the Dataset to use or 'all'\n";
		std::cout << "   AtomName: AtomName or 'all'\n";
		std::cout << " Example: ChiSq(all;all)\n";
		std::cout << "          ChiSq(600_Trosy;all)\n";
		std::cout << "          ChiSq(600_Trosy;15N)\n";
	} else if (func == "CorrMatrix") {
		std::cout << " CorrMatrix()\n";
	} else if (func == "SetInternalFixThres") {
		std::cout << " SetInternalFixThres(Thres)\n";
		std::cout << "   Thres : Threshold\n\n";
		std::cout
				<< " This function sets the threshold for internal fixing of parameters,\n";
		std::cout
				<< " thus if RMSD( (dy/dp)/y )_{y \\in datapoints} < Thres then the parameter p \n";
		std::cout << " will be fixed for future iterations.\n\n";
		std::cout << " Example: SetInternalFixThres(1E-9)\n";
	} else if (func == "echo") {
		std::cout << " echo(text)\n";
		std::cout << " Example: echo(03N)\n";
	} else if (func == "FreeLocalParam") {
		std::cout << " FreeLocalParam([AtomName];[Param];[false/true])\n";
		std::cout << " Example: FreeLocalParam(all;Omega;false)\n";
		std::cout << "          FreeLocalParam(15N;DeltaO;false)\n";
	} else if (func == "FreeGlobalParam") {
		std::cout << " FreeGlobalParam([Param];[false/true])\n";
		std::cout << " Example: FreeGlobalParam(kex;false)\n";
	} else if (func == "Minimize") {
		std::cout << " Minimize([OptionLine])\n";
		std::cout << "   [OptionLine]: Opt[0]=Val[0];Opt[1]=Val[1];....\n";
		std::cout << "       Opt: print=y/n;    [ChiSq during fitting] \n";
		std::cout << "       Opt: tol=(double)  [Tolerance] \n";
		std::cout
				<< "       Opt: maxiter=(int) [maximum number of iterations] \n";
	} else if (func == "PrintData") {
		std::cout << " PrintData([Directory/STDOUT])\n";
		std::cout
				<< "   Directory: Directory in which the data files are saved\n";
		std::cout << "      STDOUT: Print the output to the screen\n\n";
		std::cerr << " PrintData([Directory/STDOUT];[AtomName])\n";
		std::cerr
				<< "   AtomName : Only print out the data for atom AtomName\n";
	} else if (func == "PrintParam") {
		std::cout << " PrintParam([FileName];[Param])\n";
		std::cout << "   FileName: Name of output file (can also be STDOUT)\n";
		std::cout << "   Param   : Parameter to print out\n";
		std::cout << " Example: PrintParam(STDOUT;DeltaO)\n";
		std::cout << "          PrintParam(12N.out;12N)\n";
		std::cout << "          PrintParam(STDOUT;global)\n";
	} else if (func == "ReadDataset") {
		std::cout << " ReadDataset([DatasetFile])\n";
	} else if (func == "ReadParam") {
		std::cout << " ReadParam([param];[file];[NameCol];[ValCol])\n";
		std::cout << " Example: ReadParam(Omega;peak.list;0;1)\n";
	} else if (func == "ReadParam_Global") {
		std::cout << " ReadParam_Global([GlobalParamFile])\n";
	} else if (func == "SetRateType") {
		std::cout << " SetRateType([rateType])" << std::endl;
		std::cout << "        rateType could be \"standard\" or \"arrhenius\"" << std::endl;
	} else if (func == "ReadParam_Local") {
		std::cout << " ReadParam_Local([LocalParamFile])\n";
	} else if (func == "SetLocalParam") {
		std::cout << " SetLocalParam([AtomName];[Param];[value])\n";
		std::cout << " Example: SetLocalParam(all;Omega;55.3)\n";
		std::cout << "          SetLocalParam(15N;DeltaO;5.2)\n";
	} else {
		std::cout << "\b\b\b\b\b\b\b\b\b\b The function: " << func
				<< " is not recognised" << std::endl;
	}
	std::cout << std::endl;
}

std::string strip(std::string is) {
	char line[MAX_STRING_LENGTH];
	char oline[MAX_STRING_LENGTH];
	ClearBuf(line, sizeof(line));
	ClearBuf(oline, sizeof(line));
	sprintf(line, "%s", is.c_str());
	int sc = 0;
	for (unsigned int i = 0; i < is.length(); i++) {
		if (!(line[i] == ' ' || line[i] == '\t')) {
			oline[sc++] = line[i];
		}
	}
	return std::string(oline);
}

void WriteHeader() {
	std::cout << " * \n";
	std::cout << " *                         --- CATIA --- \n";
	std::cout
			<< " *         (Cpmg, Antitrosy, and Trosy Intelligent Analysis)\n";
	std::cout << " *                          Version " << CATIA_VERSION << "p"
			<< std::endl;
	std::cout << " * \n";
	std::cout << " * Copyright 2007-2008 by D. Flemming Hansen (DFH)\n";
	std::cout << " * All rights reserved.\n *" << std::endl;
	std::cout
			<< " * This software and  its related software and  documentation is provided \n";
	std::cout
			<< " * 'as is' without express or implied warrenty. The author (DFH) makes no \n";
	std::cout
			<< " * warrenties as to any matter  whatsoever  with  respect to the  program \n";
	std::cout
			<< " * and the related software and documentation. In particular, any and all \n";
	std::cout
			<< " * warranties of merchantability  and  fitness for any particular purpose \n";
	std::cout << " * are expressly excluded.\n *\n";
	std::cout
			<< " *                                         flemming@pound.med.utoronto.ca\n";
	std::cout << " *\n *\n";
	std::cout
			<< " * If you use this  program to conduct academic  research you should cite \n";
	std::cout
			<< " * the program as:\n * \n * CATIA (Cpmg, Anti-trosy, and Trosy Intelligent Analysis)  version "
			<< CATIA_VERSION << "p\n * (2008); D. Flemming Hansen\n";
	//  std::cout<<" * \n * NOT FOR COMMERCIAL USE\n";
	std::cout << " *\n";
}


int main(int argc, char** argv) {

  
  WriteHeader();
	Catia cpmg;
	char str[BUFSIZ];
	bool done = false;
	fprintf(stdout, "CATIA>");
	
	while ((!done) && (fgets(str, BUFSIZ, stdin))) {

	  
	  std::string com(str);
		com = strip(com);
		com = com.substr(0, com.find('\n'));
		if (com.length() < 2) {
			fprintf(stdout, "CATIA>");
			continue;
		}
		if (com.substr(0, 1) == "#") {
			fprintf(stdout, "CATIA>");
			continue;
		}
		std::string Func = com.substr(0, com.find('('));
		std::string Arg = com.substr(com.find('(') + 1, com.find(')')
				- com.find('(') - 1);
		std::vector<std::string> Args = split(Arg, ";");
		//

		
		if (lowercase(com.substr(0, 2)) == "ex" || lowercase(com.substr(0, 2))
				== "qu") {
			if (!(com.find('(') < com.npos) || Args[0].length() == 0) {
				returnval = 0;
			} else {
				returnval = (int) (atof(Args[0].c_str()) + 0.5);
			}
			done = true;
			continue;
		} else if (lowercase(Func.substr(0, 2)) == "he") {
			if (!(com.find('(') < com.npos) || Args[0].length() == 0) {
				WriteHelp("all");
			} else {
				WriteHelp(Args[0]);
			}
		} else if (Func == "CorrMatrix") {
			//
			// only the global parameters.
			std::vector<std::string> Param4Covar;

			std::map<std::string, int>::iterator itSI2;
			unsigned int dim = 0;
			for (itSI2 = cpmg.LastFitParamName.begin(); itSI2
					!= cpmg.LastFitParamName.end(); ++itSI2) {
				std::string ThisName = itSI2->first;
				if (ThisName.substr(ThisName.find('_'), ThisName.npos
						- ThisName.find('_')) == "_global") {
					dim++;
				}
			}
			if (dim > 0) {
				Param4Covar.resize(dim);
				for (itSI2 = cpmg.LastFitParamName.begin(); itSI2
						!= cpmg.LastFitParamName.end(); ++itSI2) {
					std::string ThisName = itSI2->first;
					if (ThisName.substr(ThisName.find('_'), ThisName.npos
							- ThisName.find('_')) == "_global") {
						Param4Covar[itSI2->second] = itSI2->first;
					}
				}
			}
			if (Param4Covar.size() > 0) {
				std::cout << std::endl;
				std::cout << "# *** Correlation Matrix ***" << std::endl;
				ublas::matrix<double> Covar(cpmg.CorrelationM(Param4Covar));
				for (unsigned int i = 0; i < Param4Covar.size(); i++) {
					std::cout.width(10);
					std::cout << Param4Covar[i].substr(0, Param4Covar[i].find(
							"_global"));
					for (unsigned int j = 0; j < Param4Covar.size(); j++) {
						char luder[20];
						sprintf(luder, "%7.3f", Covar(i, j));
						std::cout << std::string(luder);
					}
					std::cout << std::endl;
				}
			} else {
				std::cout << std::endl;
			}
		} else if (Func == "ChiSq") {
			if (Args.size() != 2) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				std::map<std::string, int>::iterator itSI;
				std::map<std::string, double> ChiSqA;
				//Initialize ChiSq Array
				for (itSI = cpmg.Atoms.begin(); itSI != cpmg.Atoms.end(); ++itSI) {
					ChiSqA[itSI->first] = 0.;
				}
				double chisq = 0.;
				double dd = 0.;
				if (lowercase(Args[0]) == "all") {
					for (unsigned int d = 0; d < cpmg.Datasets.size(); d++) {
						if (lowercase(Args[1]) == "all") {
							for (itSI = cpmg.Atoms.begin(); itSI
									!= cpmg.Atoms.end(); ++itSI) {
								//Check that this atom is in the Dataset.
								if (cpmg.Datasets[d]._atomNameToLocalAtomNumber.find(
										cpmg.AtomNumber2AtomName(itSI->second))
										== cpmg.Datasets[d]._atomNameToLocalAtomNumber.end()) {
									continue;
								} else {
								  dd = cpmg.Enorm(cpmg.Datasets[d],
								  			itSI->second);
									chisq += dd * dd;
									ChiSqA[itSI->first] += dd * dd;
								}
								//
							}
						} else {
							// is this atom here?
							if (!(cpmg.Datasets[d]._atomNameToLocalAtomNumber.find(
									Args[1])
									== cpmg.Datasets[d]._atomNameToLocalAtomNumber.end())) {

							  dd = cpmg.Enorm(cpmg.Datasets[d],
							  			cpmg.AtomName2AtomNumber(Args[1]));
								chisq += dd * dd;
							}
						}
					}
				} else {
					char line[MAX_STRING_LENGTH];
					ClearBuf(line, sizeof(line));
					sprintf(line, "%s", Args[0].c_str());
					for (unsigned int i = 0; i < Args[0].size(); i++) {
						if (line[i] == ' ' || line[i] == '\t') {
							line[i] = '_';
						}
					}
					Args[0] = std::string(line);
					signed int dset = -1;
					//
					//Find what Dataset we are looking for
					for (unsigned int d = 0; d < cpmg.Datasets.size(); d++) {
						if (cpmg.Datasets[d]._id == Args[0]) {
							dset = d;
							break;
						}
					}
					if (dset == -1) {
						std::cerr << " No Dataset matched name: " << Args[0]
								<< std::endl;
						done = false;
						continue;
					} else {
						if (lowercase(Args[1]) == "all") {
							for (itSI = cpmg.Atoms.begin(); itSI
									!= cpmg.Atoms.end(); ++itSI) {
							  dd = cpmg.Enorm(cpmg.Datasets[dset],
							  			itSI->second);
								chisq += dd * dd;
								ChiSqA[itSI->first] += dd * dd;
							}
						} else {
						  dd = cpmg.Enorm(cpmg.Datasets[dset],
						  			cpmg.AtomName2AtomNumber(Args[1]));
							chisq += dd * dd;
						}
					}
				}
				if (lowercase(Args[1]) == "all") {
					for (itSI = cpmg.Atoms.begin(); itSI != cpmg.Atoms.end(); ++itSI) {
						fprintf(stdout, "%10s %10g\n", (itSI->first).c_str(),
								ChiSqA[itSI->first]);
					}
					std::cout << "Total:" << chisq << std::endl;
				} else {
					std::cout << chisq << std::endl;
				}
				}
		}  else if (Func == "echo") {
			std::cout << Convert(Arg);
		} else if (Func == "SetInternalFixThres") {
			if (Args.size() != 1) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				cpmg.SetInternalFixThres(atof(Args[0].c_str()));
			}
		} else if (Func == "FreeLocalParam") {
			if (Args.size() != 3) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				if (lowercase(Args[2].substr(0, 1)) == "f") {
					cpmg.FreeLocalParam(Args[0], Args[1], false);
				} else if (lowercase(Args[2].substr(0, 1)) == "t") {
					cpmg.FreeLocalParam(Args[0], Args[1], true);
				} else {
					WriteHelp(Func);
					done = false;
					ErrorLine(com);
				}
			}
			}  else if (Func == "FreeGlobalParam") {
			if (Args.size() != 2) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				if (lowercase(Args[1].substr(0, 1)) == "f") {
					cpmg.FreeGlobalParam(Args[0], false);
				} else if (lowercase(Args[1].substr(0, 1)) == "t") {
					cpmg.FreeGlobalParam(Args[0], true);
				} else {
					WriteHelp(Func);
					done = false;
					ErrorLine(com);
				}
			}
			}  else if (Func == "Minimize") {
			bool conv;
			Arg = Arg + ";min=lm";
			cpmg.Minimize("all", "all", Arg, conv); //NOTE
			if (!conv) {
				std::cerr << " Minimization did not converge!\n";
				done = false;
			}
			}  else if (Func == "PrintData") {
			if (!(Args.size() == 1 || Args.size() == 2)) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else if (Args.size() == 1) {
				std::map<std::string, int>::iterator itSI;
				for (itSI = cpmg.Atoms.begin(); itSI != cpmg.Atoms.end(); ++itSI) {
					if (Args[0] == "STDOUT") {
						std::cout << "#AtomName: " << itSI->first << "\n";
						std::cout << "#\n";
						std::cout << cpmg.PrintParam(itSI->first);
						cpmg.PrintData(itSI->first, std::cout);
					} else {
						std::ofstream ofs(
								(Args[0] + "/" + itSI->first + ".dat").c_str());
						ofs << "#AtomName: " << itSI->first << "\n";
						ofs << "#\n";
						ofs << cpmg.PrintParam(itSI->first);
						cpmg.PrintData(itSI->first, ofs);
						ofs.close();
					}
				}
				} else if (Args.size() == 2) {
				if (cpmg.AtomName2AtomNumber(Args[1]) == -1) {
					std::cerr << " The residue name -->" << Args[1]
							<< "<-- is not known\n";
					WriteHelp(Func);
					ErrorLine(com);
					done = false;
				} else {
					if (Args[0] == "STDOUT") {
						std::cout << "#AtomName: " << Args[1] << "\n";
						std::cout << "#\n";
						std::cout << cpmg.PrintParam(Args[1]);
						cpmg.PrintData(Args[1], std::cout);
					} else {
						std::ofstream ofs(
								(Args[0] + "/" + Args[1] + ".dat").c_str());
						ofs << "#AtomName: " << Args[1] << "\n";
						ofs << "#\n";
						ofs << cpmg.PrintParam(Args[1]);
						cpmg.PrintData(Args[1], ofs);
						ofs.close();
					}
				}
			}
			} else if (Func == "PrintParam") {
			if (Args.size() != 2) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				if (Args[0] == "STDOUT") {
					std::cout << cpmg.PrintParam(Args[1]);
				} else {
					std::ofstream ofs(Args[0].c_str());
					ofs << cpmg.PrintParam(Args[1]);
				}

			}
		} else if (Func == "ReadDataset") {
			if (Args.size() != 1) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				std::string file = Args[0];
				cpmg.ReadDataset(file);
			}
		} else if (Func == "ReadParam") {
			if (Args.size() != 4) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				cpmg.ReadParam(Args[0], Args[1], (int) (atof(Args[2].c_str())
						+ 0.5), (int) (atof(Args[3].c_str()) + 0.5));
			}
		} else if (Func == "ReadParam_Global") {
			if (Args.size() != 1) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				std::string file = Args[0];
				cpmg.ReadParam_Global(file);
			}
		} else if (Func == "ReadParam_Local") {
			if (Args.size() != 1) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				std::string file = Args[0];
				cpmg.ReadParam_Local(file);
			}
		} else if (Func == "SetRateType") {
			if (Args.size() != 1) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
				std::string rateType = Args[0];
				cpmg.SetRateType(rateType);
			}
			//
			}  else if (Func == "SetLocalParam") {
			if (Args.size() != 3) {
				WriteHelp(Func);
				done = false;
				ErrorLine(com);
			} else {
			  cpmg.SetLocalParam(Args[0], Args[1], atof(Args[2].c_str()));
			}
			} else {
			std::cerr << " Function: " << Func << " is not recognised"
					<< std::endl;
			std::cerr << " Type 'help' to get information\n" << std::endl;
		}
		//
		//       -------   DONE!
		//
		fprintf(stdout, "CATIA>");
	  
	}
	if (!done) {
		std::cout << " EOF detected -- GoodBye!\n" << std::endl;
	} else {
		if (returnval == 0) {
			std::cout << " GoodBye!\n" << std::endl;
		} else {
			std::cout << " PROGRAM ABORTED \n" << std::endl;
		}
	}
	return (returnval);
  
}

