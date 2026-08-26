/**************************************************/
/* parser and setup     */
/* A.Baldwin            */
/* 4th Jan 2017         */
/**************************************************/
#ifndef PARSE_CPP
#define PARSE_CPP

#include "deconMain.hpp"

inline bool exists_test1 (const std::string& name) 
{
    if (FILE *file = fopen(name.c_str(), "r")) 
      {
        fclose(file);
        return true;
      } 
    else 
      {
        return false;
      }   
}



//read vector of vectors to get info
void decon::parse(string file_name)
{
  cout << "Parsing: " << file_name << endl;
  //read file in to a 2D vector of strings.
  vector<vector<string> > inputfile=MakeFileVec(file_name);

  baseFile = "False";

  for (size_t i = 0; i < inputfile.size(); ++i) {
    const vector<string> &inny = inputfile[i];

    if (inny.size() < 2 || inny[0].empty() || inny[0][0] == '#') {
      continue;
    }

    const string key = inny[0];
    const string value = inny[1];

    if (key == "dim") {
      dim = stoi(value);
      cout << key << " : " << dim << endl;
    } else if (key == "dmax") {
      noiseVal = stof(value);
      cout << key << " : " << noiseVal << endl;
    } else if (key == "fac") {
      fac = stof(value);
      cout << key << " : " << fac << endl;
    } else if (key == "voigt1") {
      voigt1 = stof(value);
      cout << key << " : " << voigt1 << endl;
    } else if (key == "voigt2") {
      voigt2 = stof(value);
      cout << key << " : " << voigt2 << endl;
    } else if (key == "voigt3") {
      voigt3 = stof(value);
      cout << key << " : " << voigt3 << endl;
    } else if (key == "voigt4") {
      voigt4 = stof(value);
      cout << key << " : " << voigt4 << endl;
    } else if (key == "lor1") {
      lor1 = stof(value);
      cout << key << " : " << lor1 << endl;
    } else if (key == "lor2") {
      lor2 = stof(value);
      cout << key << " : " << lor2 << endl;
    } else if (key == "lor3") {
      lor3 = stof(value);
      cout << key << " : " << lor3 << endl;
    } else if (key == "lor4") {
      lor4 = stof(value);
      cout << key << " : " << lor4 << endl;
    } else if (key == "sig1") {
      sig1 = stof(value);
      cout << key << " : " << sig1 << endl;
    } else if (key == "sig2") {
      sig2 = stof(value);
      cout << key << " : " << sig2 << endl;
    } else if (key == "sig3") {
      sig3 = stof(value);
      cout << key << " : " << sig3 << endl;
    } else if (key == "sig4") {
      sig4 = stof(value);
      cout << key << " : " << sig4 << endl;
    } else if (key == "squash") {
      squash = stof(value);
      cout << key << " : " << squash << endl;
    } else if (key == "dec3d") {
      bore = stof(value);
      cout << key << " : " << bore << endl;
    } else if (key == "recon") {
      recon = (stoi(value) != 0);
      cout << key << " : " << recon << endl;
    } else if (key == "pseudo3D") {
      pseudo3D = (stoi(value) != 0);
      cout << key << " : " << pseudo3D << endl;
    } else if (key == "pseudo2DFit") {
      pseudo2DFit = (stoi(value) != 0);
      cout << key << " : " << pseudo2DFit << endl;
    } else if (key == "pseudo2DOutput") {
      pseudo2DOutput = (stoi(value) != 0);
      cout << key << " : " << pseudo2DOutput << endl;
    } else if (key == "symmy") {
      symmode = stof(value);
      cout << key << " : " << symmode << endl;
    } else if (key == "rand") {
      rand = stof(value);
      cout << key << " : " << rand << endl;
    } else if (key == "ncpus") {
      ncpus = stof(value);
      cout << key << " : " << ncpus << endl;
    } else if (key == "maxIter") {
      maxIter = stoi(value);
      cout << key << " : " << maxIter << endl;
    } else if (key == "maxIter3D") {
      maxIter3D = stoi(value);
      cout << key << " : " << maxIter3D << endl;
    } else if (key == "conv") {
      convVal = stof(value);
      cout << key << " : " << convVal << endl;
    } else if (key == "infile") {
      infile = value.c_str();
      cout << key << " : " << infile << endl;
    } else if (key == "peakList") {
      peakfile = value; 
      cout << key << " : " << peakfile << endl;
    } else if (key == "baseFile") {
      baseFile = value;
      cout << key << " : " << baseFile << endl;
    } else if (key == "initFile") {
      initFile = value;
      cout << key << " : " << initFile << endl;
    } else if (key == "FIT") {
      FIT=true;
      cout <<key << " : " << FIT << endl;
    } else if (key == "FitPhase") {
      const string v=value;
      FitPhase = !(v == "0" || v == "false" || v == "False" || v == "FALSE" || v == "N" || v == "n");
      cout << key << " : " << FitPhase << endl;
    } else if (key == "enhance") {
      enhance = (stoi(value) != 0);
      cout << key << " : " << enhance << endl;
    }
    else if (key == "FitRad") {
      FitRad= stof(value);
      cout <<key << " : " << FitRad << endl;
    }
    else if (key == "FitF1") {
      FitF1 = stof(value);
      FitRadFix = true;
      cout << key << " : " << FitF1 << endl;
    }
    else if (key == "FitF2") {
      FitF2 = stof(value);
      FitRadFix = true;
      cout << key << " : " << FitF2 << endl;
    }
    else if (key == "FitWidthRestrict") {
      const string v = value;
      FitWidthRestrict = !(v == "0" || v == "false" || v == "False" || v == "FALSE" || v == "N" || v == "n");
      cout << key << " : " << FitWidthRestrict << endl;
    }
  }

}



#endif
