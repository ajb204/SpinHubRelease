#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <StringMethods.h>
#include <Abort.h>

using boost::trim_copy;

void Catia::ReadParam_Global(std::string infile){
  std::ifstream ifs(infile.c_str());
  if(!ifs){
    std::cerr<<" Could not open the inputfile "<<infile<<"\n";
    std::cerr<<" Function .ReadParam_Global()\n";
    std::cerr<<std::endl;
    Abort(1);
  };
  char line[MAX_STRING_LENGTH];
  ClearBuf(line,sizeof(line));
  while(!ifs.eof()){
    ifs.getline(line,sizeof(line));
    std::string l(line);
    std::vector<std::string> its=split(l,"=");
    if(l.length()<1||its.size()<2){
      continue;
    };
    its[0]=trim_copy(its[0]);
    GlobalParam[its[0]]=atof(its[1].c_str());
    GlobalParamF[its[0]]=1;  // by default the parameter is free
    GlobalParamE[its[0]]=-1; // not fitted yet
  };
  return;
};
