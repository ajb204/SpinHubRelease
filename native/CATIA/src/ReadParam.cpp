/*
  Flemming September 20 2007
  Read parameters from a file.

  AtomName must be in column number NameCol
  Value of the parameter must be in column number ValCol

  column are number as convention from 0,1,2,3,4,...

*/

#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <StringMethods.h>
#include <Abort.h>

using boost::trim_copy;

void Catia::ReadParam(std::string param,std::string infile,int NameCol,int ValCol){
  int GlobalAtom;
  //
  std::ifstream ifs(infile.c_str());
  if(!ifs){
    std::cerr<<" Could not open the inputfile "<<infile<<"\n";
    std::cerr<<" Function .ReadParam()\n";
    std::cerr<<std::endl;
    Abort(1);
  };
  char line[MAX_STRING_LENGTH];
  //
  while(!ifs.eof()){
    ClearBuf(line,sizeof(line));
    ifs.getline(line,sizeof(line));
    Tab2Space(line,sizeof(line));
    std::string l(line);
    std::vector<std::string> its=split(l," ");
    // do we have a comment line ?
    if(l.length()<1||its.size()<2||line[0]=='#'){
      continue;
    };
    if(!(NameCol<its.size()&&ValCol<its.size())){
      std::cerr<<" You are trying to access column number: "<<NameCol<<" and "<<ValCol<<"\n";
      std::cerr<<" of the file "<<infile<<", but only "<<its.size()<<" columns are available\n";
      std::cerr<<" Function .ReadParam();\n";
      std::cerr<<std::endl;
      Abort(1);
    };
    its[NameCol]=trim_copy(its[NameCol]);
    //see if we have an atom with name its[NameCol]
    if(AtomName2AtomNumber(its[NameCol])==-1){
      /*
	std::cerr<<" Warning: The atom: "<<its[NameCol];
	std::cerr<<" does not exists in the object\n";
	std::cerr<<" Function .ReadParam();\n";
	std::cerr<<std::endl;
      */
      continue;
    } else {
      GlobalAtom=AtomName2AtomNumber(its[NameCol]);
    };
    // Modified on September 24 2007:
    // if we dont have the parameters - just create it!
    /*
    //Do we have the parameter 'param' for this atom?
    if (LocalParam[GlobalAtom].find(param)==LocalParam[GlobalAtom].end()){
      std::cerr<<" The parameter: "<<param<<" is not declared for atom: ";
      std::cerr<<its[NameCol]<<"\n";
      std::cerr<<" Function ReadParam();\n";
      std::cerr<<std::endl;
      Abort(1);
    };
    */
    LocalParam[GlobalAtom][param]=atof(its[ValCol].c_str());
  };
  ifs.close();
  return;
};
