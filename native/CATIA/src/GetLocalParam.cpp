#include<Catia.h>
#include <Abort.h>

double Catia::GetLocalParam(std::string AtomName,std::string ParamName){
  //
  int GlobalAtomNumber=AtomName2AtomNumber(AtomName);
  //
  // Check is atom name is available in the object
  if(GlobalAtomNumber==-1){
    std::cerr<<" The atom: "<<AtomName<<" is not available in the object \n";
    std::cerr<<" Function .GetLocalParam();\n";
    std::cerr<<std::endl;
    Abort(1);
  };
  if(BResolveParam(LocalParam[GlobalAtomNumber],ParamName)){
    return LocalParam[GlobalAtomNumber][ParamName];
  } else {
    std::cerr<<" The parameter "<<ParamName<<" is not available for atom: "<<AtomName<<" \n";
    std::cerr<<" Function .GetLocalParam();\n";
    std::cerr<<std::endl;
    Abort(1);
  };
};
