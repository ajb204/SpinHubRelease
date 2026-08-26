#include<Catia.h>
#include <Abort.h>

double Catia::GetGlobalParam(std::string name){
  if(BResolveParam(GlobalParam,name)){
    return GlobalParam[name];
  } else {
    std::cerr<<" "<<name<<" is not available in the set of Global Parameters\n";
    std::cerr<<" Function .GetGlobalParam();\n";
    std::cerr<<std::endl;
    Abort(1);
  };
};
