#include <Catia.h>
#include <Abort.h>

void Catia::AddAtom(std::string Name){
  //First check that the Name is not already in the Atoms map
  if(!(AtomName2AtomNumber(Name)==-1)){
    std::cerr<<" Your are trying to add an AtomName, which already\n";
    std::cerr<<" exists in the map<string,int> Atoms"<<std::endl;
    Abort(1);
  };
  int nta=Atoms.size();
  Atoms[Name]=nta;
  return;
};
