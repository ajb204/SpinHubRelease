#include<Catia.h>

int Catia::AtomName2AtomNumber(std::string Name){
  std::map<std::string,int>::iterator it;

  it=Atoms.find(Name);
  if ( it==Atoms.end() ){
    return -1;
  } else {
    return it->second;
  };
};
