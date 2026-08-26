#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <Dataset.h>
#include <StringMethods.h>

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

using boost::to_lower_copy;
using boost::trim_copy;
//
void Catia::PrintData(std::string atom,std::ostream& OS){
  std::vector<int> Atoms2Print;
  char line[MAX_STRING_LENGTH];
  //
  if(to_lower_copy(trim_copy(atom))=="all" ||to_lower_copy(trim_copy(atom))=="*" ){
    std::map<std::string,int>::iterator itIS; //iterator:double<-string
    for(itIS=Atoms.begin();itIS!=Atoms.end();++itIS){
      Atoms2Print.push_back(itIS->second);
    };
  } else {
    Atoms2Print.push_back(AtomName2AtomNumber(atom));
  };
  //
  double swap;
  bool print;
  for(unsigned int a=0;a<Atoms2Print.size();a++){
    for(unsigned int i=0;i<Datasets.size();i++){
      /* Now updated
      for(j=0;j<Datasets[i].Atom.size();j++){
	if(Datasets[i].Atom[j]==Atoms2Print[a]){
	  print=true;
	  break;
	};
      };
      */
      //check that we have this atom
      if(Datasets[i]._atomNameToLocalAtomNumber.find(AtomNumber2AtomName(Atoms2Print[a]))==
	 Datasets[i]._atomNameToLocalAtomNumber.end() ){
	print=false;
      } else {
	print=true;
      };
      int j=-1;
      if(print){
	j = (Datasets[i]._atomNameToLocalAtomNumber[AtomNumber2AtomName(Atoms2Print[a])]);
	OS<<"#\n";
	OS<<"#Atom:    "<<AtomNumber2AtomName(Atoms2Print[a])<<std::endl;
        OS<<"#Field:   "<<Datasets[i]._sfrq<<std::endl;
        OS<<"#Temperature:   "<<Datasets[i]._temperature<<std::endl;
	OS<<"#DataSet: "<<Datasets[i]._id<<std::endl;
	ClearBuf(line,sizeof(line));
	sprintf(line, "#%12s %13s %13s %13s\n","nu_cpmg","R2_exp","Esd(R2_exp)","R2_calc");
	OS<<line;
	//lets do some sorting .. just for fun!
	for (unsigned int k=0;k<Datasets[i].ncyc[j].size();k++){
	  for (unsigned int l=k;l<Datasets[i].ncyc[j].size();l++){
	    if(Datasets[i].ncyc[j][k]>Datasets[i].ncyc[j][l]){
	      SWAP(Datasets[i].ncyc[j][k],Datasets[i].ncyc[j][l]);
	      SWAP(Datasets[i].R2_exp[j][k],Datasets[i].R2_exp[j][l]);
	      SWAP(Datasets[i].R2_esd[j][k],Datasets[i].R2_esd[j][l]);
	      SWAP(Datasets[i].R2_calc[j][k],Datasets[i].R2_calc[j][l]);
	    };
	  };
	};
	for (unsigned int k=0;k<Datasets[i].ncyc[j].size();k++){
	  ClearBuf(line,sizeof(line));
	  sprintf(line, "%13.6e %13.6e %13.6e %13.6e\n",
		  Datasets[i].ncyc[j][k],
		  Datasets[i].R2_exp[j][k],
		  Datasets[i].R2_esd[j][k],
		  Datasets[i].R2_calc[j][k]);
	  OS<<line;
	};
	OS<<"\n   \n   \n";
      };
    };
  };
};
#undef SWAP
