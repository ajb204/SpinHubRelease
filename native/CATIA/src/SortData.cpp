#include <Catia.h>
#include <Dataset.h>

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
//
void Catia::SortData(){
  double swap;
  for(unsigned int i=0;i<Datasets.size();i++){
    for(unsigned int j=0;j<Datasets[i]._localToGlobalAtomIndex.size();j++){
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
    };
  };
};
