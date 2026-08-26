#include "deconMain.hpp"

/***********************************************************************/
/* Bayesian Deconvolver for NMR data                                   */
/***********************************************************************/
// Based on an algorithm from Michael Marty
// 3rd March 2014. A Baldwin
//
// Specify:
// 1) file name
// 2) sigma for peak shape in 2 dimensions
// 3) ppm values to define a square for noise determination
// 4) a convergence value for the algorithm
//
// File will create out/output.txt, the deconvolved intensities
// and out/reconv.txt, containing raw data and back calcuated
// fit to data
//
// TO INCLUDE:
// Currently all ppm values above twice the noise are considered as
// possible targets for giving rise to peaks. This is overfitting.
// An approach such as Ftests need to be incorporated to trim the
// list of possible peaks going in.
//
// The peak shape is currently gaussian. This is not correct. Something
// better needs to be incorporated. Data will have to be analysed perhaps
// with a gaus-to-lorentz transform to make this easy?
//
// Then we're ready for the big time :)
//



/*

//Performs an F-test. I have left it the crude way, which is to simply test if the reduced chi-squared is lower. With GSL, it performs an actual F-test.
double ftest(double chisq_simple,int p1,double chisq_complex,int p2,int n)
{
  double v1=n-p1;  //datapoints minus parameter number for model 1 (simple) dof simple
  double v2=n-p2;  //datapoints minus parameter number for model 2 (complex) dof complex
  double Z=v1-v2;        // (difference in degrees of freedom
  double F=(chisq_simple-chisq_complex)/((v1-v2)*(chisq_complex/v2));   // the 'F' statistic. useless on its own.
  //double ex=v2/(v2+Z*F);

  cout << (chisq_simple-chisq_complex)/chisq_complex*(v2/(v1-v2)) << " " << v2 << " " << v1 << " " << Z << endl;
  //cout << chisq_complex/v2 << " " << chisq_simple/v1 << endl;*/
/*
  double plvl=1.0;
  if(F>1E3)
    plvl=0.0;
  else
    plvl = gsl_sf_beta_inc(v2/2.,Z/2.,ex);

  cout << plvl << endl;

  return plvl;
*/
/*return 0;
  }*/





/*
//work ou the max value in the noise box
double GetNoiseMax(double *dataDH,double *dataDC,double *dataInt,int lines,double noiseYmax,double noiseYmin,double noiseXmax,double noiseXmin)
{
  double noiseVal=0.0;
  for(int i=0;i<lines;i++)
    if(dataDH[i]>noiseXmin && dataDH[i]<noiseXmax && dataDC[i]>noiseYmin && dataDC[i]<noiseYmax)
      if(dataInt[i]>noiseVal)
	noiseVal=dataInt[i];
  return noiseVal;
}
void makefabs(double *peak,int *Size)
{
  for(int i=0;i<Size[0];i++)
    for(int j=0;j<Size[1];j++)
	peak[i+j*Size[0]]=fabs(peak[i+j*Size[0]]);
}*/


//void GetReady(decon &dec, parse &input,int argc,char *inputfile);




int main(int argc,char *argv[])
{

//    if (argv[0] == 'Decon_windows'){
//        cout << "Windows" << endl;
//    }

  decon dec;
  dec.splash(argc,argv,false);


  switch(dec.dim){
  case 1: //1D unidec
    if (dec.STD ==1){
      dec.Protocol1D(argc,argv);
    }
    else{
      dec.Protocol1D(argc,argv);
    }
    break;
  case 2: //2D unidec / physical pseudo2D restrained FIT
    if (dec.pseudo2DFit) dec.Protocol2PFit(argc,argv);
    else dec.Protocol2D(argc,argv);
    break;
  case 3: //3D using 1D slices
    if (dec.pseudo3D) dec.Protocol3P(argc,argv);
    else dec.Protocol3D(argc,argv);
    break;
  case 4: //4D via 2D slices
    dec.Protocol4D(argc,argv);
    break;}

    //case 5: //3D complete
    // dec.Protocol3Dfull(argc,argv);
    //break;}
  cout << "exiting cleanly." << endl;

  return 0;
}
