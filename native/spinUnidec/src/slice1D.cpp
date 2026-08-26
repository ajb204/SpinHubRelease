/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE1D_H
#define SLICE1D_H

#include "slice1D.hpp"

  double slice1D::ApplyIter()
  {
    double tack=0.0;

    if(SPARSE)
      {
	//CalcSpec(); has just rebuilt sparseDB from the current non-zero DB entries.
        for(size_t p=0; p<sparseDB.size(); ++p)
          {
            const int ii=sparseDB[p].i;
            DB[ii]=DB[ii]*fabs(DI[ii]/DS[ii]);
            tack+=double(DB[ii]);
          }
      }
    else
      {
        for(int i=0;i<si;++i)
          {
            if(fabs(DB[i])>0.0)
              {
                DB[i]=DB[i]*fabs(DI[i]/DS[i]);
                tack+=double(DB[i]);
              }
            else
              DB[i]=0;
          }
      }

    return tack;
  }
    double slice1D::ApplyIterBase()
  {
    //bef 1.03333e+11
    //aft 1.30966e+10 1.24056e+09 9.94294e+09 0.126741

    //py
    //bef [1.05365806e+11]
    //aft [8.80716735e+10] [1.26496562e+09] 0.835865801136026
    //cout << "bef " << DBb[bval] << endl;
    double tack=0.0;
    float winny=0.0;
    for(int i=0;i<windowBvals.size();++i) //for each blur
      {
	winny+=(DI[windowBvals[i]]/DS[windowBvals[i]]);
	//cout << "winny: " << ivals[windowBvals[i]] << endl;
      }
    winny/=windowBvals.size(); //complete average
    DBb[bval]=DBb[bval]*winny;
    //cout << windowBvals.size() << endl;
    //cout << "aft " << DBb[bval] << " " << DI[bval] << " " << DS[bval] << " " << winny << endl;
    //exit(100);
    return tack;
  }

  //bef 3.65169e+11
  //aft 3.47453e+10 1.24056e+09 1.30864e+10 0.0951486

  double slice1D::GetChi2()
  {
    double chi2=0.0;
    for(int i=0;i<size;++i)
      chi2+=pow((DS[i]-DI[i]),2.);
    return chi2;
  }

  void slice1D::SetBinit(string initFile)
  {
    vector<vector<string> > binit= MakeFileVec(initFile);
    cout << binit.size() << " " << si << endl;
    vector<int> loc;
    cout << "Cross peaks above theshold before Binit: " << CountElements() << endl;
    cout << binit.size() << endl;
    for(int ii=0;ii<binit.size();++ii)
      {
	if(binit[ii].size()==4) // Reading in from a correlate.3 file
	  {

	    float test=atof(binit[ii][2].c_str()); // Reading in from a correlate.3 file - 3rd column
	    if(test!=0)
	    {
		        int imin=0;
		        float vmin=fabs(ivals[0]-test);
		        for(int i=1;i<si;++i)
		        {
		              float tast=fabs(ivals[i]-test);
	                if(tast<vmin)
                  {
	                    imin=i;
                      vmin=tast;
                  }
		        }
		        loc.push_back(imin);
	    }
	  }
      }
    for(int i=0;i<si;++i)
      {
	int tig=0;
	for(int j=0;j<loc.size();++j)
	  {
	    if(loc[j]==i)
	      {
		tig=1;
		break;
	      }
	  }
	if(tig)
	  DB[i]=DI[i];
	else
	  DB[i]=0.0;

      }
    cout << binit.size() << " " << loc.size() << endl;
    cout << "Cross peaks above theshold after Binit: " << CountElements() << endl;


  }

  void slice1D::Read()
  {

    string newy="out/slice2d/"+refn+".proj.out";
    //cout << "reading: " <<  newy << endl;

    vector<vector<string> > newv=MakeFileVec(newy);

    for(int j=0;j<newv.size();++j)
      {
	if(newv[j].size()==2)
	  {
	    raw1D raw;
	    raw.x=atof(newv[j][0].c_str());
	    raw.y=atof(newv[j][1].c_str());
	    spec.push_back(raw); //store row vector for 2D
	  }
      }
    //reshape

    size=spec.size();


    imin=spec[0].x;
    si=size;
    imax=spec[si-1].x;


    //cout << "Hello! I am " << refi << " " << refx << " " << refy << endl;
    //cout << "spectrum size: " << size << endl;
    //cout << " x dimensions: " << imin << " " << imax << " " << si << endl; //carbon

#ifdef DOUBLE1D
    ivals=new double[si]; //1D ppm values
#else
    ivals=new float[si]; //1D ppm values
#endif

    for(int j=0;j<si;++j)
      ivals[j]=imin+(j/(si*1.-1.))*(imax-imin);

    double mintmp;
    double reffo;
    switch(indirect){
    case 0:
      reffo=refx;
      break;
    case 1:
      reffo=refy;
      break;}
    refi=0;
    mintmp=fabs(ivals[0]-reffo);
    for(int j=0;j<si;++j)
      if(fabs(ivals[j]-reffo)<mintmp)
	{
	  refi=j;
	  mintmp=fabs(ivals[j]-reffo);
	}


    /*refj=0;
    mintmp=fabs(jvals[0]-refy);
    for(int j=0;j<sj;++j)
      {
	if(fabs(jvals[j]-refy)<mintmp)
	  {
	  refj=j;
	  mintmp=fabs(jvals[j]-refy);
	  }
	  }*/
    //cout << " refi " << refi << " " << ivals[refi] << " " << refx << endl;
    //cout << " refj " << refj << " " << jvals[refj] << " " << refy << endl;

    //cout << " making square " << endl;
    MakeSquare();
    for(int i=0;i<si;++i)
	DI[i]=spec[i].y;


    //cout << "doine " << endl;

    //Get maximum intensity in slice
    maxInt=fabs(DI[0]);
    for(int i=1;i<si;++i)
      if(fabs(DI[i])>maxInt)
	maxInt=fabs(DI[i]);


    //cout << "done" << endl;


    SetIndex(); //go through peak list and figure out all positions
  }


  void slice1D::ReadPipe(pipe &pipefile)
  {

    // pipey=1;
    MakeSquare();
    //set memory for peak FT

    for(int i=0;i<si;++i)
      ivals[i]=imin+(i/(si*1.-1.))*(imax-imin);

    // cout << " x(i) dimension: " << imin << " " << imax << " " << si << "sig " << sig1 << endl; //carbon
    //exit(100);
    maxInt=0;
    for (int i = 0; i < si; ++i)
	{
	  DI[i]=pipefile.GetVal1D(i);
	  if(DI[i]>maxInt)
	    maxInt=DI[i];
	}
    //InitBlur();

    //cout << "Maximum intensity: " << maxInt << endl;

    //cout << "Initialising deltas" << endl;
    InitBlur(); //set blur matrix

}


void slice1D::ReadPipeFrom3D(pipe &pipefile,int j,int k)
  {

    // pipey=1;
    MakeSquare();
    //set memory for peak FT

    for(int i=0;i<si;++i)
      ivals[i]=imin+(i/(si*1.-1.))*(imax-imin);

    //copy up intensity information from 3D.
    //need to replace this and do this from pipe file.
    //pipe.SetPos3D(0,j,k);
    //fread(&DI,sizeof(float),size,pipe.A);

    //could be done as a single fread command.
    maxInt=0.0;
    for(int i=0;i<size;++i) //get value from pipefile
      {
	DI[i]=pipefile.GetVal3D(i,j,k);
	//DI[i] = sliceLib3D[0].DI[i + j * size + k * size * sj];
	if(fabs(DI[i])>maxInt)
	  maxInt=fabs(DI[i]);
      }
    //cout << "Maximum intensity: " << maxInt << endl;
    //cout << "Initialising deltas" << endl;
    InitBlur(); //set blur matrix

    //cout << "maxint" << " " << inst.maxInt << endl;
  }






void slice1D::ReadPipe(){
#ifdef DOUBLE1D
    ivals=new double[si]; //1D ppm values
#else
    ivals=new float[si]; //1D ppm values
#endif
    for(int i=0;i<si;++i)
      ivals[i]=imin+(i/(si*1.-1.))*(imax-imin);

      double mintmp;
      double reffo;
      switch(indirect){
      case 0:
        reffo=refx;
        break;
      case 1:
        reffo=refy;
        break;}
      refi=0;
      mintmp=fabs(ivals[0]-reffo);
      for(int j=0;j<si;++j)
      if(fabs(ivals[j]-reffo)<mintmp)
    	{
    	  refi=j;
    	  mintmp=fabs(ivals[j]-reffo);
    	}


    // cout << "Hello! I am " << refi << " " << refx << " " << refy << endl;
    // cout << "spectrum size: " << size << endl;
    // cout << " x dimensions: " << imin << " " << imax << " " << si << endl; //carbon


    //cout << "done" << endl;

    //set memory for peak FT
    //P=new fftw_complex[si/2+1];
    //C=new fftw_complex [si/2+1];

    //p1   = fftw_plan_dft_r2c_1d(si,&DB[0],&C[0],FFTW_ESTIMATE); //put DB FT into C
    //pinv = fftw_plan_dft_c2r_1d(si,&C[0],&DS[0],FFTW_ESTIMATE); //put C FT into DS


    //cout << "done" << endl;


    SetIndex(); //go through peak list and figure out all positions

    /*
    double mintmp;
    double reffo;
    switch(indirect){
    case 0:
      reffo=refx;
      break;
    case 1:
      reffo=refy;
      break;}
    refi=0;
    mintmp=fabs(ivals[0]-reffo);
    for(int j=0;j<si;++j)
      if(fabs(ivals[j]-reffo)<mintmp)
	{
	  refi=j;
	  mintmp=fabs(ivals[j]-reffo);
	}
    */

    /*refj=0;
    mintmp=fabs(jvals[0]-refy);
    for(int j=0;j<sj;++j)
      {
	if(fabs(jvals[j]-refy)<mintmp)
	  {
	  refj=j;
	  mintmp=fabs(jvals[j]-refy);
	  }
	  }*/
    //cout << " refi " << refi << " " << ivals[refi] << " " << refx << endl;
    //cout << " refj " << refj << " " << jvals[refj] << " " << refy << endl;

    //cout << " making square " << endl;

  }



  void slice1D::ReadBase(pipe &pipefile)
  {
    BASE=true;
    // pipey=1;
    MakeSquareBase();
    //set memory for peak FT

    //for(int i=0;i<si;++i)
    //  ivals[i]=imin+(i/(si*1.-1.))*(imax-imin);

    cout << " x(i) dimension: " << imin << " " << imax << " " << si << "sig " << sig1 << endl; //carbon

    maxIntBase=0;
    // float centre_offset =

  //   int imin2=0;
  //   float vmin=abs(baseCentre-ivals[0]);
  //   for(int i=1;i<si;++i) //set all frequencies to nonzero
  //     {
	// float test=abs(baseCentre-ivals[i]);
	// if(test<vmin)
	//   {
	//     vmin=test;
	//     imin2=i;
	//   }
  //     }
  //

    float res = (fabs(imax-imin)/(float)si);
    // cout <<si-imin2 << " " << ivals[si-imin2] << endl;
    int roll = (int)((baseCentre-imax)/res);
    // cout <<si-imin2 << " " << (baseCentre-imax)/res << " " << imax << " " << ivals[imin2]-ivals[0] << endl;
    // exit(19);
    for (int i = 0; i < si; ++i)
	{
    if (i +roll < si){
      // cout << i << " " << i+roll << " " << ivals[i] << " " << ivals[i+roll] << endl;
      DIb[i+roll] = pipefile.GetVal1D(i);
    }
    else{

      // cout << i << " " << i-si+roll << " " << ivals[i] << " " << ivals[i-si+roll] << endl;


      DIb[i-si+roll] = pipefile.GetVal1D(i);
    }
      // DIb[i]=pipefile.GetVal1D(i);


	  if(DIb[i]>maxIntBase)
	    maxIntBase=DIb[i];
	}
    //InitBlur();

    cout << "Maximum intensity: " << maxIntBase << endl;

    cout << "Initialising deltas" << endl;

    GetPeakBase(); //copy DIb into DBb, FT, save peak shape function.
    InitBlurBase(); //set blur matrix
    cout << "X cross peaks above theshold: " << CountElements() << endl;

  }





  void slice1D::MakeSquare()
  {

    // cout << "declaring main arrays" << endl;
    // cout << si << endl;
#ifdef DOUBLE1D
   DI=new double[si];
    DB=new double[si];
    DS=new double[si];
    ivals=new double[si];
    if(symmode)
      {
	DBA=new double[peaks];
	DBR=new double[peaks];
      }
    // cout << "setting values " << sizeof(float)*si << endl;
    // cout << "setting values " << endl;
    //DP=new double[si];
    //memset(DP, 0, sizeof(double)*si);
    //set memory for peak FT
    P=new fftw_complex[si/2+1];
    C=new fftw_complex[si/2+1];
#else
    DI=new float[si];
    DB=new float[si];
    DS=new float[si];
    ivals=new float[si];
    if(symmode)
      {
	DBA=new float[peaks];
	DBR=new float[peaks];
      }
    // cout << "setting values " << sizeof(float)*si << endl;
    // cout << "setting values " << endl;
    //DP=new double[si];
    //memset(DP, 0, sizeof(double)*si);
    //set memory for peak FT
    P=new fftwf_complex[si/2+1];
    C=new fftwf_complex[si/2+1];
#endif

    memset(DB, 0, SIZEMEM*si);
    memset(DS, 0, SIZEMEM*si);



    // cout << "setting values " << sizeof(C) << " "<<sizeof(DB) << endl;

#ifdef DOUBLE1D
    p1   = fftw_plan_dft_r2c_1d(si,&DB[0],&C[0],FFTW_MEASURE); //put DB FT into C
    pinv = fftw_plan_dft_c2r_1d(si,&C[0],&DS[0],FFTW_MEASURE); //put C FT into DS
    //p1   = fftwf_plan_dft_r2c_1d(si,&DB[0],&C[0],FFTW_PATIENT); //put DB FT into C
    // cout << "setting values " << endl;
    //pinv = fftwf_plan_dft_c2r_1d(si,&C[0],&DS[0],FFTW_PATIENT); //put C FT into DS
#else
    //p1   = fftwf_plan_dft_r2c_1d(si,&DB[0],&C[0],FFTW_ESTIMATE); //put DB FT into C
    //pinv = fftwf_plan_dft_c2r_1d(si,&C[0],&DS[0],FFTW_ESTIMATE); //put C FT into DS
    p1   = fftwf_plan_dft_r2c_1d(si,&DB[0],&C[0],FFTW_MEASURE|FFTW_DESTROY_INPUT); //put DB FT into C
    pinv = fftwf_plan_dft_c2r_1d(si,&C[0],&DS[0],FFTW_MEASURE|FFTW_DESTROY_INPUT); //put C FT into DS

    //p1   = fftwf_plan_dft_r2c_1d(si,&DB[0],&C[0],FFTW_PATIENT); //put DB FT into C
    // cout << "setting values " << endl;
    //pinv = fftwf_plan_dft_c2r_1d(si,&C[0],&DS[0],FFTW_PATIENT); //put C FT into DS
#endif    
    // cout << "setting values " << endl;

  }


  void slice1D::MakeSquareBase()
  {

    cout << "declaring arrays for baselinine" << endl;

#ifdef DOUBLE1D
    DIb=new double[si];
    DBb=new double[si];
    DSb=new double[si];
    DIb2=new double[si];
    //set memory for peak FT
    Pb=new fftw_complex[si/2+1];
    Cb=new fftw_complex [si/2+1];
    p1b   = fftw_plan_dft_r2c_1d(si,&DBb[0],&Cb[0],FFTW_MEASURE); //put DB FT into C
    pinvb = fftw_plan_dft_c2r_1d(si,&Cb[0],&DSb[0],FFTW_MEASURE); //put C FT into DS
#else
    DIb=new float[si];
    DBb=new float[si];
    DSb=new float[si];
    DIb2=new float[si];
    //set memory for peak FT
    Pb=new fftwf_complex[si/2+1];
    Cb=new fftwf_complex [si/2+1];
    p1b   = fftwf_plan_dft_r2c_1d(si,&DBb[0],&Cb[0],FFTW_MEASURE); //put DB FT into C
    pinvb = fftwf_plan_dft_c2r_1d(si,&Cb[0],&DSb[0],FFTW_MEASURE); //put C FT into DS
#endif
    //cout << "setting values " << endl;
    memset(DSb, 0, SIZEMEM*si);

    //DP=new double[si];
    //memset(DP, 0, sizeof(double)*si);

  }



#ifdef DOUBLE1D
int slice1D::DoIndex(double ref,double *vals,int ii)
#else
  int slice1D::DoIndex(float ref,float *vals,int ii)
#endif
  {
    //Is reference within the chemical shift range?
    double maxy=max(vals[0],vals[ii-1]);
    double miny=min(vals[0],vals[ii-1]);

    double incr=(maxy-miny);
    double dd=vals[1]-vals[0];

    //alias if required
    int cnt=0;
    while(ref>=maxy)
      {
	ref=ref-incr-dd;
	if(cnt==100)
	  {
	    cout << "to much folding. check the peaklist" << endl;
	    exit(100);
	  }
	cnt++;
      }
    cnt=0;
    while(ref<=miny)
      {
      ref=ref+incr+dd;
      if(cnt==100)
	{
	  cout << "to much folding. check the peaklist" << endl;
	  exit(100);
	}
      cnt++;
      }


    int imin=0;
    double minV=fabs(ivals[0]-ref);
    for(int i=1;i<ii;++i) //loop over frequencencies of slice.
      if(fabs(vals[i]-ref)<minV)
	{ //found. store in array.
	  minV=fabs(vals[i]-ref);
	  imin=i;
	}
    return imin;
  }

 //get integer index of closest point in ppm
  void slice1D::SetIndex()
  {

    for(int p=0;p<peaks;p++) //For each peak
      {
	peakList[p].indexI=DoIndex(peakList[p].y,ivals,si);

      }
    //check to make sure indexing is okay
    //for(int i=0;i<line;i++)
    //    cout << indexC[i] << " " << dataC[indexL[i]] << endl;
    return;
  }

  void slice1D::BlurFrom3D(float *ref){

  }

  void slice1D::MapBlur()
  {
    memset(DB, 0, SIZEMEM*size); //reset blur...
    for(int p=0;p<peaks;++p) //add blura values to it
      DB[peakList[p].indexI]+=DBA[p];
  }

  void slice1D::UnMapBlur()
  {
    for(int p=0;p<peaks;++p) //reset blur...
      if(DBR[p]==1) //if row is still alive...
	{
	  if(DBA[p]==0) //if no intensity is left here...
	    DBR[p]=0; //set rows to zero if we don't wany anything here.
	  if(DBR[p]==1) //if the row is still alive...
	    DBA[p]=DB[peakList[p].indexI];
	}
  }


  void slice1D::InitBlur()
  {
    if(symmode) //restrict choices further
      { //only allow specific locations to have deltas
	/*
	vector <int> county;
	for(int p=0;p<peaks;++p) //set DB to nonzero if in index
	  {
	    int cnt=0;
	    for(int q=0;q<peaks;++q) //set DB to nonzero if in index
	      if(peakList[p].index==peakList[q].index)
		cnt+=1;
	    county.push_back(cnt);

	    }*/

	for(int p=0;p<peaks;++p) //set DB to nonzero if in index
	  {
	    //cout << peakList[refi].x << " " << refx << " " << refy << endl;
	    if(fabs(DI[peakList[p].indexI])>noiseVal)
	      {
		DBR[p]=1; //place on this position
		DBA[p]=DI[peakList[p].indexI];

		  /*
		if(i==refi)
		  DBA[i]=DI[peakList[i].index];
		else
		  if(peakList[i].index!=peakList[refi].index)
		    DBA[i]=DI[peakList[i].index]/county[i];
		  */
	      }
	    else
	      {
		DBR[p]=0; //never place on this position
		DBA[p]=0;
	      }
	  }
	MapBlur();//DBA->DB
      }
    else
      { //allow all spectral frequencies to have deltas
	//int cnt=0;
	// cout << noiseVal << endl;
	for(int i=0;i<si;++i) //set all frequencies to nonzero

	  if(fabs(DI[i])>noiseVal)
	    {
	    DB[i]=DI[i];
	    //cnt++;
	    }
	  else
	    DB[i]=0;
	//cout << noiseVal << " " << cnt << " " << si << endl;
      }
  }


  void slice1D::InitBlurBase()
  {
    //int cnt=0;
    // cout << noiseVal << endl;

    int imin2=0;
    float vmin=abs(baseCentre-ivals[0]);
    for(int i=1;i<si;++i) //set all frequencies to nonzero
      {
	float test=abs(baseCentre-ivals[i]);
	if(test<vmin)
	  {
	    vmin=test;
	    imin2=i;
	  }
      }





    DBb[imin2]=DI[imin2]/DIb[0];

    //cout << DBb[imin] << " " << DI[imin] << " " << DIb[0] << endl;
    //cout << DBb[imin] << " " << DI[imin]/DIb[imin] << " " << DIb[imin] << endl;

    //exit(100);
    bval=imin2;
    //exit(100);
    for(int i=windowB;i>0;--i)
      windowBvals.push_back(bval-i);
    for(int i=0;i<=windowB;++i)
      windowBvals.push_back(bval+i);

    //for(int i=0;i<windowBvals.size();++i)
    //  cout << windowBvals[i] << endl;




    if(windowBvals.size()==0)
      	windowBvals.push_back(bval);
  }



  void slice1D::SetBlur(double *ref)
  {
    for(int i=0;i<si;++i)
      DB[i]=ref[i];
  }

  //function to return a gaussian at a specific point


  void slice1D::SetPeak(const slice1D &inst)
  {
    for(int i=0;i<si/2+1;i++) //goes with dim1 (core)
      {
	P[i][0]=inst.P[i][0];
	P[i][1]=inst.P[i][1];
      }
    /*for(int i=0;i<si;i++) //goes with dim1 (core)
      {
	DP[i]=inst.DP[i];
	DP[i]=inst.DP[i];
	}*/
    peakShape=inst.peakShape;
  }

  void slice1D::GetPeak()
  {
    peakShape.resize(si);
#ifdef DOUBLE1D
    double dC=ivals[1]-ivals[0];
#else
    float dC=ivals[1]-ivals[0];
#endif
    for(int i=0;i<si;++i)
      peakShape[i]=Peak(dC*i,0,sig1,lor1,voigt1)
                  +Peak(dC*i,dC*si,sig1,lor1,voigt1);

#ifdef DOUBLE1D
    double *DBtmp=new double[si];
#else
    float *DBtmp=new float[si];
#endif
    memcpy(DBtmp,DB,si*SIZEMEM);
    for(int i=0;i<si;++i)
      DB[i]=peakShape[i];
#ifdef DOUBLE1D
    fftw_execute(p1);
#else
    fftwf_execute(p1);
#endif
    for(int i=0;i<si/2+1;++i)
      {
        P[i][0]=C[i][0];
        P[i][1]=C[i][1];
      }
    memcpy(DB,DBtmp,si*SIZEMEM);
    delete [] DBtmp;
  }

  void slice1D::GetPeakBase()
  {
    //p2   = fftw_plan_dft_r2c_1d(si,&DP[0],&P[0],FFTW_ESTIMATE);   //put DP FT into P
    //float *DBtmp;
    //DBtmp=new float[si];
    //memset(DBtmp, 0, sizeof(float)*si); //set
    //cout << "blah cross peaks above theshold: " << CountElements() << endl;
    //memcpy(DBtmp,DBb,si*sizeof(float)); //populate DBtmp
    //memcpy(DBb,DBtmp,si*sizeof(float));
    cout << "Yh cross peaks above theshold: " << CountElements() << endl;
    //float dC=ivals[1]-ivals[0]; //carbon delta (d1)
    for(int i=0;i<si;i++)
      DBb[i]=DIb[i];



    //cout << "blah cross peaks above theshold: " << CountElements() << endl;
 #ifdef DOUBLE1D
    fftw_execute(p1b); //fourier transform DP into P
    #else
   fftwf_execute(p1b); //fourier transform DP into P
    #endif
    //cout << "blah cross peaks above theshold: " << CountElements() << endl;
    for(int i=0;i<si/2+1;i++) //goes with dim2 (C)
      {
	Pb[i][0]=Cb[i][0];
	Pb[i][1]=Cb[i][1];
      }
    memset(DBb, 0, SIZEMEM*si); //blank it...

    // memset(DB, 0, sizeof(float)*size);
    //memcpy(DB,DBtmp,si*sizeof(float));//restore
    //cout << "blah cross peaks above theshold: " << CountElements() << endl;
    //memset(DB, 0, sizeof(double)*size); //reset blur...
  }


  void slice1D::CalcSpec()
  {
    if(SPARSE)
      {
        BuildSparseDB(0.0);
        memset(DS,0,SIZEMEM*si);
        const size_t n=sparseDB.size();
        //#pragma omp parallel for schedule(static)
        for(long long out=0; out<(long long)n; ++out)
          {
            const SparsePt1D &o=sparseDB[out];
            double acc=0.0;
            for(size_t src=0;src<n;++src)
              {
                const SparsePt1D &sp=sparseDB[src];
                int d=abs(o.i-sp.i);
                if(d>si/2) d=si-d;
                acc += sp.val*peakShape[d];
              }
            DS[o.i]=acc;
          }
	
      }
    else
      {
        double scale=1.0/si;
#ifdef DOUBLE1D
        fftw_execute(p1);
#else
        fftwf_execute(p1);
#endif
        for(int i=0;i<si/2+1;++i)
          {
            double a=C[i][0], b=C[i][1], c=P[i][0], d=P[i][1];
            C[i][0]=(a*c-b*d)*scale;
            C[i][1]=(b*c+a*d)*scale;
          }
#ifdef DOUBLE1D
        fftw_execute(pinv);
#else
        fftwf_execute(pinv);
#endif
      }

    if(BASE)
      {
        double scale=1.0/si;
#ifdef DOUBLE1D
        fftw_execute(p1b);
#else
        fftwf_execute(p1b);
#endif
        for(int i=0;i<si/2+1;++i)
          {
            double a=Cb[i][0], b=Cb[i][1], c=Pb[i][0], d=Pb[i][1];
            Cb[i][0]=(a*c-b*d)*scale;
            Cb[i][1]=(b*c+a*d)*scale;
          }
#ifdef DOUBLE1D
        fftw_execute(pinvb);
#else
        fftwf_execute(pinvb);
#endif
#ifdef DOUBLE1D
        if(p1b) fftw_destroy_plan(p1b);
        if(pinvb) fftw_destroy_plan(pinvb);
#else
        if(p1b) fftwf_destroy_plan(p1b);
        if(pinvb) fftwf_destroy_plan(pinvb);
#endif
        for(int i=0;i<si;++i)
          DS[i]+=DSb[i];
      }
  }

  void slice1D::BuildSparseDB(double cutoff)
  {
    sparseDB.clear();
    for(int i=0;i<si;++i)
      {
        double v=DB[i];
        if(fabs(v)>cutoff)
          sparseDB.push_back({i,v});
      }
  }

  void slice1D::PrintSpec()
  {
    string line="out/slice2d/"+refn+".proj.decon";
    //cout << "opening: " << line<<endl;
    FILE *fp;
    fp=fopen(line.c_str(),"w");
    for (int i=0;i<si;++i) //carbon
      {
	int ii=i;
	fprintf(fp,"%f\t%e\t%e\t%e\n",ivals[i],DB[ii],DS[ii],DI[ii]);
      }
    fclose(fp);
  }

  void slice1D::PrintSpecPure()
  {
    string line="out/slice2d/"+refn+".proj.decon";
    //cout << "opening: " << line<<endl;
    FILE *fp;
    fp=fopen(line.c_str(),"w");
    for (int i=0;i<si;++i) //carbon
      {
	int ii=i;
	fprintf(fp,"%f\t%e\t%e\t%e\n",ivals[i],DB[ii],DS[ii],DI[ii]);
      }
    fclose(fp);
  }

  void slice1D::Squash()
  {
    // Keep the 1D clustering behaviour consistent with the multidimensional
    // modes: use a linewidth-sized neighbourhood, process the strongest
    // maxima first, and retain the integrated intensity at the strongest
    // point within each claimed neighbourhood.
#ifdef DOUBLE1D
    int window_limits_i = max(1, int(double(squash_window_i)/2.0));
    vector<double> max_values;
    vector<double> maximae_sorted;
#else
    int window_limits_i = max(1, int(float(squash_window_i)/2.0));
    vector<float> max_values;
    vector<float> maximae_sorted;
#endif
    vector<int> maximae;
    int max_counter = 0;
    bool maximum = true;

    for(int i=0;i<si;++i)
      {
        if (DB[i] != 0)
          {
            maximum = true;
            int bottom_edgei = max(i-1, 0);
            int top_edgei = min(i+2, si);
            for (int i2 = bottom_edgei; i2 < top_edgei; ++i2)
              {
                if (fabs(DB[i2]) > fabs(DB[i]))
                  maximum = false;
              }

            if (maximum == true)
              {
                maximae.push_back(i);
                // Maxima detection is magnitude based, so ranking should be too.
                max_values.push_back(fabs(DB[i]));
              }
          }
      }

    max_counter=maximae.size();
    while(maximae_sorted.size()!=max_counter)
      {
        int ind=max_element(max_values.begin(),max_values.end())-max_values.begin();
        maximae_sorted.push_back(maximae[ind]);
        max_values.erase(max_values.begin()+ind);
        maximae.erase(maximae.begin()+ind);
      }

    // Process strongest first.  Once a strong maximum claims a local region,
    // weaker maxima in that region see zero intensity and cannot create a
    // second peak from the same physical resonance.
    for (int i_count=0; i_count < max_counter; ++i_count)
      {
        int ii = maximae_sorted[i_count];

        // This maximum may already have been absorbed by a stronger neighbour.
        if (DB[ii] == 0.0)
          continue;

        int i = ii;
        i = max(i, window_limits_i);
        i = min(i, si-window_limits_i-1);

        double sumo = 0.0;
        double maxy = DB[ii];
        int maxii = ii;
        for (int i2 = i-window_limits_i; i2 < i+window_limits_i+1; ++i2)
          {
            sumo += DB[i2];
            if (fabs(DB[i2]) > fabs(maxy))
              {
                maxii = i2;
                maxy = DB[i2];
              }
            DB[i2]=0.0;
          }
        DB[maxii] = sumo;
      }
  }

  void slice1D::Cull(float frac)
  {
    for(int i=0;i<si;++i)
      if(fabs(DB[i])<noiseVal*frac)
	     DB[i]=0.0;
  }



  void slice1D::CullCentral(double Lim)
  {
    for(int p=0;p<peaks;++p) //for each 1D data slice
      if(p!=refp) //for all peaks other than diagonals...
	{
	  if( fabs(ivals[peakList[refp].indexI] - ivals[peakList[p].indexI]) < Lim)
	    { //if too close to diagonal....
	      DBA[p]=0; //crush
	      DBR[p]=0; //crush
	    }
	}
  }


  int slice1D::CountElements()
    {
      int cnt=0;
      for(int i=0;i<si;++i)
	if(fabs(DB[i])>0.0){
	  cnt+=1;
  }
      return cnt;
    }


  int slice1D::CountElementsSym()
  {
    int cnt=0;
    for(int p=0;p<peaks;++p)
      if(fabs(DBA[p])>0)
	cnt+=1;
    return cnt;
  }

  int slice1D::correlate(FILE *out_pt)
  {
    int cnt=0;
    for(int j=0;j<si;++j)
      if(fabs(DB[j])>0.0)
	{
	  cnt++;
	  fprintf(out_pt,"%s_%i\t%f\t%f\t%f\t%e\n",refn.c_str(), cnt,refx,refy,ivals[j],DB[j]);
	}
    return cnt;
  }

  // Standalone/pseudo2D projection output.  Keep correlate() unchanged because
  // 3D bore/slice workflows require refn, refx and refy in their five-column list.
  int slice1D::correlatePure1D(FILE *out_pt)
  {
    int cnt=0;
    for(int j=0;j<si;++j)
      if(fabs(DB[j])>0.0)
        {
          cnt++;
          fprintf(out_pt,"%i\t%f\t%e\n",cnt,ivals[j],DB[j]);
        }
    return cnt;
  }

  int slice1D::correlateBase(FILE *out_pt)
  {
    int cnt=0;
    for(int j=0;j<si;++j)
      if(fabs(DBb[j])>0.0)
	{
	  cnt++;
	  fprintf(out_pt,"%s\t%f\t%f\t%f\t%e\n",refn.c_str(),refx,refy,ivals[j],DBb[j]);
	}
    return cnt;
  }


#endif
