/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE1DUSTA_C
#define SLICE1DUSTA_C


#include "slice1DUSTA.hpp"






  double slice1Dusta::ApplyIter()
  {
    double tack=0.0;
    for(int i=0;i<si;++i) //for each blur
      {
      	if(fabs(DB[i])>0.0 )
        { //if value is sensible,
          // cout << DB[i] << " " << DI[i] << " " << DS[i] << endl;
      	  DB[i]=DB[i]*fabs(DI[i]/DS[i]);   //increment (retain sign with fabs)
          //tack+=fabs(DB[i]);
          tack+=double (DB[i]);
        }
        else
      	  DB[i]=0;

      }

     //cout << tack << endl;
    return tack;
  }
  double slice1Dusta::ApplyIterBase()
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

  double slice1Dusta::GetChi2()
  {
    double chi2=0.0;
    for(int i=0;i<size;++i)
      chi2+=pow((DS[i]-DI[i]),2.);
    return chi2;
  }

  void slice1Dusta::SetBinit(string initFile)
  {
    vector<vector<string> > binit= MakeFileVec(initFile);
    cout << binit.size() << " " << si << endl;
    vector<int> loc;
    cout << "Cross peaks above theshold before Binit: " << CountElements() << endl;
    cout << binit.size() << endl;
    for(int ii=0;ii<binit.size();++ii)
      {
	//if(binit[ii].size()==4) // Reading in from a correlate.3 file
	if(binit[ii].size()==1) // Reading in from a correlate.3 file //MARCH 31ST 2025
	  {

	    //float test=atof(binit[ii][2].c_str()); // Reading in from a correlate.3 file - 3rd column
	    float test=atof(binit[ii][0].c_str()); //march 31st 2025
	    if(fabs(test) > 0.0)
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

  void slice1Dusta::Read()
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


  void slice1Dusta::ReadPipe(pipe &pipefile)
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

    cout << "Maximum intensity: " << maxInt << endl;

    cout << "Initialising deltas" << endl;
    InitBlur(); //set blur matrix

}

void slice1Dusta::ReadPipe(){
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



  void slice1Dusta::ReadBase(pipe &pipefile)
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





  void slice1Dusta::MakeSquare()
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

    //31st March 2025: these are the statements from original
    //p1   = fftwf_plan_dft_r2c_1d(si,&DB[0],&C[0],FFTW_MEASURE); //put DB FT into C
    //pinv = fftwf_plan_dft_c2r_1d(si,&C[0],&DS[0],FFTW_MEASURE); //put C FT into DS


      p1   = fftwf_plan_dft_r2c_1d(si,&DB[0],&C[0],FFTW_MEASURE|FFTW_DESTROY_INPUT); //put DB FT into C
      pinv = fftwf_plan_dft_c2r_1d(si,&C[0],&DS[0],FFTW_MEASURE|FFTW_DESTROY_INPUT); //put C FT into DS

    //p1   = fftwf_plan_dft_r2c_1d(si,&DB[0],&C[0],FFTW_PATIENT); //put DB FT into C
    // cout << "setting values " << endl;
    //pinv = fftwf_plan_dft_c2r_1d(si,&C[0],&DS[0],FFTW_PATIENT); //put C FT into DS
#endif    
    // cout << "setting values " << endl;

  }


  void slice1Dusta::MakeSquareBase()
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
  int slice1Dusta::DoIndex(double ref,double *vals,int ii)
#else
  int slice1Dusta::DoIndex(float ref,float *vals,int ii)
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
  void slice1Dusta::SetIndex()
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

  void slice1Dusta::BlurFrom3D(float *ref){

  }

  void slice1Dusta::MapBlur()
  {
    memset(DB, 0, SIZEMEM*size); //reset blur...
    for(int p=0;p<peaks;++p) //add blura values to it
      DB[peakList[p].indexI]+=DBA[p];
  }

  void slice1Dusta::UnMapBlur()
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


  void slice1Dusta::InitBlur()
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


  void slice1Dusta::InitBlurBase()
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



  void slice1Dusta::SetBlur(double *ref)
  {
    for(int i=0;i<si;++i)
      DB[i]=ref[i];
  }

  //function to return a gaussian at a specific point


  void slice1Dusta::SetPeak(const slice1D &inst)
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
  }

  void slice1Dusta::GetPeak()
  {
    //p2   = fftw_plan_dft_r2c_1d(si,&DP[0],&P[0],FFTW_ESTIMATE);   //put DP FT into P
#ifdef DOUBLE1D
    double *DBtmp;
    DBtmp=new double[si];
    memset(DBtmp, 0, SIZEMEM*si);
    // cout << "Cross peaks above theshold: " << CountElements() << endl;
    memcpy(DBtmp,DB,si*SIZEMEM);
    memcpy(DB,DBtmp,si*SIZEMEM);
    // cout << "Cross peaks above theshold: " << CountElements() << endl;
    double dC=ivals[1]-ivals[0]; //carbon delta (d1)
    for(int i=0;i<si;i++){ //goes with dim2 (C)
      DB[i]=Peak(dC*i,0,sig1, lor1, voigt1)+ Peak(dC*i,dC*(si),sig1, lor1, voigt1); //put the peakj in the four corners
      //makes C the same as python - note shifting of maximum to far right.
      // DB[i]=PseudoVoigt(dC*i,-dC,lor1,voigt1)+PseudoVoigt(dC*i,dC*si-dC,lor1,voigt1);
      //DB[i]=PseudoVoigt(dC*i,0,lor1,voigt1)+PseudoVoigt(dC*i,dC*si,lor1,voigt1);
      // cout << sig1 << " " << lor1 << " " << voigt1 << " " << DB[i] << endl;
    }
    //FILE *fp=fopen("test.out","w");
    //correlate(fp);
    //exit(100);
    //cout << "blah cross peaks above theshold: " << CountElements() << endl;
    fftw_execute(p1); //fourier transform DP into P
#else
    float *DBtmp;
    DBtmp=new float[si];
    memset(DBtmp, 0, SIZEMEM*si);
    // cout << "Cross peaks above theshold: " << CountElements() << endl;
    memcpy(DBtmp,DB,si*SIZEMEM);
    memcpy(DB,DBtmp,si*SIZEMEM);
    // cout << "Cross peaks above theshold: " << CountElements() << endl;
    float dC=ivals[1]-ivals[0]; //carbon delta (d1)
    for(int i=0;i<si;i++){ //goes with dim2 (C)
      DB[i]=Peak(dC*i,0,sig1, lor1, voigt1)+ Peak(dC*i,dC*(si),sig1, lor1, voigt1); //put the peakj in the four corners
      //makes C the same as python - note shifting of maximum to far right.
      // DB[i]=PseudoVoigt(dC*i,-dC,lor1,voigt1)+PseudoVoigt(dC*i,dC*si-dC,lor1,voigt1);
      //DB[i]=PseudoVoigt(dC*i,0,lor1,voigt1)+PseudoVoigt(dC*i,dC*si,lor1,voigt1);
      // cout << sig1 << " " << lor1 << " " << voigt1 << " " << DB[i] << endl;
    }

    //FILE *fp=fopen("test.out","w");
    //correlate(fp);
    //exit(100);
    //cout << "blah cross peaks above theshold: " << CountElements() << endl;
    fftwf_execute(p1); //fourier transform DP into P
#endif


    //cout << "blah cross peaks above theshold: " << CountElements() << endl;
    for(int i=0;i<si/2+1;i++) //goes with dim2 (C)
      {
	P[i][0]=C[i][0];
	P[i][1]=C[i][1];
      }
    // memset(DB, 0, sizeof(float)*size);
    memcpy(DB,DBtmp,si*SIZEMEM);
    //cout << "blah cross peaks above theshold: " << CountElements() << endl;
    //memset(DB, 0, sizeof(double)*size); //reset blur...
  }


  void slice1Dusta::GetPeakBase()
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


  void slice1Dusta::CalcSpec()
  {
    double scale = 1.0 / (si);

#ifdef DOUBLE1D
    fftw_execute(p1); //fourier transform blur
#else
    fftwf_execute(p1); //fourier transform blur
#endif

    //cout << "done blkur FT " << endl;
    for (int i = 0; i < si/2+1; ++i)
      {
	int ii=i;
	double a=C[ii][0]; //blur FT real
	double b=C[ii][1]; //blur FT imag
	double c=P[ii][0]; //peak FT real
	double d=P[ii][1]; //peak FT imag
	C[ii][0]=(a*c-b*d)*scale;
	C[ii][1]=(b*c+a*d)*scale;
      }

    /*inverseFT C, and place in DS*/
#ifdef DOUBLE1D
    fftw_execute(pinv);
#else
    fftwf_execute(pinv);
#endif

    if(BASE)
      { //if we have a baseline, then add this to the calculated spectrum.
	//it's possible we can do away with Cb: C is a memory place holder and can probalby do this on its own.

#ifdef DOUBLE1D
	fftw_execute(p1b); //fourier transform blur
  #else
  	fftwf_execute(p1b); //fourier transform blur
  #endif

	//cout << "done blkur FT " << endl;
	for (int i = 0; i < si/2+1; ++i)
	  {
	    int ii=i;
	    double a=Cb[ii][0]; //blur FT real
	    double b=Cb[ii][1]; //blur FT imag
	    double c=Pb[ii][0]; //peak FT real
	    double d=Pb[ii][1]; //peak FT imag
	    Cb[ii][0]=(a*c-b*d)*scale;
	    Cb[ii][1]=(b*c+a*d)*scale;
	  }

	/*inverseFT Cb, and place in DSb*/
#ifdef DOUBLE1D
	fftw_execute(pinvb);
#else
	fftwf_execute(pinvb);
#endif


	for (int i = 0; i < si; ++i)
	  DS[i]+=DSb[i];
      }
  }

  void slice1Dusta::PrintSpec()
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

  void slice1Dusta::PrintSpecPure()
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

  void slice1Dusta::Squash()
  {
    //cout << squash_window_i << endl;
    //int window_limits_i = int(float(squash_window_i)/2.0-0.5); //turn into integers
#ifdef DOUBLE1D
    int window_limits_i = round(double(squash_window_i)/2.0); //turn into integers
    vector<double> max_values; //log value of maxima
    vector<double> maximae_sorted; //vector, to be sorted by max intensity 

#else
   int window_limits_i = round(float(squash_window_i)/2.0); //turn into integers
    vector<float> max_values; //log value of maxima
    vector<float> maximae_sorted; //vector, to be sorted by max intensity
#endif
    //cout << "window_limits_i: " << window_limits_i << endl;
    vector<int> maximae;      //log location of maxima
    int max_counter = 0;
    bool maximum = true;
    for(int i=0;i<si;++i)
	{
	  if (fabs(DB[i]) > 0.0) //if there is intensity here...
	    {
	      maximum = true;
	      int bottom_edgei = max(i-1, 0);
	      int top_edgei = min(i+2, si);
	      //cout <<bottom_edgei << " " << ii << " " << top_edgei << " " << si<<endl;
	      //python
	      //int x1=max(1,ii);
	      //int x=min(x1,si-2);
	      //bottom_edgei=x-1;
	      //top_edgei=x+2;
	      for (int i2 = bottom_edgei; i2 < top_edgei; ++i2)   //walk around the peak
		{
		  //int ii2 = i2;
		  if (fabs(DI[i2]) > fabs(DI[i])) //NOTE: this is needed for a correct benchmark (7/4/25)
		  //if (fabs(DB[i2]) > fabs(DB[i])) //even though this seems better?
		    {
		    //if ((DB[ii2]) > (DB[ii]))
		      maximum = false;
		    }
		}
	      //cout << endl;
	      // cout << i << " " << maximum << " " << ivals[i] << endl;


	      if (maximum == true) {
		// cout << ii << endl;
		maximae.push_back(i); //save array value
		max_values.push_back(DB[i]); //save intensity value
	      }
	    }
	}
    // exit(100);
    //cout << DB[14316] << endl;
    //for (int i=0;i<maximae.size();++i)
    //  cout << maximae[i] << endl;
    //cout << "Maxima: " << maximae.size() << endl;
    //exit(100);

    max_counter=maximae.size();  //here's how many maxima we're dealing with
    while(maximae_sorted.size()!=max_counter) //get indicies of sorted intensities
      {
	//maxval: *max_element(balls.begin(),balls.end())
	//maxind: max_element(balls.begin(),balls.end())-balls.begin()
	int ind= max_element(max_values.begin(),max_values.end())-max_values.begin(); //index of max
	maximae_sorted.push_back(maximae[ind] ); //save reference value
	max_values.erase(max_values.begin()+ind); //remove entry from max_values
	maximae.erase(maximae.begin()+ind);  //remove entry from maximae
      }

    for (int i_count=0; i_count < max_counter; ++i_count){
      //int ii = maximae[i_count];  //I think we were previously taking unsorted intensities?
      //int ii = maximae_sorted[i_count]; //take intensities in order
      int ii = maximae_sorted[max_counter-i_count-1]; //take intensities in order
      int i = ii; //recalc indicies
      //cout << maximae_sorted[i_count] << " " << DB[ii] << endl;
      //add up intensity within region
      double sumo = 0.0;  //soak up intensity
      //double maxy=DB[ii]; //get local maximum (negatively affets benchmark: keep it off (7/4/25))
        i = max(i, window_limits_i);
        i = min(i, si-window_limits_i);
	  for (int i2 = i - window_limits_i; i2 < i + window_limits_i + 1; ++i2)
	    {
	      int ii2 = i2;
	      sumo += DB[ii2];
	      //if (fabs(DB[ii2]) > fabs(maxy)) //make sure we track local maximum
	      //{ //if we have a local maximum...
	      //ii = ii2; //move current maxima
	      //maxy=DB[ii2]; //update local maximum
	      //	      	 }
	      DB[ii2]=0.0; //set the soaked up signal to zero
	    }
	  DB[ii] = sumo; //dump all signal into local maximum
    }
    return;
  }

  void slice1Dusta::Cull(float frac)
  {
    for(int i=0;i<si;++i)
      if(fabs(DB[i])<noiseVal*frac)
	     DB[i]=0.0;
  }



  void slice1Dusta::CullCentral(double Lim)
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

  
  int slice1Dusta::CountElements()
  {
    int cnt=0;
    for(int i=0;i<si;++i)
      {
	if(fabs(DB[i])> 0.0)
	  {
	    //cout << i << " " << DB[i] << endl;
	    cnt+=1;
	  }
	
      }
    
    
    return cnt;
  }
  

  int slice1Dusta::CountElementsSym()
  {
    int cnt=0;
    for(int p=0;p<peaks;++p)
      if(fabs(DBA[p])>0)
	cnt+=1;
    return cnt;
  }

  int slice1Dusta::correlate(FILE *out_pt)
  {
    int cnt=0;
    for(int j=0;j<si;++j)
      if(fabs(DB[j])>0.0)
	{
	  cnt++;
	  //fprintf(out_pt,"%s_%i\t%f\t%f\t%f\t%e\n",refn.c_str(), cnt,refx,refy,ivals[j],DB[j]);
	  //removed march31st as uSTA project does not want this.
	  fprintf(out_pt,"%s\t%f\t%f\t%f\t%e\n",refn.c_str(),refx,refy,ivals[j],DB[j]);
	}
    return cnt;
  }

  int slice1Dusta::correlateBase(FILE *out_pt)
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
