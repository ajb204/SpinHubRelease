#include <chrono>
/*************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE2D_H
#define SLICE2D_H


#include "slice2D.hpp"





  void slice2D::LocalMax(double &hm,double &cm,int j)
  {
    int ic=peakList[j].indexI;
    int jc=peakList[j].indexJ;

    int go=1;
    while(go)
      {

	double vals[4];
	int    li[4];
	int    lj[4];

	if(ic+1<si)
	  {li[0]=ic+1;lj[0]=jc;}
	else
	  {li[0]=ic;lj[0]=jc;}

	if(ic-1>=0)
	  {li[1]=ic-1;lj[1]=jc;}
	else
	  {li[1]=ic;lj[1]=jc;}

	if(jc+1<sj)
	  {li[2]=ic;  lj[2]=jc+1;}
	else
	  {li[2]=ic;  lj[2]=jc;}

	if(jc-1>=0)
	  {li[3]=ic;  lj[3]=jc-1;}
	else
	  {li[3]=ic;  lj[3]=jc;}

	double curr =fabs(DI[ic+jc*si]);
	for(int i=0;i<4;++i)
	  vals[i]=fabs(DI[li[i]+lj[i]*si])-curr;
	int imax=0;
	double max=vals[0];
	for(int i=1;i<4;++i)
	  if(vals[i]>max)
	    {
	      imax=i;
	      max=vals[i];
	    }
	if(max>0)
	  {
	    ic=li[imax];
	    jc=lj[imax];
	  }
	else
	  {
	    hm=fabs(ivals[ic]-peakList[j].x);
	    cm=fabs(jvals[jc]-peakList[j].y);
	    return;
	  }
      }
  }

//reads in human readable data file. outdated.
  void slice2D::Read()
  {

    string newy="out/slice2d/"+refn+".dat.out";
    //cout << "reading: " <<  newy << endl;
    vector<vector<string> > newv=MakeFileVec(newy);
    for(int j=0;j<newv.size();++j)
      {
	if(newv[j].size()==3)
	  {
	    raw2D raw;
	    raw.x=atof(newv[j][0].c_str());
	    raw.y=atof(newv[j][1].c_str());
	    raw.z=atof(newv[j][2].c_str());

	    spec.push_back(raw); //store row vector for 2D
	  }
      }
   //reshape


    //cout << "Hello! I am "  << refn << " " << refx << " " << refy << endl;
    size=spec.size();
    // cout << "spectrum size: " << size << endl;

    imin=spec[0].x;
    jmin=spec[0].y;

    for(int j=1;j<size;++j)
      if(spec[j].x==imin)
	{
	  si=j;
	  break;
	}
    sj=size/si;

    imax=spec[si-1].x;
    jmax=spec[spec.size()-1].y;

    //cout << " x dimensions: " << imin << " " << imax << " " << si << endl; //proton
    //cout << " y dimensions: " << jmin << " " << jmax << " " << sj << endl; //carbon


    //cout << " refi " << refi << " " << ivals[refi] << " " << refx << endl;
    //cout << " refj " << refj << " " << jvals[refj] << " " << refy << endl;

    MakeSquare();

    //cout << "setting values " << endl;
    for(int i=0;i<si;++i)
      for(int j=0;j<sj;++j)
	DI[i+j*si]=spec[i+j*si].z;

    maxInt=fabs(DI[0]);
    for(int i=1;i<size;++i)
      if(fabs(DI[i])>maxInt)
	maxInt=fabs(DI[i]);
    //cout << "doine " << endl;


    for(int j=0;j<si;++j)
      ivals[j]=imin+(j/(si*1.-1.))*(imax-imin);
    for(int j=0;j<sj;++j)
      jvals[j]=jmin+(j/(sj*1.-1.))*(jmax-jmin);

    refi=DoIndex(refx,ivals,si);
    refj=DoIndex(refy,jvals,sj);

    //cout << refn << endl;
    //cout << "INIT: " << DI[refi+refj*si] << " " << noiseVal << endl;

    /*double mintmp;
      refi=0;
      mintmp=fabs(ivals[0]-refx);
      for(int j=0;j<si;++j)
      if(fabs(ivals[j]-refx)<mintmp)
      {
      refi=j;
      mintmp=fabs(ivals[j]-refx);
      }*/
    /*
    refj=0;
    mintmp=fabs(jvals[0]-refy);
    for(int j=0;j<sj;++j)
    {
	if(fabs(jvals[j]-refy)<mintmp)
	{
	refj=j;
	mintmp=fabs(jvals[j]-refy);
	}
	}*/
    SetIndex();


  }

  void slice2D::MapBlur()
  {
  
    memset(DB, 0, SIZEMEM*size); //reset blur...
    for(int p=0;p<peaks;++p) //add blura values to it
      {
	if(DBR[p]==1)
	  {
      
	    int ii=peakList[p].indexI+peakList[p].indexJ*si;
	    DB[ii]+=DBA[p];
	    //DB[ii]=(DB[ii]+DBA[p])*0.5;
	  }
      }
    
  }

  void slice2D::UnMapBlur()
  {
    for(int p=0;p<peaks;++p) //reset blur...
      if(DBR[p]==1) //if row is still alive...
	{
	  if(DBA[p]==0) //if no intensity is left here...
	    DBR[p]=0; //set rows to zero if we don't wany anything here.
	  if(DBR[p]==1) //if the row is still alive...
	    {
	      int ii=peakList[p].indexI+peakList[p].indexJ*si;
	      DBA[p]=DB[ii];
	    }
	  else
	    DBA[p]=0.0;
	}
  }

  void slice2D::ReadPipe(pipe &pipefile)
  {
    pipey=1;
    size=si*sj;
    MakeSquare();
    //set memory for peak FT

    for(int i=0;i<si;++i)
      ivals[i]=imin+(i/(si*1.-1.))*(imax-imin);
    for(int j=0;j<sj;++j)
      jvals[j]=jmin+(j/(sj*1.-1.))*(jmax-jmin);

    cout << " x(i) dimension: " << imin << " " << imax << " " << si << "sig " << sig1 << endl; //carbon
    cout << " y(j) dimension: " << jmin << " " << jmax << " " << sj << "sig " << sig2 << endl; //proton

   maxInt=0;
    for (int i = 0; i < si; ++i)
      for (int j = 0; j < sj; ++j)
	{
	  int ii=i+j*si;
	  DI[ii]=pipefile.GetVal2D(i,j);
	  if(DI[ii]>maxInt)
	    maxInt=DI[ii];
	}
    //InitBlur();

    cout << "Maximum intensity: " << maxInt << endl;

    cout << "Initialising deltas" << endl;
    InitBlur(); //set blur matrix

  }

  void slice2D::InitBlurPeak()
  {

    BORE=true;
#ifdef DOUBLE2D
    DBA=new double[peaks];
    DBR=new double[peaks];
#else
    DBA=new float[peaks];
    DBR=new float[peaks];
#endif
    
    SetIndex(); //map the peak list to positions in the spectrum
    memset(DB, 0, SIZEMEM*size); //reset blur...
    
    for (int p=0; p<peaks; p++){ //initialise blur at intensities found in spectrum
      //peakEntry current_peak = sliceLib2D[0].peakList[p];
      int ii=peakList[p].indexI+peakList[p].indexJ*si;
      //cout << ii << endl;
      DBR[p]=1; //place on this position
      DBA[p]=DI[ii];
    }


    MapBlur(); //DBA->DB

    
    
  }

  
  void slice2D::InitBlur()
  {
    if(symmode) //restrict choices further
      { //only allow specific locations to have deltas
	for(int p=0;p<peaks;++p) //set DB to nonzero if in index
	  {

	    int ii=peakList[p].indexI+peakList[p].indexJ*si;

	    //cout << peakList[refi].x << " " << refx << " " << refy << endl;
	    // cout << refn << " " << peakList[p].name << " " << DI[ii] << " " << peakList[p].x << " " << peakList[p].y << endl;
	    if(fabs(DI[ii])>noiseVal)
	      {
		DBR[p]=1; //place on this position
		DBA[p]=DI[ii];
	      }
	    else
	      {
		DBR[p]=0; //never place on this position
		DBA[p]=0;
	      }
	  }

	//cout <<"CRRentral: " << DBA[refp] << " " << DBR[refp] << endl;

	MapBlur();//DBA->DB
      }
    else
      { //allow all spectral frequencies to have deltas
	//int blurs=0;
	for(int i=0;i<si;++i) //set all frequencies to nonzero
	  for(int j=0;j<sj;++j) //set all frequencies to nonzero
	    {
	      int ii=i+j*si;
	      if(fabs(DI[ii])>noiseVal)
		{
		  //blurs+=1;
		  DB[ii]=DI[ii];
		}
	      else
		DB[ii]=0;
	    }
      }
  }

  double slice2D::ApplyIter()
  {
    double tack=0.0;

    // CalcSpec() has just rebuilt sparseDB from the current non-zero DB values.
    // Update only those active entries instead of scanning the full si*sj array.
    for(size_t p=0; p<sparseDB.size(); ++p)
      {
        const int ii=sparseDB[p].ii;
        DB[ii]=DB[ii]*fabs(DI[ii]/DS[ii]); // increment (retain sign with fabs)
        tack+=fabs(DB[ii]);
      }

    return tack;
  }

  double slice2D::ApplyIter2(){
    CalcSpec();
    //double tack=0;
    DB_sum = 0;
    for(int i=0;i<si;++i) //for each blur
      for(int j=0;j<sj;++j) //for each blur
	{
	  int ii = i + j * sj;
	  if (fabs(DB[ii]) > 0.0) //if value is sensible,
	    {
	      //tack += DB[ii];
	      DB[ii] = DB[ii] * (DI[ii] / DS[ii]);   //increment (retain sign with fabs)
	      //tack-=DB[ii];
	      DB_sum += DB[ii];
	    }
	  //else
	  //  DB[ii] = 0;
	}
    return DB_sum;
  }


#ifdef DOUBLE2D
  int slice2D::DoIndex(double ref,double *vals,int ii)
  #else
 int slice2D::DoIndex(double ref,float *vals,int ii)
  #endif
  {
    //Is reference within the chemical shift range?
    double maxy=max(vals[0],vals[ii-1]);
    double miny=min(vals[0],vals[ii-1]);

    double incr=fabs(maxy-miny);
    double dd=fabs(vals[1]-vals[0]);
    // cout << dd << " " << incr << " " << maxy << " " << miny << endl;
    //alias if required
    int cnt=0;
    while(ref>=maxy)
      {
	ref=ref-incr-dd;
	if(cnt==100)
	  {
	    cout << "too much folding. check the peaklist" << endl;
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
	  cout << "too much folding. check the peaklist" << endl;
	  exit(100);
	}
      cnt++;
      }

    int imin=0;
    double minV=fabs(vals[0]-ref);
    for(int i=1;i<ii;++i) //loop over frequencencies of slice.
      if(fabs(vals[i]-ref)<minV)
	{ //found. store in array.
	  minV=fabs(vals[i]-ref);
	  imin=i;
	}
    return imin;
  }



 //get integer index of closest point in ppm
  void slice2D::SetIndex()
  {
      cout << peaks << endl;
//      for(int p=0;p<peaks;p++) //For each peak
//      {
//
//          cout << peakList[p].x << " " << peakList[p].y << endl;
//      }
    for(int p=0;p<peaks;p++) //For each peak
      {
//        cout << p << " " << ivals[0] << " " << si << endl;
//          cout << peakList[20].x << " " << ivals[0] << " " << si << endl;
    // cout << peakList[p].x << " " << ivals[0] << endl;
	peakList[p].indexI=DoIndex(peakList[p].x,ivals,si);
	peakList[p].indexJ=DoIndex(peakList[p].y,jvals,sj);
  
//    cout << peakList[p].indexI << " " << peakList[p].indexJ << endl;
	/*
	{
	  double pk=peakList[p].y;//get carbon reference
	  peakList[p].index=0; //search to find nearest index
	  double minV=fabs(jvals[0]-pk);
	  for(int j=1;j<sj;++j) //loop over frequencencies of slice.
	    if(fabs(jvals[j]-pk)<minV)
	      { //found. store in array.
		minV=fabs(jvals[j]-pk);
		peakList[p].index=j;
	      }
	}
	{
	  double pk=peakList[p].x;//get proton reference
	  peakList[p].indexH=0; //search to find nearest index
	  double minV=fabs(ivals[0]-pk);
	  for(int i=1;i<si;++i) //loop over frequencencies of slice.
	    if(fabs(ivals[i]-pk)<minV)
	      { //found. store in array.
		minV=fabs(ivals[i]-pk);
		peakList[p].indexH=i;
	      }
	      }*/
      }
    return;
  }


  void slice2D::MakeSquare()
  {
 
#ifdef DOUBLE2D
   //cout << "declaring main arrays" << endl;
    DI=new double[size];
    DB=new double[size];
    DS=new double[size];

    ivals=new double[si];
    jvals=new double[sj];
   
    //DP=new double[size];
    //memset(DP, 0, sizeof(double)*size);

    if(symmode)
      {
	  DBA=new double[peaks];
	  DBR=new double[peaks];
    //memset(DBA, 0, SIZEMEM*peaks);
    //memset(DBR, 0, SIZEMEM*peaks);
      }

    P=new fftw_complex[(si/2+1)*sj];
    C=new fftw_complex[(si/2+1)*sj];
#else
   //cout << "declaring main arrays" << endl;
    DI=new float[size];
    DB=new float[size];
    DS=new float[size];

    ivals=new float[si];
    jvals=new float[sj];

    //DP=new double[size];
    //memset(DP, 0, sizeof(double)*size);

    if(symmode)
      {
	  DBA=new float[peaks];
	  DBR=new float[peaks];
    //memset(DBA, 0, SIZEMEM*peaks);
    //memset(DBR, 0, SIZEMEM*peaks);
      }


    P=new fftwf_complex[(si/2+1)*sj];
    C=new fftwf_complex[(si/2+1)*sj];
#endif

  memset(DB, 0, SIZEMEM*size);
    memset(DS, 0, SIZEMEM*size);


    //p1   = fftw_plan_dft_r2c_2d(sj,si,&DB[0],&C[0],FFTW_ESTIMATE); //put DB FT into C

    //pinv = fftw_plan_dft_c2r_2d(sj,si,&C[0],&DS[0],FFTW_ESTIMATE); //put C FT into DS


#ifdef DOUBLE2D
    p1   = fftw_plan_dft_r2c_2d(sj,si,&DB[0],&C[0],FFTW_MEASURE); //put DB FT into C
    pinv = fftw_plan_dft_c2r_2d(sj,si,&C[0],&DS[0],FFTW_MEASURE); //put C FT into DS
#else
    //p1   = fftwf_plan_dft_r2c_2d(sj,si,&DB[0],&C[0],FFTW_ESTIMATE|FFTW_DESTROY_INPUT); //put DB FT into C
    //pinv = fftwf_plan_dft_c2r_2d(sj,si,&C[0],&DS[0],FFTW_ESTIMATE|FFTW_DESTROY_INPUT); //put C FT into DS
    p1   = fftwf_plan_dft_r2c_2d(sj,si,&DB[0],&C[0],FFTW_MEASURE|FFTW_DESTROY_INPUT); //put DB FT into C
    pinv = fftwf_plan_dft_c2r_2d(sj,si,&C[0],&DS[0],FFTW_MEASURE|FFTW_DESTROY_INPUT); //put C FT into DS
    //p1   = fftwf_plan_dft_r2c_2d(sj,si,&DB[0],&C[0],FFTW_PATIENT|FFTW_DESTROY_INPUT); //put DB FT into C
    //pinv = fftwf_plan_dft_c2r_2d(sj,si,&C[0],&DS[0],FFTW_PATIENT|FFTW_DESTROY_INPUT); //put C FT into DS


#endif

    //p1   = fftwf_plan_dft_r2c_2d(sj,si,&DB[0],&C[0],FFTW_PATIENT); //put DB FT into C
    //pinv = fftwf_plan_dft_c2r_2d(sj,si,&C[0],&DS[0],FFTW_PATIENT); //put C FT into DS



  }

  void slice2D::SetBlur(double *ref)
  {
    for(int i=0;i<si;++i)
      DB[i+refj*si]=ref[i];
  }

  //function to return a gaussian at a specific point
//  double Gaus(double x,double x0,double sig){
//    return exp(-pow(x-x0,2.)/(2.*sig*sig));}
//  //function to return a gaussian at a specific point
//  double Lorentz(double x,double x0,double sig){
//    return pow(sig/2,2)/(pow(x-x0,2)+pow(sig/2,2));}
//  //function to return a 2D peak value at a given point
//  double Peak(double x1,double x2,double sig1)
//  {
//    switch(peaky){
//    case 0:
//      return Gaus(x1,x2,sig1);
//      break;
//    case 1:
//      return Lorentz(x1,x2,sig1);
//      break;}
//  }

  void slice2D::SetPeak(const slice2D &inst)
  {
    for(int ii=0;ii<(si/2+1)*sj;++ii)
      //for(int i=0;i<si;i++) //goes with dim2 (C)
      //for(int j=0;j<sj/2+1;j++) //goes with dim1 (H)
	  {
	    //int ii=i+j*si;
	    P[ii][0]=inst.P[ii][0];
	    P[ii][1]=inst.P[ii][1];
	  }
  }

  void slice2D::GetPeak()
  {
    cout << "Making peak" << endl;
    pki.resize(si);
    pkj.resize(sj);
#ifdef DOUBLE2D
    double di=ivals[1]-ivals[0];
    double dj=jvals[1]-jvals[0];
#else
    float di=ivals[1]-ivals[0];
    float dj=jvals[1]-jvals[0];
#endif
    for(int i=0;i<si;++i)
      pki[i]=Peak(di*i,0,sig1,lor1,voigt1)+Peak(di*i,di*si,sig1,lor1,voigt1);
    for(int j=0;j<sj;++j)
      pkj[j]=Peak(dj*j,0,sig2,lor2,voigt2)+Peak(dj*j,dj*sj,sig2,lor2,voigt2);

#ifdef DOUBLE2D
    double *DBcpy=new double[size];
#else
    float *DBcpy=new float[size];
#endif
    memcpy(DBcpy,DB,SIZEMEM*size);
    for(int i=0;i<si;++i)
      for(int j=0;j<sj;++j)
        DB[i+j*si]=pki[i]*pkj[j];
#ifdef DOUBLE2D
    fftw_execute(p1);
#else
    fftwf_execute(p1);
#endif
    int si2=si/2+1;
    for(int i=0;i<si2;++i)
      for(int j=0;j<sj;++j)
        {
          int ii=i+j*si2;
          P[ii][0]=C[ii][0];
          P[ii][1]=C[ii][1];
        }
    memcpy(DB,DBcpy,SIZEMEM*size);
    delete [] DBcpy;
    cout << "Peak complete." << endl;
  }

  // Set to 1 to benchmark the separable sparse convolution.
  // Set to 0 to use the original sparse O(n^2) convolution.
#ifndef SEPARABLE_SPARSE_2D
#define SEPARABLE_SPARSE_2D 1
#endif

  void slice2D::CalcSpec()
  {
    if(SPARSE)
      {
        BuildSparseDB(0.0);
        memset(DS,0,SIZEMEM*size);
        const size_t n=sparseDB.size();

#if SEPARABLE_SPARSE_2D
        // sparseDB is ordered by j then i.  Build the compact set of occupied
        // source rows and their start/end offsets.  The first 1D convolution
        // produces only occupied-j rows, so the intermediate is unique_j * si.
        std::vector<int> occupiedJ;
        std::vector<size_t> rowBegin;
        std::vector<size_t> rowEnd;
        occupiedJ.reserve(n);
        rowBegin.reserve(n);
        rowEnd.reserve(n);

        size_t start=0;
        while(start<n)
          {
            const int j=sparseDB[start].j;
            size_t end=start+1;
            while(end<n && sparseDB[end].j==j) ++end;
            occupiedJ.push_back(j);
            rowBegin.push_back(start);
            rowEnd.push_back(end);
            start=end;
          }

        const size_t nj=occupiedJ.size();
        std::vector<double> intermediate(nj*(size_t)si, 0.0);

#ifdef PROFILE_SEPARABLE_2D
        const auto sep_t0=std::chrono::high_resolution_clock::now();
#endif

        // First 1D convolution: along i, independently for each occupied j.
        #pragma omp parallel for schedule(static)
        for(long long r=0; r<(long long)nj; ++r)
          {
            const size_t begin=rowBegin[r];
            const size_t end=rowEnd[r];
            double *dst=&intermediate[(size_t)r*(size_t)si];
            for(int oi=0; oi<si; ++oi)
              {
                double acc=0.0;
                for(size_t src=begin; src<end; ++src)
                  {
                    const SparsePt2D &sp=sparseDB[src];
                    int di=abs(oi-sp.i);
                    if(di>si/2) di=si-di;
                    acc += sp.val*pki[di];
                  }
                dst[oi]=acc;
              }
          }

#ifdef PROFILE_SEPARABLE_2D
        const auto sep_t1=std::chrono::high_resolution_clock::now();
#endif

        // Second 1D convolution: for each active output, combine only the
        // occupied source-j rows from the intermediate.
        #pragma omp parallel for schedule(static)
        for(long long out=0; out<(long long)n; ++out)
          {
            const SparsePt2D &o=sparseDB[(size_t)out];
            double acc=0.0;
            for(size_t r=0; r<nj; ++r)
              {
                int dj=abs(o.j-occupiedJ[r]);
                if(dj>sj/2) dj=sj-dj;
                acc += intermediate[r*(size_t)si+(size_t)o.i]*pkj[dj];
              }
            DS[o.ii]=acc;
          }

#ifdef PROFILE_SEPARABLE_2D
        const auto sep_t2=std::chrono::high_resolution_clock::now();
        const double first_ms=std::chrono::duration<double,std::milli>(sep_t1-sep_t0).count();
        const double second_ms=std::chrono::duration<double,std::milli>(sep_t2-sep_t1).count();
        std::cout << "[PROFILE_SEPARABLE2D] sparse_n=" << n
                  << " unique_j=" << nj
                  << " intermediate=" << intermediate.size()
                  << " first_i_ms=" << first_ms
                  << " second_j_ms=" << second_ms
                  << " total_ms=" << (first_ms+second_ms) << std::endl;
#endif

#else
        #pragma omp parallel for schedule(static)
        for(long long out=0; out<(long long)n; ++out)
          {
            const SparsePt2D &o=sparseDB[(size_t)out];
            double acc=0.0;
            for(size_t src=0;src<n;++src)
              {
                const SparsePt2D &sp=sparseDB[src];
                int di=abs(o.i-sp.i);
                int dj=abs(o.j-sp.j);
                if(di>si/2) di=si-di;
                if(dj>sj/2) dj=sj-dj;
                acc += sp.val*pki[di]*pkj[dj];
              }
            DS[o.ii]=acc;
          }
#endif
      }
    else
      {
        double scale=1.0/size;
#ifdef DOUBLE2D
        fftw_execute(p1);
#else
        fftwf_execute(p1);
#endif
        int si2=si/2+1;
        for(int i=0;i<si2;++i)
          for(int j=0;j<sj;++j)
            {
              int ii=i+j*si2;
              double a=C[ii][0], b=C[ii][1], c=P[ii][0], d=P[ii][1];
              C[ii][0]=(a*c-b*d)*scale;
              C[ii][1]=(b*c+a*d)*scale;
            }
#ifdef DOUBLE2D
        fftw_execute(pinv);
#else
        fftwf_execute(pinv);
#endif
      }
  }

  void slice2D::BuildSparseDB(double cutoff)
  {
    sparseDB.clear();
    for(int j=0;j<sj;++j)
      for(int i=0;i<si;++i)
        {
          int ii=i+j*si;
          double v=DB[ii];
          if(fabs(v)>cutoff)
            sparseDB.push_back({i,j,ii,v});
        }
  }

  void slice2D::PrintSpec()
  {
    string line ="out/slice2d/"+refn+".dat.decon";
    FILE *fp;
    fp=fopen(line.c_str(),"w");

    for (int j=0;j<sj;++j) //proton
      {
	for (int i=0;i<si;++i) //carbon
	  {
	    int ii=i+j*si;

	    fprintf(fp,"%f\t%f\t%e\t%e\t%e\n",ivals[i],jvals[j],DB[ii],DS[ii],DI[ii]);
	  }
	fprintf(fp,"\n");
      }
    fclose(fp);

    string lone ="out/slice2d/"+refn+".dat.blura";
    //cout << "opening: " << line<<endl;
    fp=fopen(lone.c_str(),"w");

    //fprintf(out_pt,"%f\t%f\t%f\t%f\t%s\n" , d1,d3,val,vol,newv[n][2].c_str());

    for (int p=0;p<peaks;++p) //proton
      {
	int ii=peakList[p].indexI+peakList[p].indexJ*si;
	fprintf(fp,"%s\t%f\t%f\t%e\t%e\t%e\t%e\t%i\n",peakList[p].name.c_str(),peakList[p].x,peakList[p].y,DBA[p],DB[ii],DS[ii],DI[ii],ii);
      }
    fprintf(fp,"\n");

    fclose(fp);


  }


  void slice2D::PrintSpecPure()
  {
    printf("hi");
    string line ="out/2d.decon";
    FILE *fp;
    fp=fopen(line.c_str(),"w");

    for (int j=0;j<sj;++j) //proton
      {
	for (int i=0;i<si;++i) //carbon
	  {
	    int ii=i+j*si;

	    fprintf(fp,"%f\t%f\t%e\t%e\t%e\n",ivals[i],jvals[j],DS[ii],DB[ii],DI[ii]);
	  }
	fprintf(fp,"\n");
      }
    fclose(fp);

    /*
    string lone ="out/slice2d/"+refn+".dat.blura";
    //cout << "opening: " << line<<endl;
    fp=fopen(lone.c_str(),"w");

    //fprintf(out_pt,"%f\t%f\t%f\t%f\t%s\n" , d1,d3,val,vol,newv[n][2].c_str());

    for (int p=0;p<peaks;++p) //proton
      {
	int ii=peakList[p].indexI+peakList[p].indexJ*si;
	fprintf(fp,"%s\t%f\t%f\t%e\t%e\t%e\t%e\t%i\n",peakList[p].name.c_str(),peakList[p].x,peakList[p].y,DBA[p],DB[ii],DS[ii],DI[ii],ii);
      }
    fprintf(fp,"\n");

    fclose(fp);
    */

  }



  double slice2D::GetChi2()
  {
    double chi2=0.0;
    for(int i=0;i<size;++i)
      chi2+=pow((DS[i]-DI[i]),2.);
    return chi2;
  }



  void slice2D::Squash()
  {
#ifdef DOUBLE2D
      int window_limits_i = int(double(squash_window_i)/2.0-0.5); //turn into integers
      int window_limits_j = int(double(squash_window_j)/2.0-0.5); //turn into integers
      vector<double> max_values; //log value of maxima
      vector<double> maximae_sorted; //vector, to be sorted by max intensity
#else
      int window_limits_i = int(float(squash_window_i)/2.0-0.5); //turn into integers
      int window_limits_j = int(float(squash_window_j)/2.0-0.5); //turn into integers
      vector<float> max_values; //log value of maxima
      vector<float> maximae_sorted; //vector, to be sorted by max intensity
#endif


      vector<int> maximae;      //log location of maxima
      int max_counter = 0;
      bool maximum = true;
      for(int i=0;i<si;++i)
          for(int j=0;j<sj;++j)
          {
              int ii=i+si*j;
              if (fabs(DB[ii])> 0.0) //if there is intensity here...
              {
                  maximum = true;
                  int bottom_edgei = max(i-1, 0);
                  int bottom_edgej = max(j-1, 0);

                  int top_edgei = min(i+1, si);
                  int top_edgej = min(j+1, sj);


                  for (int i2 = bottom_edgei; i2 < top_edgei; ++i2)   //walk around the peak
                      for (int j2 = bottom_edgej; j2 < top_edgej; ++j2)  //walk around the peak
                      {

                          int ii2 = i2 + si * j2;
                          if (fabs(DB[ii2]) > fabs(DB[ii]))
                              maximum = false;
                      }

                  if (maximum == true) {
                      maximae.push_back(ii); //save array value
                      max_values.push_back(DB[ii]); //save intensity value
                  }
              }
          }

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
          int ii = maximae_sorted[i_count]; //take intensities in order
          int i = ii % si; //recalc indicies
          int j = ii / si; //recalc indicies

          i = max(i, window_limits_i);
          i = min(i, si-window_limits_i-1);
          j = max(j, window_limits_j);
          j = min(j, sj-window_limits_j-1);

          //add up intensity within region
          double sumo = 0.0;  //soak up intensity
          double maxy=DB[ii]; //get local maximum
          for (int i2 = i - window_limits_i; i2 < i + window_limits_i + 1; ++i2)
              for (int j2 = j - window_limits_j; j2 < j + window_limits_j + 1; ++j2)
              {
                  // if (pow(((i-i2)*1.0/(window_limits_i+1*1.0)),2.0)+pow(((j-j2)*1.0/(window_limits_j+1*1.0)),2.0)<=1.0)
                  // //(pow((i-i2),2)/pow(window_limits_i, 2) + pow((j-j2),2)/pow(window_limits_j, 2)<=1.0)
                  {
                    // cout << (pow(((i-i2)*1.0/(window_limits_i*1.0)),2.0)+pow(((j-j2)*1.0/(window_limits_j*1.0)),2.0)) << " " << i-i2 << " " << j-j2 << " " << window_limits_i << " " << window_limits_j << endl;
                    int ii2 = i2 + si * j2;
                    sumo += DB[ii2];
                    if (fabs(DB[ii2]) > fabs(maxy)) //make sure we track local maximum
                    { //if we have a local maximum...
                        ii = ii2; //move current maxima
                        maxy=DB[ii2]; //update local maximum
                    }
                    DB[ii2]=0.0; //set the soaked up signal to zero
                }
              }
          DB[ii] = sumo; //dump all signal into local maximum
          // exit(10);
      }
      return;
  }

    void slice2D::Cull(float frac)
  {
    for(int ii=0;ii<si*sj;++ii)
      if (fabs(DB[ii]) < noiseVal*frac)
	DB[ii] = 0.0;
  }


  void slice2D::CullCentral(double Hlim,double Clim)
  {
    for(int p=0;p<peaks;++p) //for each 1D data slice
      if(p!=refp) //for all peaks other than diagonals...
	{
	  if( fabs(ivals[peakList[refp].indexI] - ivals[peakList[p].indexI]) < Hlim)
	    if(fabs(jvals[peakList[refp].indexJ] - jvals[peakList[p].indexJ])< Clim )
	      { //if too close to diagonal....
		DBA[p]=0; //crush
		DBR[p]=0; //crush
	      }
	}
    //cout <<"Central: " << DBA[refp] << " " << DBR[refp] << endl;

  }



  int slice2D::CountElements()
  {
    int cnt=0;
    switch(symmode){
    case 0:
      for(int i=0;i<size;++i)
	if(fabs(DB[i])>0)
	  cnt+=1;
      break;
    case 1:
      for(int i=0;i<size;++i)
	if(fabs(DB[i])>0)
	  cnt+=1;
      break;}
    return cnt;

  }

  
  
  int slice2D::correlate(FILE* out_pt,string infile)
  {
    int cnt=0;
    if(BORE) //if we are doing a restricted peak list calculation...
      {
	cout << "peaks: " << peaks << endl;
	for (int p=0; p<peaks; p++){ //initialise blur at intensities found in spectrum
	  //peakEntry current_peak = sliceLib2D[0].peakList[p];
	  //int ii=peakList[p].indexI+peakList[p].indexJ*si;
	  //cout << ii << endl;
	  //DBR[p]=1; //place on this position
	  //DBA[p]=DI[ii];
	  //cout << peakList[p].name << " " << ii << " " << peakList[p].x << " " << peakList[p].y << " " << DBA[p] << " " << DB[ii] << endl;
	  cnt++;
	  fprintf(out_pt,"%s\t%f\t%f\t%e\n",peakList[p].name.c_str(),peakList[p].x,peakList[p].y,DBA[p]);
	}
      }
    else
      {
	//test this: take file name, and figure out nuclei from the file name. should get this from labb?
	string nuc1;
	string nuc2;
	
	{//split the projection file name into two to get nuclei IDs
	  cout << "  splitting projection file name into two" << endl;
	  string sentence= infile;
	  cout << "  file to split: " << sentence << endl;
	  istringstream iss(sentence);
	  std::vector<std::string> tokens;
	  std::string token;
	  while (std::getline(iss, token, '/')) {
	    if (!token.empty())
	      tokens.push_back(token);
	  }
	  //for (int i=0;i<tokens.size();++i)
	  //  cout << tokens[i] << endl;
	  //cout << tokens[tokens.size()-1] << endl;
	  istringstream is2(tokens[tokens.size()-1]);
	  std::vector<std::string> tokens2;
	  std::string token2;
	  while (std::getline(is2, token2, '.')) {
	    if (!token2.empty())
	      tokens2.push_back(token2);
	  }
	  //cout << tokens2.size() << endl;
	  //for (int i=0;i<tokens2.size();++i)
	  //  cout << tokens2[i] << endl;
	  nuc1=tokens2[0];
	  nuc2=tokens2[1];
	  //  exit(100);
	}
	
	for(int i=0;i<si;++i)
	  for(int j=0;j<sj;++j)
	    {
	      int ii=i+si*j;
	      if(fabs(DB[ii])>0)
		{
		  cnt++;
		  fprintf(out_pt,"%i\t%f\t%f\t%e\n",cnt,ivals[i],jvals[j],DB[ii]);
		}
	    }
      }
    
    return cnt;
  }
  


#endif
