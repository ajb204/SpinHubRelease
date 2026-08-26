/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE4D_H
#define SLICE4D_H

#include "slice4D.hpp"


/*4D FTs running at floating point precision to preserve memory*/
//1st April 2019: DS removed to save memory: DB now does all the work
//works because we don't need to save DBs in this mode.



  void slice4D::Read()
  {
    MakeSquare();
    cout << " x(i) dimensions: " << imin << " " << imax << " " << si << " " << sig1 << endl; //proton //direct
    cout << " y(j) dimensions: " << jmin << " " << jmax << " " << sj << " " << sig2 << endl; //carbon
    cout << " z(k) dimensions: " << kmin << " " << kmax << " " << sk << " " << sig3 << endl; //proton
    cout << " a(l) dimensions: " << lmin << " " << lmax << " " << sl << " " << sig4 << endl; //carbon

    //max max max max
    //for(int i=0;i<si;++i)
    //  ivals[i]=imin+(i/(si*1.-1.))*(imax-imin);
    for(int i=0;i<si;++i)
      ivals[i]=imax-(i/(si*1.-1.))*(imax-imin);
    //for(int j=0;j<sj;++j)
    //  jvals[j]=jmin+(j/(sj*1.-1.))*(jmax-jmin);
    for(int j=0;j<sj;++j)
      jvals[j]=jmax-(j/(sj*1.-1.))*(jmax-jmin);

    for(int k=0;k<sk;++k)
      kvals[k]=kmin+(k/(sk*1.-1.))*(kmax-kmin);
    for(int l=0;l<sl;++l)
      lvals[l]=lmin+(l/(sl*1.-1.))*(lmax-lmin);

    //exit(100);
    SetIndex(); //go through peak list and figure out all positions
  }

  void slice4D::MakeSquare()
  {
    cout << "Declaring memory " << endl;
    size=si*sj*sk*sl;
    size2=(si/2+1)*sj*sk*sl;
    //cout << "spectrum size: " << size << endl;
    //get chemical shift arrays (1D)
#ifdef DOUBLE4D
    ivals=new double[si];
    jvals=new double[sj];
    kvals=new double[sk];
    lvals=new double[sl];
    //cout << "declaring main arrays" << endl;
    //DB=new double[size];
    //DS=new double[size];
    DB=new double[size];
    //DS=new float[size];
    //memset(DS, 0, sizeof(float)*size);
    //memset(DB, 0, sizeof(float)*size); //will blank after getpeak

    //DI=new double[size]; //don't need to store 4D data
    //memset(DI, 0, sizeof(double)*size);

    //P=new fftwf_complex[(si/2+1)*sj*sk*sl]; //peak shape function
    //C=new fftwf_complex[(si/2+1)*sj*sk*sl]; //will hold FT of DB

    P=fftw_alloc_complex(size2); //peak shape function

    //technically: I do not need to declare this. If I can do this in place.
    //then I do not need both DB and C.
    C=fftw_alloc_complex(size2); //will hold FT of DB

#else
    ivals=new float[si];
    jvals=new float[sj];
    kvals=new float[sk];
    lvals=new float[sl];
    //cout << "declaring main arrays" << endl;
    //DB=new double[size];
    //DS=new double[size];
    DB=new float[size];
    //DS=new float[size];
    //memset(DS, 0, sizeof(float)*size);
    //memset(DB, 0, sizeof(float)*size); //will blank after getpeak

    //DI=new double[size]; //don't need to store 4D data
    //memset(DI, 0, sizeof(double)*size);

    //P=new fftwf_complex[(si/2+1)*sj*sk*sl]; //peak shape function
    //C=new fftwf_complex[(si/2+1)*sj*sk*sl]; //will hold FT of DB

    P=fftwf_alloc_complex(size2); //peak shape function

    //technically: I do not need to declare this. If I can do this in place.
    //then I do not need both DB and C.
    C=fftwf_alloc_complex(size2); //will hold FT of DB


#endif


    n=new int[4];
    n[0]=sl;
    n[1]=sk;
    n[2]=sj;
    n[3]=si;

#ifdef DOUBLE4D
    p1   = fftw_plan_dft_r2c(4,n,&DB[0],&C[0],FFTW_ESTIMATE|FFTW_DESTROY_INPUT); //put DB FT intoC
    //p1   = fftwf_plan_dft_r2c(4,n,&DB[0],&C[0],FFTW_ESTIMATE); //put DB FT intoC
    //pinv = fftwf_plan_dft_c2r(4,n,&C[0],&DS[0],FFTW_ESTIMATE); //put C FT into DS
    pinv = fftw_plan_dft_c2r(4,n,&C[0],&DB[0],FFTW_ESTIMATE|FFTW_DESTROY_INPUT ); //put C FT into DS
    //p1   = fftw_plan_dft_r2c(4,n,&DB[0],&C[0],FFTW_MEASURE); //put DB FT into C
    //pinv = fftw_plan_dft_c2r(4,n,&C[0],&DS[0],FFTW_MEASURE); //put C FT into DS
#else
    p1   = fftwf_plan_dft_r2c(4,n,&DB[0],&C[0],FFTW_ESTIMATE|FFTW_DESTROY_INPUT); //put DB FT intoC
    //p1   = fftwf_plan_dft_r2c(4,n,&DB[0],&C[0],FFTW_ESTIMATE); //put DB FT intoC
    //pinv = fftwf_plan_dft_c2r(4,n,&C[0],&DS[0],FFTW_ESTIMATE); //put C FT into DS
    pinv = fftwf_plan_dft_c2r(4,n,&C[0],&DB[0],FFTW_ESTIMATE|FFTW_DESTROY_INPUT ); //put C FT into DS
    //p1   = fftw_plan_dft_r2c(4,n,&DB[0],&C[0],FFTW_MEASURE); //put DB FT into C
    //pinv = fftw_plan_dft_c2r(4,n,&C[0],&DS[0],FFTW_MEASURE); //put C FT into DS

#endif

    cout << "Done" << endl;
  }


  //int DoIndex(double ref,doule *vals,int ii)
 #ifdef DOUBLE4D
   int slice4D::DoIndex(double ref,double *vals,int ii)
 #else
   int slice4D::DoIndex(double ref,float *vals,int ii)
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
	if(cnt==10)
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
      if(cnt==10)
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
  void slice4D::SetIndex()
  {
    for(int p=0;p<peaks;p++) //For each peak
      {
	peakList[p].indexI=DoIndex(peakList[p].x,ivals,si);//c
	peakList[p].indexJ=DoIndex(peakList[p].y,jvals,sj);//h
	peakList[p].indexK=DoIndex(peakList[p].x,kvals,sk);//c
	peakList[p].indexL=DoIndex(peakList[p].y,lvals,sl);//h

	//cout << peakList[p].y << " " << lvals[peakList[p].indexL] << endl;
	//cout << peakList[p].x << " " << kvals[peakList[p].indexK] << endl;
	//cout << peakList[p].y << " " << jvals[peakList[p].indexJ] << endl;
	//cout << peakList[p].x << " " << ivals[peakList[p].indexI] << endl;
	//exit(100);

	/*{
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
	/*	{
	  double pk=peakList[p].y;//get carbon reference
	  peakList[p].index=0; //search to find nearest index
	  double minV=fabs(jvals[0]-pk);
	  for(int j=1;j<sj;++j) //loop over frequencencies of slice.
	    if(fabs(jvals[j]-pk)<minV)
	      { //found. store in array.
		minV=fabs(jvals[j]-pk);
		peakList[p].index=j;
	      }
	      }*/
	/*	{
	  double pk=peakList[p].x;//get proton reference
	  peakList[p].indexA=0; //search to find nearest index
	  double minV=fabs(kvals[0]-pk);
	  for(int k=1;k<sk;++k) //loop over frequencencies of slice.
	    if(fabs(kvals[k]-pk)<minV)
	      { //found. store in array.
		minV=fabs(kvals[k]-pk);
		peakList[p].indexA=k;
	      }
	}
	{
	  double pk=peakList[p].y;//get carbon reference
	  peakList[p].indexB=0; //search to find nearest index
	  double minV=fabs(lvals[0]-pk);
	  for(int l=1;l<sl;++l) //loop over frequencencies of slice.
	    if(fabs(lvals[l]-pk)<minV)
	      { //found. store in array.
		minV=fabs(lvals[l]-pk);
		peakList[p].indexB=l;
	      }
	      }*/

      }
    return;
  }

  void slice4D::BlankDBFull() //completely blank DB
  {
    memset(DB, 0, SIZEMEM*size); //blank DB
  }

  void slice4D::BlankDB() //blank just the 2D slice
  {
    for(int p=0;p<peaks;++p)
      {
	int kslice=peakList[p].indexK;
	int lslice=peakList[p].indexL;
	int loc=kslice*si*sj+lslice*si*sj*sk;
	memset(DB+loc, 0, SIZEMEM*si*sj);
	//for(int i=0;i<si;++i)
	//  for(int j=0;j<sj;++j)
	//    {
	//int ii=i+j*si+kslice*si*sj+lslice*si*sj*sk;
	//      DB[ii]=0; //2D is proton,carbon
	//    }
      }
  }


  int slice4D::CountElements()
  {
    int cnt=0;
    for(int i=0;i<size;++i)
      {
	if(fabs(DB[i]) > 0.0)
	  cnt+=1;
      }
    return cnt;

  }

#ifdef DOUBLE2D
  void slice4D::SetBlur(int p,double *DB_2D) //copy in 2D blurs to 4D
  #else
  void slice4D::SetBlur(int p,float *DB_2D) //copy in 2D blurs to 4D
 #endif
  {
    //inst.SetBlur(i,sliceLib1D[i].DB); //set blur matrix
    //if(peakList[p].name=="A50C-H")
    //  cout << peakList[p].name << endl;
    int kslice=peakList[p].indexK;
    int lslice=peakList[p].indexL;
    //cout << kvals[kslice] << " " << lvals[lslice] << " " <<  peakList[p].x << " " << peakList[p].y << " " << endl;
    int slicey=kslice*si*sj+lslice*si*sj*sk; //location of slice
    for(int i=0;i<si*sj;++i)
	{
	  int ii=i+slicey;
	  //float testy=ref[i];
	  //if(fabs(testy)>noiseLim)
	  DB[ii]+=float(DB_2D[i]); //2D is proton,carbon
	}

  }

#ifdef DOUBLE2D
  void slice4D::ReadBlur(int p,double *DS_2D) //copy blurs back to 2D DS
  #else
  void slice4D::ReadBlur(int p,float *DS_2D) //copy blurs back to 2D DS
  #endif
  {
    int kslice=peakList[p].indexK;
    int lslice=peakList[p].indexL;
    int slicey=kslice*si*sj+lslice*si*sj*sk; //location of slice
    for(int i=0;i<si*sj;++i)
      {
	  int ii=i+slicey;
	  DS_2D[i]=DB[ii]; //DB currently DS
      }

    //for(int i=0;i<sliceLib2D[p].si;++i) //update slice2
    //	  for(int j=0;j<sliceLib2D[p].sj;++j) //update slice2
    //	  sliceLib2D[p].DS[i+j*sliceLib2D[p].si]=sliceLib4D[0].DB[i+j*sliceLib4D[0].si+kslice*sliceLib4D[0].si*sliceLib4D[0].sj+lslice*sliceLib4D[0].si*sliceLib4D[0].sj*sliceLib4D[0].sk];
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

  void slice4D::SetPeak(slice4D inst)
  {
    //int si2=si/2+1;
    //for(int i=0;i<si2;i++) //goes with dim1 (C)
    //  for(int j=0;j<sj;j++) //goes with dim2 (H)
    //	for(int k=0;k<sk;k++) //goes with dim3 (H)
    //	  for(int l=0;l<sl;l++) //goes with dim4 (C)
    for(int ii=0;ii<size2;++ii)
      {
	//int ii=i+j*si2+k*si2*sj+l*si2*sj*sk;
	P[ii][0]=inst.P[ii][0];
	P[ii][1]=inst.P[ii][1];
      }
  }


  void slice4D::GetPeak()
  {
    cout << "Making peak" << endl;
    pki.resize(si); pkj.resize(sj); pkk.resize(sk); pkl.resize(sl);
#ifdef DOUBLE4D
    double di=ivals[1]-ivals[0], dj=jvals[1]-jvals[0], dk=kvals[1]-kvals[0], dl=lvals[1]-lvals[0];
#else
    float di=ivals[1]-ivals[0], dj=jvals[1]-jvals[0], dk=kvals[1]-kvals[0], dl=lvals[1]-lvals[0];
#endif
    for(int i=0;i<si;++i)
      pki[i]=Peak(di*i,0,sig1,lor1,voigt1)+Peak(di*i,di*si,sig1,lor1,voigt1);
    for(int j=0;j<sj;++j)
      pkj[j]=Peak(dj*j,0,sig2,lor2,voigt2)+Peak(dj*j,dj*sj,sig2,lor2,voigt2);
    for(int k=0;k<sk;++k)
      pkk[k]=Peak(dk*k,0,sig3,lor3,voigt3)+Peak(dk*k,dk*sk,sig3,lor3,voigt3);
    for(int l=0;l<sl;++l)
      pkl[l]=Peak(dl*l,0,sig4,lor4,voigt4)+Peak(dl*l,dl*sl,sig4,lor4,voigt4);

    for(int l=0;l<sl;++l)
      for(int k=0;k<sk;++k)
        for(int j=0;j<sj;++j)
          for(int i=0;i<si;++i)
            {
              int ii=i+j*si+k*si*sj+l*si*sj*sk;
              DB[ii]=pki[i]*pkj[j]*pkk[k]*pkl[l];
            }
#ifdef DOUBLE4D
    fftw_execute(p1);
#else
    fftwf_execute(p1);
#endif
    for(int ii=0;ii<size2;++ii)
      {
        P[ii][0]=C[ii][0];
        P[ii][1]=C[ii][1];
      }
    cout << "Peak complete" << endl;
  }

  void slice4D::RunFFT_1()
  {
    //if fft2 has been run, save it current wisdom
    //fftw_export_to_filename("FFT2");
    //fftw_forget_wisdom(); //does this save me memory and hence time?
 #ifdef DOUBLE4D
    fftw_execute(p1); //fourier transform DB
#else
  fftwf_execute(p1); //fourier transform DB
#endif 
  }

  void slice4D::RunFFT_inv()
  {
    //if fft2 has been run, save it current wisdom
    //fftw_export_to_filename("FFT2");
    //fftw_forget_wisdom(); //does this save me memory and hence time?
#ifdef DOUBLE4D
    fftw_execute(p1); //fourier transform DB
#else
    fftwf_execute(p1); //fourier transform DB
#endif
  }

  void slice4D::CalcSpec()
  {
    cout << "calcspec4D" << endl;
    if(SPARSE)
      {
        BuildSparseDB(0.0);
        const size_t n=sparseDB.size();
        memset(DB,0,SIZEMEM*size);
        #pragma omp parallel for schedule(static)
        for(long long out=0; out<(long long)n; ++out)
          {
            const SparsePt4D &o=sparseDB[out];
            double acc=0.0;
            for(size_t src=0;src<n;++src)
              {
                const SparsePt4D &sp=sparseDB[src];
                int di=abs(o.i-sp.i), dj=abs(o.j-sp.j), dk=abs(o.k-sp.k), dl=abs(o.l-sp.l);
                if(di>si/2) di=si-di;
                if(dj>sj/2) dj=sj-dj;
                if(dk>sk/2) dk=sk-dk;
                if(dl>sl/2) dl=sl-dl;
                acc += sp.val*pki[di]*pkj[dj]*pkk[dk]*pkl[dl];
              }
            DB[o.ii]=acc;
          }
      }
    else
      {
        double scale=1.0/size;
#ifdef DOUBLE4D
        fftw_execute(p1);
#else
        fftwf_execute(p1);
#endif
        for(int ii=0;ii<size2;++ii)
          {
            double a=C[ii][0], b=C[ii][1], c=P[ii][0], d=P[ii][1];
            C[ii][0]=(a*c-b*d)*scale;
            C[ii][1]=(b*c+a*d)*scale;
          }
#ifdef DOUBLE4D
        fftw_execute(pinv);
#else
        fftwf_execute(pinv);
#endif
      }
  }

  void slice4D::BuildSparseDB(double cutoff)
  {
    sparseDB.clear();
    for(int l=0;l<sl;++l)
      for(int k=0;k<sk;++k)
        for(int j=0;j<sj;++j)
          for(int i=0;i<si;++i)
            {
              int ii=i+j*si+k*si*sj+l*si*sj*sk;
              double v=DB[ii];
              if(fabs(v)>cutoff)
                sparseDB.push_back({i,j,k,l,ii,v});
            }
  }

  void slice4D::PrintSpec()
  {
    for(int p=0;p<peaks;++p)
      {
	string line ="out/slice2d/"+peakList[p].name+".dat.decon";
	//cout << "opening: " << line<<endl;
	FILE *fp;
	fp=fopen(line.c_str(),"w");

	//cout << peakList[p].name << " " << peakList[p].x << " " << kvals[peakList[p].indexH] << " " << peakList[p].y << " " << jvals[peakList[p].index] << endl;

	int kslice=peakList[p].indexK; //proton in d3
	int lslice=peakList[p].indexL; //carbon in d4

	for(int i=0;i<si;++i) //proton
	  {
	    for (int j=0;j<sj;++j) //carbon
	      {
		int ii=i+j*si+kslice*si*sj+lslice*si*sj*sk;
		fprintf(fp,"%f\t%f\t%e\t%e\t%f\t%f\n",ivals[i],jvals[j],DB[ii],DB[ii],kvals[kslice],lvals[lslice]);
	      }
	    fprintf(fp,"\n");
	  }
	fclose(fp);
      }


    FILE *fp;
    {
      cout << "Making simulated projections" << endl;
      string line ="out/xy.decon";
      fp=fopen(line.c_str(),"w");
      for (int i=0;i<si;++i)
	{
	  for (int j=0;j<sj;++j)
	    {
	      double sum=0.0;
	      for (int k=0;k<sk;++k)
		for (int l=0;l<sl;++l)
		  {
		    int ii=i+j*si+k*si*sj+l*si*sj*sk;
		    sum+=DB[ii];
		  }
	      if(fabs(sum)<=noiseVal)
		sum=0;
	      fprintf(fp,"%f\t%f\t%e\n",ivals[i],jvals[j],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }

   {
      string line ="out/za.decon";
      fp=fopen(line.c_str(),"w");
      for (int k=0;k<sk;++k)
	{
	  for (int l=0;l<sl;++l)
	    {
	      double sum=0.0;
	      for (int i=0;i<si;++i)
		for (int j=0;j<sj;++j)
		  {
		    int ii=i+j*si+k*si*sj+l*si*sj*sk;
		    sum+=DB[ii];
		  }
	      if(fabs(sum)<=noiseVal)
		sum=0;
	      fprintf(fp,"%f\t%f\t%e\n",kvals[k],lvals[l],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }

    {
      string line ="out/xz.decon";
      fp=fopen(line.c_str(),"w");
      for (int k=0;k<sk;++k)

	{
	  for (int i=0;i<si;++i)
	    {
	      double sum=0.0;
	      for (int j=0;j<sj;++j)
		for (int l=0;l<sl;++l)
		  {
		    int ii=i+j*si+k*si*sj+l*si*sj*sk;
		    sum+=DB[ii];
		  }
	      if(fabs(sum)<=noiseVal)
		sum=0;
	      fprintf(fp,"%f\t%f\t%e\n",kvals[k],ivals[i],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }

    {
      string line ="out/yz.decon";
      fp=fopen(line.c_str(),"w");

      for (int l=0;l<sl;++l)
	{
	  for (int j=0;j<sj;++j)
	    {
	      double sum=0.0;
	      for (int i=0;i<si;++i)
		for (int k=0;k<sk;++k)
		  {
		    int ii=i+j*si+k*si*sj+l*si*sj*sk;
		    sum+=DB[ii];
		  }
	      if(fabs(sum)<=noiseVal)
		sum=0;
	      fprintf(fp,"%f\t%f\t%e\n",lvals[l],jvals[j],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }
  }//end printspec



#endif
