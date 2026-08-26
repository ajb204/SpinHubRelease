/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef SLICE3D_H
#define SLICE3D_H

#include "slice3D.hpp"




  void slice3D::SetMem()
  {

    size=si*sj*sk;
    size2=(si/2+1)*sj*sk;
    //cout << "spectrum size: " << size << endl;

    //get chemical shift arrays (1D)
#ifdef DOUBLE3D
    ivals=new double[si];
    jvals=new double[sj];
    kvals=new double[sk];
#else
    ivals=new float[si];
    jvals=new float[sj];
    kvals=new float[sk];
#endif

    for(int i=0;i<si;++i)
      ivals[i]=imin+(i/(si*1.-1.))*(imax-imin);
    for(int j=0;j<sj;++j)
      jvals[j]=jmin+(j/(sj*1.-1.))*(jmax-jmin);
    for(int k=0;k<sk;++k)
      kvals[k]=kmin+(k/(sk*1.-1.))*(kmax-kmin);

    if(sig1==0)
      sig1=abs(ivals[1]-ivals[0]);
    if(sig2==0)
      sig2=abs(jvals[1]-jvals[0]);
    if(sig3==0)
      sig3=abs(kvals[1]-kvals[0]);

    cout << " x(i) dimension: " << imin << " " << imax << " " << si << "sig " << sig1 << endl; //carbon
    cout << " y(j) dimension: " << jmin << " " << jmax << " " << sj << "sig " << sig2 << endl; //proton
    cout << " z(k) dimension: " << kmin << " " << kmax << " " << sk << "sig " << sig3 << endl; //proton


    cout << "Making square " << endl;
    MakeSquare();
    // cout << "done" << endl;


    cout << "Setting indices" << endl;
    SetIndex(); //go through peak list and figure out all positions

    cout << "Setting peak" << endl;
    GetPeak();  //define peak shape function
    cout << "done" << endl;
  }

  void slice3D::ReadPipe(pipe &pipefile)
  {
    pipey=1;
    //declare memory only if reading in full 3D
#ifdef DOUBLE3D
    DI=new double[size];
#else
    DI=new float[size];    
#endif

    maxInt=0;
    int c1=pipefile.size[pipefile.dimord[0]];
    int c2=pipefile.size[pipefile.dimord[1]]*c1;
    // cout << c1 << " " << c2/c1 << endl;
    //
    fseek(pipefile.A, 512*sizeof(float), 0);
    float *Reffy; 
    Reffy=new float[size];
    fread(Reffy,sizeof(float),size,pipefile.A);

    for (int i = 0; i < si; ++i)
      for (int j = 0; j < sj; ++j)
	for (int k = 0; k < sk; ++k) {
	  int ri= k + j * c1 + i * c2;
	  int ii=i+j*si+k*si*sj;
	  // cout << ii << " " << ri << endl;
	  DI[ii]=Reffy[ri];
	  if(DI[ii]>maxInt)
	    maxInt=DI[ii];
	}
    
    
    
    //      for (int i = 0; i < si; ++i)
    //      for (int j = 0; j < sj; ++j)
    //	for (int k = 0; k < sk; ++k)
    //	{
    //	  int ii=i+j*si+k*si*sj;
    //	  DI[ii]=pipefile.GetVal3D(i,j,k);
    //	  if(((100*ii)/(si*sj*sk))
    //	  if(DI[ii]>maxInt)
    //	    maxInt=DI[ii];
    //	}
    delete [] Reffy;
    InitBlur();
    
    cout << "Maximum intensity: " << maxInt << endl;
  }

  void slice3D::InitBlurPeakRestricted()
  {
    SetIndexRestricted();
    memset(DB, 0, SIZEMEM*size);

    // Seed only the user-supplied locations.  Multiple restrained peaks that
    // alias to the same voxel intentionally share that voxel rather than
    // creating a new unconstrained location.
    for(int p=0; p<peaks; ++p)
      {
        const int ii=peakList[p].indexI
                   + peakList[p].indexJ*si
                   + peakList[p].indexK*si*sj;
        DB[ii]=DI[ii];
      }
  }

  void slice3D::InitBlur()
  {
    //int blurs=0;
    for(int i=0;i<si;++i) //set all frequencies to nonzero
      for(int j=0;j<sj;++j) //set all frequencies to nonzero
	for(int k=0;k<sk;++k) //set all frequencies to nonzero
	  {
	    int ii=i+j*si+k*si*sj;
	    if(fabs(DI[ii])>noiseVal)
	      {
	      DB[ii]=DI[ii];
	      //blurs++;
	      }
	    else
	      DB[ii]=0;
	  }
    // cout << "StartingDeltas: " << blurs << endl;
  }


  void slice3D::MakeSquare()
  {

    //cout << "declaring main arrays" << endl;
    //DI=new double[size];
#ifdef DOUBLE3D
    DB=new double[size];
    DS=new double[size];
    P=new fftw_complex[size2];
    C=new fftw_complex[size2];
#else
    DB=new float[size];
    DS=new float[size];
     P=new fftwf_complex[size2];
    C=new fftwf_complex[size2];
#endif

    //cout << "setting values " << endl;
    memset(DS, 0, SIZEMEM*size);
    //memset(DB, 0, sizeof(double)*size);

    //DP=new double[size];
    //memset(DP, 0, sizeof(double)*size);


    //memset(DI, 0, sizeof(double)*size);
    //set memory for peak FT
   
    cout << "planning 3D FT" << endl;

#ifdef DOUBLE3D
    p1   = fftw_plan_dft_r2c_3d(sk,sj,si,&DB[0],&C[0],FFTW_ESTIMATE); //put DB FT into C
    //p2   = fftw_plan_dft_r2c_3d(sk,sj,si,&DP[0],&P[0],FFTW_ESTIMATE); //put DB FT into C
    pinv = fftw_plan_dft_c2r_3d(sk,sj,si,&C[0],&DS[0],FFTW_ESTIMATE); //put C FT into DS
    //p1   = fftwf_plan_dft_r2c_3d(sk,sj,si,&DB[0],&C[0],FFTW_MEASURE); //put DB FT into C
    //pinv = fftwf_plan_dft_c2r_3d(sk,sj,si,&C[0],&DS[0],FFTW_MEASURE); //put C FT into DS
#else
    p1   = fftwf_plan_dft_r2c_3d(sk,sj,si,&DB[0],&C[0],FFTW_ESTIMATE|FFTW_DESTROY_INPUT); //put DB FT into C
    //p2   = fftw_plan_dft_r2c_3d(sk,sj,si,&DP[0],&P[0],FFTW_ESTIMATE); //put DB FT into C
    pinv = fftwf_plan_dft_c2r_3d(sk,sj,si,&C[0],&DS[0],FFTW_ESTIMATE|FFTW_DESTROY_INPUT); //put C FT into DS
    //p1   = fftwf_plan_dft_r2c_3d(sk,sj,si,&DB[0],&C[0],FFTW_MEASURE); //put DB FT into C
    //pinv = fftwf_plan_dft_c2r_3d(sk,sj,si,&C[0],&DS[0],FFTW_MEASURE); //put C FT into DS

#endif


    cout << "plans complete" << endl;

  }


#ifdef DOUBLE3D
  int slice3D::DoIndex(double ref,double *vals,int ii,int p)
#else
  int slice3D::DoIndex(double ref,float *vals,int ii,int p)
#endif
  {
    //Is reference within the chemical shift range?
    double maxy=max(vals[0],vals[ii-1]);
    double miny=min(vals[0],vals[ii-1]);

    double incr=fabs(maxy-miny);
    double dd=fabs(vals[1]-vals[0]);


    //alias if required
    int cnt=0;
    while(ref>maxy)
      {
//          cout << ref << " " << maxy << " " << miny << endl;
	ref=ref-incr-dd;
//	cout << ref << endl;
//        cout << ref << " " << maxy << " " << miny << endl;

	if(cnt==100)
	  {
	    cout << "too much folding. check the peaklist" << endl;
	    exit(100);
	  }
	cnt++;
      }
    cnt=0;
    while(ref<miny)
      {
      ref=ref+incr+dd;
      if(cnt==100)
	{
	  cout << "to much folding. check the peaklist" << endl;
	  exit(100);
	}
      cnt++;
      }
//          exit(100);


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

  // Map a true full 3D peak list (f1,f2,f3) to the internal i,j,k axes.
  // Full-list f1 is k (fastest), f2 is j, and f3 is i (slowest).
  void slice3D::SetIndexRestricted()
  {
    cout << "Indexing " << peaks << " restrained 3D peaks" << endl;
    for(int p=0; p<peaks; ++p)
      {
        peakList[p].indexI=DoIndex(peakList[p].z,ivals,si,p);
        peakList[p].indexJ=DoIndex(peakList[p].y,jvals,sj,p);
        peakList[p].indexK=DoIndex(peakList[p].x,kvals,sk,p);
      }
  }

  //get integer index of closest point in ppm
  void slice3D::SetIndex()
  {
    cout << "Size of peaks: " << peaks << endl;
    for(int p=0;p<peaks;p++) //For each peak
      {

	peakList[p].indexI=DoIndex(peakList[p].y,ivals,si,p);
	peakList[p].indexJ=DoIndex(peakList[p].y,jvals,sj,p);
	peakList[p].indexK=DoIndex(peakList[p].x,kvals,sk,p);

	/*
	  double pk=peakList[p].y;//get carbon reference
	  peakList[p].index=0; //search to find nearest index
	  double minV=fabs(jvals[0]-pk);
	  for(int j=1;j<sj;++j) //loop over frequencencies of slice.
	    if(fabs(jvals[j]-pk)<minV)
	      { //found. store in array.
		minV=fabs(jvals[j]-pk);
		peakList[p].index=j;
	      }
	  //cout << peakList[p].name << endl;
	  //cout << "C " << jvals[peakList[p].index] << " " << pk << " " << peakList[p].index << endl;
	}
	{
	  double pk=peakList[p].x;//get proton reference
	  peakList[p].indexH=0; //search to find nearest index
	  double minV=fabs(kvals[0]-pk);
	  for(int k=1;k<sk;++k) //loop over frequencencies of slice.
	    if(fabs(kvals[k]-pk)<minV)
	      { //found. store in array.
		minV=fabs(kvals[k]-pk);
		peakList[p].indexH=k;
	      }
	  //cout << "H " << kvals[peakList[p].indexH] << " " << pk << " " << peakList[p].indexH << endl;
	}
	{




	  double pk=peakList[p].y;//get proton reference
	  peakList[p].indexI=0; //search to find nearest index
	  double minV=fabs(ivals[0]-pk);
	  for(int i=1;k<si;++i) //loop over frequencencies of slice.
	    if(fabs(ivals[i]-pk)<minV)
	      { //found. store in array.
		minV=fabs(ivals[i]-pk);
		peakList[p].indexI=k;
	      }
	  //cout << "H " << kvals[peakList[p].indexH] << " " << pk << " " << peakList[p].indexH << endl;
	  }*/

      }

    //check to make sure indexing is okay
    //for(int i=0;i<line;i++)
    //    cout << indexC[i] << " " << dataC[indexL[i]] << endl;
    return;
  }

#ifdef DOUBLE1D
  void slice3D::SetBlur(int p,double *ref)
#else
  void slice3D::SetBlur(int p,float *ref)
#endif

  {
    //inst.SetBlur(i,sliceLib1D[i].DB); //set blur matrix
    refj=peakList[p].indexJ; //index of carbon
    refk=peakList[p].indexK; //index of proton
    for(int i=0;i<si;++i)
      if(fabs(ref[i]) > 0.0)
	DB[i+refj*si+refk*si*sj]+=ref[i];
  }


  double slice3D::ApplyIter()
  {
    double tack=0.0;

    if(SPARSE)
      {
        // CalcSpec() has just rebuilt sparseDB from the active DB entries.
        for(size_t p=0; p<sparseDB.size(); ++p)
          {
            const int ii=sparseDB[p].ii;
            DB[ii]=DB[ii]*fabs(DI[ii]/DS[ii]);
            tack+=double(DB[ii]);
          }
      }
    else
      {
        for(int i=0;i<si;++i)
          for(int j=0;j<sj;++j)
            for(int k=0;k<sk;++k)
              {
                const int ii=i+j*si+k*si*sj;
                if(fabs(DB[ii])>0.0)
                  {
                    DB[ii]=DB[ii]*fabs(DI[ii]/DS[ii]);
                    tack+=double(DB[ii]);
                  }
                else
                  DB[ii]=0;
              }
      }

    return tack;
  }

  //function to return a gaussian at a specific point
//  double Gaus(double x,double x0,double sig){
//    return exp(-pow(x-x0,2.)/(2.*sig*sig));}
//  //function to return a gaussian at a specific point
//  double Lorentz(double x,double x0,double sig){
//    return pow(sig/2,2)/(pow(x-x0,2)+pow(sig/2,2));}
//  //function to return a 2D peak value at a given point
//  double Peak(double x1,double x2,double sig)
//  {
//    switch(peaky){
//    case 0:
//      return Gaus(x1,x2,sig);
//      break;
//    case 1:
//      return Lorentz(x1,x2,sig);
//      break;}
//  }

  /*  void SetPeak(slice3D inst)
  {
  int si2=si/2+1;
    for(int i=0;i<si2;i++) //goes with dim2 (C)
      for(int j=0;j<sj;j++) //goes with dim1 (H)
	for(int k=0;k<sk;k++) //goes with dim1 (H)
	  {
	    int ii=i+j*si2+k*si2*sj;
	    P[ii][0]=inst.P[ii][0];
	    P[ii][1]=inst.P[ii][1];
	  }

    //unneccessary
    for(int i=0;i<si;i++) //goes with dim2 (C)
      for(int j=0;j<sj;j++) //goes with dim1 (H)
	for(int k=0;k<sk;k++) //goes with dim1 (H)
	  {
	    int ii=i+j*si+k*si*sj;
	    DP[ii]=inst.DP[ii];
	  }
	  }*/


/*
void slice3D::GetPeak()
{
  cout << "Setting 3D peak: " << sig1 << " " << sig2 << " " << sig3 << endl;

  peakI.assign(si, 0.0);
  peakJ.assign(sj, 0.0);
  peakK.assign(sk, 0.0);

#ifdef DOUBLE3D
  double di = ivals[1] - ivals[0];
  double dj = jvals[1] - jvals[0];
  double dk = kvals[1] - kvals[0];
  double *DBtmp = new double[size];
#else
  float di = ivals[1] - ivals[0];
  float dj = jvals[1] - jvals[0];
  float dk = kvals[1] - kvals[0];
  float *DBtmp = new float[size];
#endif

  // build reusable 1D wrapped peak tables
  for (int i = 0; i < si; ++i)
    peakI[i] = Peak(di * i, 0, sig1, lor1, voigt1) +
               Peak(di * i, di * si, sig1, lor1, voigt1);

  for (int j = 0; j < sj; ++j)
    peakJ[j] = Peak(dj * j, 0, sig2, lor2, voigt2) +
               Peak(dj * j, dj * sj, sig2, lor2, voigt2);

  for (int k = 0; k < sk; ++k)
    peakK[k] = Peak(dk * k, 0, sig3, lor3, voigt3) +
               Peak(dk * k, dk * sk, sig3, lor3, voigt3);

  if (sparseSpec)
  {
    // sparse mode only needs the 1D tables
    delete [] DBtmp;
    return;
  }

  // full FFT path stays available
  memcpy(DBtmp, DB, size * SIZEMEM);

  for (int i = 0; i < si; ++i)
    for (int j = 0; j < sj; ++j)
      for (int k = 0; k < sk; ++k)
        DB[i + j * si + k * si * sj] = peakI[i] * peakJ[j] * peakK[k];

#ifdef DOUBLE3D
  fftw_execute(p1);
#else
  fftwf_execute(p1);
#endif

  for (int ii = 0; ii < size2; ++ii)
  {
    P[ii][0] = C[ii][0];
    P[ii][1] = C[ii][1];
  }

  memcpy(DB, DBtmp, size * SIZEMEM);
  delete [] DBtmp;
}
*/



  void slice3D::GetPeak()
  {
    cout << "Setting 3D peak: " << sig1 << " " << sig2 << " " << sig3 << endl;
#ifdef DOUBLE3D
    double di=ivals[1]-ivals[0]; //carbon delta (d1)
    double dj=jvals[1]-jvals[0]; //proton delta (d2)
    double dk=kvals[1]-kvals[0]; //proton delta (d2)
    pki = new double[si];
    pkj = new double[sj];
    pkk = new double[sk];
    double *DBtmp;
    DBtmp=new double[size];
#else
   float di=ivals[1]-ivals[0]; //carbon delta (d1)
    float dj=jvals[1]-jvals[0]; //proton delta (d2)
    float dk=kvals[1]-kvals[0]; //proton delta (d2)
    pki = new float[si];
    pkj = new float[sj];
    pkk = new float[sk];
    float *DBtmp;
    DBtmp=new float[size];
#endif
    // cout << '1' << endl;

    memcpy(DBtmp,DB,size*SIZEMEM);
    cout << sig1 << " " << sig2 << " " << sig3 << " " << lor1 << " " << lor2 << " " << lor3 << endl;
    cout << voigt1 << " " << voigt2 << " " << voigt3 << endl;

    for(int i=0;i<si;++i)
	pki[i]=Peak(di*i,0,sig1,lor1,voigt1)+Peak(di*i,di*si,sig1, lor1, voigt1);
    for(int j=0;j<sj;++j)
	pkj[j]=Peak(dj*j,0,sig2,lor2, voigt2)+Peak(dj*j,dj*sj,sig2,lor2, voigt2);
    for(int k=0;k<sk;++k)
	pkk[k]=Peak(dk*k,0,sig3,lor3, voigt3)+Peak(dk*k,dk*sk,sig3,lor3, voigt3);


    for(int i=0;i<si;++i) //goes with dim2 (C)
      for(int j=0;j<sj;++j) //goes with dim1 (H)
	for(int k=0;k<sk;++k) //goes with dim1 (H)
	  DB[i+j*si+k*si*sj]=pki[i]*pkj[j]*pkk[k]; //store peak in complete peak list

#ifdef DOUBLE3D
    fftw_execute(p1); //fourier transform DP into P
#else
    fftwf_execute(p1); //fourier transform DP into P
#endif
    for(int ii=0;ii<size2;++ii)
      {
	//int ii=i+j*si2+k*si2*sj;
	P[ii][0]=C[ii][0];
	P[ii][1]=C[ii][1];
      }
    memcpy(DB,DBtmp,size*SIZEMEM);
    delete [] DBtmp;
  }

  void slice3D::BlankDB()
  {
    for(int p=0;p<peaks;++p)
      {
	int jslice=peakList[p].indexJ;
	int kslice=peakList[p].indexK;
	int loc=jslice*si+kslice*si*sj;
	memset(DB+loc, 0, SIZEMEM*si);
      }
  }

  void slice3D::BlankDBfull()
  {
    memset(DB, 0, SIZEMEM*size);
  }

static inline int wrapDist(int a, int b, int n)
{
  int d = abs(a - b);
  return (d <= n / 2) ? d : (n - d);
}


  void slice3D::CalcSpec()
  {
    if(SPARSE)
      {
	//const double scale = 1.0;// / double(size);
	// refresh sparse source list from current DB values
	BuildSparseDB(0.0); //take all values in DB bigger than zero.
	//could add extra checks to stop this when the peak shape function gets small.
	// Only clear the entries we are going to write, or zero all DS if you
	// want a fully consistent dense DS array.
	// For maximum safety:
	memset(DS, 0, SIZEMEM * size);
	const size_t n = sparseDB.size();
	if (n == 0) return;

	#pragma omp parallel for schedule(static)
	for (long long out = 0; out < (long long)n; ++out)
	  {
	    const SparsePt3D &o = sparseDB[out];
	    double acc = 0.0;
	    
	    for (size_t src = 0; src < n; ++src)
	      {
		const SparsePt3D &s = sparseDB[src];
		
		int di = wrapDist(o.i, s.i, si);
		int dj = wrapDist(o.j, s.j, sj);
		int dk = wrapDist(o.k, s.k, sk);
		
		acc += s.val * pki[di] * pkj[dj] * pkk[dk];
	      }
	    
	    DS[o.ii] = acc ;
	  }
      }
    else
      {
	
	double scale = 1.0 / (size);
	
	
	//cout << "transforming 3D blur...." << endl;
#ifdef DOUBLE3D
	fftw_execute(p1); //fourier transform blur
#else
	fftwf_execute(p1); //fourier transform blur
#endif
	//int si2=si/2+1;
	//cout << "done blur FT " << endl;
	//for (int i = 0; i < si2; ++i)
	//for (int j = 0; j < sj; ++j)
	//	for (int k = 0; k < sk; ++k)
	for(int ii=0;ii<size2;++ii)
	  {
	    double a=C[ii][0]; //blur FT real
	    double b=C[ii][1]; //blur FT imag
	    double c=P[ii][0]; //peak FT real
	    double d=P[ii][1]; //peak FT imag
	    C[ii][0]=(a*c-b*d)*scale;
	    C[ii][1]=(b*c+a*d)*scale;
	  }
	//cout << "doing inverse... " << endl;
#ifdef DOUBLE3D
	fftw_execute(pinv);
#else
	fftwf_execute(pinv);
#endif
      }
  }



void slice3D::BuildSparseDB(double cutoff)
{
  sparseDB.clear();

  // A restrained 3D reconstruction has one allowed delta voxel per input
  // peak.  Build the sparse support from those exact voxels, once each.
  //
  // The historical implementation scanned the complete i-line for every
  // peak's (j,k) pair.  That is appropriate to the old 2D-reference/bore
  // representation, but is wrong for a full 3D restrained list: several
  // peaks commonly share the same (j,k) parent position, so the same active
  // voxels were inserted repeatedly.  CalcSpec then counted those sources
  // multiple times on every iteration, driving all reconstructed amplitudes
  // rapidly to zero.
  if(peaks > 0)
    {
      sparseDB.reserve(peaks);
      vector<int> seen;
      seen.reserve(peaks);

      for(int p=0; p<peaks; ++p)
        {
          const int i=peakList[p].indexI;
          const int j=peakList[p].indexJ;
          const int k=peakList[p].indexK;
          const int ii=i+j*si+k*si*sj;
          const double v=DB[ii];

          if(fabs(v) <= cutoff)
            continue;

          bool duplicate=false;
          for(size_t q=0; q<seen.size(); ++q)
            if(seen[q] == ii)
              {
                duplicate=true;
                break;
              }
          if(duplicate)
            continue;

          seen.push_back(ii);
          SparsePt3D pt;
          pt.i=i;
          pt.j=j;
          pt.k=k;
          pt.ii=ii;
          pt.val=v;
          sparseDB.push_back(pt);
        }
      return;
    }

  // Defensive fallback for a 3D calculation without a peak list.
  for(int ii=0; ii<size; ++ii)
    {
      const double v=DB[ii];
      if(fabs(v) <= cutoff)
        continue;
      const int k=ii/(si*sj);
      const int rem=ii-k*si*sj;
      const int j=rem/si;
      const int i=rem-j*si;
      SparsePt3D pt;
      pt.i=i; pt.j=j; pt.k=k; pt.ii=ii; pt.val=v;
      sparseDB.push_back(pt);
    }
}
  void slice3D::PrintProj()
  {
    FILE *fp;
    {
      cout << "Making simulated projections" << endl;
      string line ="out/slice2d/xy.dat.out";
      fp=fopen(line.c_str(),"w");
      for (int i=0;i<si;++i)
	{
	  for (int j=0;j<sj;++j)
	    {
	      double sum=0.0;
	      for (int k=0;k<sk;++k)
		{
		  int ii=i+j*si+k*si*sj;
		  if(fabs(DI[ii])>noiseVal)
		    sum+=DI[ii];
		}
	      fprintf(fp,"%f\t%f\t%e\n",ivals[i],jvals[j],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }

    {
      string line ="out/slice2d/xz.dat.out";
      fp=fopen(line.c_str(),"w");
      for (int i=0;i<si;++i)
	{
	  for (int k=0;k<sk;++k)
	    {
	      double sum=0.0;
	      for (int j=0;j<sj;++j)
		{
		  int ii=i+j*si+k*si*sj;
		  if(fabs(DI[ii])>noiseVal)
		    sum+=DI[ii];
		}
	      fprintf(fp,"%f\t%f\t%e\n",ivals[i],kvals[k],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }

    {
      string line ="out/slice2d/yz.dat.out";
      fp=fopen(line.c_str(),"w");
      for (int j=0;j<sj;++j)

	{
	  for (int k=0;k<sk;++k)
	    {
	      double sum=0.0;
	      for (int i=0;i<si;++i)
		{
		  int ii=i+j*si+k*si*sj;
		  if(fabs(DI[ii])>noiseVal)
		    sum+=DI[ii];
		}
	      //double sum = accumulate(&DI[j*si+k*si*sj], &DI[j*si+k*si*sj]+si, 0,[](int a, int b){return b>noiseVal? a+b: a;});
	      //if(fabs(sum)<=noiseVal)
	      //	sum=0;
	      fprintf(fp,"%f\t%f\t%e\n",kvals[k],jvals[j],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }
  }


  void slice3D::PrintSpec()
  {
    /*cout << "printing spectrum " << endl;

    for(int p=0;p<peaks;++p) //for each slice
      {
	string line ="out/slice2d/"+peakList[p].name+".dat.3d.decon";
	//cout << "opening: " << line<<endl;
	FILE *fp;
	fp=fopen(line.c_str(),"w");

	int jslice=peakList[p].indexJ; //get index of second dimension

	double wid = 0.4/2.;//need to get this from input

	for (int k=0;k<sk;++k) //proton (f3)
	  {
	    if( kvals[k]>(peakList[p].x-wid) && kvals[k]<(peakList[p].x+wid))
	      {
		for(int i=0;i<si;++i) //core dimension
		  {
		    int ii=i+jslice*si+k*si*sj;
		    if(pipey==0)
		      fprintf(fp,"%f\t%f\t%f\t%e\t%e\t%e\n",ivals[i],jvals[jslice],kvals[k],DB[ii],DS[ii],DP[ii]);
		    else
		      fprintf(fp,"%f\t%f\t%f\t%e\t%e\t%e\t%e\n",ivals[i],jvals[jslice],kvals[k],DB[ii],DS[ii],DP[ii],DI[ii]);
		  }
	      }
	    fprintf(fp,"\n");
	  }
	fclose(fp);
      }
    */
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
	      //double sum2=0.0;
	      for (int k=0;k<sk;++k)
		{
		  int ii=i+j*si+k*si*sj;
		  sum+=DS[ii];
		  //sum2+=DP[ii];
		}
	      if(fabs(sum)<=noiseVal)
		sum=0;
	      fprintf(fp,"%f\t%f\t%e\n",jvals[j],ivals[i],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }

    {
      string line ="out/xz.decon";
      fp=fopen(line.c_str(),"w");
      for (int i=0;i<si;++i)
	{
	  for (int k=0;k<sk;++k)
	    {
	      double sum=0.0;
	      //double sum2=0.0;
	      for (int j=0;j<sj;++j)
		{
		  int ii=i+j*si+k*si*sj;
		  sum+=DS[ii];
		  //sum2+=DP[ii];
		}
	      if(fabs(sum)<=noiseVal)
		sum=0;
	      fprintf(fp,"%f\t%f\t%e\n",ivals[i],kvals[k],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }

    {
      string line ="out/yz.decon";
      fp=fopen(line.c_str(),"w");
      for (int k=0;k<sk;++k)
	{
	  for (int j=0;j<sj;++j)
	    {
	      double sum=0.0;
	      //double sum2=0.0;
	      for (int i=0;i<si;++i)
		{
		  int ii=i+j*si+k*si*sj;
		  sum+=DS[ii];
		  //sum2+=DP[ii];
		}
	      if(fabs(sum)<=noiseVal)
		sum=0;
	      fprintf(fp,"%f\t%f\t%e\n",kvals[k],jvals[j],sum);
	    }
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }
  }//end printspec

  double slice3D::GetChi2()
  {
    double chi2=0.0;
    for(int i=0;i<size;++i)
      chi2+=pow((DS[i]-DI[i]),2.);
    return chi2;
  }

    void slice3D::Squash(){
#ifdef DOUBLE3D
        int window_limits_i = int(double(squash_window_i)/2.0); //turn into integers
        int window_limits_j = int(double(squash_window_j)/2.0); //turn into integers
        int window_limits_k = int(double(squash_window_k)/2.0);
        vector<double> max_values; //log value of maxima

#else
       int window_limits_i = int(float(squash_window_i)/2.0); //turn into integers
        int window_limits_j = int(float(squash_window_j)/2.0); //turn into integers
        int window_limits_k = int(float(squash_window_k)/2.0);
        vector<float> max_values; //log value of maxima
#endif

        cout << "squashing " << window_limits_i << endl;

        vector<vector<int> > maximae;      //log location of maxima

        int max_counter = 0;
        bool maximum = true;
        for(int i=0;i<si;++i)
            for(int j=0;j<sj;++j)
                for (int k=0; k<sk; ++k)
                {
                    int ii=i+si*j+k*si*sj;
                    if (fabs(DB[ii]) > 0.0) //if there is intensity here...
                    {
                        maximum = true;

                        //               if (i != 0 && j != 0  && k != 0 && i != si - 1 && j != sj - 1  && k != sk-1)
                        //              {
                        int bottom_edgei = max(i-1, 0);
                        int bottom_edgej = max(j-1, 0);
                        int bottom_edgek = max(k-1, 0);


                        int top_edgei = min(i+2, si);
                        int top_edgej = min(j+2, sj);
                        int top_edgek = min(k+2, sk);
                        for (int i2 = bottom_edgei; i2 < top_edgei; ++i2)   //walk around the peak
                            for (int j2 = bottom_edgej; j2 < top_edgej; ++j2)  //walk around the peak
                                for (int k2= bottom_edgek; k2 < top_edgek; ++k2)
                                {
                                    int ii2 = i2 + si * j2 + k2*si*sj;

                                    if (fabs(DB[ii2]) > fabs(DB[ii]))
                                        maximum = false;
                                }
                        //             } else {
                        //                 maximum = false;
                        //             }
                        if (maximum == true) {
                            vector<int> coords;
                            coords.push_back(i);
                            coords.push_back(j);
                            coords.push_back(k);
                            maximae.push_back(coords); //save array value
                            max_values.push_back(DB[ii]); //save intensity value
                        }
                    }
                }

        max_counter=maximae.size();  //here's how many maxima we're dealing with


        vector<vector<int> > maximae_sorted; //vector, to be sorted by max intensity
        cout << max_counter << " : max counter" << endl;
        while(maximae_sorted.size()!=max_counter) //get indicies of sorted intensities
        {
            //maxval: *max_element(balls.begin(),balls.end())
            //maxind: max_element(balls.begin(),balls.end())-balls.begin()
            int ind= min_element(max_values.begin(),max_values.end())-max_values.begin(); //index of max
            maximae_sorted.push_back(maximae[ind] ); //save reference value
            max_values.erase(max_values.begin()+ind); //remove entry from max_values
            maximae.erase(maximae.begin()+ind);  //remove entry from maximae
        }
        int charlie123 = 0;
        for (int i_count=0; i_count < max_counter; ++i_count){
            //int ii = maximae[i_count];  //I think we were previously taking unsorted intensities?
            vector<int> max_values = maximae_sorted[i_count]; //take intensities in order
            int i = max_values[0];
            int j= max_values[1];
            int k = max_values[2];
            int ii = i+si*j+k*si*sj;

            double sumo = 0.0;  //soak up intensity

            double maxy=DB[ii]; //get local maximum

            i = max(i, window_limits_i);
            i = min(i, si-window_limits_i);
            j = max(j, window_limits_j);
            j = min(j, sj-window_limits_j);
            k = max(k, window_limits_k);
            k = min(k, sk-window_limits_k);

            // float total = 0;
            // for(int i2=0;i2<si;++i2)
            //     for(int j2=0;j2<sj;++j2)
            //         for (int k2=0; k2<sk; ++k2)
            //         {
            //             int ii2 = i2 + si * j2 + k2 * si * sj;
            //             if (DB[ii2] > 0.0)
            //               total += 1;
            //         }
            // cout << "total: " << total << endl;

            for (int i2 = i - window_limits_i; i2 < i + window_limits_i + 1; ++i2)
                for (int j2 = j - window_limits_j; j2 < j + window_limits_j + 1; ++j2)
                    for (int k2 = k - window_limits_k; k2 < k + window_limits_k + 1; ++k2)
                    {
                        int ii2 = i2 + si * j2 + k2 * si * sj;
                        sumo += DB[ii2];

                        if (fabs(DB[ii2]) > fabs(maxy)) //make sure we track local maximum
                        { //if we have a local maximum...

                            ii = ii2; //move current maxima
                            maxy=DB[ii2]; //update local maximum
                            charlie123 += 1;
                            cout << charlie123 << endl;
                        }
                        DB[ii2]=0.0; //set the soaked up signal to zero

                    }
            DB[ii] = sumo; //dump all signal into local maximum


        }
        return;


    }

  void slice3D::Cull(float frac)
  {
    for(int i=0;i<size;++i)
      if(fabs(DB[i])<noiseVal*frac)
	 DB[i]=0.0;
  }

  int slice3D::CountElements()
  {
    int cnt=0;
    for(int i=0;i<size;++i)
      if(fabs(DB[i])> 0.0)
	cnt+=1;
    return cnt;
  }

void slice3D::AddIfNew(vector<vector <int> > &peakSum,int j,int k)
  {
    int tag=0;
    for(int ii=0;ii<peakSum.size();++ii)
      {
	if(peakSum[ii][0]==j)
	  if(peakSum[ii][1]==k)
	    tag=1;
      }
    if(tag==0)
      {
	vector<int> newy;
	newy.push_back(j);
	newy.push_back(k);
	peakSum.push_back(newy);
      }
    return;
  }

  int slice3D::correlateRestricted(FILE *out_pt)
  {
    // Preserve the authoritative input names and coordinates.  The intensity
    // remains column dim+1, matching the GUI Full nD list reader.
    for(int p=0; p<peaks; ++p)
      {
        const int ii=peakList[p].indexI
                   + peakList[p].indexJ*si
                   + peakList[p].indexK*si*sj;
        fprintf(out_pt,"%s\t%f\t%f\t%f\t%e\t%e\t%e\n",
                peakList[p].name.c_str(),
                peakList[p].x, peakList[p].y, peakList[p].z,
                DB[ii], DS[ii], DI[ii]);
      }
    return peaks;
  }

  int slice3D::correlate(string projListFile,FILE *out_pt)
  {
    FILE *out_pk;
    // The projected 2D peak list belongs alongside the main 3D peak list.
    // Its name is supplied by Protocol3D via decon::correlate().
    string inpk=projListFile;

    out_pk=fopen(inpk.c_str(),"w");

    vector<vector<int> > peakSum;
    vector<int > peakLog;
    for(int i=0;i<si;i++)
      for(int j=0;j<sj;++j)
	for(int k=0;k<sk;++k)
	  {
	    int ii=i+j*si+k*si*sj;
	    if(fabs(DB[ii]) > 0.0)
	      { //collect only unique j,k combinations for 2D peaklist
		AddIfNew(peakSum,j,k);
		peakLog.push_back(peakSum.size());
		//cout << peakSum.size() << endl;
	      }
	  }
    cout << "2D PEAKS DETECTED: " << peakSum.size() << endl;

    for(int pk=0;pk<peakSum.size();++pk)
      {
	int j=peakSum[pk][0];
	int k=peakSum[pk][1];
	//cout << pk+1 << " " << j << " " << k << endl;
	fprintf(out_pk,"%i\t%f\t%f\n",pk+1,jvals[j],kvals[k]);
      }
    fclose(out_pk);

    //second - print out 3D peaks, but indexed by first peak
    int cnt=0;
    for(int i=0;i<si;i++)
      for(int j=0;j<sj;++j)
	for(int k=0;k<sk;++k)
	  {
	    int ii=i+j*si+k*si*sj;
	    if(fabs(DB[ii]) > 0.0)
	      {
		fprintf(out_pt,"%i\t%f\t%f\t%f\t%e\t%e\t%e\n",peakLog[cnt],kvals[k],jvals[j],ivals[i],DB[ii],DS[ii],DI[ii]);
		cnt++;

	      }

	  }
    return cnt;
  }



#endif
