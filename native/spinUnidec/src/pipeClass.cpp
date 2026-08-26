/**************************************************/
/* nmrPipe class         */
/* A.Baldwin     */
/* 19th March 2018*/
/**************************************************/

#ifndef PIPE_CPP
#define PIPE_CPP


#include "deconMain.hpp"
#include "pipeClass.hpp"


int pipe::LetterTest(char test)
  {
    string newy;
    newy+=test;
    string t1="abcdefghijklmnopqrstuvwxyz";
    string t2="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    //cout << test << endl;
    //cout << t1.compare(newy) << endl;
    //cout << t2.compare(newy) << endl;

    if(t1.compare(newy)<0)
      return 1;
    if(t2.compare(newy)<0)
      return 1;

    return 0;
  }

string pipe::GetLabel(int posy)
  {
    string output;
    int pos;
    const int sz=sizeof(char);
    bitset<sz> bot;
    //string t1="abcdefghijklmnopqrstuvwxyz";
    //for (int i=0;i<t1.length();++i)
    //  LetterTest(t1[i]);
    for(int i=0;i<4;++i)
      {
	pos=posy*sizeof(float)+i*sz;
	fseek(A,pos,0);
	fread(&bot,sz,1,A);
	char ne=char (bot.to_ulong());
  // cout << "bot "<< ne << endl;
        // if(LetterTest(ne)==1) //verify if this is a character a-Z A-Z or not
	  output+= ne;
      }
   //cout << output << " " << output.length() << endl;

   return output;
  }


void pipe::readHeader()
  {

    cout << iName << endl;
    A   = fopen(iName.c_str(),"rb");
    fread(prehead,sizeof(float),512,A);

    // print all header arguments to the screen
  //   for(int i=0;i<512;++i)
  //    {
	// //if(abs(prehead[i]-26.871)<1)
	// //if(prehead[i]==2048)
  //     std::cout << i << " " << prehead[i] << std::endl;
  //     }

    dim=prehead[9]; //FDDIMCOUNT

    size =new int[4];
    dimord=new int[4];
    sw=new float[4];
    car=new float[4];
    obs=new float[4];
    orig=new float[4];
    x1=new int[4];
    xn=new int[4];
    center=new int[4];

    /*
    for(int i=0;i<512;++i)
      {
      if(int(prehead[i])==100)
	cout << prehead[i] << " " << i << endl;
      if(int(prehead[i])==200)
	cout << prehead[i] << " " << i << endl;
      if(int(prehead[i])==684)
	cout << prehead[i] << " " << i << endl;
	}*/

    //size[FDDIMORDERx-1]=int(FD{dim}SIZE
    //PRETTY STRANGE: either index is wrong or is not being used systematically
    dimord[0]=int(prehead[24])-1;
    dimord[1]=int(prehead[25])-1;
    dimord[2]=int(prehead[26])-1;
    dimord[3]=int(prehead[27])-1;


    size[dimord[0]]=int(prehead[99]); //FDSIZE; //99 //128       //FDF1FTSIZE // 98  //'FDF1TDSIZE': '387'
    size[dimord[1]]=int(prehead[219]); //FDF3FTSIZE//200  //64   //FDF2FTSIZE // 96   'FDF2TDSIZE': '386',
    size[dimord[2]]=int(prehead[15]);//FDF3SIZE //15 128       //FDF3FTSIZE // 128  'FDF3TDSIZE': '388',
    size[dimord[3]]=int(prehead[32]);//FDF4SIZE //32  236       //FDF4FTSIZE // 201 'FDF4TDSIZE': '389'
    // 'FDDIMORDER4': '27', 'FDDIMORDER3': '26',
    //'FDDIMORDER2': '25', 'FDDIMORDER1': '24',

    vector <string> labTmp;
    labTmp.push_back(GetLabel(16)); //FDF1LABEL
    labTmp.push_back(GetLabel(18));//FDF2LABEL
    labTmp.push_back(GetLabel(20));//FDF3LABEL
    labTmp.push_back(GetLabel(22));//FDF4LABEL
    for(int i=0;i<4;++i)
      labels.push_back(labTmp[dimord[i]]);

     cout << "dimorder:" << endl;
    for(int i=0;i<4;++i)
      cout << dimord[i]+1 << " " << labels[i] << endl;
    cout << endl;


    obs[0]=prehead[218]; //FDF1OBS
    obs[1]=prehead[119]; //FDF2OBS
    obs[2]=prehead[10];  //FDF3OBS
    obs[3]=prehead[28]; //FDF4OBS

    sw[0]=prehead[229]; //FDF1SW
    sw[1]=prehead[100]; //FDF2SW
    sw[2]=prehead[11];  //FDF3SW
    sw[3]=prehead[29]; //FDF4SW

    car[0]=prehead[67]; //FDF1CAR
    car[1]=prehead[66]; //FDF2CAR
    car[2]=prehead[68]; //FDF3CAR
    car[3]=prehead[69]; //FDF4CAR

    orig[0]=prehead[249]; //FDF1ORIG
    orig[1]=prehead[191]; //FDF2ORIG
    orig[2]=prehead[12]; //FDF3ORIG
    orig[3]=prehead[30]; //FDF4ORIG

    x1[0]=prehead[259]; //FDF1X1
    x1[1]=prehead[257]; //FDF2X1
    x1[2]=prehead[261]; //FDF3X1
    x1[3]=prehead[263]; //FDF4X1

    xn[0]=prehead[260]; //FDF1XN
    xn[1]=prehead[258]; //FDF2XN
    xn[2]=prehead[262]; //FDF3XN
    xn[3]=prehead[264]; //FDF4XN

    center[0]=prehead[80]; //FDF1CENTER
    center[1]=prehead[79]; //FDF2CENTER
    center[2]=prehead[81]; //FDF3CENTER
    center[3]=prehead[82]; //FDF4CENTER

    for(int i=0;i<4;++i) //adjust center if needed (if extraction used)
      car[i] = car[i] + sw[i] / size[i] * ( center[i] - 1. - size[i]/2.)/obs[i];


    xvals=new float*[4];
    for(int i=0;i<4;++i)
      xvals[i]=new float[size[i]];

    for(int i=0;i<4;++i) //for each dimension
      for(int j=0;j<size[i];++j) //calculate ppm values
	xvals[i][j]=car[i]+sw[i]/2./obs[i]-(j)*(sw[i]*1.)/size[i]/obs[i];


    //'FDF1FTFLAG': '222',
    //'FDF2FTFLAG': '220'
    //FDF3FTFLAG': '13
    //'FDF4FTFLAG': '31'
    //std::cout << prehead[222] << " "<<  prehead[220] << " " << prehead[13] << " " << prehead[31] << std::endl;

    std::cout << "dimensions: " << dim << std::endl;
    std::cout << "Label   " << "Size   " << "sw      " << "car    " << "obs   " << "min(ppm)   " << "max(ppm)   " << std::endl;


    //    for (int i=0;i<dim;++i) //print value to screen
        for (int i=0;i<4;++i) //print value to screen
	  {
	    //string lab=labels[0];
	    //for(int j=0;j<
      
	    std::cout << dimord[i] << " " << labels[i] << " " << size[i] << " " << sw[i] << " " << car[i] << " " << obs[i] << " " << car[i]+sw[i]/obs[i]/2. << " " << car[i]-sw[i]/obs[i]/2.+sw[i]/obs[i]/(size[i])  << std::endl;
	  }

    std::cout << std::endl;


  }

int pipe::GetIndex(float val,int col)
  {
    //int fold=0;
    //double vol=val;
    while(val>xvals[col][0])
      {
	val-=(sw[col]/obs[col]+fabs(xvals[col][1]-xvals[col][0]));
	//fold=1;
      }
    while(val<xvals[col][size[col]-1])
      {
	val+=(sw[col]/obs[col]+fabs(xvals[col][1]-xvals[col][0]));
	//fold=1;
      }
    //if(fold)
    //  std::cout << "Aliased peak: New: " << val << " Old: " << vol << std::endl;

    double min=abs(xvals[col][0]-val);
    int imin=0;
    for(int i=1;i<size[col];++i)
      {
	double test=abs(xvals[col][i]-val);
	if(test<min)
	  {
	    min=test;
	    imin=i;
	  }
      }
    return imin;
  }

void pipe::SetPos(int &i,int &j,int &k,int &l)
  {
    //int pos=(512+l+k*size[3]+(size[3]*size[2])*(i+j*size[0]))*sizeof(float);
    //int pos=(512+l+k*size[dimord[0]]+(size[dimord[0]]*size[dimord[1]])*(i+j*size[dimord[2]]))*sizeof(float);
    int pos=(512+l+k*size[dimord[0]]+(size[dimord[0]]*size[dimord[1]])*(j+i*size[dimord[2]]))*sizeof(float);

    //int pos=(512+k+j*size[dimord[0]]+i*size[dimord[0]]*size[dimord[1]])*sizeof(float);

    fseek(A,pos,0);
  }

void pipe::SetPosInt(int *loc)
  {
    long int pos=(512+loc[dimord[0]]+loc[dimord[1]]*size[dimord[0]]+(size[dimord[0]]*size[dimord[1]])*(loc[dimord[2]]+loc[dimord[3]]*size[dimord[2]]))*sizeof(float);
    /* if(pos> (size[0]*size[1]*size[2]*size[3]+512)*sizeof(float) )
    //if(cnt*sizeof(float)>lim*sizeof(float))
      {
	for(int i=0;i<4;++i)
	  cout << loc[dimord[i]] << " " << size[dimord[i]] << endl;
	int cnt=loc[dimord[0]]+loc[dimord[1]]*size[dimord[0]]+(size[dimord[0]]*size[dimord[1]])*(loc[dimord[2]]+loc[dimord[3]]*size[dimord[2]]);
	//int lim=size[0]*size[1]*size[2]*size[3]+512;
	//cout << cnt << " " << lim << " " << cnt/(lim*1.) << endl;
      cout << "shit" << endl;
      exit(100);
      }*/
    fseek(A,pos,0);
  }


void pipe::SetPos3D(int &i,int &j,int &k)
  {
    //128 615 44
    //13C Hn N15
    //i j k
    //13C 15N Hn
    //i goes with 0 13C
    //j goes with 2 15N
    //k goes with 1 HN
    //    //int pos=(512+k+j*size[1]+i*size[1]*size[2])*sizeof(float);
    //file goes with

    //i goes with C
    //j goes with N
    //k goes with H


    //int pos=(512+k+j*size[1]+i*size[1]*size[0])*sizeof(float);
    long int pos=(512+k+j*size[dimord[0]]+i*size[dimord[0]]*size[dimord[1]])*sizeof(float);
    //if(pos>size[dimord[0]]*size[dimord[1]]*size[dimord[2]]*sizeof(float)+512*sizeof(float))

    fseek(A,pos,0);
  }

void pipe::SetPos2D(int &i,int &j)
  {
    //128 615 44
    //13C Hn N15
    //i j k
    //13C 15N Hn
    //i goes with 0 13C
    //j goes with 2 15N
    //k goes with 1 HN
    //    //int pos=(512+k+j*size[1]+i*size[1]*size[2])*sizeof(float);
    //file goes with

    //i goes with C
    //j goes with N
    //k goes with H


    //int pos=(512+k+j*size[1]+i*size[1]*size[0])*sizeof(float);
    //cout << i << " " << j << " " << +j+i*size[dimord[0]] << " " << size[dimord[0]]*size[dimord[1]] << endl;
    int pos=(512+j+i*size[dimord[0]])*sizeof(float);
    fseek(A,pos,0);
  }
void pipe::SetPos1D(int &i)
  {
    //128 615 44
    //13C Hn N15
    //i j k
    //13C 15N Hn
    //i goes with 0 13C
    //j goes with 2 15N
    //k goes with 1 HN
    //    //int pos=(512+k+j*size[1]+i*size[1]*size[2])*sizeof(float);
    //file goes with

    //i goes with C
    //j goes with N
    //k goes with H


    //int pos=(512+k+j*size[1]+i*size[1]*size[0])*sizeof(float);
    //cout << i << " " << j << " " << +j+i*size[dimord[0]] << " " << size[dimord[0]]*size[dimord[1]] << endl;
    int pos=(512+i)*sizeof(float);
    fseek(A,pos,0);
  }

  //read a single datapoint specified from an ijkl vector
double pipe::GetData(int *loc)
  {
    //SetPos(indy[0],indy[1],indy[2],indy[3]);
    SetPosInt(loc);
    float Data;
    fread(&Data,sizeof(float),1,A);
    return Data;
  }

  //in 4 dimensions, find the maximum of the auto peak
void pipe::LocalMax(int p)
  {
    //int kk=GetIndex(peakList[p].x,2); //proton
    //int ll=GetIndex(peakList[p].y,3); //carbon
    //int jj=GetIndex(peakList[p].x,1); //proton
    //int ii=GetIndex(peakList[p].y,0); //carbon

    int loc[dim],lac[dim];
    //the GetData function, and this bit are the only dimension dependent parts
    loc[dimord[0]]=GetIndex(peakList[p].y,dimord[0]); //carbon
    loc[dimord[1]]=GetIndex(peakList[p].x,dimord[1]); //proton
    loc[dimord[2]]=GetIndex(peakList[p].y,dimord[2]); //carbon
    loc[dimord[3]]=GetIndex(peakList[p].x,dimord[3]); //proton

    double curr=GetData(loc);
    int go=0;
    int steps=0;
    while(go==0) //take steps uphill until we reach a summit
      {
	//vector<int> indy; //reset indy
	vector<float> inty;        //rest inty
	for(int c=0;c<2*dim;++c) //go plus and minus in 4 directions
	  {
	    for(int i=0;i<dim;++i) //reset lac
	      lac[i]=loc[i];
	    if(c<dim) //increment lac
	      lac[c%dim]+=1;
	    else
	      lac[c%dim]-=1;

	    //std::cout << c%4 << std::endl;
	    //std::cout << row[0] << " " << row[1] << " " <<row[2] << " " <<row[3] << std::endl;
	    //indy.push_back(c);
	    inty.push_back(fabs(GetData(lac))-fabs(curr) );
	    //std::cout << GetData(row) << std::endl;
	  }
	double max=inty[0];
	int cmax=0;
	for(int c=1;c<2*dim;++c)
	  if(inty[c]>max)
	    {
	      cmax=c;
	      max=inty[c];
	    }
	if(max>0)
	  { //update loc
	    if(cmax<dim)
	      loc[cmax%dim]+=1;
	    else
	      loc[cmax%dim]-=1;
	    curr=GetData(loc);
	    steps+=1;
	  }
	else
	  {
	    break;
	  }
      }

    peakList[p].xOld=peakList[p].x;
    peakList[p].yOld=peakList[p].y;

    if(steps>0) //if we've made any adjustments...
      {
	//std::cout << ref[0] << " " << ref[1] << " " << ref[2] << " " << ref[3] << std::endl;
	//std::cout << ii << " " << jj << " " << kk << " " << ll << std::endl;
	//std::cout << currRef << std::endl;
	//std::cout << curr << std::endl;
	//std::cout << steps <<std::endl;

	/*std::cout << peakList[p].x << std::endl;
	std::cout << peakList[p].y << std::endl;
	std::cout << xvals[2][ref[2]] << std::endl;
	std::cout << xvals[3][ref[3]] << std::endl;*/

	//std::cout << xvals[0][ii] << std::endl;
	//std::cout << xvals[1][jj] << std::endl;
	//std::cout << xvals[2][kk] << std::endl;
	//std::cout << xvals[3][ll] << std::endl;


	//std::cout << "moved in proton: " << abs(xvals[2][kk]-xvals[2][ref[2]]) << std::endl;
	//std::cout << "moved in carbon: " << abs(xvals[3][ll]-xvals[3][ref[3]]) << std::endl;

	//std::cout << "protonDiff: " << abs(xvals[2][kk]-xvals[1][jj]) << std::endl;
	//std::cout << "carbonDiff: " << abs(xvals[3][ll]-xvals[0][ii]) << std::endl;


	//peakList[p].x=xvals[dimord[1]][loc[dimord[1]]];
	//peakList[p].y=xvals[dimord[0]][loc[dimord[0]]];
	peakList[p].x=xvals[dimord[3]][loc[dimord[3]]];
	peakList[p].y=xvals[dimord[2]][loc[dimord[2]]];

      }
    //cout << peakList[p].name << " " << peakList[p].x << " " << peakList[p].y << " " << peakList[p].xOld << " " << peakList[p].yOld << endl;
  }

  //write a peaklist file
void pipe::WritePeak()
  {
    FILE *fp;
    fp=fopen("out\\peakNew.list","w");
    for(int p=0;p<peaks;++p)
      fprintf(fp,"%s\t%f\t%f\n",peakList[p].name.c_str(),peakList[p].y,peakList[p].x);
    fclose(fp);
  }


void pipe::WritePipe4D(float *arr,int si,int sj,int sk, int sl)
  {
    string oName =iName+".decon";
    cout << "Writing to " << oName << endl;
    cout << si << " " << sj << " " << sk << endl;
    cout << sizeof(float) << endl;
    B   = fopen(oName.c_str(),"wb");
    fwrite(prehead,sizeof(float),512,B);
    cout << ftell(B) << " hello!" << endl;

    for(int i=0;i<si;++i)
      for(int j=0;j<sj;++j)
  for(int k=0;k<sk;++k)
  for (int l=0; l<sl; ++l)
    {
      //sliceLib3D[0].DS[i+j*sliceLib3D[0].si+k*sliceLib3D[0].si*sliceLib3D[0].sj];
      //float val=slicey.DS[i+j*size[0]+k*size[0]*size[2]];
      //float val=arr[i+j*size[0]+k*size[0]*size[1]];
      //int pos=(512+k+j*size[1]+i*size[1]*size[0])*sizeof(float);


      // float val=arr[i+j*si+k*sj*si];
      float val=arr[i+j*si+k*si*sj+l*si*sj*sk];
//	    float val=float(GetVal3D(i,j,k)); //write original data out (test)
//	    float val = 42;
      //cout << i << " " << j << " " << k << " " << val <<  endl;
//    float val=42;
//        int a1=ftell(B);

      fwrite(&val,sizeof(float),1,B);
//         int a2=ftell(B);
//         if( (a2-a1)!=4) {
//cout <<"shitballs " << a1 << " " << a2 << endl;
//exit(100);
//         }

   //   cout << ftell(B) << endl;
    }
//    cout << i << " " << j << " " << k << endl;
      fclose(B);
    cout << "done" << endl;
  }

#ifdef DOUBLE3D
void pipe::WritePipe3D(double *arr,int si,int sj,int sk)
  #else
  void pipe::WritePipe3D(float *arr,int si,int sj,int sk)
   #endif
  {
    string oName =iName+".decon";
    // cout << "Writing to " << oName << endl;
    // cout << si << " " << sj << " " << sk << endl;
    // cout << sizeof(float) << endl;
    B   = fopen(oName.c_str(),"wb");
    fwrite(prehead,sizeof(float),512,B);
    // cout << ftell(B) << " hello!" << endl;

    for(int i=0;i<si;++i)
      for(int j=0;j<sj;++j)
	      for(int k=0;k<sk;++k)
	  {
	    //sliceLib3D[0].DS[i+j*sliceLib3D[0].si+k*sliceLib3D[0].si*sliceLib3D[0].sj];
	    //float val=slicey.DS[i+j*size[0]+k*size[0]*size[2]];
	    //float val=arr[i+j*size[0]+k*size[0]*size[1]];
	    //int pos=(512+k+j*size[1]+i*size[1]*size[0])*sizeof(float);


	    float val=float(arr[i+j*si+k*sj*si]);
//	    float val=float(GetVal3D(i,j,k)); //write original data out (test)
//	    float val = 42;
	    //cout << i << " " << j << " " << k << " " << val <<  endl;
//    float val=42;
//        int a1=ftell(B);

	    fwrite(&val,sizeof(float),1,B);
//         int a2=ftell(B);
//         if( (a2-a1)!=4) {
//cout <<"shitballs " << a1 << " " << a2 << endl;
//exit(100);
//         }

	 //   cout << ftell(B) << endl;
	  }
//    cout << i << " " << j << " " << k << endl;
      fclose(B);
    cout << "done" << endl;
  }

void pipe::WritePipe2D(double *arr,int si,int sj,const std::string &tag)
  {
    string oName =iName+"."+tag;
    cout << "Writing to " << oName << endl;
    B   = fopen(oName.c_str(),"wb");
    fwrite(prehead,sizeof(float),512,B);
    for(int i=0;i<si;++i)
      for(int j=0;j<sj;++j)
	  {
	    //sliceLib3D[0].DS[i+j*sliceLib3D[0].si+k*sliceLib3D[0].si*sliceLib3D[0].sj];
	    //float val=slicey.DS[i+j*size[0]+k*size[0]*size[2]];
	    //float val=arr[i+j*size[0]+k*size[0]*size[1]];
	    //int pos=(512+k+j*size[1]+i*size[1]*size[0])*sizeof(float);

	    float val=float(arr[i+j*si]);
	    //float val=GetVal3D(i,j,k); //write original data out (test)
	    //cout << i << " " << j << " " << k << " " << val <<  endl;
	    fwrite(&val,sizeof(float),1,B);
	  }
      fclose(B);
    cout << "done" << endl;
  }

  
void pipe::WritePipe2D(float *arr,int si,int sj,const std::string &tag)
  {
    string oName =iName+"."+tag;
    cout << "Writing to " << oName << endl;
    B   = fopen(oName.c_str(),"wb");
    fwrite(prehead,sizeof(float),512,B);
    for(int i=0;i<si;++i)
      for(int j=0;j<sj;++j)
	  {
	    //sliceLib3D[0].DS[i+j*sliceLib3D[0].si+k*sliceLib3D[0].si*sliceLib3D[0].sj];
	    //float val=slicey.DS[i+j*size[0]+k*size[0]*size[2]];
	    //float val=arr[i+j*size[0]+k*size[0]*size[1]];
	    //int pos=(512+k+j*size[1]+i*size[1]*size[0])*sizeof(float);

	    float val=float(arr[i+j*si]);
	    //float val=GetVal3D(i,j,k); //write original data out (test)
	    //cout << i << " " << j << " " << k << " " << val <<  endl;
	    fwrite(&val,sizeof(float),1,B);
	  }
      fclose(B);
    cout << "done" << endl;
  }




#ifdef DOUBLE1D
void pipe::WritePipe1D(double *arr,int si)
  #else
  void pipe::WritePipe1D(float *arr,int si)
  #endif
  {
    string oName =iName+".decon";
    cout << "Writing to " << oName << endl;
    B   = fopen(oName.c_str(),"wb");
    fwrite(prehead,sizeof(float),512,B);
    for(int i=0;i<si;++i)
    {
      float val=float (arr[i]);
      fwrite(&val,sizeof(float),1,B);
    }
      fclose(B);
    cout << "done" << endl;
  }

  //optimise peak position based on auto peak maximum
void pipe::optPeak4D()
  {
    std::cout << "Locally optimising peaks... " << std::endl;
    for(int p=0;p<peaks;++p)
      {
	//int lslice=GetIndex(peakList[p].x,dimord[1]); //proton
	//int kslice=GetIndex(peakList[p].y,dimord[0]); //carbon
	//loc[dimord[1]]=lslice;
	//loc[dimord[0]]=kslice;
	//int kslice=GetIndex(peakList[p].x,2); //proton
	//int lslice=GetIndex(peakList[p].y,3); //carbon
	LocalMax(p);
      }

    cout << "Done local optimisation" << endl;
    for(int p=0;p<peaks;++p)
      {
	for(int q=p+1;q<peaks;++q)
	  {
	    if(peakList[p].x==peakList[q].x && peakList[p].y==peakList[q].y)
	      {
		std::cout << "Overlap in peaklist: " ;
		std::cout << peakList[p].name << " " <<peakList[q].name << ". Resetting. " << std::endl;
		peakList[p].x=peakList[p].xOld;
		peakList[p].y=peakList[p].yOld;
		peakList[q].x=peakList[q].xOld;
		peakList[q].y=peakList[q].yOld;
	      }
	  }
      }
    WritePeak();
  }

float pipe::GetVal3D(int &i,int &j, int &k)
  {
    SetPos3D(i,j,k);
    float Data2;
//      std::cout << i << " " << j << " " << k << " " << Data << endl;
      fread(&Data2,sizeof(float),1,A);
      float Data = float(Data2);
//        Data = 42;
//      cout << sizeof(Data) << endl;
//    if(==0)
//      {
//        std::cout << i << " " << j << " " << k << " " << Data << endl;
//	std::cout << "problem reading in file ... 3D ... " << std::endl;
//	exit(100);
//      }

    //fprintf(B,"%f\t%f\t%e\n",xvals[1][j],xvals[0][i],abs(Data));
    return Data;
  }




float pipe::GetVal2D(int &i,int &j)
  {
    SetPos2D(i,j);
    float Data;

    if(fread(&Data,sizeof(float),1,A)==0)
      {

  std::cout << i << "  " << j << endl;
	std::cout << "problem reading in file ... 2D ... " << std::endl;
	exit(100);
      }
    return Data;
  }

float pipe::GetVal1D(int &i)
  {
    SetPos1D(i);
    float Data;

    if(fread(&Data,sizeof(float),1,A)==0)
      {


	std::cout << "problem reading in file ... 1D ... " << std::endl;
	exit(100);
      }
    return Data;
  }


  //extract slice based on peaklist
void pipe::extract4D()
  {
    cout <<"Extracting..." << endl;
    int loc[4];
    for(int p=0;p<peaks;++p)
      {
	//0,1,2,3
	//c,p,c,p
	std::cout << "Extracting " << peakList[p].name << std::endl;
	loc[dimord[1]]=GetIndex(peakList[p].x,dimord[1]); //proton
	loc[dimord[0]]=GetIndex(peakList[p].y,dimord[0]); //carbon


	cout << p << " " << peakList[p].name << endl;
	cout << peakList[p].x << " " << xvals[dimord[1]][loc[dimord[1]]] << endl;
	cout << peakList[p].y << " " << xvals[dimord[0]][loc[dimord[0]]] << endl;

	/*if(p==48){
	  //A100 is 48
	  cout << p << " " << peakList[p].name << endl;
	  cout << peakList[p].x << " " << xvals[dimord[1]][loc[dimord[1]]] << endl;
	  cout << peakList[p].y << " " << xvals[dimord[0]][loc[dimord[0]]] << endl;
	  }*/

	//int kslice=GetIndex(peakList[p].x,2); //proton
	//int lslice=GetIndex(peakList[p].y,3); //carbon
	//         x   y    z   a
	//         i   j    k   l
	//         3   2    1   0
	//         C H(direct) H C
	//        128 236   64 128
	//so l is the 'fast moving' direcion on the inside.
	//j is the 'slowest moving' direction on the outside.
	//shape of data is given by DIMORDER.
	// dimOrd: 4   3     1   2
	//         l   k     i   j
	//         a   z     x   y

	string line ="out/slice2d/"+peakList[p].name+".dat.out";
	B = fopen(line.c_str(),"w");
	for(int i=0;i<size[dimord[2]];++i)
	  {
	    //int islice=GetIndex(xvals[0][i],0); //carbon
	    loc[dimord[2]]=GetIndex(xvals[dimord[2]][i],dimord[2]); //carbon ;
	    for(int j=0;j<size[dimord[3]];++j)
	      {
		loc[dimord[3]]=GetIndex(xvals[dimord[3]][j],dimord[3]); //proton ;
		//cout << xvals[dimord[3]][j] << " " << xvals[dimord[3]][loc[dimord[3]]] << endl; //proton
		//cout << dimord[2] << " " << xvals[dimord[2]][i] << " " << xvals[dimord[2]][loc[dimord[2]]] << " " <<  xvals[dimord[3]][j] << " " << xvals[dimord[3]][loc[dimord[3]]] << endl; //carbon
		//cout << "max: " << xvals[dimord[2]][0] << endl;
		//cout << "min: " << xvals[dimord[2]][size[dimord[2]]-1] << endl;
		//    while(val<xvals[col][size[col]-1])

		//int jslice=GetIndex(xvals[1][j],1); //proton
		//SetPos(islice,jslice,kslice,lslice);
		SetPosInt(loc);
		float Data;
		if(fread(&Data,sizeof(float),1,A)==0)
		  {
		  std::cout << "problem reading in file... " << std::endl;
		  //exit(100);
		  }
		fprintf(B,"%f\t%f\t%e\n",xvals[dimord[3]][j],xvals[dimord[2]][i],Data);
	      }
	    fprintf(B,"\n");
	  }
	fclose(B);
      }
  }


  //extract slice based on peaklist
void pipe::diag4D()
  {
    float proj2D[size[dimord[0]]*size[dimord[1]]];

    cout <<"Projecting diagonal..." << endl;
    int loc[4];
    string line ="out/diag.out";
    //B = fopen(line.c_str(),"w");
    for(int i=0;i<size[dimord[0]];++i)
      {
	for(int j=0;j<size[dimord[1]];++j)
	  {
	    loc[dimord[0]]=GetIndex(xvals[dimord[0]][i],dimord[0]); //proton ;
	    loc[dimord[1]]=GetIndex(xvals[dimord[1]][j],dimord[1]); //proton ;
	    loc[dimord[2]]=GetIndex(xvals[dimord[0]][i],dimord[2]); //proton ;
	    loc[dimord[3]]=GetIndex(xvals[dimord[1]][j],dimord[3]); //proton ;
	    SetPosInt(loc);
	    float Data;
	    if(fread(&Data,sizeof(float),1,A)==0)
	      {
		std::cout << "problem reading in file... " << std::endl;
		//exit(100);
	      }
	    proj2D[i+size[dimord[0]]*j]=Data;
	      //fprintf(B,"%f\t%f\t%e\n",xvals[dimord[0]][i],xvals[dimord[1]][j],Data);
	  }
	//fprintf(B,"\n");
      }
    //fclose(B);

    //Get input directory:
    string sentence= iName;
    //cout << sentence << endl;
    istringstream iss(sentence);
    std::vector<std::string> tokens;
    std::string token;
    while (std::getline(iss, token, '/')) {
      if (!token.empty())
	tokens.push_back(token);
    }
    string path;
    for (int i=0;i<tokens.size()-1;++i)
      path+=tokens[i]+'/';

    string i2dname=path+"projections/"+labels[dimord[0]]+"."+labels[dimord[1]]+".dat";
    // string i2dname=path+"projections/"+labels[dimord[0]]+"."+"C13n"+".dat";
    cout << i2dname << endl;

    pipe pipefile2D;
    pipefile2D.iName=i2dname;
    cout << 1 << endl;
    pipefile2D.readHeader();
    cout << 2 << endl;
    pipefile2D.WritePipe2D(proj2D,size[dimord[0]],size[dimord[1]]);
    cout << 3 << endl;


  }

  //make projections xy,yz, xz.
  //reading the inner 'a' dimension is always fast.
void pipe::project4D()
  {
    //x y z a
    //i j k l
    //1 2 3 4
    //H C H C
    std::cout <<"Making projections..." << std::endl;
    std::string p1="out/slice2d/xy.out";
    std::string p2="out/slice2d/xz.out";
    std::string p3="out/slice2d/yz.out";
    std::cout << "xy" << std::endl;
    B= fopen(p1.c_str(),"w");
    for(int j=0;j<size[1];++j)
      {
	int jslice=GetIndex(xvals[1][j],1); //proton
	for(int i=0;i<size[0];++i)
	  {
	    int islice=GetIndex(xvals[0][i],0); //carbon
	    double val=0;
	    for(int k=0;k<size[2];++k)
	      {
		int kslice=GetIndex(xvals[2][k],2); //carbon
		int lslice=0;
		//4 3 1 2
		//l k i j
		SetPos(islice,jslice,kslice,lslice);
		float Data[size[3]];
		if(fread(Data,sizeof(float),size[3],A)==0)
		  std::cout << "problem reading in file..." << std::endl;
		else
		  for(int l=0;l<size[3];++l)
		    if(abs(Data[l])>noise)
		      val+=Data[l];
	      }
	    fprintf(B,"%f\t%f\t%e\n",xvals[1][j],xvals[0][i],val);
	  }
	fprintf(B,"\n");
      }
    fclose(B);

    std::cout << "xz" << std::endl;
    B= fopen(p2.c_str(),"w");
    for(int k=0;k<size[2];++k)//z
      {
	int kslice=GetIndex(xvals[2][k],2); //carbon
	for(int i=0;i<size[0];++i)//x
	  {
	    int islice=GetIndex(xvals[0][i],0); //carbon
	    double val=0;
	    for(int j=0;j<size[1];++j)//y
	      {
		int jslice=GetIndex(xvals[1][j],1); //proton
		int lslice=0;
		//4 3 1 2
		//l k i j
		SetPos(islice,jslice,kslice,lslice);
		float Data[size[3]];
		if(fread(Data,sizeof(float),size[3],A)==0)
		    std::cout << "problem reading file..." << std::endl;
		else
		  for(int l=0;l<size[3];++l)
		    if(abs(Data[l])>noise)
		      val+=Data[l];
	      }
	    fprintf(B,"%f\t%f\t%e\n",xvals[2][k],xvals[0][i],val);
	  }
	fprintf(B,"\n");
      }
    fclose(B);


    std::cout <<"yz" << std::endl;
    B= fopen(p3.c_str(),"w");
    for(int j=0;j<size[1];++j)//y
      {
	int jslice=GetIndex(xvals[1][j],1); //proton
	for(int k=0;k<size[2];++k)//z
	  {
	    int kslice=GetIndex(xvals[2][k],2); //carbon
	    double val=0;
	    for(int i=0;i<size[0];++i)//x
	      {
		int islice=GetIndex(xvals[0][i],0); //carbon
		int lslice=0;
		//4 3 1 2
		//l k i j
		SetPos(islice,jslice,kslice,lslice);
		float Data[size[3]];
		if(fread(Data,sizeof(float),size[3],A)==0)
		  std::cout << "problem reading file..." << std::endl;
		else
		  for(int l=0;l<size[3];++l)
		    if(abs(Data[l])>noise)
		      val+=Data[l];
	      }
	    fprintf(B,"%f\t%f\t%e\n",xvals[1][j],xvals[2][k],val);
	  }
	fprintf(B,"\n");
      }
    fclose(B);


  }


#endif
