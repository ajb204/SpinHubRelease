/**************************************************/
/* decon         */
/* A.Baldwin     */
/* 3rd March 2014*/
/**************************************************/
#ifndef GENERAL_CPP
#define GENERAL_CPP


#include "general.hpp"

//stick a file in a vector of vectors
vector<vector<string> > MakeFileVec(string file_name)
{
  
  ifstream in(file_name.c_str());
  vector<vector<string> > infile;
  string line;
  while(getline(in,line)){
    istringstream iss(line);
    vector<string> tokens;
    copy(istream_iterator<string>(iss),
	 istream_iterator<string>(),
	 back_inserter(tokens));
    infile.push_back(tokens);
  }
  in.close();
  return infile;
  
}

#endif
