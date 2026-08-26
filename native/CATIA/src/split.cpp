#include<StringMethods.h>
#include <Abort.h>
/*
  splits a text string according to the deliminator
  returns the items in an vector array strings

*/

std::vector<std::string> split(std::string inp, const char* deliminator ) {

  /*
  if ( inp.length() > MAX_STRING_LENGTH ){
    std::cerr<<" Length of string exceeded MAX_STRING_LENGTH"<<std::endl;
    std::cerr<<"-->"<<inp<<"<--"<<std::endl;
    std::cerr<<" please recompile "<<std::endl;
    std::cerr<<" PROGRAM ABORTED "<<std::endl;
    Abort(1);
  };
  */
  bool go=true;
  int i=0, j=0, count=0;
  int start[inp.length()];
  int stop[inp.length()];
  int pos[inp.length()+5];

  std::string del(deliminator);
  std::vector<std::string> out;

  // Find the number of delimators in the input
  for (int i =0; i< inp.size(); i++) {
    if ( (inp.c_str()[i] == deliminator[0]) && (inp.c_str()[i] != inp.c_str()[i+1]) ) {
      count++;
    };
  };

  if (inp.find(del,0) < inp.size() ) { // Found at least one deliminator

    int check_pos=0;

    //Find the start of the first part.
    if ( inp.c_str()[0] != deliminator[0] ) {//Check if the first character is a deliminator
      start[0]=0;
    } else { //Find start of first part
      start[0]=0;
      while ( go ) {
	check_pos=inp.find(del,start[0]+1);
	if (  check_pos != start[0]+1 ) {
	  go=false;
	};
	start[0]++;
      };
    };
    int partcount=0;
    bool go2=true;
    while( go2 ) {
      stop[partcount]=inp.find(del,start[partcount]+1); //End of part
      if  ( !( stop[partcount] < inp.npos) ) {
	go2=false;
      };
      out.push_back(
	     inp.substr(start[partcount],stop[partcount]-start[partcount])
	     );
      // Find start of next part
      partcount++;
      start[partcount]=stop[partcount-1]+1;
      go=true;
      start[partcount]--; // Go one step back in the string
      while ( go ) {
	check_pos=inp.find(del,start[partcount]+1);
	if (  check_pos != start[partcount]+1 ) {
	  go=false;
	};
	start[partcount]++;
      };
    };
  }

  /* check whether the last part is empty:
     if so: Delete and resize vector
  */
  std::vector<std::string> final_out;

  if ( count == 0) {
    final_out.push_back(inp);
  } else {


    if ( out[out.size()-1] == "" ) {
      if ( out.size() == 1 ) {
	final_out.push_back(out[0]);
      } else {
	for (int i=0; i < out.size()-1 ; i++ ) {
	  final_out.push_back(out[i]);
	};
      }
    } else {
      final_out=out;
    };
  };

  return final_out;
};
