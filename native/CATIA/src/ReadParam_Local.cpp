#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <Catia.h>
#include <StringMethods.h>
#include <Abort.h>

using boost::to_lower_copy;
using boost::trim_copy;

typedef std::pair<std::string,int> pairSI_t;

void Catia::ReadParam_Local(std::string infile){
    std::vector<std::string> Format;
    //The first line must contain the format string
    std::ifstream ifs(infile.c_str());
    if(!ifs){
        std::cerr<<" Could not open the inputfile "<<infile<<"\n";
        std::cerr<<" Function .ReadParam_Local()\n";
        std::cerr<<std::endl;
        Abort(1);
    }
    char line[MAX_STRING_LENGTH];
    bool commentline=true;
    while(commentline){
        ClearBuf(line,sizeof(line));
        ifs.getline(line,sizeof(line));
        if(line[0]=='#'){
            commentline=true;
        } else {
            commentline=false;
        }
    }
    std::string firstl(line);
    std::vector<std::string> fits=split(firstl,"=");
    fits[0]=trim_copy(to_lower_copy(fits[0]));
    if(!(fits[0]=="format")){
        std::cerr<<" The first line of the ReadParam_Local input file must\n";
        std::cerr<<" contain the format string parameters - e.g.\n";
        std::cerr<<" format=(DeltaO;DeltaJ;R2CW500)\n";
        std::cerr<<" Functione: ReadParam_Local()\n";
        std::cerr<<std::endl;
        Abort(1);
    } else {
        fits[1]=fits[1].substr(fits[1].find_first_of('(')+1,
                            fits[1].find_last_of(')')-1-fits[1].find_first_of('('));
        fits=split(fits[1],";");
        for(unsigned int i=0;i<fits.size();i++){
            Format.push_back(trim_copy(fits[i]));
        }
    }
    //Allocate space on the LocalParam vector array
    LocalParam.resize(Atoms.size());
    LocalParamF.resize(Atoms.size());
    LocalParamE.resize(Atoms.size());
    LocalNotes.resize(Atoms.size());
    //
    while(!ifs.eof()){
        ClearBuf(line,sizeof(line));
        ifs.getline(line,sizeof(line));
        std::string l(line);
        std::vector<std::string> its=split(l,"=");
        if(l.length()<1||its.size()<2||line[0]=='#'){
            continue;
        }


        its[0] = trim_copy(its[0]);
        if(its[0]=="*") {
            its[0]=".+";
        }
        std::vector<int> AtomNumber;

        const boost::regex regularExpression(its[0]);
         int flag = 0;

        BOOST_FOREACH(pairSI_t atom, Atoms) {
             if ( regex_match(atom.first, regularExpression)) {
        //if ( regex_match(atom.first, its[0])) {
                AtomNumber.push_back(atom.second);
                flag = 1;
            }
        }
        its[1]=its[1].substr(its[1].find_first_of('(')+1,
                             its[1].find_last_of(')')-1-its[1].find_first_of('('));
        its=split(its[1],";");
        //Lets check that the right number of parameters are provided
        if(!(its.size()==Format.size())){
            std::cerr<<" The line "<<l<<"\n";
            std::cerr<<" contains "<< its.size() <<" elements, which is different from "<<Format.size()<<" parameters\n";
            std::cerr<<" that was specified in the format string:\n";
            std::cerr<<"->"<<firstl<<"<-"<<std::endl;
            std::cerr<<" Function: ReadParam_Local();"<<std::endl;
            Abort(1);
        }

        for(unsigned int i=0;i<AtomNumber.size();++i) {
            for(unsigned int j=0;j<its.size();++j) {
                LocalParam[AtomNumber[i]][Format[j]]=atof(its[j].c_str());
                LocalParamF[AtomNumber[i]][Format[j]]=1;
                LocalParamE[AtomNumber[i]][Format[j]]=-1;
            }
        }
    }
    return;
}
