/*Legal_notice:

Copyright 1996-2001 by Soren M. Kristensen, Department of Chemistry,
University of Copenhagen, DK-2100 Copenhagen O, Denmark.  All rights
reserved.

This software and its related software and documentation is provided 'as
is' without express or implied warrenty. The Departmet of Chemistry, The
University of Copenhagen and the author make no warrenties as to any matter
whatsoever with respect to the program and the related software and
documentation. In particular, any and all warranties of merchantability and
fitness for any particular purpose are expressly excluded.

Legal_notice*/
#include <cmath>
#include <string>
#include <list>
#include <vector>
#include <complex>
#include "fuda_classes.H"
#include "fuda_utils.H"
#include "fuda_cftypes.H"

// Declaration of the cftype declaration functions.
namespace FUDA 
{
  void declare_cftypes_base(Fuda *fuda, std::string tag);
  void declare_cftypes_relax(Fuda *fuda, std::string tag);
  void declare_cftypes_user(Fuda *fuda, std::string tag);
}

// Here we declare all the cftypes.
void FUDA::declare_cftypes(Fuda *fuda, std::string tag)
{
  // We catch all exceptions here if not caught in the declare routines.
  try {

    // Call each of the cftype declaration routines.
    declare_cftypes_base(fuda, tag);
    declare_cftypes_relax(fuda, tag);
    declare_cftypes_user(fuda, tag);

  }
  catch(Uferr::FtypeNameInvalid& e) {
    std::cout << "declare_cftypes : ftype name invalid: "
	      << e.name << "\n";
    exit(1);
  }
  catch(Uferr::FtypeEmptyFtypeVec) {
    std::cout << "declare_cftypes : Empty ftype vec\n";
    exit(1);
  }
  catch(Uferr::FtypeRefInvalid) {
    std::cout << "declare_cftypes : FtypeRefInvalid\n";
    exit(1);
  }
  catch(Uferr::ParamNameInvalid& e) {
    std::cout << "declare_cftypes : param name invalid: "
	      << e.name << "\n";
    exit(1);
  }
  catch(Uferr::FtypeNameAlreadyUsed& e) {
    std::cout << "declare_cftypes : ftype name already used: "
	      << e.name << "\n";
    exit(1);
  }
  catch(Uferr::NumParamDescrInvalid& e) {
    std::cout << "declare_cftypes : number of parameter "
	      << "descriptors invalid: "
	      << e.n << "\n";
    exit(1);
  }
  catch(Uferr::NumVarInvalid& e) {
    std::cout << "declare_cftypes : number of variables invalid: "
	      << e.n << "\n";
    exit(1);
  }
  catch(...) {
    std::cout << "declare_cftypes : unspecified exception\n";
    exit(1);
  }
}






