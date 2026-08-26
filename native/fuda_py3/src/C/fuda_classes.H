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

#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <cstdio>
#include <cmath>

#define FUDA_VERSION_MAJOR (0)
#define FUDA_VERSION_MINOR (5)

// Pre-declarations of classes.
class Ptype;
class Param;
class Dtype;
class Data; class DataPtr_less;
class Ftype;
class Func;
class Fuda;
class FuncCall;
class IntIndex;
class DataDtype;


namespace FUDA
{
  // The different kinds a parameter may be.
  enum Param_kind {PKIND_EXPL, PKIND_PARAM, PKIND_CONST};

  // The different types of Ftypes.
  enum Ftype_type {FTYPE_CFUNC, FTYPE_PROD, FTYPE_COMP, FTYPE_SUM,
		   FTYPE_PYTHON};

  /* CFUNC fitting function typedef. The contents of p may only be read, never
     modified. 

     The cmd specifies the action to take by the function:
     1 : calculate value, do not address dp and dp_flg.
     2 : calculate value and dp according to dp_flg array.
     3 : calculate dp according dp_flg array, do not address value. 

     The return value specifies a status:
     0 : success.
     1 : error. */
  typedef int (*CFUNC)(void *fs,
		       double p[],
		       int dp_flg[],
		       double dp[],
		       double *value);
  
}

/************************************************************************/
/* Ptype: Parameter-type structure. The parameter type is a type
   declaration for later declaration of many parameters of the same
   type. Since all parameters has an associated paramter type, later in a
   calculation, all parameters of a certain type can be manipulated in a
   particular way. Also the parameter type has default values for many of
   the settings for a paramter */
class Ptype
{
  Fuda *fuda;         /* Pointer to owner fuda structure*/
  FUDA::Param_kind kind;      /* Parameter kind */
  std::string name;           /* Name of paramter type */

  double init;           /* Initial value of paramter */
  bool free_flg;         /* Is the paramter free or fixed in the */
  unsigned int bounds;   /* Apply bounds when minimizing */
  double lower;          /* Lower allowed value of parameter */
  double upper;          /* Upper allowed value */
  double delta;          /* delta_param defining approximate order of
                           magnitude of parameter, used in simplex. */
  double esd;           /* e.s.d. of parameter */

  // Constructor.
  Ptype(Fuda *a_fuda, std::string& name, FUDA::Param_kind a_kind);

  // Destructor.
  ~Ptype();

 public:
  friend class Fuda;

  Fuda *get_fuda();
  void get_name(std::string& name);
  FUDA::Param_kind get_kind();
  bool is_free();
  void set_free(bool a_free);
  void free();
  void fix();
  bool is_referenced();
  double get_init();
  void set_init(double a_init);
  double get_value();
  void set_value(double a_init);
  unsigned int get_bounds();
  void set_bounds(unsigned int a_bounds);
  double get_lower();
  void set_lower(double a_lower);
  double get_upper();
  void set_upper(double a_upper);
  double get_delta();
  void set_delta(double a_delta);
  double get_esd();
  void set_esd(double a_esd);

  void print();
};

typedef std::list<Ptype*>::iterator Ptype_iterator;

/************************************************************************/
/* Param: Paramter structure. This structure holds information about
   a dependent parameter or an independent variable depending on hte value
   of depend. A dependent parameter is a fitting parameter if free and a
   constant if not free. A paramter is associated with a parameter type and
   alsow holds information about fitting specific settings */
class Param
{
  Fuda *fuda;         /* Pointer to owner fuda structure*/
  FUDA::Param_kind kind;  /* Defines the kind of the parameter. */
  std::string name;     /* Name of paramter */
  double value;         /* Actual value of paramter */

  Ptype *ptype;         /* Ref. to the type of the parameter. The ptype
			   may be NULL. */
  double init;          /* Initial value of paramter */
  bool free_flg;         /* Is the paramter free or fixed in the
                           calculation */
  unsigned int bounds;  /* Apply bounds when minimizing */
  double lower;         /* Lower allowed value if bounds */
  double upper;         /* Upper allowed value if bounds */
  double delta;         /* delta_param defining approximate order of
                           magnitude of parameter, used in simplex. */
  double esd;           /* e.s.d. of parameter. Can be set by user
                           or by update */

  /* Evaluation specific data */
  struct Eval
  {
    bool active;        /* Set if parameter is an explanatory,
                           dependent or constant parameter in
                           minimization. */
    unsigned int index; /* If active: index is the index of the
                           paramter in the fitting run. The index
                           refers to a position in the explanatory
                           array, the dependent variable array or the
                           constant array depending on the kind and
                           free */
  } eval;

  // Constructors.
  Param(Fuda *a_fuda, std::string& name, Ptype *type);
  Param(Fuda *a_fuda, std::string& name, FUDA::Param_kind a_kind);

  // Destructor.
  ~Param();

 public:
  friend class Fuda;

  Fuda *get_fuda();
  void get_name(std::string& name);  
  void get_ptype_name(std::string& name);  
  Ptype *get_ptype();
  FUDA::Param_kind get_kind();
  bool is_free();
  void set_free(bool a_free);
  void free();
  void fix();
  bool is_referenced();
  double get_init();
  void set_init(double a_init);
  void init_value();
  double get_value();
  void set_value(double a_value);
  unsigned int get_bounds();
  void set_bounds(unsigned int a_bounds);
  double get_lower();
  void set_lower(double a_lower);
  double get_upper();
  void set_upper(double a_upper);
  double get_esd();
  void set_esd(double a_esd);
  double get_delta();
  void set_delta(double a_delta);
  bool is_eval_active();
  unsigned int get_eval_index();

  void print();
};

typedef std::list<Param*>::iterator Param_iterator;


/************************************************************************/
/* Dtype: Data type structure. Every data point and every function is
   associated with a data type. Functions with a particular data type
   only applies to data with the same data type. This provides a way
   of having different kinds of data for different functions.

   pvec must be a vector of references to parameters marked as
   independent variables and helps identify the data points of the
   corresponding type. */
class Dtype
{
  Fuda *fuda;             /* Pointer to owner fuda structure*/
  std::string name;         /* Data type name */
  std::vector<Param*> pvec; /* references to parameters */
  bool purge_flg;         /* Apply purging of data */
  std::vector<double> purge_radius; /* Purging radius for data x[] */

  /* Evaluation structure is set up prior to calling the minimizer */
  struct Eval
  {
    std::list<Func*> funcs; /* List of references to Func records */
    std::list<Data*> data; /* List of references to Data
                                           of this dtype */
  } eval;

  // Constructor.
  Dtype(Fuda *a_fuda, std::string& a_name, unsigned int dim);

  // Destructor.
  ~Dtype();

 public:
  friend class Fuda;

  Fuda *get_fuda();
  void get_name(std::string& a_name);
  unsigned int get_dim();
  void set_ip(unsigned int i, Param *p);
  Param *get_ip(unsigned int i);
  void set_purge(bool flg);
  bool is_purge();
  bool is_referenced();
  void set_purge_radius(unsigned int ip, double radius);
  double get_purge_radius(unsigned int ip);
  void purge();
    
  void print();
};


typedef std::list<Dtype*>::iterator Dtype_iterator;


// Comprare two pointers to Data objects.
class DataPtr_less
{
public:
  bool operator()(Data*,Data*);
};


/************************************************************************/
/* Data: Data point structure. A data point record holds information
   about a single experimental data point, including independent variable
   vector x, data point value and uncertainty.

   Data points are associated with a Dtype which holds information about the
   dimension of the x vector. Data records belonging to different Dtypes can
   never be mixed mixed by convention. 

   The data record is minimal to minimize storage requirements */
class Data
{
  Func *func;              /* Ref. to the func the data record is
			      associated with */
  std::vector<double> x;   /* Independent variable vector */
  double init;             /* Initial value set for data point */
  double value;            /* Current experimental value used in
                              subsequent calculations. For usual
                              fitting calculations this value will be
                              identical to init_value, but for
                              Montecarlo calculations this value is
                              set to init plus some random error
                              contribution */
  double u;                /* uncertainty of data value */

  Data(Func *a_func, unsigned int dim);

 public:
  friend class Fuda;
  friend class Func;
  friend void Dtype::purge();
  friend bool DataPtr_less::operator()(Data*, Data*);

  ~Data();
  Func *get_func();
  Dtype *get_dtype();
  unsigned int get_dim();
  double get_ix(unsigned int i);
  void set_ix(unsigned int i, double a_x);
  double get_init();
  void set_init(double a_value);
  void init_value();
  void montecarlo_value(double random_factor);
  double get_value();
  void set_value(double a_value);
  double get_u();
  void  set_u(double a_u);
  void print();
};

typedef std::list<Data>::iterator Data_iterator;


/************************************************************************/
/* Ftype: Function-type classes. */


// Base (abstract) Ftype class.
class Ftype
{
  /* Pointer to owner fuda structure */
  Fuda *fuda;

  /* Function type name */
  std::string name;

  /* Function description string */
  std::string descr;

  /* Number of parameters */
  unsigned int nparam;

  /* Parameter names */
  std::vector<std::string> p_name;

  /* Parameter descriptor strings */
  std::vector<std::string> p_descr;

  /* True if the associated function supports derivative calculation */
  std::vector<bool> p_deriv;

  /* true if parameter is a variable */
  std::vector<bool> p_var;

  /* Number of variable parameters */
  unsigned int nvar;

  /* Vector with indeces of the variable parameters (kind is PARAM or
     EXPL). p_index must reference each parameter at most once. The
     order of p_index defines the order of the parameters of functions
     of this ftype. Pararmeters not present in p_index become CONST
     kind pararmeters named according to the names in p_names. */
  std::vector<unsigned int> var_index;

  void set_var_index();

 public:
  Ftype(Fuda *a_fuda, std::string a_name);
  virtual ~Ftype();

  Fuda *get_fuda();
  void set_nparam(unsigned int a_nparam);
  void set_descr(std::string a_descr);
  void set_p_name(unsigned int i, std::string a_name);
  void set_p_descr(unsigned int i, std::string a_descr);
  void set_p_var(unsigned int i, bool var_flg);
  void set_p_deriv(unsigned int i, bool deriv_flg);
  void set_param(unsigned int i, bool var_flg, bool deriv_flg,
                 std::string a_name, std::string a_descr);
  void get_name(std::string& a_name);
  unsigned int get_nparam();
  unsigned int get_nvar();
  void get_p_name(unsigned int i, std::string& a_name);
  void get_p_descr(unsigned int i, std::string& a_descr);
  bool is_p_deriv(unsigned int i);
  bool is_p_var(unsigned int i);
  unsigned int get_var_index(unsigned int i);
  void get_descr(std::string& a_descr);
  virtual double call(std::vector<double>& p_vec) = 0;
  virtual void call(double p[], int dp_flg[], double dp[], double *value) = 0;
  virtual FUDA::Ftype_type get_type() = 0;
  virtual void print();
};

typedef std::list<Ftype*>::iterator Ftype_iterator;


// FtypeCFUNC.
class FtypeCFUNC : public Ftype
{
  FUDA::CFUNC func;       /* Pointer to compiled function */
  void *fs;             /* Pointer to function record which is given
			   to the func and so provides a way to supply
			   extra information to func */
  unsigned int fs_size; /* Size of fs record in bytes */

  FtypeCFUNC (Fuda *a_fuda, std::string a_name, unsigned int nparam,
	      FUDA::CFUNC a_func_ptr, void *a_fs, unsigned int a_fs_size);
  ~FtypeCFUNC();

 public:
  friend class Fuda;

  FUDA::CFUNC get_func();
  void *get_fs();
  unsigned int get_fs_size();
  double call(std::vector<double>& p_vec);
  void call(double p[], int dp_flg[], double dp[], double *value);
  FUDA::Ftype_type get_type();
  void print();
};


// FtypePROD.
class FtypePROD : public Ftype
{
  bool scale_flg;             /* If true, a scale parameter is included. */
  std::vector<Ftype*> ftypes; /* References to the ftypes forming the product
				 ftype */
  std::vector<double> norm;   /* Normalization vector used by int call
				 function. */
  
  FtypePROD (Fuda *a_fuda, std::string a_name,
	     bool a_scale_flg, std::vector<Ftype*>& a_ftypes);
  ~FtypePROD();

 public:
  friend class Fuda;

  Ftype *get_ftype(unsigned int i);
  bool is_scaled();
  unsigned int get_nftype();
  double call(std::vector<double>& p_vec);
  void call(double p[], int dp_flg[], double dp[], double *value);
  FUDA::Ftype_type get_type();
  void print();
};


// FtypeCOMP defines a composite f(g()) ftype.
class FtypeCOMP : public Ftype
{
  bool scale_flg;  /* If true, a scale parameter is included. */
  Ftype* f; /* A reference to the outer function type f() */
  Ftype* g; /* A reference to the inner function type g() */
  unsigned int p_index; /* The index of the paraemter in f() that
			   g() replaces */

  FuncCall *fcall; /* Workspace for f function call */
  FuncCall *gcall; /* Workspace for g function call */
  
  FtypeCOMP (Fuda *a_fuda, std::string a_name,
	     bool a_scale_flg, unsigned int p_index,
	     Ftype *a_f, Ftype *a_g);
  ~FtypeCOMP();

 public:
  friend class Fuda;

  Ftype *get_f();
  Ftype *get_g();
  unsigned int get_p_index();
  bool is_scaled();
  double call(std::vector<double>& p_vec);
  void call(double p[], int dp_flg[], double dp[], double *value);
  FUDA::Ftype_type get_type();
  void print();
};


// FtypeSUM.
class FtypeSUM : public Ftype
{
  bool scale_flg;             /* If true, a scale parameter is included. */
  std::vector<Ftype*> ftypes; /* References to the ftypes forming the product
				 ftype */
  
  FtypeSUM (Fuda *a_fuda, std::string a_name,
	     bool a_scale_flg, std::vector<Ftype*>& a_ftypes);
  ~FtypeSUM();

 public:
  friend class Fuda;

  Ftype *get_ftype(unsigned int i);
  bool is_scaled();
  unsigned int get_nftype();
  double call(std::vector<double>& p_vec);
  void call(double p[], int dp_flg[], double dp[], double *value);
  FUDA::Ftype_type get_type();
  void print();
};


/************************************************************************/
/* Some structures needed by Func. */
class IntIndex
{
 public:
  unsigned int iv;
  unsigned int val;
};

typedef std::list<IntIndex>::iterator IntIndex_iterator;


/************************************************************************/
/* Func: Function structure. */
class Func
{
  Fuda *fuda;         /* Pointer to owner fuda structure*/
  std::string name;     /* Function name */
  bool use;             /* Use function in fitting run */

  Ftype *ftype;         /* Function type ref. */
  std::vector<Param*> params; /* Array of refs. to parameters */

  Dtype *dtype;         /* Data type reference */
  std::list<Data> data; /* data points  associated with function */

  /* Evaluation structure is set up by Func::eval_x_init(),
     Func::eval_pc_init() and Fuda::eval_init() */
  struct Eval
  {
    std::vector<IntIndex> x; /* Vector of (func. arg. #, expl. var #) */
    std::vector<IntIndex> p; /* Vector of (func. arg. #, dep. param #) */
    std::vector<IntIndex> c; /* Vector of (func. arg. #, const. #) */
  } eval;

  Func(Fuda *a_fuda, std::string& a_name,
       Ftype *a_ftype, Dtype *a_dtype);
  ~Func();
  void set_param(unsigned int i, Param *p); // Set a parameter
  void set_var(unsigned int i, Param *p); // Set a variable parameter
  void eval_x_init(); /* Setup eval.x vector completely. */
  void eval_pc_init(); /* Setup eval.p[].iv and eval.c[].iv values. */

 public:
  friend class Fuda;

  Fuda *get_fuda();
  void get_name(std::string& name);
  bool get_use();
  void set_use(bool a_use);
  Param *get_param(unsigned int i);
  Param *get_var(unsigned int i);
  Ftype *get_ftype();
  Dtype *get_dtype();
  double call(std::vector<double>& expl_vec);
  double call_by_expl(std::vector<double>& expl_vec);
  double call_by_var(std::vector<double>& v_vec);
  double call_by_var(std::vector<double>& v_vec,
		     std::vector<double>& dv_vec);
  double call_by_param(std::vector<double>& p_vec);

  unsigned int get_ndata();
  unsigned int get_nexpl();
  Data *add_data();
  void delete_data();
  void print_data();
  Data_iterator data_begin();
  Data_iterator data_end();

  void print();
};

typedef std::list<Func*>::iterator Func_iterator;


/************************************************************************/
/* Fuda. The Fuda structure contains lists of all the relevant fuda
   record types and provides an interface to these lists which ensures
   that the lists are consistent with each other */
class Fuda
{
  // Lists.
  std::list<Ptype*> ptypes; // List of parameter types.
  std::list<Param*> params; // List of paramters.
  std::list<Dtype*> dtypes; // List of data types.
  std::list<Ftype*> ftypes; // List of function types.
  std::list<Func*>  funcs;  // List of functions.

  /* Sync flags related to lists. The difference between the sync and
     rsync flags is that e.g. if param_sync=0 and param_rsync=1
     then some param has been modified since the last eval_init()
     call, but all references to parameters are still valid. This is
     the case if some param has been modified. In this case fuctions
     relying on the eval pointer lists, like eval_get_free(i), are
     still available. If param_rsync=0, some parameter has changed
     reference or has been deleted and eval_get_free(i) is not
     allowed, because a pointer in the eval structure may reference a
     non-existing parameter which is illegal.
   */
  bool param_sync;
  bool param_rsync;
  bool dtype_sync;
  bool dtype_rsync;
  bool func_sync;
  bool func_rsync;

  // Data sync flags.
  bool data_sync;
  bool data_rsync;

  // Fuda workspace.
  unsigned int wsize;
  bool wlocked;
  double *wspace;
  
  // Eval structure.
  struct Eval
  {
    bool sync; /* Indicates that eval struct is in sync with fundamental
		  Fuda lists, incl. ptypes, params, dtypes,... . */

    unsigned int mod_count; /* The counter is incremented when the
                               eval structure is cleared and the value
                               can be used by other classes to check
                               if eval has been modified */
    
    std::vector<Param*> x; /* Vector of ptrs to explanatory parms. */
    std::vector<Param*> p; /* Vector of ptrs to free params. */
    std::vector<Param*> c; /* vector of ptrs to constant (fixed) params. */
    std::vector<Data*> d; /* Data record pointer vector. */

  } eval;

  // Function call structure.
  FuncCall *fcall;

  // Private functions.
  void param_set_sync();
  void param_set_rsync();
  void dtype_set_sync();
  void dtype_set_rsync();
  void func_set_sync();
  void func_set_rsync();
  void data_set_sync();
  void data_set_rsync();
  void eval_set_sync();
  void eval_incr_mod_count();
  void eval_set_mod_count(unsigned int value);
  void wspace_deallocate();
  void wspace_allocate(unsigned int minsize);
  unsigned int wspace_size();  
  double *wspace_ref(unsigned int minsize);
  bool wspace_lock();
  void wspace_unlock();

 public:
  Fuda();

  // Ptype stuff.
  Ptype *ptype_add(std::string name, FUDA::Param_kind a_kind);
  void ptype_del(std::string& name);
  void ptype_del(const char name[]);
  void ptype_del_all();
  Ptype *ptype_find(std::string& name);
  Ptype *ptype_find(const char name[]);
  Ptype_iterator ptype_begin();
  Ptype_iterator ptype_end();
  unsigned int get_nptype();

  // Param stuff.
  Param *param_add(std::string name, std::string ptype_name);
  Param *param_add(std::string name, FUDA::Param_kind a_kind);
  Param *param_find(std::string& name);
  Param *param_find(const char name[]);
  void param_del(std::string& name);
  void param_del(const char name[]);
  void param_del_all();
  Param_iterator param_begin();
  Param_iterator param_end();
  unsigned int get_nparam();
  void param_clear_sync();
  void param_clear_rsync();
  bool param_is_sync();
  bool param_is_rsync();
  

  // Dtype stuff.
  Dtype *dtype_add(std::string name, std::list<std::string>& p_names);
  void dtype_del(std::string& name);
  void dtype_del(const char name[]);
  void dtype_del_all();
  Dtype *dtype_find(std::string& name);
  Dtype *dtype_find(const char name[]);
  Dtype_iterator dtype_begin();
  Dtype_iterator dtype_end();
  unsigned int get_ndtype();
  void dtype_clear_sync();
  void dtype_clear_rsync();
  bool dtype_is_sync();
  bool dtype_is_rsync();

  // Data stuff.
  void data_clear_sync();
  void data_clear_rsync();
  bool data_is_sync();
  bool data_is_rsync();

  // Ftype stuff.
  Ftype *ftype_add(Ftype *ftype);
  Ftype *ftype_add_cfunc(std::string name,
			 unsigned int nparam,
			 FUDA::CFUNC func,
			 void *fs=0,
			 unsigned int fs_size=0);
  Ftype *ftype_add_prod(std::string name,
			bool a_scale_flg,
			std::vector<Ftype*>& ftvec);
  Ftype *ftype_add_sum(std::string name,
		       bool a_scale_flg,
		       std::vector<Ftype*>& ftvec);
  Ftype *ftype_add_comp(std::string name,
			bool a_scale_flg,
			unsigned int p_index,
			Ftype *f, Ftype *g);
  Ftype *ftype_find(std::string& name);
  Ftype *ftype_find(const char name[]);
  Ftype_iterator ftype_begin();
  Ftype_iterator ftype_end();
  unsigned int get_nftype();

  // Func stuff.
  Func *func_add(std::string name,
                 std::string ftype_name,
		 std::string dtype_name,
                 std::list<std::string>& param_names);
  void func_del(std::string& name);
  void func_del(const char name[]);
  void func_del_all();
  Func *func_find(std::string& name);
  Func *func_find(const char name[]);
  void func_set_param(std::string fname,
		      unsigned int pindex,
		      std::string pname);
  Func_iterator func_begin();
  Func_iterator func_end();
  unsigned int get_nfunc();
  void func_clear_sync();
  void func_clear_rsync();
  bool func_is_sync();
  bool func_is_rsync();

  // eval stuff.
  bool eval_is_sync();
  bool eval_is_rsync();
  unsigned int eval_get_mod_count();
  void eval_clear(); // Clear the eval sync flag and structures.
  void eval_init(); // Initialize the eval structure.
  double eval_call(Dtype& dtype, std::vector<double>& xvec,
                   double pvec[], double dpvec[], bool dp_flg);
  double eval_call(unsigned int d_index, double pvec[],
                   double dpvec[], bool dp_flg);
  double eval_obs_calc(unsigned int d_index, double pvec[],
		       double dpvec[], bool dp_flg);
  void eval_data_init();
  void eval_data_recalc();
  void eval_data_random();
  double eval_enorm(double pvec[]);
  double eval_enorm();
  unsigned int eval_get_nfree(); // Get number of free parameters in fit.
  Param *eval_get_free(unsigned int i); // Return Pointer to free param.
  unsigned int eval_get_nexpl(); // Get number of expl. parameters in fit.
  Param *eval_get_expl(unsigned int i); // Return Pointer to explanatory param.
  unsigned int eval_get_nconst(); // Get # of const. (fixed) param. in fit.
  Param *eval_get_const(unsigned int i); // Return Pointer to constant param.
  unsigned int eval_get_ndata(); // Return the number of records in eval.d.
  Data *eval_get_data(unsigned int i); // Return Pointer to Data record.
  Dtype *eval_get_dtype(unsigned int i); // Return Pointer to Dtype record.

  // Misc.
  void get_version(std::string& version_string);
  unsigned int get_version_major();
  unsigned int get_version_minor();
  FuncCall *fcall_get_ref();
  void print();
};


class FuncCall
{
  bool defined;      // Set if p, dp,... has been allocated.
  double *p;         // double array of psize.
  double *dp;        // double array of psize.
  int *dp_flg;       // int array of psize.
  int *dp_flg_false; // int array of psize with zeros.
  unsigned int psize; /* >= # parmeters for ftype with
			largest # of params. */ 

  FuncCall();
  ~FuncCall();
  void set_size(unsigned int size);
  unsigned int get_size();

 public:
  friend class Fuda;
  friend class Func;
  friend class Ftype;
  friend class FtypeCFUNC;
  friend class FtypePROD;
  friend class FtypeCOMP;
  friend class FtypeSUM;
};



/************************************************************************/
/* Uferr. Uf error signalling namespace holding all exceptions created by
   uf routines. */
namespace Uferr
{
  struct PtypeNameAlreadyUsed {
    std::string name;
    PtypeNameAlreadyUsed(std::string n) : name(n) {}
  };
  struct PtypeReferencedByParam {
    std::string name;
    PtypeReferencedByParam(std::string n) : name(n) {}
  };
  struct ParamNameAlreadyUsed {
    std::string name;
    ParamNameAlreadyUsed(std::string n) : name(n) {}
  };
  struct PtypeNameInvalid {
    std::string name;
    PtypeNameInvalid(std::string n) : name(n) {}
  };
  struct DtypeNameAlreadyUsed {
    std::string name;
    DtypeNameAlreadyUsed(std::string n) : name(n) {}
  };
  struct DtypeReferencedByFunc {
    std::string name;
    DtypeReferencedByFunc(std::string n) : name(n) {}
  };
  struct FtypeNameAlreadyUsed {
    std::string name;
    FtypeNameAlreadyUsed(std::string n) : name(n) {}
  };
  struct FuncNameAlreadyUsed {
    std::string name;
    FuncNameAlreadyUsed(std::string n) : name(n) {}
  };
  struct ParamNameInvalid {
    std::string name;
    ParamNameInvalid(std::string n) : name(n) {}
  };
  struct ParamReferencedByFunc {
    std::string name;
    ParamReferencedByFunc(std::string n) : name(n) {}
  };
  struct FuncNameInvalid {
    std::string name;
    FuncNameInvalid(std::string n) : name(n) {}
  };
  struct FtypeNameInvalid {
    std::string name;
    FtypeNameInvalid(std::string n) : name(n) {}
  };
  struct FtypeFudaRefMismatch {};
  struct FtypeEmptyFtypeVec {};
  struct FtypeRefInvalid {};
  struct DtypeNameInvalid {
    std::string name;
    DtypeNameInvalid(std::string n) : name(n) {}
  };
  struct FtypeCallError {
    std::string msg;
    FtypeCallError(std::string m) : msg(m) {}
  };
  struct NparamInvalid {
    unsigned int nparam;
    NparamInvalid(unsigned int n) : nparam(n) {}
  };
  struct DtypeExplVarMissing {
    std::string name;
    DtypeExplVarMissing(std::string n) : name(n) {}
  };
  struct ExplVarInvalid {
    std::string name;
    ExplVarInvalid(std::string n) : name(n) {}
  };
  struct ParamNotExplanatory {
    std::string name;
    ParamNotExplanatory(std::string n) : name(n) {}
  };
  struct NoData {};
  struct NoFreeParam {};
  struct NparamGtNdata {
    unsigned int nparam, ndata;
    NparamGtNdata(unsigned int np, unsigned int nd) : nparam(np), ndata(nd) {}
  };
  struct ParamIndexInvalid {
    unsigned int i;
    ParamIndexInvalid(unsigned int pi) : i(pi) {}
  };
  struct DtypeExplIndexInvalid {
    unsigned int i;
    DtypeExplIndexInvalid(unsigned int pi) : i(pi) {}
  };
  struct ParamNotConst {
    unsigned int i;
    ParamNotConst(unsigned int pi) : i(pi) {}
    ParamNotConst() : i(0) {}
  };
  struct ParamNotVar {
    unsigned int i;
    ParamNotVar(unsigned int pi) : i(pi) {}
    ParamNotVar() : i(0) {}
  };
  struct DataUInvalid {
    double u;
    DataUInvalid(double datau) : u(datau) {}
  };
  struct EvalNotSync {};
  struct ParamKindInvalid {};
  struct NumParamDescrInvalid {
    unsigned int n;
    NumParamDescrInvalid(unsigned int a_n) : n(a_n) {}
  };
  struct VarIndexInvalid {
    unsigned int i;
    VarIndexInvalid(unsigned int a_i) : i(a_i) {}
  };
  struct VarIndexSizeInvalid {
    unsigned int size;
    VarIndexSizeInvalid(unsigned int a_size) : size(a_size) {}
  };
  struct NumVarInvalid {
    unsigned int n;
    NumVarInvalid(unsigned int a_n) : n(a_n) {}
  };
  struct NumParamInvalid {
    unsigned int n;
    NumParamInvalid(unsigned int a_n) : n(a_n) {}
  };
  struct EvalFreeIndexInvalid {
    unsigned int i;
    EvalFreeIndexInvalid(unsigned int a_i) : i(a_i) {}
  };
  struct EvalExplIndexInvalid {
    unsigned int i;
    EvalExplIndexInvalid(unsigned int a_i) : i(a_i) {}
  };
  struct EvalConstIndexInvalid {
    unsigned int i;
    EvalConstIndexInvalid(unsigned int a_i) : i(a_i) {}
  };
  struct EvalDataIndexInvalid {
    unsigned int i;
    EvalDataIndexInvalid(unsigned int a_i) : i(a_i) {}
  };
  struct DataIndexInvalid {
    unsigned int i;
    DataIndexInvalid(unsigned int a_i) : i(a_i) {}
  };
  struct EvalCallFuncError {
    std::string name;
    unsigned int value;
    EvalCallFuncError(std::string n, unsigned int i) : name(n), value(i) {}
  };
  struct FuncCallInvalSize {
    unsigned int i;
    FuncCallInvalSize(unsigned int a_i) : i(a_i) {}
  };
  struct NumExplInvalid {
    unsigned int n;
    NumExplInvalid(unsigned int a_n) : n(a_n) {}
  };
  struct DataDimMismatch {};
}

