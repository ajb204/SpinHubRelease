# fuda module: FUnction and Data Analysis tool.

# fuda [FUnction and Data Analysis tool] is a user interface to the
# fudalib module. It wraps all the functions in the fudalib module and
# supplement with some new functions.

import fudalib

######################################################################
#
# Wrapper functions for fudalib.
#
######################################################################

# Param functions.
def param(*arguments,**keywords):
    return fudalib.param(*arguments, **keywords)

def param_del(*arguments,**keywords):
    return fudalib.param_del(*arguments, **keywords)

def param_del_all(*arguments,**keywords):
    return fudalib.param_del_all(*arguments, **keywords)

def param_get(*arguments,**keywords):
    return fudalib.param_get(*arguments, **keywords)

def param_init_value(*arguments,**keywords):
    return fudalib.param_init_value(*arguments, **keywords)

def param_get_all(*arguments,**keywords):
    return fudalib.param_get_all(*arguments, **keywords)

def param_exists(*arguments,**keywords):
    return fudalib.param_exists(*arguments, **keywords)

def param_is_referenced(*arguments,**keywords):
    return fudalib.param_is_referenced(*arguments, **keywords)

# Ptype functions.
def ptype(*arguments,**keywords):
    return fudalib.ptype(*arguments, **keywords)

def ptype_get_all(*arguments,**keywords):
    return fudalib.ptype_get_all(*arguments, **keywords)

def ptype_del(*arguments,**keywords):
    return fudalib.ptype_del(*arguments, **keywords)

def ptype_del_all(*arguments,**keywords):
    return fudalib.ptype_del_all(*arguments, **keywords)

def ptype_exists(*arguments,**keywords):
    return fudalib.ptype_exists(*arguments, **keywords)

def ptype_is_referenced(*arguments,**keywords):
    return fudalib.ptype_is_referenced(*arguments, **keywords)

# Ftype functions.
def ftype_python(*arguments,**keywords):
    return fudalib.ftype_python(*arguments, **keywords)

def ftype_product(*arguments,**keywords):
    return fudalib.ftype_product(*arguments, **keywords)

def ftype_sum(*arguments,**keywords):
    return fudalib.ftype_sum(*arguments, **keywords)

def ftype_composite(*arguments,**keywords):
    return fudalib.ftype_composite(*arguments, **keywords)

def ftype_get_all(*arguments,**keywords):
    return fudalib.ftype_get_all(*arguments, **keywords)

def ftype_exists(*arguments,**keywords):
    return fudalib.ftype_exists(*arguments, **keywords)

def ftype_call(*arguments,**keywords):
    return fudalib.ftype_call(*arguments, **keywords)

def ftype_get_param(*arguments,**keywords):
    return fudalib.ftype_get_param(*arguments, **keywords)

def ftype_get_param_descr(*arguments,**keywords):
    return fudalib.ftype_get_param_descr(*arguments, **keywords)

def ftype_get_var(*arguments,**keywords):
    return fudalib.ftype_get_var(*arguments, **keywords)

def ftype_get_var_index(*arguments,**keywords):
    return fudalib.ftype_get_var_index(*arguments, **keywords)

def ftype_get(*arguments,**keywords):
    return fudalib.ftype_get(*arguments, **keywords)

# Func functions.
def func(*arguments,**keywords):
    return fudalib.func(*arguments, **keywords)

def func_del(*arguments,**keywords):
    return fudalib.func_del(*arguments, **keywords)

def func_del_all(*arguments,**keywords):
    return fudalib.func_del_all(*arguments, **keywords)

def func_get(*arguments,**keywords):
    return fudalib.func_get(*arguments, **keywords)

def func_exists(*arguments,**keywords):
    return fudalib.func_exists(*arguments, **keywords)

def func_get_all(*arguments,**keywords):
    return fudalib.func_get_all(*arguments, **keywords)

def func_set_param(*arguments,**keywords):
    return fudalib.func_set_param(*arguments, **keywords)

def func_get_param(*arguments,**keywords):
    return fudalib.func_get_param(*arguments, **keywords)

def func_get_var(*arguments,**keywords):
    return fudalib.func_get_var(*arguments, **keywords)

def func_call(*arguments,**keywords):
    return fudalib.func_call(*arguments, **keywords)

def func_call_by_expl(*arguments,**keywords):
    return fudalib.func_call_by_expl(*arguments, **keywords)

def func_call_by_var(*arguments,**keywords):
    return fudalib.func_call_by_var(*arguments, **keywords)

def func_deriv_by_var(*arguments,**keywords):
    return fudalib.func_deriv_by_var(*arguments, **keywords)

def func_call_by_param(*arguments,**keywords):
    return fudalib.func_call_by_param(*arguments, **keywords)

def func_get_current(*arguments,**keywords):
    return fudalib.func_get_current(*arguments, **keywords)

def func_current(*arguments,**keywords):
    return fudalib.func_current(*arguments, **keywords)

# Dtype functions.
def dtype(*arguments,**keywords):
    return fudalib.dtype(*arguments, **keywords)

def dtype_del(*arguments,**keywords):
    return fudalib.dtype_del(*arguments, **keywords)

def dtype_del_all(*arguments,**keywords):
    return fudalib.dtype_del_all(*arguments, **keywords)

def dtype_get(*arguments,**keywords):
    return fudalib.dtype_get(*arguments, **keywords)

def dtype_exists(*arguments,**keywords):
    return fudalib.dtype_exists(*arguments, **keywords)

def dtype_is_referenced(*arguments,**keywords):
    return fudalib.dtype_is_referenced(*arguments, **keywords)

def dtype_get_all(*arguments,**keywords):
    return fudalib.dtype_get_all(*arguments, **keywords)

def dtype_set_purge(*arguments,**keywords):
    return fudalib.dtype_set_purge(*arguments, **keywords)

def dtype_set_purge_radius(*arguments,**keywords):
    return fudalib.dtype_set_purge_radius(*arguments, **keywords)

def dtype_get_purge_radius(*arguments,**keywords):
    return fudalib.dtype_get_purge_radius(*arguments, **keywords)

# Data functions.
def data(*arguments,**keywords):
    return fudalib.data(*arguments, **keywords)

def data_del(*arguments,**keywords):
    return fudalib.data_del(*arguments, **keywords)

def data_del_all(*arguments,**keywords):
    return fudalib.data_del_all(*arguments, **keywords)

# Eval functions.
def eval_init(*arguments,**keywords):
    return fudalib.eval_init(*arguments, **keywords)

def eval_get(*arguments,**keywords):
    return fudalib.eval_get(*arguments, **keywords)

def eval_get_free(*arguments,**keywords):
    return fudalib.eval_get_free(*arguments, **keywords)

def eval_get_expl(*arguments,**keywords):
    return fudalib.eval_get_expl(*arguments, **keywords)

def eval_get_const(*arguments,**keywords):
    return fudalib.eval_get_const(*arguments, **keywords)

def eval_get_data(*arguments,**keywords):
    return fudalib.eval_get_data(*arguments, **keywords)

def eval_get_dtype(*arguments,**keywords):
    return fudalib.eval_get_dtype(*arguments, **keywords)

def eval_get_func(*arguments,**keywords):
    return fudalib.eval_get_func(*arguments, **keywords)

def eval_enorm(*arguments,**keywords):
    return fudalib.eval_enorm(*arguments, **keywords)

def eval_data_recalc(*arguments,**keywords):
    return fudalib.eval_data_recalc(*arguments, **keywords)

def eval_data_random(*arguments,**keywords):
    return fudalib.eval_data_random(*arguments, **keywords)

def eval_call(*arguments,**keywords):
    return fudalib.eval_call(*arguments, **keywords)

def eval_deriv(*arguments,**keywords):
    return fudalib.eval_deriv(*arguments, **keywords)

# lm functions.
def lm_minimize(*arguments,**keywords):
    return fudalib.lm_minimize(*arguments, **keywords)

def lm(*arguments,**keywords):
    return fudalib.lm(*arguments, **keywords)

def lm_get(*arguments,**keywords):
    return fudalib.lm_get(*arguments, **keywords)

def lm_update_param(*arguments,**keywords):
    return fudalib.lm_update_param(*arguments, **keywords)

def lm_get_value(*arguments,**keywords):
    return fudalib.lm_get_value(*arguments, **keywords)

def lm_get_esd(*arguments,**keywords):
    return fudalib.lm_get_esd(*arguments, **keywords)

def lm_get_covar(*arguments,**keywords):
    return fudalib.lm_get_covar(*arguments, **keywords)

# Misc functions.
def rand_set_seed(*arguments,**keywords):
    return fudalib.rand_set_seed(*arguments, **keywords)

def rand_get_seed(*arguments,**keywords):
    return fudalib.rand_get_seed(*arguments, **keywords)

def rand_gauss(*arguments,**keywords):
    return fudalib.rand_gauss(*arguments, **keywords)

def rand_uniform(*arguments,**keywords):
    return fudalib.rand_uniform(*arguments, **keywords)

def z_distrib_p(*arguments,**keywords):
    return fudalib.z_distrib_p(*arguments, **keywords)

def z_distrib_crit(*arguments,**keywords):
    return fudalib.z_distrib_crit(*arguments, **keywords)

def chi2_distrib_p(*arguments,**keywords):
    return fudalib.chi2_distrib_p(*arguments, **keywords)

def chi2_distrib_crit(*arguments,**keywords):
    return fudalib.chi2_distrib_crit(*arguments, **keywords)

def f_distrib_p(*arguments,**keywords):
    return fudalib.f_distrib_p(*arguments, **keywords)

def f_distrib_crit(*arguments,**keywords):
    return fudalib.f_distrib_crit(*arguments, **keywords)

def dump(*arguments,**keywords):
    return fudalib.fuda_print(*arguments, **keywords)

######################################################################
#
# Extension functions.
#
######################################################################

import sys, string, re
from types import *

# Fuda exception string raised on error.
error='fuda.error'


def func_use(name,use=1):
    return fudalib.func_use(name,use)

def func_use_all(use=1):
    return fudalib.func_use_all(use)

def func_disuse(name):
    return fudalib.func_use(name,0)

def func_disuse_all():
    return fudalib.func_use_all(0)

def func_dump():
    for fname in func_get_all():
        ftype,dtype,nparam,nvar,use,ndata=func_get(fname,'ftype',\
                                                   'dtype','nparam',\
                                                   'nvar','use','ndata')
        print()
        print('Function:   %s' % (fname,))
        print()
        print('ftype:      %s' % (ftype,))
        print('dtype:      %s' % (dtype,))
        print('nparm:      %6d' % (nparam,))
        print('nvar:       %6d' % (nvar,))
        print('use:        %6d' % (use,))
        print('ndata:      %6d' % (ndata,))
        print()
        print('Function parameters:')
        for ip in range(nparam): 
            p_name=func_get_param(fname,ip)
            p_free=param_get(p_name,'free')
            print('%3d  %-16s  %1d' % (ip,p_name,p_free))

def _lm_report_cprint_f(*args):
    for keyword in args:
        print('%-12s : %12.6g' % (keyword,lm_get(keyword)))

def _lm_report_cprint_d(*args):
    for keyword in args:
        print('%-12s : %12d' % (keyword,lm_get(keyword)))

def lm_report(type='result',headings=1,verbose=1):

    if type=='control' or type=='all':
        print()
        if headings:
            print('Levenberg-Marquardt (lm) control variables:')
        _lm_report_cprint_d('maxfev','nprint','numderiv')
        _lm_report_cprint_f('numderiv_eps','tol','ftol','xtol','gtol',\
                            'ctol','factor')

    if type=='status' or type=='all':
        print()
        if headings:
            print('Levenberg-Marquardt (lm) status flags:')
        _lm_report_cprint_f('minimized','fit_ok','fit_converged','sync')

    if type=='fit' or type=='all':
        print()
        if headings:
            print('Levenberg-Marquardt (lm) fit output variables:')
        if lm_get('fit_ok'):
            _lm_report_cprint_d('info','nfev','njev')
            _lm_report_cprint_f('enorm','sd')
        else:
            print('(No fit)')

    if type=='param' or type=='all':
        print()
        if headings:
            print('Levenberg-Marquardt (lm) estimated parameters:')
        if lm_get('fit_ok') and lm_get('sync'):
            print("%-10s %16s %16s" % ('Parameter','value','esd'))
            if lm_get('fit_converged'):
                for i in range(eval_get('nfree')):
                    print("%-16s %16.6g %16.6g" % \
                          (eval_get_free(i), lm_get_value(i), lm_get_esd(i)))
            else:
                # When not converged, we dont have esd's.
                for i in range(eval_get('nfree')):
                    print("%-16s %16.6g %16s" % \
                          (eval_get_free(i), lm_get_value(i), '*'*16))
                
        else:
            print('(No fit)')
                    
    if type=='result':
        print()
        print('Levenberg-Marquardt (lm) fit result:')
        print('------------------------------------')
        if lm_get('fit_ok'):

            print('Standard deviation (sd)                : %12.6g' %\
                  (lm_get('sd'),))
            print('Euclidian norm (enorm)                 : %12.6g' %\
                  (lm_get('enorm'),))
            print('Number of data points                  : %d' %\
                  (eval_get('ndata'),))
            print('Number of free parameters:             : %d' %\
                  (eval_get('nfree'),))
            print('Number of function evaluations (nfev)  : %d' %\
                  (lm_get('nfev'),))
            print('Number of derivative evaluations (njev): %d' %\
                  (lm_get('njev'),))
            if lm_get('fit_converged'):
                info=lm_get('info')
                if info==1:
                    print('Fit converged according to ftol')
                elif info==2:
                    print('Fit converged according to xtol')
                elif info==3:
                    print('fit converged according to ftol and xtol')
                elif info==4:
                    print('Fit converged according to gtol')
                else:
                    print('fuda.lm_report - fuda in panic')
                    sys.exit(1)
            else:
                print('Fit did not converge')
            lm_report('param')
        else:
            print('(No fit)')


def data_read_ascii(fname,cols):
    # Read data from file fname from the columns specified in cols.

    # Check that cols is a tuple of integers and save largest col in maxcol.
    if not (type(cols) != TupleType or type(cols) != ListType):
        raise error('2nd argument must be a sequence of integers')

    maxcol=0
    for col in cols:
        if type(col) != IntType:
            raise error('2nd argument must be a sequence of integers')
        if col>maxcol:
            col=maxcol

    # Get dimension of the current functions dtype and check with cols.
    dt=func_get(func_get_current(),'dtype')
    dt_dim=dtype_get(dt,'dim')

    # We need the explanatory parameters + value + uncertainty.
    if dt_dim+2!=len(cols):
        raise error('number of columns (%d) invalid. Should be %d' %\
              (len(cols),dt_dim+2))

    # Open the file.
    if fname=='-':
        # Read from standard input.
        f=sys.stdin
        file_flg=0
    else:
        # Read from file.
        f=open(fname,"r")
        file_flg=1

    # Loop over lines.
    comment=re.compile(r'^\s*[#]')
    trim_left=re.compile(r'^\s+')
    trim_right=re.compile(r'\s+$')
    separator=re.compile(r'\s+')
    lcount=0
    dcount=0
    line=f.readline()
    while line!='':
        lcount=lcount+1

        # We skip comments.
        if comment.search(line):
            pass
        else:
            # Trim leading and trailing white space.
            line=trim_left.sub('',line)
            line=trim_right.sub('',line)

            # Split in words.
            w=separator.split(line)

            # Do we have enough words.
            if len(w)-1<maxcol:
                if file_flg:
                    f.close()
                raise error('insufficient number of words on input line %d' %\
                      lcount)

            # Setup tuple of values to pass to data.
            data_args = []
            for col in cols:
                try:
                    data_args.append(float(w[col]))
                except ValueError:
                    raise error('invalid float in col %d input line %d: %s' %\
                          (col,lcount,w[col]))

            # Define data.
            data(*tuple(data_args))
            dcount=dcount+1

        # Read next line.
        line=f.readline()

    # Close file.
    if file_flg:
        f.close()

    # Return number of data points read.
    return dcount



# Declare default product ftypes.
ftype_product('poly',1,('norm_poly',))
ftype_product('exp',1,('norm_exp',))
ftype_product('pow',1,('norm_pow',))
ftype_product('exp_decay',1,('norm_exp_decay',))
ftype_product('cos',1,('norm_cos',))
ftype_product('lore1d',1,('norm_lore',))
ftype_product('lore2d',1,('norm_lore','norm_lore'))
ftype_product('gausslore1d',1,('norm_gausslore',))
ftype_product('gausslore2d',1,('norm_gausslore','norm_gausslore'))
ftype_product('dlore1d',1,('norm_dlore',))
ftype_product('dlore2d',1,('norm_dlore','norm_dlore'))

