/* toms611.f downloaded from http://people.sc.fsu.edu/~jburkardt/f77_src/toms611/toms611.f
 * converted using f2c -C++ toms611.f
 *
 * turned all static variable to thread-static
 */
/* toms611.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static __thread integer c__2 = 2;
static __thread integer c__1 = 1;
static __thread integer c__4 = 4;
static __thread integer c_n1 = -1;
static __thread integer c__3 = 3;
static __thread integer c__6 = 6;
static __thread integer c__5 = 5;

integer imdcon_(integer *k)
{
    /* Initialized data */

    static __thread integer mdperm[3] = { 2,4,1 };

    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer i1mach_(integer *);



/*  ***  return integer machine-dependent constants  *** */

/*     ***  k = 1 means return standard output unit number.   *** */
/*     ***  k = 2 means return alternate output unit number.  *** */
/*     ***  k = 3 means return  input unit number.            *** */
/*          (note -- k = 2, 3 are used only by test programs.) */

/*  +++  port version follows... */

    ret_val = i1mach_(&mdperm[(0 + (0 + (*k - 1 << 2))) / 4]);
/*  +++  end of port version  +++ */

/*  +++  non-port version follows... */
/*     integer mdcon(3) */
/*     data mdcon(1)/6/, mdcon(2)/8/, mdcon(3)/5/ */
/*     imdcon = mdcon(k) */
/*  +++  end of non-port version  +++ */

/* L999: */
    return ret_val;
/*  ***  last card of imdcon follows  *** */
} /* imdcon_ */

doublereal rmdcon_(integer *k)
{
    /* Initialized data */

    static __thread doublereal zero = 0.;
    static __thread struct {
	doublereal e_1;
	doublereal fill_2[1];
	} equiv_0 = { 0. };

    static __thread struct {
	doublereal e_1;
	doublereal fill_2[1];
	} equiv_1 = { 0. };

    static __thread struct {
	doublereal e_1;
	doublereal fill_2[1];
	} equiv_2 = { 0. };


    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
#define big ((doublereal *)&equiv_0)
#define eta ((doublereal *)&equiv_1)
#define bigi ((integer *)&equiv_0)
#define etai ((integer *)&equiv_1)
    extern doublereal d1mach_(integer *);
#define machei ((integer *)&equiv_2)
#define machep ((doublereal *)&equiv_2)


/*  ***  return machine dependent constants used by nl2sol  *** */

/* +++  comments below contain data statements for various machines.  +++ */
/* +++  to convert to another machine, place a c in column 1 of the   +++ */
/* +++  data statement line(s) that correspond to the current machine +++ */
/* +++  and remove the c from column 1 of the data statement line(s)  +++ */
/* +++  that correspond to the new machine.                           +++ */


/*  ***  the constant returned depends on k... */

/*  ***        k = 1... smallest pos. eta such that -eta exists. */
/*  ***        k = 2... square root of eta. */
/*  ***        k = 3... unit roundoff = smallest pos. no. machep such */
/*  ***                 that 1 + machep .gt. 1 .and. 1 - machep .lt. 1. */
/*  ***        k = 4... square root of machep. */
/*  ***        k = 5... square root of big (see k = 6). */
/*  ***        k = 6... largest machine no. big such that -big exists. */

/* /+ */
/* / */

/*  +++  ibm 360, ibm 370, or xerox  +++ */

/*     data big/z7fffffffffffffff/, eta/z0010000000000000/, */
/*    1     machep/z3410000000000000/ */

/*  +++  data general  +++ */

/*     data big/0.7237005577d+76/, eta/0.5397605347d-78/, */
/*    1     machep/2.22044605d-16/ */

/*  +++  dec 11  +++ */

/*     data big/1.7d+38/, eta/2.938735878d-39/, machep/2.775557562d-17/ */

/*  +++  hp3000  +++ */

/*     data big/1.157920892d+77/, eta/8.636168556d-78/, */
/*    1     machep/5.551115124d-17/ */

/*  +++  honeywell  +++ */

/*     data big/1.69d+38/, eta/5.9d-39/, machep/2.1680435d-19/ */

/*  +++  dec10  +++ */

/*     data big/"377777100000000000000000/, */
/*    1     eta/"002400400000000000000000/, */
/*    2     machep/"104400000000000000000000/ */

/*  +++  burroughs  +++ */

/*     data big/o0777777777777777,o7777777777777777/, */
/*    1     eta/o1771000000000000,o7770000000000000/, */
/*    2     machep/o1451000000000000,o0000000000000000/ */

/*  +++  control data  +++ */

/*     data big/37767777777777777777b,37167777777777777777b/, */
/*    1     eta/00014000000000000000b,00000000000000000000b/, */
/*    2     machep/15614000000000000000b,15010000000000000000b/ */

/*  +++  prime  +++ */

/*     data big/1.0d+9786/, eta/1.0d-9860/, machep/1.4210855d-14/ */

/*  +++  univac  +++ */

/*     data big/8.988d+307/, eta/1.2d-308/, machep/1.734723476d-18/ */

/*  +++  vax  +++ */

/*     data big/1.7d+38/, eta/2.939d-39/, machep/1.3877788d-17/ */

/*  +++  cray 1  +++ */

/*     data bigi(1)/577767777777777777777b/, */
/*    1     bigi(2)/000007777777777777776b/, */
/*    2     etai(1)/200004000000000000000b/, */
/*    3     etai(2)/000000000000000000000b/, */
/*    4     machei(1)/377224000000000000000b/, */
/*    5     machei(2)/000000000000000000000b/ */

/*  +++  port library -- requires more than just a data statement... +++ */

    if (*big > zero) {
	goto L1;
    }
    *big = d1mach_(&c__2);
    *eta = d1mach_(&c__1);
    *machep = d1mach_(&c__4);
L1:

/*  +++ end of port +++ */

/* -------------------------------  body  -------------------------------- */

    switch (*k) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	case 5:  goto L50;
	case 6:  goto L60;
    }

L10:
    ret_val = *eta;
    goto L999;

L20:
    ret_val = sqrt(*eta * 256.) / 16.;
    goto L999;

L30:
    ret_val = *machep;
    goto L999;

L40:
    ret_val = sqrt(*machep);
    goto L999;

L50:
    ret_val = sqrt(*big / 256.) * 16.;
    goto L999;

L60:
    ret_val = *big;

L999:
    return ret_val;
/*  ***  last card of rmdcon follows  *** */
} /* rmdcon_ */

#undef machep
#undef machei
#undef etai
#undef bigi
#undef eta
#undef big


/* Subroutine */ int sumsl_(integer *n, doublereal *d__, doublereal *x, S_fp 
	calcf, S_fp calcg, integer *iv, integer *liv, integer *lv, doublereal 
	*v, integer *uiparm, doublereal *urparm, U_fp ufparm)
{
    /* Initialized data */

    static __thread integer nfcall = 6;
    static __thread integer nfgcal = 7;
    static __thread integer g = 28;
    static __thread integer toobig = 2;
    static __thread integer vneed = 4;
    static __thread integer nextv = 47;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static __thread doublereal f;
    static __thread integer g1, nf, iv1;
    extern /* Subroutine */ int deflt_(integer *, integer *, integer *, 
	    integer *, doublereal *), sumit_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *);


/*  ***  minimize general unconstrained objective function using   *** */
/*  ***  analytic gradient and hessian approx. from secant update  *** */

/*     dimension v(71 + n*(n+15)/2), uiparm(*), urparm(*) */

/*  ***  purpose  *** */

/*        this routine interacts with subroutine  sumit  in an attempt */
/*     to find an n-vector  x*  that minimizes the (unconstrained) */
/*     objective function computed by  calcf.  (often the  x*  found is */
/*     a local minimizer rather than a global one.) */

/* --------------------------  parameter usage  -------------------------- */

/* n........ (input) the number of variables on which  f  depends, i.e., */
/*                  the number of components in  x. */
/* d........ (input/output) a scale vector such that  d(i)*x(i), */
/*                  i = 1,2,...,n,  are all in comparable units. */
/*                  d can strongly affect the behavior of sumsl. */
/*                  finding the best choice of d is generally a trial- */
/*                  and-error process.  choosing d so that d(i)*x(i) */
/*                  has about the same value for all i often works well. */
/*                  the defaults provided by subroutine deflt (see iv */
/*                  below) require the caller to supply d. */
/* x........ (input/output) before (initially) calling sumsl, the call- */
/*                  er should set  x  to an initial guess at  x*.  when */
/*                  sumsl returns,  x  contains the best point so far */
/*                  found, i.e., the one that gives the least value so */
/*                  far seen for  f(x). */
/* calcf.... (input) a subroutine that, given x, computes f(x).  calcf */
/*                  must be declared external in the calling program. */
/*                  it is invoked by */
/*                       call calcf(n, x, nf, f, uiparm, urparm, ufparm) */
/*                  when calcf is called, nf is the invocation */
/*                  count for calcf.  nf is included for possible use */
/*                  with calcg.  if x is out of bounds (e.g., if it */
/*                  would cause overflow in computing f(x)), then calcf */
/*                  should set nf to 0.  this will cause a shorter step */
/*                  to be attempted.  (if x is in bounds, then calcf */
/*                  should not change nf.)  the other parameters are as */
/*                  described above and below.  calcf should not change */
/*                  n, p, or x. */
/* calcg.... (input) a subroutine that, given x, computes g(x), the gra- */
/*                  dient of f at x.  calcg must be declared external in */
/*                  the calling program.  it is invoked by */
/*                       call calcg(n, x, nf, g, uiparm, urparm, ufaprm) */
/*                  when calcg is called, nf is the invocation */
/*                  count for calcf at the time f(x) was evaluated.  the */
/*                  x passed to calcg is usually the one passed to calcf */
/*                  on either its most recent invocation or the one */
/*                  prior to it.  if calcf saves intermediate results */
/*                  for use by calcg, then it is possible to tell from */
/*                  nf whether they are valid for the current x (or */
/*                  which copy is valid if two copies are kept).  if g */
/*                  cannot be computed at x, then calcg should set nf to */
/*                  0.  in this case, sumsl will return with iv(1) = 65. */
/*                  (if g can be computed at x, then calcg should not */
/*                  changed nf.)  the other parameters to calcg are as */
/*                  described above and below.  calcg should not change */
/*                  n or x. */
/* iv....... (input/output) an integer value array of length liv (see */
/*                  below) that helps control the sumsl algorithm and */
/*                  that is used to store various intermediate quanti- */
/*                  ties.  of particular interest are the initialization/ */
/*                  return code iv(1) and the entries in iv that control */
/*                  printing and limit the number of iterations and func- */
/*                  tion evaluations.  see the section on iv input */
/*                  values below. */
/* liv...... (input) length of iv array.  must be at least 60.  if liv */
/*                  is too small, then sumsl returns with iv(1) = 15. */
/*                  when sumsl returns, the smallest allowed value of */
/*                  liv is stored in iv(lastiv) -- see the section on */
/*                  iv output values below.  (this is intended for use */
/*                  with extensions of sumsl that handle constraints.) */
/* lv....... (input) length of v array.  must be at least 71+n*(n+15)/2. */
/*                  (at least 77+n*(n+17)/2 for smsno, at least */
/*                  78+n*(n+12) for humsl).  if lv is too small, then */
/*                  sumsl returns with iv(1) = 16.  when sumsl returns, */
/*                  the smallest allowed value of lv is stored in */
/*                  iv(lastv) -- see the section on iv output values */
/*                  below. */
/* v........ (input/output) a floating-point value array of length lv */
/*                  (see below) that helps control the sumsl algorithm */
/*                  and that is used to store various intermediate */
/*                  quantities.  of particular interest are the entries */
/*                  in v that limit the length of the first step */
/*                  attempted (lmax0) and specify convergence tolerances */
/*                  (afctol, lmaxs, rfctol, sctol, xctol, xftol). */
/* uiparm... (input) user integer parameter array passed without change */
/*                  to calcf and calcg. */
/* urparm... (input) user floating-point parameter array passed without */
/*                  change to calcf and calcg. */
/* ufparm... (input) user external subroutine or function passed without */
/*                  change to calcf and calcg. */

/*  ***  iv input values (from subroutine deflt)  *** */

/* iv(1)...  on input, iv(1) should have a value between 0 and 14...... */
/*             0 and 12 mean this is a fresh start.  0 means that */
/*                  deflt(2, iv, liv, lv, v) */
/*             is to be called to provide all default values to iv and */
/*             v.  12 (the value that deflt assigns to iv(1)) means the */
/*             caller has already called deflt and has possibly changed */
/*             some iv and/or v entries to non-default values. */
/*             13 means deflt has been called and that sumsl (and */
/*             sumit) should only do their storage allocation.  that is, */
/*             they should set the output components of iv that tell */
/*             where various subarrays arrays of v begin, such as iv(g) */
/*             (and, for humsl and humit only, iv(dtol)), and return. */
/*             14 means that a storage has been allocated (by a call */
/*             with iv(1) = 13) and that the algorithm should be */
/*             started.  when called with iv(1) = 13, sumsl returns */
/*             iv(1) = 14 unless liv or lv is too small (or n is not */
/*             positive).  default = 12. */
/* iv(inith).... iv(25) tells whether the hessian approximation h should */
/*             be initialized.  1 (the default) means sumit should */
/*             initialize h to the diagonal matrix whose i-th diagonal */
/*             element is d(i)**2.  0 means the caller has supplied a */
/*             cholesky factor  l  of the initial hessian approximation */
/*             h = l*(l**t)  in v, starting at v(iv(lmat)) = v(iv(42)) */
/*             (and stored compactly by rows).  note that iv(lmat) may */
/*             be initialized by calling sumsl with iv(1) = 13 (see */
/*             the iv(1) discussion above).  default = 1. */
/* iv(mxfcal)... iv(17) gives the maximum number of function evaluations */
/*             (calls on calcf) allowed.  if this number does not suf- */
/*             fice, then sumsl returns with iv(1) = 9.  default = 200. */
/* iv(mxiter)... iv(18) gives the maximum number of iterations allowed. */
/*             it also indirectly limits the number of gradient evalua- */
/*             tions (calls on calcg) to iv(mxiter) + 1.  if iv(mxiter) */
/*             iterations do not suffice, then sumsl returns with */
/*             iv(1) = 10.  default = 150. */
/* iv(outlev)... iv(19) controls the number and length of iteration sum- */
/*             mary lines printed (by itsum).  iv(outlev) = 0 means do */
/*             not print any summary lines.  otherwise, print a summary */
/*             line after each abs(iv(outlev)) iterations.  if iv(outlev) */
/*             is positive, then summary lines of length 78 (plus carri- */
/*             age control) are printed, including the following...  the */
/*             iteration and function evaluation counts, f = the current */
/*             function value, relative difference in function values */
/*             achieved by the latest step (i.e., reldf = (f0-v(f))/f01, */
/*             where f01 is the maximum of abs(v(f)) and abs(v(f0)) and */
/*             v(f0) is the function value from the previous itera- */
/*             tion), the relative function reduction predicted for the */
/*             step just taken (i.e., preldf = v(preduc) / f01, where */
/*             v(preduc) is described below), the scaled relative change */
/*             in x (see v(reldx) below), the step parameter for the */
/*             step just taken (stppar = 0 means a full newton step, */
/*             between 0 and 1 means a relaxed newton step, between 1 */
/*             and 2 means a double dogleg step, greater than 2 means */
/*             a scaled down cauchy step -- see subroutine dbldog), the */
/*             2-norm of the scale vector d times the step just taken */
/*             (see v(dstnrm) below), and npreldf, i.e., */
/*             v(nreduc)/f01, where v(nreduc) is described below -- if */
/*             npreldf is positive, then it is the relative function */
/*             reduction predicted for a newton step (one with */
/*             stppar = 0).  if npreldf is negative, then it is the */
/*             negative of the relative function reduction predicted */
/*             for a step computed with step bound v(lmaxs) for use in */
/*             testing for singular convergence. */
/*                  if iv(outlev) is negative, then lines of length 50 */
/*             are printed, including only the first 6 items listed */
/*             above (through reldx). */
/*             default = 1. */
/* iv(parprt)... iv(20) = 1 means print any nondefault v values on a */
/*             fresh start or any changed v values on a restart. */
/*             iv(parprt) = 0 means skip this printing.  default = 1. */
/* iv(prunit)... iv(21) is the output unit number on which all printing */
/*             is done.  iv(prunit) = 0 means suppress all printing. */
/*             default = standard output unit (unit 6 on most systems). */
/* iv(solprt)... iv(22) = 1 means print out the value of x returned (as */
/*             well as the gradient and the scale vector d). */
/*             iv(solprt) = 0 means skip this printing.  default = 1. */
/* iv(statpr)... iv(23) = 1 means print summary statistics upon return- */
/*             ing.  these consist of the function value, the scaled */
/*             relative change in x caused by the most recent step (see */
/*             v(reldx) below), the number of function and gradient */
/*             evaluations (calls on calcf and calcg), and the relative */
/*             function reductions predicted for the last step taken and */
/*             for a newton step (or perhaps a step bounded by v(lmaxs) */
/*             -- see the descriptions of preldf and npreldf under */
/*             iv(outlev) above). */
/*             iv(statpr) = 0 means skip this printing. */
/*             iv(statpr) = -1 means skip this printing as well as that */
/*             of the one-line termination reason message.  default = 1. */
/* iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d */
/*             (on a fresh start only).  iv(x0prt) = 0 means skip this */
/*             printing.  default = 1. */

/*  ***  (selected) iv output values  *** */

/* iv(1)........ on output, iv(1) is a return code.... */
/*             3 = x-convergence.  the scaled relative difference (see */
/*                  v(reldx)) between the current parameter vector x and */
/*                  a locally optimal parameter vector is very likely at */
/*                  most v(xctol). */
/*             4 = relative function convergence.  the relative differ- */
/*                  ence between the current function value and its lo- */
/*                  cally optimal value is very likely at most v(rfctol). */
/*             5 = both x- and relative function convergence (i.e., the */
/*                  conditions for iv(1) = 3 and iv(1) = 4 both hold). */
/*             6 = absolute function convergence.  the current function */
/*                  value is at most v(afctol) in absolute value. */
/*             7 = singular convergence.  the hessian near the current */
/*                  iterate appears to be singular or nearly so, and a */
/*                  step of length at most v(lmaxs) is unlikely to yield */
/*                  a relative function decrease of more than v(sctol). */
/*             8 = false convergence.  the iterates appear to be converg- */
/*                  ing to a noncritical point.  this may mean that the */
/*                  convergence tolerances (v(afctol), v(rfctol), */
/*                  v(xctol)) are too small for the accuracy to which */
/*                  the function and gradient are being computed, that */
/*                  there is an error in computing the gradient, or that */
/*                  the function or gradient is discontinuous near x. */
/*             9 = function evaluation limit reached without other con- */
/*                  vergence (see iv(mxfcal)). */
/*            10 = iteration limit reached without other convergence */
/*                  (see iv(mxiter)). */
/*            11 = stopx returned .true. (external interrupt).  see the */
/*                  usage notes below. */
/*            14 = storage has been allocated (after a call with */
/*                  iv(1) = 13). */
/*            17 = restart attempted with n changed. */
/*            18 = d has a negative component and iv(dtype) .le. 0. */
/*            19...43 = v(iv(1)) is out of range. */
/*            63 = f(x) cannot be computed at the initial x. */
/*            64 = bad parameters passed to assess (which should not */
/*                  occur). */
/*            65 = the gradient could not be computed at x (see calcg */
/*                  above). */
/*            67 = bad first parameter to deflt. */
/*            80 = iv(1) was out of range. */
/*            81 = n is not positive. */
/* iv(g)........ iv(28) is the starting subscript in v of the current */
/*             gradient vector (the one corresponding to x). */
/* iv(lastiv)... iv(44) is the least acceptable value of liv.  (it is */
/*             only set if liv is at least 44.) */
/* iv(lastv).... iv(45) is the least acceptable value of lv.  (it is */
/*             only set if liv is large enough, at least iv(lastiv).) */
/* iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e., */
/*             function evaluations). */
/* iv(ngcall)... iv(30) is the number of gradient evaluations (calls on */
/*             calcg). */
/* iv(niter).... iv(31) is the number of iterations performed. */

/*  ***  (selected) v input values (from subroutine deflt)  *** */

/* v(bias)..... v(43) is the bias parameter used in subroutine dbldog -- */
/*             see that subroutine for details.  default = 0.8. */
/* v(afctol)... v(31) is the absolute function convergence tolerance. */
/*             if sumsl finds a point where the function value is less */
/*             than v(afctol) in absolute value, and if sumsl does not */
/*             return with iv(1) = 3, 4, or 5, then it returns with */
/*             iv(1) = 6.  this test can be turned off by setting */
/*             v(afctol) to zero.  default = max(10**-20, machep**2), */
/*             where machep is the unit roundoff. */
/* v(dinit).... v(38), if nonnegative, is the value to which the scale */
/*             vector d is initialized.  default = -1. */
/* v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the */
/*             very first step that sumsl attempts.  this parameter can */
/*             markedly affect the performance of sumsl. */
/* v(lmaxs).... v(36) is used in testing for singular convergence -- if */
/*             the function reduction predicted for a step of length */
/*             bounded by v(lmaxs) is at most v(sctol) * abs(f0), where */
/*             f0  is the function value at the start of the current */
/*             iteration, and if sumsl does not return with iv(1) = 3, */
/*             4, 5, or 6, then it returns with iv(1) = 7.  default = 1. */
/* v(rfctol)... v(32) is the relative function convergence tolerance. */
/*             if the current model predicts a maximum possible function */
/*             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0) */
/*             at the start of the current iteration, where  f0  is the */
/*             then current function value, and if the last step attempt- */
/*             ed achieved no more than twice the predicted function */
/*             decrease, then sumsl returns with iv(1) = 4 (or 5). */
/*             default = max(10**-10, machep**(2/3)), where machep is */
/*             the unit roundoff. */
/* v(sctol).... v(37) is the singular convergence tolerance -- see the */
/*             description of v(lmaxs) above. */
/* v(tuner1)... v(26) helps decide when to check for false convergence. */
/*             this is done if the actual function decrease from the */
/*             current step is no more than v(tuner1) times its predict- */
/*             ed value.  default = 0.1. */
/* v(xctol).... v(33) is the x-convergence tolerance.  if a newton step */
/*             (see v(nreduc)) is tried that has v(reldx) .le. v(xctol) */
/*             and if this step yields at most twice the predicted func- */
/*             tion decrease, then sumsl returns with iv(1) = 3 (or 5). */
/*             (see the description of v(reldx) below.) */
/*             default = machep**0.5, where machep is the unit roundoff. */
/* v(xftol).... v(34) is the false convergence tolerance.  if a step is */
/*             tried that gives no more than v(tuner1) times the predict- */
/*             ed function decrease and that has v(reldx) .le. v(xftol), */
/*             and if sumsl does not return with iv(1) = 3, 4, 5, 6, or */
/*             7, then it returns with iv(1) = 8.  (see the description */
/*             of v(reldx) below.)  default = 100*machep, where */
/*             machep is the unit roundoff. */
/* v(*)........ deflt supplies to v a number of tuning constants, with */
/*             which it should ordinarily be unnecessary to tinker.  see */
/*             section 17 of version 2.2 of the nl2sol usage summary */
/*             (i.e., the appendix to ref. 1) for details on v(i), */
/*             i = decfac, incfac, phmnfc, phmxfc, rdfcmn, rdfcmx, */
/*             tuner2, tuner3, tuner4, tuner5. */

/*  ***  (selected) v output values  *** */

/* v(dgnorm)... v(1) is the 2-norm of (diag(d)**-1)*g, where g is the */
/*             most recently computed gradient. */
/* v(dstnrm)... v(2) is the 2-norm of diag(d)*step, where step is the */
/*             current step. */
/* v(f)........ v(10) is the current function value. */
/* v(f0)....... v(13) is the function value at the start of the current */
/*             iteration. */
/* v(nreduc)... v(6), if positive, is the maximum function reduction */
/*             possible according to the current model, i.e., the func- */
/*             tion reduction predicted for a newton step (i.e., */
/*             step = -h**-1 * g,  where  g  is the current gradient and */
/*             h is the current hessian approximation). */
/*                  if v(nreduc) is negative, then it is the negative of */
/*             the function reduction predicted for a step computed with */
/*             a step bound of v(lmaxs) for use in testing for singular */
/*             convergence. */
/* v(preduc)... v(7) is the function reduction predicted (by the current */
/*             quadratic model) for the current step.  this (divided by */
/*             v(f0)) is used in testing for relative function */
/*             convergence. */
/* v(reldx).... v(17) is the scaled relative change in x caused by the */
/*             current step, computed as */
/*                  max(abs(d(i)*(x(i)-x0(i)), 1 .le. i .le. p) / */
/*                     max(d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p), */
/*             where x = x0 + step. */

/* -------------------------------  notes  ------------------------------- */

/*  ***  algorithm notes  *** */

/*        this routine uses a hessian approximation computed from the */
/*     bfgs update (see ref 3).  only a cholesky factor of the hessian */
/*     approximation is stored, and this is updated using ideas from */
/*     ref. 4.  steps are computed by the double dogleg scheme described */
/*     in ref. 2.  the steps are assessed as in ref. 1. */

/*  ***  usage notes  *** */

/*        after a return with iv(1) .le. 11, it is possible to restart, */
/*     i.e., to change some of the iv and v input values described above */
/*     and continue the algorithm from the point where it was interrupt- */
/*     ed.  iv(1) should not be changed, nor should any entries of iv */
/*     and v other than the input values (those supplied by deflt). */
/*        those who do not wish to write a calcg which computes the */
/*     gradient analytically should call smsno rather than sumsl. */
/*     smsno uses finite differences to compute an approximate gradient. */
/*        those who would prefer to provide f and g (the function and */
/*     gradient) by reverse communication rather than by writing subrou- */
/*     tines calcf and calcg may call on sumit directly.  see the com- */
/*     ments at the beginning of sumit. */
/*        those who use sumsl interactively may wish to supply their */
/*     own stopx function, which should return .true. if the break key */
/*     has been pressed since stopx was last invoked.  this makes it */
/*     possible to externally interrupt sumsl (which will return with */
/*     iv(1) = 11 if stopx returns .true.). */
/*        storage for g is allocated at the end of v.  thus the caller */
/*     may make v longer than specified above and may allow calcg to use */
/*     elements of g beyond the first n as scratch storage. */

/*  ***  portability notes  *** */

/*        the sumsl distribution tape contains both single- and double- */
/*     precision versions of the sumsl source code, so it should be un- */
/*     necessary to change precisions. */
/*        only the functions imdcon and rmdcon contain machine-dependent */
/*     constants.  to change from one machine to another, it should */
/*     suffice to change the (few) relevant lines in these functions. */
/*        intrinsic functions are explicitly declared.  on certain com- */
/*     puters (e.g. univac), it may be necessary to comment out these */
/*     declarations.  so that this may be done automatically by a simple */
/*     program, such declarations are preceded by a comment having c/+ */
/*     in columns 1-3 and blanks in columns 4-72 and are followed by */
/*     a comment having c/ in columns 1 and 2 and blanks in columns 3-72. */
/*        the sumsl source code is expressed in 1966 ansi standard */
/*     fortran.  it may be converted to fortran 77 by commenting out all */
/*     lines that fall between a line having c/6 in columns 1-3 and a */
/*     line having c/7 in columns 1-3 and by removing (i.e., replacing */
/*     by a blank) the c in column 1 of the lines that follow the c/7 */
/*     line and precede a line having c/ in columns 1-2 and blanks in */
/*     columns 3-72.  these changes convert some data statements into */
/*     parameter statements, convert some variables from real to */
/*     character*4, and make the data statements that initialize these */
/*     variables use character strings delimited by primes instead */
/*     of hollerith constants.  (such variables and data statements */
/*     appear only in modules itsum and parck.  parameter statements */
/*     appear nearly everywhere.)  these changes also add save state- */
/*     ments for variables given machine-dependent constants by rmdcon. */

/*  ***  references  *** */

/* 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), algorithm 573 -- */
/*             an adaptive nonlinear least-squares algorithm, acm trans. */
/*             math. software 7, pp. 369-383. */

/* 2.  dennis, j.e., and mei, h.h.w. (1979), two new unconstrained opti- */
/*             mization algorithms which use function and gradient */
/*             values, j. optim. theory applic. 28, pp. 453-482. */

/* 3.  dennis, j.e., and more, j.j. (1977), quasi-newton methods, motiva- */
/*             tion and theory, siam rev. 19, pp. 46-89. */

/* 4.  goldfarb, d. (1976), factorized variable metric methods for uncon- */
/*             strained optimization, math. comput. 30, pp. 796-811. */

/*  ***  general  *** */

/*     coded by david m. gay (winter 1980).  revised summer 1982. */
/*     this subroutine was written in connection with research */
/*     supported in part by the national science foundation under */
/*     grants mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, */
/*     and mcs-7906671. */
/* . */

/* ----------------------------  declarations  --------------------------- */


/* deflt... supplies default iv and v input components. */
/* sumit... reverse-communication routine that carries out sumsl algo- */
/*             rithm. */


/*  ***  subscripts for iv   *** */


/* /6 */
    /* Parameter adjustments */
    --x;
    --d__;
    --iv;
    --v;
    --uiparm;
    --urparm;

    /* Function Body */
/* /7 */
/*     parameter (nextv=47, nfcall=6, nfgcal=7, g=28, toobig=2, vneed=4) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

    if (iv[1] == 0) {
	deflt_(&c__2, &iv[1], liv, lv, &v[1]);
    }
    iv1 = iv[1];
    if (iv1 == 12 || iv1 == 13) {
	iv[vneed] += *n;
    }
    if (iv1 == 14) {
	goto L10;
    }
    if (iv1 > 2 && iv1 < 12) {
	goto L10;
    }
    g1 = 1;
    if (iv1 == 12) {
	iv[1] = 13;
    }
    goto L20;

L10:
    g1 = iv[g];

L20:
    sumit_(&d__[1], &f, &v[g1], &iv[1], liv, lv, n, &v[1], &x[1]);
    if ((i__1 = iv[1] - 2) < 0) {
	goto L30;
    } else if (i__1 == 0) {
	goto L40;
    } else {
	goto L50;
    }

L30:
    nf = iv[nfcall];
    (*calcf)(n, &x[1], &nf, &f, &uiparm[1], &urparm[1], (U_fp)ufparm);
    if (nf <= 0) {
	iv[toobig] = 1;
    }
    goto L20;

L40:
    (*calcg)(n, &x[1], &iv[nfgcal], &v[g1], &uiparm[1], &urparm[1], (U_fp)
	    ufparm);
    goto L20;

L50:
    if (iv[1] != 14) {
	goto L999;
    }

/*  ***  storage allocation */

    iv[g] = iv[nextv];
    iv[nextv] = iv[g] + *n;
    if (iv1 != 13) {
	goto L10;
    }

L999:
    return 0;
/*  ***  last card of sumsl follows  *** */
} /* sumsl_ */

/* Subroutine */ int smsno_(integer *n, doublereal *d__, doublereal *x, S_fp 
	calcf, integer *iv, integer *liv, integer *lv, doublereal *v, integer 
	*uiparm, doublereal *urparm, U_fp ufparm)
{
    /* Initialized data */

    static __thread integer nfcall = 6;
    static __thread integer toobig = 2;

    static __thread integer nf;
    static __thread doublereal fx;
    extern /* Subroutine */ int snoit_(doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);


/*  ***  minimize general unconstrained objective function using */
/*  ***  finite-difference gradients and secant hessian approximations. */

/*     dimension v(77 + n*(n+17)/2), uiparm(*), urparm(*) */

/*  ***  purpose  *** */

/*        this routine interacts with subroutine  snoit  in an attempt */
/*     to find an n-vector  x*  that minimizes the (unconstrained) */
/*     objective function computed by  calcf.  (often the  x*  found is */
/*     a local minimizer rather than a global one.) */

/*  ***  parameters  *** */

/*        the parameters for smsno are the same as those for sumsl */
/*     (which see), except that calcg is omitted.  instead of calling */
/*     calcg to obtain the gradient of the objective function at x, */
/*     smsno calls sgrad2, which computes an approximation to the */
/*     gradient by finite (forward and central) differences using the */
/*     method of ref. 1.  the following input component is of interest */
/*     in this regard (and is not described in sumsl). */

/* v(eta0)..... v(42) is an estimated bound on the relative error in the */
/*             objective function value computed by calcf... */
/*                  (true value) = (computed value) * (1 + e), */
/*             where abs(e) .le. v(eta0).  default = machep * 10**3, */
/*             where machep is the unit roundoff. */

/*        the output values iv(nfcall) and iv(ngcall) have different */
/*     meanings for smsno than for sumsl... */

/* iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e., */
/*             function evaluations) excluding those made only for */
/*             computing gradients.  the input value iv(mxfcal) is a */
/*             limit on iv(nfcall). */
/* iv(ngcall)... iv(30) is the number of function evaluations made only */
/*             for computing gradients.  the total number of function */
/*             evaluations is thus  iv(nfcall) + iv(ngcall). */

/*  ***  reference  *** */

/* 1. stewart, g.w. (1967), a modification of davidon*s minimization */
/*        method to accept difference approximations of derivatives, */
/*        j. assoc. comput. mach. 14, pp. 72-83. */
/* . */
/*  ***  general  *** */

/*     coded by david m. gay (winter 1980).  revised sept. 1982. */
/*     this subroutine was written in connection with research */
/*     supported in part by the national science foundation under */
/*     grants mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, */
/*     and mcs-7906671. */


/* ----------------------------  declarations  --------------------------- */


/* snoit.... oversees computation of finite-difference gradient and */
/*         calls sumit to carry out sumsl algorithm. */


/*  ***  subscripts for iv   *** */


/* /6 */
    /* Parameter adjustments */
    --x;
    --d__;
    --iv;
    --v;
    --uiparm;
    --urparm;

    /* Function Body */
/* /7 */
/*     parameter (nfcall=6, toobig=2) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

L10:
    snoit_(&d__[1], &fx, &iv[1], liv, lv, n, &v[1], &x[1]);
    if (iv[1] > 2) {
	goto L999;
    }

/*     ***  compute function  *** */

    nf = iv[nfcall];
    (*calcf)(n, &x[1], &nf, &fx, &uiparm[1], &urparm[1], (U_fp)ufparm);
    if (nf <= 0) {
	iv[toobig] = 1;
    }
    goto L10;


L999:
    return 0;
/*  ***  last card of smsno follows  *** */
} /* smsno_ */

/* Subroutine */ int sumit_(doublereal *d__, doublereal *fx, doublereal *g, 
	integer *iv, integer *liv, integer *lv, integer *n, doublereal *v, 
	doublereal *x)
{
    /* Initialized data */

    static __thread integer cnvcod = 55;
    static __thread integer dg = 37;
    static __thread integer g0 = 48;
    static __thread integer inith = 25;
    static __thread integer irc = 29;
    static __thread integer kagqt = 33;
    static __thread integer mode = 35;
    static __thread integer model = 5;
    static __thread integer mxfcal = 17;
    static __thread integer mxiter = 18;
    static __thread integer nfcall = 6;
    static __thread integer nfgcal = 7;
    static __thread integer ngcall = 30;
    static __thread integer niter = 31;
    static __thread integer nwtstp = 34;
    static __thread integer radinc = 8;
    static __thread integer restor = 9;
    static __thread integer step = 40;
    static __thread integer stglim = 11;
    static __thread integer stlstg = 41;
    static __thread integer toobig = 2;
    static __thread integer vneed = 4;
    static __thread integer xirc = 13;
    static __thread integer x0 = 43;
    static __thread integer dgnorm = 1;
    static __thread integer dinit = 38;
    static __thread integer dstnrm = 2;
    static __thread integer dst0 = 3;
    static __thread integer f = 10;
    static __thread integer f0 = 13;
    static __thread integer fdif = 11;
    static __thread integer gthg = 44;
    static __thread integer gtstep = 4;
    static __thread integer incfac = 23;
    static __thread integer lmat = 42;
    static __thread integer lmax0 = 35;
    static __thread integer lmaxs = 36;
    static __thread integer nextv = 47;
    static __thread integer nreduc = 6;
    static __thread integer preduc = 7;
    static __thread integer radfac = 16;
    static __thread integer radius = 8;
    static __thread integer rad0 = 9;
    static __thread integer reldx = 17;
    static __thread integer tuner4 = 29;
    static __thread integer tuner5 = 30;
    static __thread doublereal half = .5;
    static __thread doublereal negone = -1.;
    static __thread doublereal one = 1.;
    static __thread doublereal onep2 = 1.2;
    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static __thread integer i__, k, l;
    static __thread doublereal t;
    static __thread integer w, z__, g01, x01, dg1, temp1, step1;
    extern /* Subroutine */ int dbdog_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), deflt_(integer *, 
	    integer *, integer *, integer *, doublereal *), parck_(integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *);
    static __thread integer dummy;
    extern /* Subroutine */ int assst_(integer *, integer *, integer *, 
	    doublereal *), lvmul_(integer *, doublereal *, doublereal *, 
	    doublereal *), vcopy_(integer *, doublereal *, doublereal *), 
	    itsum_(doublereal *, doublereal *, integer *, integer *, integer *
	    , integer *, doublereal *, doublereal *), vaxpy_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern logical stopx_(integer *);
    extern doublereal v2norm_(integer *, doublereal *);
    static __thread integer nwtst1;
    extern /* Subroutine */ int lupdat_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *);
    extern doublereal dotprd_(integer *, doublereal *, doublereal *), reldst_(
	    integer *, doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int wzbfgs_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *), livmul_(integer *, 
	    doublereal *, doublereal *, doublereal *);
    static __thread integer lstgst;
    extern /* Subroutine */ int litvmu_(integer *, doublereal *, doublereal *,
	     doublereal *), ltvmul_(integer *, doublereal *, doublereal *, 
	    doublereal *), vscopy_(integer *, doublereal *, doublereal *), 
	    vvmulp_(integer *, doublereal *, doublereal *, doublereal *, 
	    integer *);


/*  ***  carry out sumsl (unconstrained minimization) iterations, using */
/*  ***  double-dogleg/bfgs steps. */

/*  ***  parameter declarations  *** */


/* --------------------------  parameter usage  -------------------------- */

/* d.... scale vector. */
/* fx... function value. */
/* g.... gradient vector. */
/* iv... integer value array. */
/* liv.. length of iv (at least 60). */
/* lv... length of v (at least 71 + n*(n+13)/2). */
/* n.... number of variables (components in x and g). */
/* v.... floating-point value array. */
/* x.... vector of parameters to be optimized. */

/*  ***  discussion  *** */

/*        parameters iv, n, v, and x are the same as the corresponding */
/*     ones to sumsl (which see), except that v can be shorter (since */
/*     the part of v that sumsl uses for storing g is not needed). */
/*     moreover, compared with sumsl, iv(1) may have the two additional */
/*     output values 1 and 2, which are explained below, as is the use */
/*     of iv(toobig) and iv(nfgcal).  the value iv(g), which is an */
/*     output value from sumsl (and smsno), is not referenced by */
/*     sumit or the subroutines it calls. */
/*        fx and g need not have been initialized when sumit is called */
/*     with iv(1) = 12, 13, or 14. */

/* iv(1) = 1 means the caller should set fx to f(x), the function value */
/*             at x, and call sumit again, having changed none of the */
/*             other parameters.  an exception occurs if f(x) cannot be */
/*             (e.g. if overflow would occur), which may happen because */
/*             of an oversized step.  in this case the caller should set */
/*             iv(toobig) = iv(2) to 1, which will cause sumit to ig- */
/*             nore fx and try a smaller step.  the parameter nf that */
/*             sumsl passes to calcf (for possible use by calcg) is a */
/*             copy of iv(nfcall) = iv(6). */
/* iv(1) = 2 means the caller should set g to g(x), the gradient vector */
/*             of f at x, and call sumit again, having changed none of */
/*             the other parameters except possibly the scale vector d */
/*             when iv(dtype) = 0.  the parameter nf that sumsl passes */
/*             to calcg is iv(nfgcal) = iv(7).  if g(x) cannot be */
/*             evaluated, then the caller may set iv(nfgcal) to 0, in */
/*             which case sumit will return with iv(1) = 65. */
/* . */
/*  ***  general  *** */

/*     coded by david m. gay (december 1979).  revised sept. 1982. */
/*     this subroutine was written in connection with research supported */
/*     in part by the national science foundation under grants */
/*     mcs-7600324 and mcs-7906671. */

/*        (see sumsl for references.) */

/* +++++++++++++++++++++++++++  declarations  ++++++++++++++++++++++++++++ */

/*  ***  local variables  *** */


/*     ***  constants  *** */


/*  ***  no intrinsic functions  *** */

/*  ***  external functions and subroutines  *** */


/* assst.... assesses candidate step. */
/* dbdog.... computes double-dogleg (candidate) step. */
/* deflt.... supplies default iv and v input components. */
/* dotprd... returns inner product of two vectors. */
/* itsum.... prints iteration summary and info on initial and final x. */
/* litvmu... multiplies inverse transpose of lower triangle times vector. */
/* livmul... multiplies inverse of lower triangle times vector. */
/* ltvmul... multiplies transpose of lower triangle times vector. */
/* lupdt.... updates cholesky factor of hessian approximation. */
/* lvmul.... multiplies lower triangle times vector. */
/* parck.... checks validity of input iv and v values. */
/* reldst... computes v(reldx) = relative step size. */
/* stopx.... returns .true. if the break key has been pressed. */
/* vaxpy.... computes scalar times one vector plus another. */
/* vcopy.... copies one vector to another. */
/* vscopy... sets all elements of a vector to a scalar. */
/* vvmulp... multiplies vector by vector raised to power (componentwise). */
/* v2norm... returns the 2-norm of a vector. */
/* wzbfgs... computes w and z for lupdat corresponding to bfgs update. */

/*  ***  subscripts for iv and v  *** */


/*  ***  iv subscript values  *** */

/* /6 */
    /* Parameter adjustments */
    --iv;
    --v;
    --x;
    --g;
    --d__;

    /* Function Body */
/* /7 */
/*     parameter (cnvcod=55, dg=37, g0=48, inith=25, irc=29, kagqt=33, */
/*    1           mode=35, model=5, mxfcal=17, mxiter=18, nfcall=6, */
/*    2           nfgcal=7, ngcall=30, niter=31, nwtstp=34, radinc=8, */
/*    3           restor=9, step=40, stglim=11, stlstg=41, toobig=2, */
/*    4           vneed=4, xirc=13, x0=43) */
/* / */

/*  ***  v subscript values  *** */

/* /6 */
/* /7 */
/*     parameter (dgnorm=1, dinit=38, dstnrm=2, dst0=3, f=10, f0=13, */
/*    1           fdif=11, gthg=44, gtstep=4, incfac=23, lmat=42, */
/*    2           lmax0=35, lmaxs=36, nextv=47, nreduc=6, preduc=7, */
/*    3           radfac=16, radius=8, rad0=9, reldx=17, tuner4=29, */
/*    4           tuner5=30) */
/* / */

/* /6 */
/* /7 */
/*     parameter (half=0.5d+0, negone=-1.d+0, one=1.d+0, onep2=1.2d+0, */
/*    1           zero=0.d+0) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

    i__ = iv[1];
    if (i__ == 1) {
	goto L50;
    }
    if (i__ == 2) {
	goto L60;
    }

/*  ***  check validity of iv and v input values  *** */

    if (iv[1] == 0) {
	deflt_(&c__2, &iv[1], liv, lv, &v[1]);
    }
    if (iv[1] == 12 || iv[1] == 13) {
	iv[vneed] += *n * (*n + 13) / 2;
    }
    parck_(&c__2, &d__[1], &iv[1], liv, lv, n, &v[1]);
    i__ = iv[1] - 2;
    if (i__ > 12) {
	goto L999;
    }
    switch (i__) {
	case 1:  goto L180;
	case 2:  goto L180;
	case 3:  goto L180;
	case 4:  goto L180;
	case 5:  goto L180;
	case 6:  goto L180;
	case 7:  goto L120;
	case 8:  goto L90;
	case 9:  goto L120;
	case 10:  goto L10;
	case 11:  goto L10;
	case 12:  goto L20;
    }

/*  ***  storage allocation  *** */

L10:
    l = iv[lmat];
    iv[x0] = l + *n * (*n + 1) / 2;
    iv[step] = iv[x0] + *n;
    iv[stlstg] = iv[step] + *n;
    iv[g0] = iv[stlstg] + *n;
    iv[nwtstp] = iv[g0] + *n;
    iv[dg] = iv[nwtstp] + *n;
    iv[nextv] = iv[dg] + *n;
    if (iv[1] != 13) {
	goto L20;
    }
    iv[1] = 14;
    goto L999;

/*  ***  initialization  *** */

L20:
    iv[niter] = 0;
    iv[nfcall] = 1;
    iv[ngcall] = 1;
    iv[nfgcal] = 1;
    iv[mode] = -1;
    iv[model] = 1;
    iv[stglim] = 1;
    iv[toobig] = 0;
    iv[cnvcod] = 0;
    iv[radinc] = 0;
    v[rad0] = zero;
    if (v[dinit] >= zero) {
	vscopy_(n, &d__[1], &v[dinit]);
    }
    if (iv[inith] != 1) {
	goto L40;
    }

/*     ***  set the initial hessian approximation to diag(d)**-2  *** */

    l = iv[lmat];
    i__1 = *n * (*n + 1) / 2;
    vscopy_(&i__1, &v[l], &zero);
    k = l - 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k += i__;
	t = d__[i__];
	if (t <= zero) {
	    t = one;
	}
	v[k] = t;
/* L30: */
    }

/*  ***  compute initial function value  *** */

L40:
    iv[1] = 1;
    goto L999;

L50:
    v[f] = *fx;
    if (iv[mode] >= 0) {
	goto L180;
    }
    iv[1] = 2;
    if (iv[toobig] == 0) {
	goto L999;
    }
    iv[1] = 63;
    goto L300;

/*  ***  make sure gradient could be computed  *** */

L60:
    if (iv[nfgcal] != 0) {
	goto L70;
    }
    iv[1] = 65;
    goto L300;

L70:
    dg1 = iv[dg];
    vvmulp_(n, &v[dg1], &g[1], &d__[1], &c_n1);
    v[dgnorm] = v2norm_(n, &v[dg1]);

    if (iv[cnvcod] != 0) {
	goto L290;
    }
    if (iv[mode] == 0) {
	goto L250;
    }

/*  ***  allow first step to have scaled 2-norm at most v(lmax0)  *** */

    v[radius] = v[lmax0];

    iv[mode] = 0;


/* -----------------------------  main loop  ----------------------------- */


/*  ***  print iteration summary, check iteration limit  *** */

L80:
    itsum_(&d__[1], &g[1], &iv[1], liv, lv, n, &v[1], &x[1]);
L90:
    k = iv[niter];
    if (k < iv[mxiter]) {
	goto L100;
    }
    iv[1] = 10;
    goto L300;

/*  ***  update radius  *** */

L100:
    iv[niter] = k + 1;
    if (k > 0) {
	v[radius] = v[radfac] * v[dstnrm];
    }

/*  ***  initialize for start of next iteration  *** */

    g01 = iv[g0];
    x01 = iv[x0];
    v[f0] = v[f];
    iv[irc] = 4;
    iv[kagqt] = -1;

/*     ***  copy x to x0, g to g0  *** */

    vcopy_(n, &v[x01], &x[1]);
    vcopy_(n, &v[g01], &g[1]);

/*  ***  check stopx and function evaluation limit  *** */

L110:
    if (! stopx_(&dummy)) {
	goto L130;
    }
    iv[1] = 11;
    goto L140;

/*     ***  come here when restarting after func. eval. limit or stopx. */

L120:
    if (v[f] >= v[f0]) {
	goto L130;
    }
    v[radfac] = one;
    k = iv[niter];
    goto L100;

L130:
    if (iv[nfcall] < iv[mxfcal]) {
	goto L150;
    }
    iv[1] = 9;
L140:
    if (v[f] >= v[f0]) {
	goto L300;
    }

/*        ***  in case of stopx or function evaluation limit with */
/*        ***  improved v(f), evaluate the gradient at x. */

    iv[cnvcod] = iv[1];
    goto L240;

/* . . . . . . . . . . . . .  compute candidate step  . . . . . . . . . . */

L150:
    step1 = iv[step];
    dg1 = iv[dg];
    nwtst1 = iv[nwtstp];
    if (iv[kagqt] >= 0) {
	goto L160;
    }
    l = iv[lmat];
    livmul_(n, &v[nwtst1], &v[l], &g[1]);
    v[nreduc] = half * dotprd_(n, &v[nwtst1], &v[nwtst1]);
    litvmu_(n, &v[nwtst1], &v[l], &v[nwtst1]);
    vvmulp_(n, &v[step1], &v[nwtst1], &d__[1], &c__1);
    v[dst0] = v2norm_(n, &v[step1]);
    vvmulp_(n, &v[dg1], &v[dg1], &d__[1], &c_n1);
    ltvmul_(n, &v[step1], &v[l], &v[dg1]);
    v[gthg] = v2norm_(n, &v[step1]);
    iv[kagqt] = 0;
L160:
    dbdog_(&v[dg1], lv, n, &v[nwtst1], &v[step1], &v[1]);
    if (iv[irc] == 6) {
	goto L180;
    }

/*  ***  check whether evaluating f(x0 + step) looks worthwhile  *** */

    if (v[dstnrm] <= zero) {
	goto L180;
    }
    if (iv[irc] != 5) {
	goto L170;
    }
    if (v[radfac] <= one) {
	goto L170;
    }
    if (v[preduc] <= onep2 * v[fdif]) {
	goto L180;
    }

/*  ***  compute f(x0 + step)  *** */

L170:
    x01 = iv[x0];
    step1 = iv[step];
    vaxpy_(n, &x[1], &one, &v[step1], &v[x01]);
    ++iv[nfcall];
    iv[1] = 1;
    iv[toobig] = 0;
    goto L999;

/* . . . . . . . . . . . . .  assess candidate step  . . . . . . . . . . . */

L180:
    x01 = iv[x0];
    v[reldx] = reldst_(n, &d__[1], &x[1], &v[x01]);
    assst_(&iv[1], liv, lv, &v[1]);
    step1 = iv[step];
    lstgst = iv[stlstg];
    if (iv[restor] == 1) {
	vcopy_(n, &x[1], &v[x01]);
    }
    if (iv[restor] == 2) {
	vcopy_(n, &v[lstgst], &v[step1]);
    }
    if (iv[restor] != 3) {
	goto L190;
    }
    vcopy_(n, &v[step1], &v[lstgst]);
    vaxpy_(n, &x[1], &one, &v[step1], &v[x01]);
    v[reldx] = reldst_(n, &d__[1], &x[1], &v[x01]);

L190:
    k = iv[irc];
    switch (k) {
	case 1:  goto L200;
	case 2:  goto L230;
	case 3:  goto L230;
	case 4:  goto L230;
	case 5:  goto L200;
	case 6:  goto L210;
	case 7:  goto L220;
	case 8:  goto L220;
	case 9:  goto L220;
	case 10:  goto L220;
	case 11:  goto L220;
	case 12:  goto L220;
	case 13:  goto L280;
	case 14:  goto L250;
    }

/*     ***  recompute step with changed radius  *** */

L200:
    v[radius] = v[radfac] * v[dstnrm];
    goto L110;

/*  ***  compute step of length v(lmaxs) for singular convergence test. */

L210:
    v[radius] = v[lmaxs];
    goto L150;

/*  ***  convergence or false convergence  *** */

L220:
    iv[cnvcod] = k - 4;
    if (v[f] >= v[f0]) {
	goto L290;
    }
    if (iv[xirc] == 14) {
	goto L290;
    }
    iv[xirc] = 14;

/* . . . . . . . . . . . .  process acceptable step  . . . . . . . . . . . */

L230:
    if (iv[irc] != 3) {
	goto L240;
    }
    step1 = iv[step];
    temp1 = iv[stlstg];

/*     ***  set  temp1 = hessian * step  for use in gradient tests  *** */

    l = iv[lmat];
    ltvmul_(n, &v[temp1], &v[l], &v[step1]);
    lvmul_(n, &v[temp1], &v[l], &v[temp1]);

/*  ***  compute gradient  *** */

L240:
    ++iv[ngcall];
    iv[1] = 2;
    goto L999;

/*  ***  initializations -- g0 = g - g0, etc.  *** */

L250:
    g01 = iv[g0];
    vaxpy_(n, &v[g01], &negone, &v[g01], &g[1]);
    step1 = iv[step];
    temp1 = iv[stlstg];
    if (iv[irc] != 3) {
	goto L270;
    }

/*  ***  set v(radfac) by gradient tests  *** */

/*     ***  set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))  *** */

    vaxpy_(n, &v[temp1], &negone, &v[g01], &v[temp1]);
    vvmulp_(n, &v[temp1], &v[temp1], &d__[1], &c_n1);

/*        ***  do gradient tests  *** */

    if (v2norm_(n, &v[temp1]) <= v[dgnorm] * v[tuner4]) {
	goto L260;
    }
    if (dotprd_(n, &g[1], &v[step1]) >= v[gtstep] * v[tuner5]) {
	goto L270;
    }
L260:
    v[radfac] = v[incfac];

/*  ***  update h, loop  *** */

L270:
    w = iv[nwtstp];
    z__ = iv[x0];
    l = iv[lmat];
    wzbfgs_(&v[l], n, &v[step1], &v[w], &v[g01], &v[z__]);

/*     ** use the n-vectors starting at v(step1) and v(g01) for scratch.. */
    lupdat_(&v[temp1], &v[step1], &v[l], &v[g01], &v[l], n, &v[w], &v[z__]);
    iv[1] = 2;
    goto L80;

/* . . . . . . . . . . . . . .  misc. details  . . . . . . . . . . . . . . */

/*  ***  bad parameters to assess  *** */

L280:
    iv[1] = 64;
    goto L300;

/*  ***  print summary of final iteration and other requested items  *** */

L290:
    iv[1] = iv[cnvcod];
    iv[cnvcod] = 0;
L300:
    itsum_(&d__[1], &g[1], &iv[1], liv, lv, n, &v[1], &x[1]);

L999:
    return 0;

/*  ***  last line of sumit follows  *** */
} /* sumit_ */

/* Subroutine */ int snoit_(doublereal *d__, doublereal *fx, integer *iv, 
	integer *liv, integer *lv, integer *n, doublereal *v, doublereal *x)
{
    /* Initialized data */

    static __thread integer vneed = 4;
    static __thread doublereal zero = 0.;
    static __thread integer eta0 = 42;
    static __thread integer f = 10;
    static __thread integer g = 28;
    static __thread integer lmat = 42;
    static __thread integer nextv = 47;
    static __thread integer nfgcal = 7;
    static __thread integer ngcall = 30;
    static __thread integer niter = 31;
    static __thread integer sgirc = 57;
    static __thread integer toobig = 2;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static __thread integer i__, j, k, w, g1, iv1, alpha;
    extern /* Subroutine */ int deflt_(integer *, integer *, integer *, 
	    integer *, doublereal *), sumit_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *), sgrad2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int vscopy_(integer *, doublereal *, doublereal *)
	    ;


/*  ***  iteration driver for smsno... */
/*  ***  minimize general unconstrained objective function using */
/*  ***  finite-difference gradients and secant hessian approximations. */

/*     dimension v(77 + n*(n+17)/2) */

/*  ***  purpose  *** */

/*        this routine interacts with subroutine  sumit  in an attempt */
/*     to find an n-vector  x*  that minimizes the (unconstrained) */
/*     objective function  fx = f(x)  computed by the caller.  (often */
/*     the  x*  found is a local minimizer rather than a global one.) */

/*  ***  parameters  *** */

/*        the parameters for snoit are the same as those for sumsl */
/*     (which see), except that calcf, calcg, uiparm, urparm, and ufparm */
/*     are omitted, and a parameter  fx  for the objective function */
/*     value at x is added.  instead of calling calcg to obtain the */
/*     gradient of the objective function at x, snoit calls sgrad2, */
/*     which computes an approximation to the gradient by finite */
/*     (forward and central) differences using the method of ref. 1. */
/*     the following input component is of interest in this regard */
/*     (and is not described in sumsl). */

/* v(eta0)..... v(42) is an estimated bound on the relative error in the */
/*             objective function value computed by calcf... */
/*                  (true value) = (computed value) * (1 + e), */
/*             where abs(e) .le. v(eta0).  default = machep * 10**3, */
/*             where machep is the unit roundoff. */

/*        the output values iv(nfcall) and iv(ngcall) have different */
/*     meanings for smsno than for sumsl... */

/* iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e., */
/*             function evaluations) excluding those made only for */
/*             computing gradients.  the input value iv(mxfcal) is a */
/*             limit on iv(nfcall). */
/* iv(ngcall)... iv(30) is the number of function evaluations made only */
/*             for computing gradients.  the total number of function */
/*             evaluations is thus  iv(nfcall) + iv(ngcall). */

/*  ***  references  *** */

/* 1. stewart, g.w. (1967), a modification of davidon*s minimization */
/*        method to accept difference approximations of derivatives, */
/*        j. assoc. comput. mach. 14, pp. 72-83. */
/* . */
/*  ***  general  *** */

/*     coded by david m. gay (august 1982). */

/* ----------------------------  declarations  --------------------------- */


/* deflt.... supplies default parameter values. */
/* dotprd... returns inner product of two vectors. */
/* sgrad2... computes finite-difference gradient approximation. */
/* sumit.... reverse-communication routine that does sumsl algorithm. */
/* vscopy... sets all elements of a vector to a scalar. */


/*  ***  subscripts for iv   *** */


/* /6 */
    /* Parameter adjustments */
    --iv;
    --v;
    --x;
    --d__;

    /* Function Body */
/* /7 */
/*     parameter (dtype=16, eta0=42, f=10, g=28, lmat=42, nextv=47, */
/*    1           nfcall=6, nfgcal=7, ngcall=30, niter=31, sgirc=57, */
/*    2           toobig=2, vneed=4) */
/* / */
/* /6 */
/* /7 */
/*     parameter (one=1.d+0, zero=0.d+0) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

    iv1 = iv[1];
    if (iv1 == 1) {
	goto L10;
    }
    if (iv1 == 2) {
	goto L50;
    }
    if (iv[1] == 0) {
	deflt_(&c__2, &iv[1], liv, lv, &v[1]);
    }
    iv1 = iv[1];
    if (iv1 == 12 || iv1 == 13) {
	iv[vneed] = iv[vneed] + (*n << 1) + 6;
    }
    if (iv1 == 14) {
	goto L10;
    }
    if (iv1 > 2 && iv1 < 12) {
	goto L10;
    }
    g1 = 1;
    if (iv1 == 12) {
	iv[1] = 13;
    }
    goto L20;

L10:
    g1 = iv[g];

L20:
    sumit_(&d__[1], fx, &v[g1], &iv[1], liv, lv, n, &v[1], &x[1]);
    if ((i__1 = iv[1] - 2) < 0) {
	goto L999;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L70;
    }

/*  ***  compute gradient  *** */

L30:
    if (iv[niter] == 0) {
	vscopy_(n, &v[g1], &zero);
    }
    j = iv[lmat];
    k = g1 - *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = dotprd_(&i__, &v[j], &v[j]);
	++k;
	j += i__;
/* L40: */
    }
/*     ***  undo increment of iv(ngcall) done by sumit  *** */
    --iv[ngcall];
/*     ***  store return code from sgrad2 in iv(sgirc)  *** */
    iv[sgirc] = 0;
/*     ***  x may have been restored, so copy back fx... *** */
    *fx = v[f];
    goto L60;

/*     ***  gradient loop  *** */

L50:
    if (iv[toobig] == 0) {
	goto L60;
    }
    iv[nfgcal] = 0;
    goto L10;

L60:
    g1 = iv[g];
    alpha = g1 - *n;
    w = alpha - 6;
    sgrad2_(&v[alpha], &d__[1], &v[eta0], fx, &v[g1], &iv[sgirc], n, &v[w], &
	    x[1]);
    if (iv[sgirc] == 0) {
	goto L10;
    }
    ++iv[ngcall];
    goto L999;

L70:
    if (iv[1] != 14) {
	goto L999;
    }

/*  ***  storage allocation  *** */

    iv[g] = iv[nextv] + *n + 6;
    iv[nextv] = iv[g] + *n;
    if (iv1 != 13) {
	goto L10;
    }

L999:
    return 0;
/*  ***  last card of snoit follows  *** */
} /* snoit_ */

/* Subroutine */ int dbdog_(doublereal *dig, integer *lv, integer *n, 
	doublereal *nwtstp, doublereal *step, doublereal *v)
{
    /* Initialized data */

    static __thread doublereal half = .5;
    static __thread doublereal one = 1.;
    static __thread doublereal two = 2.;
    static __thread doublereal zero = 0.;
    static __thread integer bias = 43;
    static __thread integer dgnorm = 1;
    static __thread integer dstnrm = 2;
    static __thread integer dst0 = 3;
    static __thread integer grdfac = 45;
    static __thread integer gthg = 44;
    static __thread integer gtstep = 4;
    static __thread integer nreduc = 6;
    static __thread integer nwtfac = 46;
    static __thread integer preduc = 7;
    static __thread integer radius = 8;
    static __thread integer stppar = 5;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static __thread integer i__;
    static __thread doublereal t, t1, t2, cfact, relax, cnorm, gnorm, rlambd, ghinvg,
	    femnsq, ctrnwt, nwtnrm;


/*  ***  compute double dogleg step  *** */

/*  ***  parameter declarations  *** */


/*  ***  purpose  *** */

/*        this subroutine computes a candidate step (for use in an uncon- */
/*     strained minimization code) by the double dogleg algorithm of */
/*     dennis and mei (ref. 1), which is a variation on powell*s dogleg */
/*     scheme (ref. 2, p. 95). */

/* --------------------------  parameter usage  -------------------------- */

/*    dig (input) diag(d)**-2 * g -- see algorithm notes. */
/*      g (input) the current gradient vector. */
/*     lv (input) length of v. */
/*      n (input) number of components in  dig, g, nwtstp,  and  step. */
/* nwtstp (input) negative newton step -- see algorithm notes. */
/*   step (output) the computed step. */
/*      v (i/o) values array, the following components of which are */
/*             used here... */
/* v(bias)   (input) bias for relaxed newton step, which is v(bias) of */
/*             the way from the full newton to the fully relaxed newton */
/*             step.  recommended value = 0.8 . */
/* v(dgnorm) (input) 2-norm of diag(d)**-1 * g -- see algorithm notes. */
/* v(dstnrm) (output) 2-norm of diag(d) * step, which is v(radius) */
/*             unless v(stppar) = 0 -- see algorithm notes. */
/* v(dst0) (input) 2-norm of diag(d) * nwtstp -- see algorithm notes. */
/* v(grdfac) (output) the coefficient of  dig  in the step returned -- */
/*             step(i) = v(grdfac)*dig(i) + v(nwtfac)*nwtstp(i). */
/* v(gthg)   (input) square-root of (dig**t) * (hessian) * dig -- see */
/*             algorithm notes. */
/* v(gtstep) (output) inner product between g and step. */
/* v(nreduc) (output) function reduction predicted for the full newton */
/*             step. */
/* v(nwtfac) (output) the coefficient of  nwtstp  in the step returned -- */
/*             see v(grdfac) above. */
/* v(preduc) (output) function reduction predicted for the step returned. */
/* v(radius) (input) the trust region radius.  d times the step returned */
/*             has 2-norm v(radius) unless v(stppar) = 0. */
/* v(stppar) (output) code telling how step was computed... 0 means a */
/*             full newton step.  between 0 and 1 means v(stppar) of the */
/*             way from the newton to the relaxed newton step.  between */
/*             1 and 2 means a true double dogleg step, v(stppar) - 1 of */
/*             the way from the relaxed newton to the cauchy step. */
/*             greater than 2 means 1 / (v(stppar) - 1) times the cauchy */
/*             step. */

/* -------------------------------  notes  ------------------------------- */

/*  ***  algorithm notes  *** */

/*        let  g  and  h  be the current gradient and hessian approxima- */
/*     tion respectively and let d be the current scale vector.  this */
/*     routine assumes dig = diag(d)**-2 * g  and  nwtstp = h**-1 * g. */
/*     the step computed is the same one would get by replacing g and h */
/*     by  diag(d)**-1 * g  and  diag(d)**-1 * h * diag(d)**-1, */
/*     computing step, and translating step back to the original */
/*     variables, i.e., premultiplying it by diag(d)**-1. */

/*  ***  references  *** */

/* 1.  dennis, j.e., and mei, h.h.w. (1979), two new unconstrained opti- */
/*             mization algorithms which use function and gradient */
/*             values, j. optim. theory applic. 28, pp. 453-482. */
/* 2. powell, m.j.d. (1970), a hybrid method for non-linear equations, */
/*             in numerical methods for non-linear equations, edited by */
/*             p. rabinowitz, gordon and breach, london. */

/*  ***  general  *** */

/*     coded by david m. gay. */
/*     this subroutine was written in connection with research supported */
/*     by the national science foundation under grants mcs-7600324 and */
/*     mcs-7906671. */

/* ------------------------  external quantities  ------------------------ */

/*  ***  functions and subroutines called  *** */


/* dotprd... returns inner product of two vectors. */
/* v2norm... returns 2-norm of a vector. */

/*  ***  intrinsic functions  *** */
/* /+ */
/* / */
/* --------------------------  local variables  -------------------------- */


/*  ***  v subscripts  *** */


/*  ***  data initializations  *** */

/* /6 */
    /* Parameter adjustments */
    --v;
    --step;
    --nwtstp;
    --dig;

    /* Function Body */
/* /7 */
/*     parameter (half=0.5d+0, one=1.d+0, two=2.d+0, zero=0.d+0) */
/* / */

/* /6 */
/* /7 */
/*     parameter (bias=43, dgnorm=1, dstnrm=2, dst0=3, grdfac=45, */
/*    1           gthg=44, gtstep=4, nreduc=6, nwtfac=46, preduc=7, */
/*    2           radius=8, stppar=5) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

    nwtnrm = v[(0 + (0 + (dst0 << 3))) / 8];
    rlambd = one;
    if (nwtnrm > zero) {
	rlambd = v[radius] / nwtnrm;
    }
    gnorm = v[dgnorm];
    ghinvg = two * v[nreduc];
    v[grdfac] = zero;
    v[nwtfac] = zero;
    if (rlambd < one) {
	goto L30;
    }

/*        ***  the newton step is inside the trust region  *** */

    v[stppar] = zero;
    v[dstnrm] = nwtnrm;
    v[gtstep] = -ghinvg;
    v[preduc] = v[nreduc];
    v[nwtfac] = -one;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	step[i__] = -nwtstp[i__];
    }
    goto L999;

L30:
    v[dstnrm] = v[radius];
/* Computing 2nd power */
    d__1 = gnorm / v[gthg];
    cfact = d__1 * d__1;
/*     ***  cauchy step = -cfact * g. */
    cnorm = gnorm * cfact;
    relax = one - v[bias] * (one - gnorm * cnorm / ghinvg);
    if (rlambd < relax) {
	goto L50;
    }

/*        ***  step is between relaxed newton and full newton steps  *** */

    v[stppar] = one - (rlambd - relax) / (one - relax);
    t = -rlambd;
    v[gtstep] = t * ghinvg;
    v[preduc] = rlambd * (one - half * rlambd) * ghinvg;
    v[nwtfac] = t;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	step[i__] = t * nwtstp[i__];
    }
    goto L999;

L50:
    if (cnorm < v[radius]) {
	goto L70;
    }

/*        ***  the cauchy step lies outside the trust region -- */
/*        ***  step = scaled cauchy step  *** */

    t = -v[radius] / gnorm;
    v[grdfac] = t;
    v[stppar] = one + cnorm / v[radius];
    v[gtstep] = -v[radius] * gnorm;
/* Computing 2nd power */
    d__1 = v[gthg] / gnorm;
    v[preduc] = v[radius] * (gnorm - half * v[radius] * (d__1 * d__1));
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	step[i__] = t * dig[i__];
    }
    goto L999;

/*     ***  compute dogleg step between cauchy and relaxed newton  *** */
/*     ***  femur = relaxed newton step minus cauchy step  *** */

L70:
    ctrnwt = cfact * relax * ghinvg / gnorm;
/*     *** ctrnwt = inner prod. of cauchy and relaxed newton steps, */
/*     *** scaled by gnorm**-1. */
/* Computing 2nd power */
    d__1 = cfact;
    t1 = ctrnwt - gnorm * (d__1 * d__1);
/*     ***  t1 = inner prod. of femur and cauchy step, scaled by */
/*     ***  gnorm**-1. */
/* Computing 2nd power */
    d__1 = cfact;
    t2 = v[radius] * (v[radius] / gnorm) - gnorm * (d__1 * d__1);
    t = relax * nwtnrm;
    femnsq = t / gnorm * t - ctrnwt - t1;
/*     ***  femnsq = square of 2-norm of femur, scaled by gnorm**-1. */
/* Computing 2nd power */
    d__1 = t1;
    t = t2 / (t1 + sqrt(d__1 * d__1 + femnsq * t2));
/*     ***  dogleg step  =  cauchy step  +  t * femur. */
    t1 = (t - one) * cfact;
    v[grdfac] = t1;
    t2 = -t * relax;
    v[nwtfac] = t2;
    v[stppar] = two - t;
/* Computing 2nd power */
    d__1 = gnorm;
    v[gtstep] = t1 * (d__1 * d__1) + t2 * ghinvg;
/* Computing 2nd power */
    d__1 = v[gthg] * t1;
    v[preduc] = -t1 * gnorm * ((t2 + one) * gnorm) - t2 * (one + half * t2) * 
	    ghinvg - half * (d__1 * d__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	step[i__] = t1 * dig[i__] + t2 * nwtstp[i__];
    }

L999:
    return 0;
/*  ***  last line of dbdog follows  *** */
} /* dbdog_ */

/* Subroutine */ int ltvmul_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* Initialized data */

    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static __thread integer i__, j, i0, ij;
    static __thread doublereal yi;


/*  ***  compute  x = (l**t)*y, where  l  is an  n x n  lower */
/*  ***  triangular matrix stored compactly by rows.  x and y may */
/*  ***  occupy the same storage.  *** */

/*     dimension l(n*(n+1)/2) */
/* /6 */
    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
/* /7 */
/*     parameter (zero=0.d+0) */
/* / */

    i0 = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yi = y[i__];
	x[i__] = zero;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    ij = i0 + j;
	    x[j] += yi * l[ij];
/* L10: */
	}
	i0 += i__;
/* L20: */
    }
/* L999: */
    return 0;
/*  ***  last card of ltvmul follows  *** */
} /* ltvmul_ */

/* Subroutine */ int lupdat_(doublereal *beta, doublereal *gamma, doublereal *
	l, doublereal *lambda, doublereal *lplus, integer *n, doublereal *w, 
	doublereal *z__)
{
    /* Initialized data */

    static __thread doublereal one = 1.;
    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static __thread doublereal a, b;
    static __thread integer i__, j, k;
    static __thread doublereal s, bj, gj;
    static __thread integer ij, jj;
    static __thread doublereal lj, wj, nu, zj;
    static __thread integer jp1, nm1, np1;
    static __thread doublereal eta, lij, ljj, theta;


/*  ***  compute lplus = secant update of l  *** */

/*  ***  parameter declarations  *** */

/*     dimension l(n*(n+1)/2), lplus(n*(n+1)/2) */

/* --------------------------  parameter usage  -------------------------- */

/*   beta = scratch vector. */
/*  gamma = scratch vector. */
/*      l (input) lower triangular matrix, stored rowwise. */
/* lambda = scratch vector. */
/*  lplus (output) lower triangular matrix, stored rowwise, which may */
/*             occupy the same storage as  l. */
/*      n (input) length of vector parameters and order of matrices. */
/*      w (input, destroyed on output) right singular vector of rank 1 */
/*             correction to  l. */
/*      z (input, destroyed on output) left singular vector of rank 1 */
/*             correction to  l. */

/* -------------------------------  notes  ------------------------------- */

/*  ***  application and usage restrictions  *** */

/*        this routine updates the cholesky factor  l  of a symmetric */
/*     positive definite matrix to which a secant update is being */
/*     applied -- it computes a cholesky factor  lplus  of */
/*     l * (i + z*w**t) * (i + w*z**t) * l**t.  it is assumed that  w */
/*     and  z  have been chosen so that the updated matrix is strictly */
/*     positive definite. */

/*  ***  algorithm notes  *** */

/*        this code uses recurrence 3 of ref. 1 (with d(j) = 1 for all j) */
/*     to compute  lplus  of the form  l * (i + z*w**t) * q,  where  q */
/*     is an orthogonal matrix that makes the result lower triangular. */
/*        lplus may have some negative diagonal elements. */

/*  ***  references  *** */

/* 1.  goldfarb, d. (1976), factorized variable metric methods for uncon- */
/*             strained optimization, math. comput. 30, pp. 796-811. */

/*  ***  general  *** */

/*     coded by david m. gay (fall 1979). */
/*     this subroutine was written in connection with research supported */
/*     by the national science foundation under grants mcs-7600324 and */
/*     mcs-7906671. */

/* ------------------------  external quantities  ------------------------ */

/*  ***  intrinsic functions  *** */
/* /+ */
/* / */
/* --------------------------  local variables  -------------------------- */


/*  ***  data initializations  *** */

/* /6 */
    /* Parameter adjustments */
    --l;
    --lplus;
    --z__;
    --w;
    --lambda;
    --gamma;
    --beta;

    /* Function Body */
/* /7 */
/*     parameter (one=1.d+0, zero=0.d+0) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

    nu = one;
    eta = zero;
    if (*n <= 1) {
	goto L30;
    }
    nm1 = *n - 1;

/*  ***  temporarily store s(j) = sum over k = j+1 to n of w(k)**2 in */
/*  ***  lambda(j). */

    s = zero;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n - i__;
/* Computing 2nd power */
	d__1 = w[j + 1];
	s += d__1 * d__1;
	lambda[j] = s;
/* L10: */
    }

/*  ***  compute lambda, gamma, and beta by goldfarb*s recurrence 3. */

    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	wj = w[j];
	a = nu * z__[j] - eta * wj;
	theta = one + a * wj;
	s = a * lambda[j];
/* Computing 2nd power */
	d__1 = theta;
	lj = sqrt(d__1 * d__1 + a * s);
	if (theta > zero) {
	    lj = -lj;
	}
	lambda[j] = lj;
	b = theta * wj + s;
	gamma[j] = b * nu / lj;
	beta[j] = (a - b * eta) / lj;
	nu = -nu / lj;
/* Computing 2nd power */
	d__1 = a;
	eta = -(eta + d__1 * d__1 / (theta - lj)) / lj;
/* L20: */
    }
L30:
    lambda[*n] = one + (nu * z__[*n] - eta * w[*n]) * w[*n];

/*  ***  update l, gradually overwriting  w  and  z  with  l*w  and  l*z. */

    np1 = *n + 1;
    jj = *n * (*n + 1) / 2;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	j = np1 - k;
	lj = lambda[j];
	ljj = l[jj];
	lplus[jj] = lj * ljj;
	wj = w[j];
	w[j] = ljj * wj;
	zj = z__[j];
	z__[j] = ljj * zj;
	if (k == 1) {
	    goto L50;
	}
	bj = beta[j];
	gj = gamma[j];
	ij = jj + j;
	jp1 = j + 1;
	i__2 = *n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    lij = l[ij];
	    lplus[ij] = lj * lij + bj * w[i__] + gj * z__[i__];
	    w[i__] += lij * wj;
	    z__[i__] += lij * zj;
	    ij += i__;
/* L40: */
	}
L50:
	jj -= j;
/* L60: */
    }

/* L999: */
    return 0;
/*  ***  last card of lupdat follows  *** */
} /* lupdat_ */

/* Subroutine */ int lvmul_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* Initialized data */

    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static __thread integer i__, j;
    static __thread doublereal t;
    static __thread integer i0, ii, ij, np1;


/*  ***  compute  x = l*y, where  l  is an  n x n  lower triangular */
/*  ***  matrix stored compactly by rows.  x and y may occupy the same */
/*  ***  storage.  *** */

/*     dimension l(n*(n+1)/2) */
/* /6 */
    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
/* /7 */
/*     parameter (zero=0.d+0) */
/* / */

    np1 = *n + 1;
    i0 = *n * (*n + 1) / 2;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = np1 - ii;
	i0 -= i__;
	t = zero;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    ij = i0 + j;
	    t += l[ij] * y[j];
/* L10: */
	}
	x[i__] = t;
/* L20: */
    }
/* L999: */
    return 0;
/*  ***  last card of lvmul follows  *** */
} /* lvmul_ */

/* Subroutine */ int sgrad2_(doublereal *alpha, doublereal *d__, doublereal *
	eta0, doublereal *fx, doublereal *g, integer *irc, integer *n, 
	doublereal *w, doublereal *x)
{
    /* Initialized data */

    static __thread doublereal c2000 = 2e3;
    static __thread doublereal four = 4.;
    static __thread doublereal hmax0 = .02;
    static __thread doublereal hmin0 = 50.;
    static __thread doublereal one = 1.;
    static __thread doublereal p002 = .002;
    static __thread doublereal three = 3.;
    static __thread doublereal two = 2.;
    static __thread doublereal zero = 0.;
    static __thread integer fh = 3;
    static __thread integer fx0 = 4;
    static __thread integer hsave = 5;
    static __thread integer xisave = 6;

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static __thread doublereal h__;
    static __thread integer i__;
    static __thread doublereal h0, gi, aai, agi, eta, afx, axi, hmin, machep, alphai,
	    axibar, afxeta, discon;
    extern doublereal rmdcon_(integer *);


/*  ***  compute finite difference gradient by stweart*s scheme  *** */

/*     ***  parameters  *** */


/* ....................................................................... */

/*     ***  purpose  *** */

/*        this subroutine uses an embellished form of the finite-differ- */
/*     ence scheme proposed by stewart (ref. 1) to approximate the */
/*     gradient of the function f(x), whose values are supplied by */
/*     reverse communication. */

/*     ***  parameter description  *** */

/*  alpha in  (approximate) diagonal elements of the hessian of f(x). */
/*      d in  scale vector such that d(i)*x(i), i = 1,...,n, are in */
/*             comparable units. */
/*   eta0 in  estimated bound on relative error in the function value... */
/*             (true value) = (computed value)*(1+e),   where */
/*             abs(e) .le. eta0. */
/*     fx i/o on input,  fx  must be the computed value of f(x).  on */
/*             output with irc = 0, fx has been restored to its original */
/*             value, the one it had when sgrad2 was last called with */
/*             irc = 0. */
/*      g i/o on input with irc = 0, g should contain an approximation */
/*             to the gradient of f near x, e.g., the gradient at the */
/*             previous iterate.  when sgrad2 returns with irc = 0, g is */
/*             the desired finite-difference approximation to the */
/*             gradient at x. */
/*    irc i/o input/return code... before the very first call on sgrad2, */
/*             the caller must set irc to 0.  whenever sgrad2 returns a */
/*             nonzero value for irc, it has perturbed some component of */
/*             x... the caller should evaluate f(x) and call sgrad2 */
/*             again with fx = f(x). */
/*      n in  the number of variables (components of x) on which f */
/*             depends. */
/*      x i/o on input with irc = 0, x is the point at which the */
/*             gradient of f is desired.  on output with irc nonzero, x */
/*             is the point at which f should be evaluated.  on output */
/*             with irc = 0, x has been restored to its original value */
/*             (the one it had when sgrad2 was last called with irc = 0) */
/*             and g contains the desired gradient approximation. */
/*      w i/o work vector of length 6 in which sgrad2 saves certain */
/*             quantities while the caller is evaluating f(x) at a */
/*             perturbed x. */

/*     ***  application and usage restrictions  *** */

/*        this routine is intended for use with quasi-newton routines */
/*     for unconstrained minimization (in which case  alpha  comes from */
/*     the diagonal of the quasi-newton hessian approximation). */

/*     ***  algorithm notes  *** */

/*        this code departs from the scheme proposed by stewart (ref. 1) */
/*     in its guarding against overly large or small step sizes and its */
/*     handling of special cases (such as zero components of alpha or g). */

/*     ***  references  *** */

/* 1. stewart, g.w. (1967), a modification of davidon*s minimization */
/*        method to accept difference approximations of derivatives, */
/*        j. assoc. comput. mach. 14, pp. 72-83. */

/*     ***  history  *** */

/*     designed and coded by david m. gay (summer 1977/summer 1980). */

/*     ***  general  *** */

/*        this routine was prepared in connection with work supported by */
/*     the national science foundation under grants mcs76-00324 and */
/*     mcs-7906671. */

/* ....................................................................... */

/*     *****  external function  ***** */

/* rmdcon... returns machine-dependent constants. */

/*     ***** intrinsic functions ***** */
/* /+ */
/* / */
/*     ***** local variables ***** */


/* /6 */
    /* Parameter adjustments */
    --x;
    --g;
    --d__;
    --alpha;
    --w;

    /* Function Body */
/* /7 */
/*     parameter (c2000=2.0d+3, four=4.0d+0, hmax0=0.02d+0, hmin0=5.0d+1, */
/*    1     one=1.0d+0, p002=0.002d+0, three=3.0d+0, */
/*    2     two=2.0d+0, zero=0.0d+0) */
/* / */
/* /6 */
/* /7 */
/*     parameter (fh=3, fx0=4, hsave=5, xisave=6) */
/* / */

/* ---------------------------------  body  ------------------------------ */

    if (*irc < 0) {
	goto L140;
    } else if (*irc == 0) {
	goto L100;
    } else {
	goto L210;
    }

/*     ***  fresh start -- get machine-dependent constants  *** */

/*     store machep in w(1) and h0 in w(2), where machep is the unit */
/*     roundoff (the smallest positive number such that */
/*     1 + machep .gt. 1  and  1 - machep .lt. 1),  and  h0 is the */
/*     square-root of machep. */

L100:
    w[1] = rmdcon_(&c__3);
    w[2] = sqrt(w[1]);

    w[fx0] = *fx;

/*     ***  increment  i  and start computing  g(i)  *** */

L110:
    i__ = abs(*irc) + 1;
    if (i__ > *n) {
	goto L300;
    }
    *irc = i__;
    afx = (d__1 = w[fx0], abs(d__1));
    machep = w[1];
    h0 = w[2];
    hmin = hmin0 * machep;
    w[xisave] = x[i__];
    axi = (d__1 = x[i__], abs(d__1));
/* Computing MAX */
    d__1 = axi, d__2 = one / d__[i__];
    axibar = max(d__1,d__2);
    gi = g[i__];
    agi = abs(gi);
    eta = abs(*eta0);
    if (afx > zero) {
/* Computing MAX */
	d__1 = eta, d__2 = agi * axi * machep / afx;
	eta = max(d__1,d__2);
    }
    alphai = alpha[i__];
    if (alphai == zero) {
	goto L170;
    }
    if (gi == zero || *fx == zero) {
	goto L180;
    }
    afxeta = afx * eta;
    aai = abs(alphai);

/*        *** compute h = stewart*s forward-difference step size. */

/* Computing 2nd power */
    d__1 = gi;
    if (d__1 * d__1 <= afxeta * aai) {
	goto L120;
    }
    h__ = two * sqrt(afxeta / aai);
    h__ *= one - aai * h__ / (three * aai * h__ + four * agi);
    goto L130;
L120:
/* Computing 2nd power */
    d__2 = aai;
    d__1 = afxeta * agi / (d__2 * d__2);
    d__3 = one / three;
    h__ = two * pow_dd(&d__1, &d__3);
    h__ *= one - two * agi / (three * aai * h__ + four * agi);

/*        ***  ensure that  h  is not insignificantly small  *** */

L130:
/* Computing MAX */
    d__1 = h__, d__2 = hmin * axibar;
    h__ = max(d__1,d__2);

/*        *** use forward difference if bound on truncation error is at */
/*        *** most 10**-3. */

    if (aai * h__ <= p002 * agi) {
	goto L160;
    }

/*        *** compute h = stewart*s step for central difference. */

    discon = c2000 * afxeta;
/* Computing 2nd power */
    d__1 = gi;
    h__ = discon / (agi + sqrt(d__1 * d__1 + aai * discon));

/*        ***  ensure that  h  is neither too small nor too big  *** */

/* Computing MAX */
    d__1 = h__, d__2 = hmin * axibar;
    h__ = max(d__1,d__2);
    if (h__ >= hmax0 * axibar) {
	d__1 = two / three;
	h__ = axibar * pow_dd(&h0, &d__1);
    }

/*        ***  compute central difference  *** */

    *irc = -i__;
    goto L200;

L140:
    h__ = -w[hsave];
    i__ = abs(*irc);
    if (h__ > zero) {
	goto L150;
    }
    w[fh] = *fx;
    goto L200;

L150:
    g[i__] = (w[fh] - *fx) / (two * h__);
    x[i__] = w[xisave];
    goto L110;

/*     ***  compute forward differences in various cases  *** */

L160:
    if (h__ >= hmax0 * axibar) {
	h__ = h0 * axibar;
    }
    if (alphai * gi < zero) {
	h__ = -h__;
    }
    goto L200;
L170:
    h__ = axibar;
    goto L200;
L180:
    h__ = h0 * axibar;

L200:
    x[i__] = w[xisave] + h__;
    w[hsave] = h__;
    goto L999;

/*     ***  compute actual forward difference  *** */

L210:
    g[*irc] = (*fx - w[fx0]) / w[hsave];
    x[*irc] = w[xisave];
    goto L110;

/*  ***  restore fx and indicate that g has been computed  *** */

L300:
    *fx = w[fx0];
    *irc = 0;

L999:
    return 0;
/*  ***  last card of sgrad2 follows  *** */
} /* sgrad2_ */

/* Subroutine */ int vvmulp_(integer *n, doublereal *x, doublereal *y, 
	doublereal *z__, integer *k)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static __thread integer i__;


/* ***  set x(i) = y(i) * z(i)**k, 1 .le. i .le. n (for k = 1 or -1)  *** */


    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    if (*k >= 0) {
	goto L20;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = y[i__] / z__[i__];
    }
    goto L999;

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	x[i__] = y[i__] * z__[i__];
    }
L999:
    return 0;
/*  ***  last card of vvmulp follows  *** */
} /* vvmulp_ */

/* Subroutine */ int wzbfgs_(doublereal *l, integer *n, doublereal *s, 
	doublereal *w, doublereal *y, doublereal *z__)
{
    /* Initialized data */

    static __thread doublereal eps = .1;
    static __thread doublereal one = 1.;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static __thread integer i__;
    static __thread doublereal cs, cy, ys, shs, theta, epsrt;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int livmul_(integer *, doublereal *, doublereal *,
	     doublereal *), ltvmul_(integer *, doublereal *, doublereal *, 
	    doublereal *);


/*  ***  compute  y  and  z  for  lupdat  corresponding to bfgs update. */

/*     dimension l(n*(n+1)/2) */

/* --------------------------  parameter usage  -------------------------- */

/* l (i/o) cholesky factor of hessian, a lower triang. matrix stored */
/*             compactly by rows. */
/* n (input) order of  l  and length of  s,  w,  y,  z. */
/* s (input) the step just taken. */
/* w (output) right singular vector of rank 1 correction to l. */
/* y (input) change in gradients corresponding to s. */
/* z (output) left singular vector of rank 1 correction to l. */

/* -------------------------------  notes  ------------------------------- */

/*  ***  algorithm notes  *** */

/*        when  s  is computed in certain ways, e.g. by  gqtstp  or */
/*     dbldog,  it is possible to save n**2/2 operations since  (l**t)*s */
/*     or  l*(l**t)*s is then known. */
/*        if the bfgs update to l*(l**t) would reduce its determinant to */
/*     less than eps times its old value, then this routine in effect */
/*     replaces  y  by  theta*y + (1 - theta)*l*(l**t)*s,  where  theta */
/*     (between 0 and 1) is chosen to make the reduction factor = eps. */

/*  ***  general  *** */

/*     coded by david m. gay (fall 1979). */
/*     this subroutine was written in connection with research supported */
/*     by the national science foundation under grants mcs-7600324 and */
/*     mcs-7906671. */

/* ------------------------  external quantities  ------------------------ */

/*  ***  functions and subroutines called  *** */

/* dotprd returns inner product of two vectors. */
/* livmul multiplies l**-1 times a vector. */
/* ltvmul multiplies l**t times a vector. */

/*  ***  intrinsic functions  *** */
/* /+ */
/* / */
/* --------------------------  local variables  -------------------------- */


/*  ***  data initializations  *** */

/* /6 */
    /* Parameter adjustments */
    --l;
    --z__;
    --y;
    --w;
    --s;

    /* Function Body */
/* /7 */
/*     parameter (eps=0.1d+0, one=1.d+0) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

    ltvmul_(n, &w[1], &l[1], &s[1]);
    shs = dotprd_(n, &w[1], &w[1]);
    ys = dotprd_(n, &y[1], &s[1]);
    if (ys >= eps * shs) {
	goto L10;
    }
    theta = (one - eps) * shs / (shs - ys);
    epsrt = sqrt(eps);
    cy = theta / (shs * epsrt);
    cs = (one + (theta - one) / epsrt) / shs;
    goto L20;
L10:
    cy = one / (sqrt(ys) * sqrt(shs));
    cs = one / shs;
L20:
    livmul_(n, &z__[1], &l[1], &y[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	z__[i__] = cy * z__[i__] - cs * w[i__];
    }

/* L999: */
    return 0;
/*  ***  last card of wzbfgs follows  *** */
} /* wzbfgs_ */

/* Subroutine */ int assst_(integer *iv, integer *liv, integer *lv, 
	doublereal *v)
{
    /* Initialized data */

    static __thread doublereal half = .5;
    static __thread doublereal one = 1.;
    static __thread doublereal onep2 = 1.2;
    static __thread doublereal two = 2.;
    static __thread doublereal zero = 0.;
    static __thread integer irc = 29;
    static __thread integer mlstgd = 32;
    static __thread integer model = 5;
    static __thread integer nfcall = 6;
    static __thread integer nfgcal = 7;
    static __thread integer radinc = 8;
    static __thread integer restor = 9;
    static __thread integer stage = 10;
    static __thread integer stglim = 11;
    static __thread integer switch__ = 12;
    static __thread integer toobig = 2;
    static __thread integer xirc = 13;
    static __thread integer afctol = 31;
    static __thread integer decfac = 22;
    static __thread integer dstnrm = 2;
    static __thread integer dst0 = 3;
    static __thread integer dstsav = 18;
    static __thread integer f = 10;
    static __thread integer fdif = 11;
    static __thread integer flstgd = 12;
    static __thread integer f0 = 13;
    static __thread integer gtslst = 14;
    static __thread integer gtstep = 4;
    static __thread integer incfac = 23;
    static __thread integer lmaxs = 36;
    static __thread integer nreduc = 6;
    static __thread integer plstgd = 15;
    static __thread integer preduc = 7;
    static __thread integer radfac = 16;
    static __thread integer rdfcmn = 24;
    static __thread integer rdfcmx = 25;
    static __thread integer reldx = 17;
    static __thread integer rfctol = 32;
    static __thread integer sctol = 37;
    static __thread integer stppar = 5;
    static __thread integer tuner1 = 26;
    static __thread integer tuner2 = 27;
    static __thread integer tuner3 = 28;
    static __thread integer xctol = 33;
    static __thread integer xftol = 34;

    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static __thread integer i__, nfc;
    static __thread doublereal gts, emax, xmax, rfac1, emaxs;
    static __thread logical goodx;


/*  ***  assess candidate step (***sol version 2.3)  *** */


/*  ***  purpose  *** */

/*        this subroutine is called by an unconstrained minimization */
/*     routine to assess the next candidate step.  it may recommend one */
/*     of several courses of action, such as accepting the step, recom- */
/*     puting it using the same or a new quadratic model, or halting due */
/*     to convergence or false convergence.  see the return code listing */
/*     below. */

/* --------------------------  parameter usage  -------------------------- */

/*  iv (i/o) integer parameter and scratch vector -- see description */
/*             below of iv values referenced. */
/* liv (in)  length of iv array. */
/*  lv (in)  length of v array. */
/*   v (i/o) real parameter and scratch vector -- see description */
/*             below of v values referenced. */

/*  ***  iv values referenced  *** */

/*    iv(irc) (i/o) on input for the first step tried in a new iteration, */
/*             iv(irc) should be set to 3 or 4 (the value to which it is */
/*             set when step is definitely to be accepted).  on input */
/*             after step has been recomputed, iv(irc) should be */
/*             unchanged since the previous return of assst. */
/*                on output, iv(irc) is a return code having one of the */
/*             following values... */
/*                  1 = switch models or try smaller step. */
/*                  2 = switch models or accept step. */
/*                  3 = accept step and determine v(radfac) by gradient */
/*                       tests. */
/*                  4 = accept step, v(radfac) has been determined. */
/*                  5 = recompute step (using the same model). */
/*                  6 = recompute step with radius = v(lmaxs) but do not */
/*                       evaulate the objective function. */
/*                  7 = x-convergence (see v(xctol)). */
/*                  8 = relative function convergence (see v(rfctol)). */
/*                  9 = both x- and relative function convergence. */
/*                 10 = absolute function convergence (see v(afctol)). */
/*                 11 = singular convergence (see v(lmaxs)). */
/*                 12 = false convergence (see v(xftol)). */
/*                 13 = iv(irc) was out of range on input. */
/*             return code i has precdence over i+1 for i = 9, 10, 11. */
/* iv(mlstgd) (i/o) saved value of iv(model). */
/*  iv(model) (i/o) on input, iv(model) should be an integer identifying */
/*             the current quadratic model of the objective function. */
/*             if a previous step yielded a better function reduction, */
/*             then iv(model) will be set to iv(mlstgd) on output. */
/* iv(nfcall) (in)  invocation count for the objective function. */
/* iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest */
/*             function reduction this iteration.  iv(nfgcal) remains */
/*             unchanged until a function reduction is obtained. */
/* iv(radinc) (i/o) the number of radius increases (or minus the number */
/*             of decreases) so far this iteration. */
/* iv(restor) (out) set to 1 if v(f) has been restored and x should be */
/*             restored to its initial value, to 2 if x should be saved, */
/*             to 3 if x should be restored from the saved value, and to */
/*             0 otherwise. */
/*  iv(stage) (i/o) count of the number of models tried so far in the */
/*             current iteration. */
/* iv(stglim) (in)  maximum number of models to consider. */
/* iv(switch) (out) set to 0 unless a new model is being tried and it */
/*             gives a smaller function value than the previous model, */
/*             in which case assst sets iv(switch) = 1. */
/* iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused */
/*             overflow). */
/*   iv(xirc) (i/o) value that iv(irc) would have in the absence of */
/*             convergence, false convergence, and oversized steps. */

/*  ***  v values referenced  *** */

/* v(afctol) (in)  absolute function convergence tolerance.  if the */
/*             absolute value of the current function value v(f) is less */
/*             than v(afctol), then assst returns with iv(irc) = 10. */
/* v(decfac) (in)  factor by which to decrease radius when iv(toobig) is */
/*             nonzero. */
/* v(dstnrm) (in)  the 2-norm of d*step. */
/* v(dstsav) (i/o) value of v(dstnrm) on saved step. */
/*   v(dst0) (in)  the 2-norm of d times the newton step (when defined, */
/*             i.e., for v(nreduc) .ge. 0). */
/*      v(f) (i/o) on both input and output, v(f) is the objective func- */
/*             tion value at x.  if x is restored to a previous value, */
/*             then v(f) is restored to the corresponding value. */
/*   v(fdif) (out) the function reduction v(f0) - v(f) (for the output */
/*             value of v(f) if an earlier step gave a bigger function */
/*             decrease, and for the input value of v(f) otherwise). */
/* v(flstgd) (i/o) saved value of v(f). */
/*     v(f0) (in)  objective function value at start of iteration. */
/* v(gtslst) (i/o) value of v(gtstep) on saved step. */
/* v(gtstep) (in)  inner product between step and gradient. */
/* v(incfac) (in)  minimum factor by which to increase radius. */
/*  v(lmaxs) (in)  maximum reasonable step size (and initial step bound). */
/*             if the actual function decrease is no more than twice */
/*             what was predicted, if a return with iv(irc) = 7, 8, 9, */
/*             or 10 does not occur, if v(dstnrm) .gt. v(lmaxs), and if */
/*             v(preduc) .le. v(sctol) * abs(v(f0)), then assst re- */
/*             turns with iv(irc) = 11.  if so doing appears worthwhile, */
/*             then assst repeats this test with v(preduc) computed for */
/*             a step of length v(lmaxs) (by a return with iv(irc) = 6). */
/* v(nreduc) (i/o)  function reduction predicted by quadratic model for */
/*             newton step.  if assst is called with iv(irc) = 6, i.e., */
/*             if v(preduc) has been computed with radius = v(lmaxs) for */
/*             use in the singular convervence test, then v(nreduc) is */
/*             set to -v(preduc) before the latter is restored. */
/* v(plstgd) (i/o) value of v(preduc) on saved step. */
/* v(preduc) (i/o) function reduction predicted by quadratic model for */
/*             current step. */
/* v(radfac) (out) factor to be used in determining the new radius, */
/*             which should be v(radfac)*dst, where  dst  is either the */
/*             output value of v(dstnrm) or the 2-norm of */
/*             diag(newd)*step  for the output value of step and the */
/*             updated version, newd, of the scale vector d.  for */
/*             iv(irc) = 3, v(radfac) = 1.0 is returned. */
/* v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input */
/*             value of v(dstnrm) -- suggested value = 0.1. */
/* v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0. */
/*  v(reldx) (in) scaled relative change in x caused by step, computed */
/*             (e.g.) by function  reldst  as */
/*                 max (d(i)*abs(x(i)-x0(i)), 1 .le. i .le. p) / */
/*                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p). */
/* v(rfctol) (in)  relative function convergence tolerance.  if the */
/*             actual function reduction is at most twice what was pre- */
/*             dicted and  v(nreduc) .le. v(rfctol)*abs(v(f0)),  then */
/*             assst returns with iv(irc) = 8 or 9. */
/* v(stppar) (in)  marquardt parameter -- 0 means full newton step. */
/* v(tuner1) (in)  tuning constant used to decide if the function */
/*             reduction was much less than expected.  suggested */
/*             value = 0.1. */
/* v(tuner2) (in)  tuning constant used to decide if the function */
/*             reduction was large enough to accept step.  suggested */
/*             value = 10**-4. */
/* v(tuner3) (in)  tuning constant used to decide if the radius */
/*             should be increased.  suggested value = 0.75. */
/*  v(xctol) (in)  x-convergence criterion.  if step is a newton step */
/*             (v(stppar) = 0) having v(reldx) .le. v(xctol) and giving */
/*             at most twice the predicted function decrease, then */
/*             assst returns iv(irc) = 7 or 9. */
/*  v(xftol) (in)  false convergence tolerance.  if step gave no or only */
/*             a small function decrease and v(reldx) .le. v(xftol), */
/*             then assst returns with iv(irc) = 12. */

/* -------------------------------  notes  ------------------------------- */

/*  ***  application and usage restrictions  *** */

/*        this routine is called as part of the nl2sol (nonlinear */
/*     least-squares) package.  it may be used in any unconstrained */
/*     minimization solver that uses dogleg, goldfeld-quandt-trotter, */
/*     or levenberg-marquardt steps. */

/*  ***  algorithm notes  *** */

/*        see (1) for further discussion of the assessing and model */
/*     switching strategies.  while nl2sol considers only two models, */
/*     assst is designed to handle any number of models. */

/*  ***  usage notes  *** */

/*        on the first call of an iteration, only the i/o variables */
/*     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and */
/*     v(preduc) need have been initialized.  between calls, no i/o */
/*     values execpt step, x, iv(model), v(f) and the stopping toler- */
/*     ances should be changed. */
/*        after a return for convergence or false convergence, one can */
/*     change the stopping tolerances and call assst again, in which */
/*     case the stopping tests will be repeated. */

/*  ***  references  *** */

/*     (1) dennis, j.e., jr., gay, d.m., and welsch, r.e. (1981), */
/*        an adaptive nonlinear least-squares algorithm, */
/*        acm trans. math. software, vol. 7, no. 3. */

/*     (2) powell, m.j.d. (1970)  a fortran subroutine for solving */
/*        systems of nonlinear algebraic equations, in numerical */
/*        methods for nonlinear algebraic equations, edited by */
/*        p. rabinowitz, gordon and breach, london. */

/*  ***  history  *** */

/*        john dennis designed much of this routine, starting with */
/*     ideas in (2). roy welsch suggested the model switching strategy. */
/*        david gay and stephen peters cast this subroutine into a more */
/*     portable form (winter 1977), and david gay cast it into its */
/*     present form (fall 1978). */

/*  ***  general  *** */

/*     this subroutine was written in connection with research */
/*     supported by the national science foundation under grants */
/*     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and */
/*     mcs-7906671. */

/* ------------------------  external quantities  ------------------------ */

/*  ***  no external functions and subroutines  *** */

/*  ***  intrinsic functions  *** */
/* /+ */
/* / */
/*  ***  no common blocks  *** */

/* --------------------------  local variables  -------------------------- */


/*  ***  subscripts for iv and v  *** */


/*  ***  data initializations  *** */

/* /6 */
    /* Parameter adjustments */
    --iv;
    --v;

    /* Function Body */
/* /7 */
/*     parameter (half=0.5d+0, one=1.d+0, onep2=1.2d+0, two=2.d+0, */
/*    1           zero=0.d+0) */
/* / */

/* /6 */
/* /7 */
/*     parameter (irc=29, mlstgd=32, model=5, nfcall=6, nfgcal=7, */
/*    1           radinc=8, restor=9, stage=10, stglim=11, switch=12, */
/*    2           toobig=2, xirc=13) */
/* / */
/* /6 */
/* /7 */
/*     parameter (afctol=31, decfac=22, dstnrm=2, dst0=3, dstsav=18, */
/*    1           f=10, fdif=11, flstgd=12, f0=13, gtslst=14, gtstep=4, */
/*    2           incfac=23, lmaxs=36, nreduc=6, plstgd=15, preduc=7, */
/*    3           radfac=16, rdfcmn=24, rdfcmx=25, reldx=17, rfctol=32, */
/*    4           sctol=37, stppar=5, tuner1=26, tuner2=27, tuner3=28, */
/*    5           xctol=33, xftol=34) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

    nfc = iv[(0 + (0 + (nfcall << 2))) / 4];
    iv[switch__] = 0;
    iv[restor] = 0;
    rfac1 = one;
    goodx = TRUE_;
    i__ = iv[irc];
    if (i__ >= 1 && i__ <= 12) {
	switch (i__) {
	    case 1:  goto L20;
	    case 2:  goto L30;
	    case 3:  goto L10;
	    case 4:  goto L10;
	    case 5:  goto L40;
	    case 6:  goto L280;
	    case 7:  goto L220;
	    case 8:  goto L220;
	    case 9:  goto L220;
	    case 10:  goto L220;
	    case 11:  goto L220;
	    case 12:  goto L170;
	}
    }
    iv[irc] = 13;
    goto L999;

/*  ***  initialize for new iteration  *** */

L10:
    iv[stage] = 1;
    iv[radinc] = 0;
    v[flstgd] = v[f0];
    if (iv[toobig] == 0) {
	goto L110;
    }
    iv[stage] = -1;
    iv[xirc] = i__;
    goto L60;

/*  ***  step was recomputed with new model or smaller radius  *** */
/*  ***  first decide which  *** */

L20:
    if (iv[model] != iv[mlstgd]) {
	goto L30;
    }
/*        ***  old model retained, smaller radius tried  *** */
/*        ***  do not consider any more new models this iteration  *** */
    iv[stage] = iv[stglim];
    iv[radinc] = -1;
    goto L110;

/*  ***  a new model is being tried.  decide whether to keep it.  *** */

L30:
    ++iv[stage];

/*     ***  now we add the possibiltiy that step was recomputed with  *** */
/*     ***  the same model, perhaps because of an oversized step.     *** */

L40:
    if (iv[stage] > 0) {
	goto L50;
    }

/*        ***  step was recomputed because it was too big.  *** */

    if (iv[toobig] != 0) {
	goto L60;
    }

/*        ***  restore iv(stage) and pick up where we left off.  *** */

    iv[stage] = -iv[stage];
    i__ = iv[xirc];
    switch (i__) {
	case 1:  goto L20;
	case 2:  goto L30;
	case 3:  goto L110;
	case 4:  goto L110;
	case 5:  goto L70;
    }

L50:
    if (iv[toobig] == 0) {
	goto L70;
    }

/*  ***  handle oversize step  *** */

    if (iv[radinc] > 0) {
	goto L80;
    }
    iv[stage] = -iv[stage];
    iv[xirc] = iv[irc];

L60:
    v[radfac] = v[decfac];
    --iv[radinc];
    iv[irc] = 5;
    iv[restor] = 1;
    goto L999;

L70:
    if (v[f] < v[flstgd]) {
	goto L110;
    }

/*     *** the new step is a loser.  restore old model.  *** */

    if (iv[model] == iv[mlstgd]) {
	goto L80;
    }
    iv[model] = iv[mlstgd];
    iv[switch__] = 1;

/*     ***  restore step, etc. only if a previous step decreased v(f). */

L80:
    if (v[flstgd] >= v[f0]) {
	goto L110;
    }
    iv[restor] = 1;
    v[f] = v[flstgd];
    v[preduc] = v[plstgd];
    v[gtstep] = v[gtslst];
    if (iv[switch__] == 0) {
	rfac1 = v[dstnrm] / v[dstsav];
    }
    v[dstnrm] = v[dstsav];
    nfc = iv[nfgcal];
    goodx = FALSE_;

L110:
    v[fdif] = v[f0] - v[f];
    if (v[fdif] > v[tuner2] * v[preduc]) {
	goto L140;
    }
    if (iv[radinc] > 0) {
	goto L140;
    }

/*        ***  no (or only a trivial) function decrease */
/*        ***  -- so try new model or smaller radius */

    if (v[f] < v[f0]) {
	goto L120;
    }
    iv[mlstgd] = iv[model];
    v[flstgd] = v[f];
    v[f] = v[f0];
    iv[restor] = 1;
    goto L130;
L120:
    iv[nfgcal] = nfc;
L130:
    iv[irc] = 1;
    if (iv[stage] < iv[stglim]) {
	goto L160;
    }
    iv[irc] = 5;
    --iv[radinc];
    goto L160;

/*  ***  nontrivial function decrease achieved  *** */

L140:
    iv[nfgcal] = nfc;
    rfac1 = one;
    v[dstsav] = v[dstnrm];
    if (v[fdif] > v[preduc] * v[tuner1]) {
	goto L190;
    }

/*  ***  decrease was much less than predicted -- either change models */
/*  ***  or accept step with decreased radius. */

    if (iv[stage] >= iv[stglim]) {
	goto L150;
    }
/*        ***  consider switching models  *** */
    iv[irc] = 2;
    goto L160;

/*     ***  accept step with decreased radius  *** */

L150:
    iv[irc] = 4;

/*  ***  set v(radfac) to fletcher*s decrease factor  *** */

L160:
    iv[xirc] = iv[irc];
    emax = v[gtstep] + v[fdif];
    v[radfac] = half * rfac1;
    if (emax < v[gtstep]) {
/* Computing MAX */
	d__1 = v[rdfcmn], d__2 = half * v[gtstep] / emax;
	v[radfac] = rfac1 * max(d__1,d__2);
    }

/*  ***  do false convergence test  *** */

L170:
    if (v[reldx] <= v[xftol]) {
	goto L180;
    }
    iv[irc] = iv[xirc];
    if (v[f] < v[f0]) {
	goto L200;
    }
    goto L230;

L180:
    iv[irc] = 12;
    goto L240;

/*  ***  handle good function decrease  *** */

L190:
    if (v[fdif] < -v[tuner3] * v[gtstep]) {
	goto L210;
    }

/*     ***  increasing radius looks worthwhile.  see if we just */
/*     ***  recomputed step with a decreased radius or restored step */
/*     ***  after recomputing it with a larger radius. */

    if (iv[radinc] < 0) {
	goto L210;
    }
    if (iv[restor] == 1) {
	goto L210;
    }

/*        ***  we did not.  try a longer step unless this was a newton */
/*        ***  step. */

    v[radfac] = v[rdfcmx];
    gts = v[gtstep];
    if (v[fdif] < (half / v[radfac] - one) * gts) {
/* Computing MAX */
	d__1 = v[incfac], d__2 = half * gts / (gts + v[fdif]);
	v[radfac] = max(d__1,d__2);
    }
    iv[irc] = 4;
    if (v[stppar] == zero) {
	goto L230;
    }
    if (v[dst0] >= zero && (v[dst0] < two * v[dstnrm] || v[nreduc] < onep2 * 
	    v[fdif])) {
	goto L230;
    }
/*             ***  step was not a newton step.  recompute it with */
/*             ***  a larger radius. */
    iv[irc] = 5;
    ++iv[radinc];

/*  ***  save values corresponding to good step  *** */

L200:
    v[flstgd] = v[f];
    iv[mlstgd] = iv[model];
    if (iv[restor] != 1) {
	iv[restor] = 2;
    }
    v[dstsav] = v[dstnrm];
    iv[nfgcal] = nfc;
    v[plstgd] = v[preduc];
    v[gtslst] = v[gtstep];
    goto L230;

/*  ***  accept step with radius unchanged  *** */

L210:
    v[radfac] = one;
    iv[irc] = 3;
    goto L230;

/*  ***  come here for a restart after convergence  *** */

L220:
    iv[irc] = iv[xirc];
    if (v[dstsav] >= zero) {
	goto L240;
    }
    iv[irc] = 12;
    goto L240;

/*  ***  perform convergence tests  *** */

L230:
    iv[xirc] = iv[irc];
L240:
    if (iv[restor] == 1 && v[flstgd] < v[f0]) {
	iv[restor] = 3;
    }
    if ((d__1 = v[f], abs(d__1)) < v[afctol]) {
	iv[irc] = 10;
    }
    if (half * v[fdif] > v[preduc]) {
	goto L999;
    }
    emax = v[rfctol] * (d__1 = v[f0], abs(d__1));
    emaxs = v[sctol] * (d__1 = v[f0], abs(d__1));
    if (v[dstnrm] > v[lmaxs] && v[preduc] <= emaxs) {
	iv[irc] = 11;
    }
    if (v[dst0] < zero) {
	goto L250;
    }
    i__ = 0;
    if (v[nreduc] > zero && v[nreduc] <= emax || v[nreduc] == zero && v[
	    preduc] == zero) {
	i__ = 2;
    }
    if (v[stppar] == zero && v[reldx] <= v[xctol] && goodx) {
	++i__;
    }
    if (i__ > 0) {
	iv[irc] = i__ + 6;
    }

/*  ***  consider recomputing step of length v(lmaxs) for singular */
/*  ***  convergence test. */

L250:
    if (iv[irc] > 5 && iv[irc] != 12) {
	goto L999;
    }
    if (v[dstnrm] > v[lmaxs]) {
	goto L260;
    }
    if (v[preduc] >= emaxs) {
	goto L999;
    }
    if (v[dst0] <= zero) {
	goto L270;
    }
    if (half * v[dst0] <= v[lmaxs]) {
	goto L999;
    }
    goto L270;
L260:
    if (half * v[dstnrm] <= v[lmaxs]) {
	goto L999;
    }
    xmax = v[lmaxs] / v[dstnrm];
    if (xmax * (two - xmax) * v[preduc] >= emaxs) {
	goto L999;
    }
L270:
    if (v[nreduc] < zero) {
	goto L290;
    }

/*  ***  recompute v(preduc) for use in singular convergence test  *** */

    v[gtslst] = v[gtstep];
    v[dstsav] = v[dstnrm];
    if (iv[irc] == 12) {
	v[dstsav] = -v[dstsav];
    }
    v[plstgd] = v[preduc];
    i__ = iv[restor];
    iv[restor] = 2;
    if (i__ == 3) {
	iv[restor] = 0;
    }
    iv[irc] = 6;
    goto L999;

/*  ***  perform singular convergence test with recomputed v(preduc)  *** */

L280:
    v[gtstep] = v[gtslst];
    v[dstnrm] = (d__1 = v[dstsav], abs(d__1));
    iv[irc] = iv[xirc];
    if (v[dstsav] <= zero) {
	iv[irc] = 12;
    }
    v[nreduc] = -v[preduc];
    v[preduc] = v[plstgd];
    iv[restor] = 3;
L290:
    if (-v[nreduc] <= v[rfctol] * (d__1 = v[f0], abs(d__1))) {
	iv[irc] = 11;
    }

L999:
    return 0;

/*  ***  last card of assst follows  *** */
} /* assst_ */

/* Subroutine */ int deflt_(integer *alg, integer *iv, integer *liv, integer *
	lv, doublereal *v)
{
    /* Initialized data */

    static __thread integer algsav = 51;
    static __thread integer covprt = 14;
    static __thread integer covreq = 15;
    static __thread integer dtype = 16;
    static __thread integer hc = 71;
    static __thread integer ierr = 75;
    static __thread integer inith = 25;
    static __thread integer inits = 25;
    static __thread integer ipivot = 76;
    static __thread integer ivneed = 3;
    static __thread integer lastiv = 44;
    static __thread integer lastv = 45;
    static __thread integer lmat = 42;
    static __thread integer mxfcal = 17;
    static __thread integer mxiter = 18;
    static __thread integer nfcov = 52;
    static __thread integer ngcov = 53;
    static __thread integer nvdflt = 50;
    static __thread integer outlev = 19;
    static __thread integer parprt = 20;
    static __thread integer parsav = 49;
    static __thread integer perm = 58;
    static __thread integer prunit = 21;
    static __thread integer qrtyp = 80;
    static __thread integer rdreq = 57;
    static __thread integer rmat = 78;
    static __thread integer solprt = 22;
    static __thread integer statpr = 23;
    static __thread integer vneed = 4;
    static __thread integer vsave = 60;
    static __thread integer x0prt = 24;
    static __thread integer miniv[2] = { 80,59 };
    static __thread integer minv[2] = { 98,71 };

    static __thread integer mv, miv;
    extern /* Subroutine */ int vdflt_(integer *, integer *, doublereal *);
    extern integer imdcon_(integer *);


/*  ***  supply ***sol (version 2.3) default values to iv and v  *** */

/*  ***  alg = 1 means regression constants. */
/*  ***  alg = 2 means general unconstrained optimization constants. */


/* imdcon... returns machine-dependent integer constants. */
/* vdflt.... provides default values to v. */


/*  ***  subscripts for iv  *** */


/*  ***  iv subscript values  *** */

/* /6 */
    /* Parameter adjustments */
    --iv;
    --v;

    /* Function Body */
/* /7 */
/*     parameter (algsav=51, covprt=14, covreq=15, dtype=16, hc=71, */
/*    1           ierr=75, inith=25, inits=25, ipivot=76, ivneed=3, */
/*    2           lastiv=44, lastv=45, lmat=42, mxfcal=17, mxiter=18, */
/*    3           nfcov=52, ngcov=53, nvdflt=50, outlev=19, parprt=20, */
/*    4           parsav=49, perm=58, prunit=21, qrtyp=80, rdreq=57, */
/*    5           rmat=78, solprt=22, statpr=23, vneed=4, vsave=60, */
/*    6           x0prt=24) */
/* / */

/* -------------------------------  body  -------------------------------- */

    if (*alg < 1 || *alg > 2) {
	goto L40;
    }
    miv = miniv[*alg - 1];
    if (*liv < miv) {
	goto L20;
    }
    mv = minv[*alg - 1];
    if (*lv < mv) {
	goto L30;
    }
    vdflt_(alg, lv, &v[1]);
    iv[1] = 12;
    iv[algsav] = *alg;
    iv[ivneed] = 0;
    iv[lastiv] = miv;
    iv[lastv] = mv;
    iv[lmat] = mv + 1;
    iv[mxfcal] = 200;
    iv[mxiter] = 150;
    iv[outlev] = 1;
    iv[parprt] = 1;
    iv[perm] = miv + 1;
    iv[prunit] = imdcon_(&c__1);
    iv[solprt] = 1;
    iv[statpr] = 1;
    iv[vneed] = 0;
    iv[x0prt] = 1;

    if (*alg >= 2) {
	goto L10;
    }

/*  ***  regression  values */

    iv[covprt] = 3;
    iv[covreq] = 1;
    iv[dtype] = 1;
    iv[hc] = 0;
    iv[ierr] = 0;
    iv[inits] = 0;
    iv[ipivot] = 0;
    iv[nvdflt] = 32;
    iv[parsav] = 67;
    iv[qrtyp] = 1;
    iv[rdreq] = 3;
    iv[rmat] = 0;
    iv[vsave] = 58;
    goto L999;

/*  ***  general optimization values */

L10:
    iv[dtype] = 0;
    iv[inith] = 1;
    iv[nfcov] = 0;
    iv[ngcov] = 0;
    iv[nvdflt] = 25;
    iv[parsav] = 47;
    goto L999;

L20:
    iv[1] = 15;
    goto L999;

L30:
    iv[1] = 16;
    goto L999;

L40:
    iv[1] = 67;

L999:
    return 0;
/*  ***  last card of deflt follows  *** */
} /* deflt_ */

doublereal dotprd_(integer *p, doublereal *x, doublereal *y)
{
    /* Initialized data */

    static __thread doublereal one = 1.;
    static __thread doublereal sqteta = 0.;
    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Local variables */
    static __thread integer i__;
    static __thread doublereal t;
    extern doublereal rmdcon_(integer *);


/*  ***  return the inner product of the p-vectors x and y.  *** */


/* /+ */
/* / */

/*  ***  rmdcon(2) returns a machine-dependent constant, sqteta, which */
/*  ***  is slightly larger than the smallest positive number that */
/*  ***  can be squared without underflowing. */

/* /6 */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
/* /7 */
/*     parameter (one=1.d+0, zero=0.d+0) */
/*     data sqteta/0.d+0/ */
/* / */

    ret_val = zero;
    if (*p <= 0) {
	goto L999;
    }
    if (sqteta == zero) {
	sqteta = rmdcon_(&c__2);
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__3 = (d__1 = x[i__], abs(d__1)), d__4 = (d__2 = y[i__], abs(d__2));
	t = max(d__3,d__4);
	if (t > one) {
	    goto L10;
	}
	if (t < sqteta) {
	    goto L20;
	}
	t = x[i__] / sqteta * y[i__];
	if (abs(t) < sqteta) {
	    goto L20;
	}
L10:
	ret_val += x[i__] * y[i__];
L20:
	;
    }

L999:
    return ret_val;
/*  ***  last card of dotprd follows  *** */
} /* dotprd_ */

/* Subroutine */ int itsum_(doublereal *d__, doublereal *g, integer *iv, 
	integer *liv, integer *lv, integer *p, doublereal *v, doublereal *x)
{
    /* Initialized data */

    static __thread integer algsav = 51;
    static __thread integer needhd = 36;
    static __thread integer nfcall = 6;
    static __thread integer nfcov = 52;
    static __thread integer ngcall = 30;
    static __thread integer ngcov = 53;
    static __thread integer niter = 31;
    static __thread integer outlev = 19;
    static __thread integer prntit = 39;
    static __thread integer prunit = 21;
    static __thread integer solprt = 22;
    static __thread integer statpr = 23;
    static __thread integer sused = 64;
    static __thread integer x0prt = 24;
    static __thread integer dstnrm = 2;
    static __thread integer f = 10;
    static __thread integer f0 = 13;
    static __thread integer fdif = 11;
    static __thread integer nreduc = 6;
    static __thread integer preduc = 7;
    static __thread integer reldx = 17;
    static __thread integer stppar = 5;
    static __thread doublereal zero = 0.;
    static __thread struct {
	char e_1[24];
	real e_2;
	} equiv_378 = { {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'g', ' ', ' ', ' ', 
		's', ' '}, (float)0. };

#define model1 ((real *)&equiv_378)

    static __thread struct {
	char e_1[24];
	real e_2;
	} equiv_379 = { {' ', 'g', ' ', ' ', ' ', 's', ' ', ' ', 'g', '-', 
		's', ' ', 's', '-', 'g', ' ', '-', 's', '-', 'g', '-', 'g', 
		'-', 's'}, (float)0. };

#define model2 ((real *)&equiv_379)


    /* Format strings */
    static __thread char fmt_30[] = "(/\002   it   nf\002,6x,\002f\002,7x,\002reld\
f\002,3x,\002preldf\002,3x,\002reldx\002,2x,\002model  stppar\002)";
    static __thread char fmt_40[] = "(/\002    it   nf\002,7x,\002f\002,8x,\002reld\
f\002,4x,\002preldf\002,4x,\002reldx\002,3x,\002stppar\002)";
    static __thread char fmt_100[] = "(i6,i5,d10.3,2d9.2,d8.1,a3,a4,2d8.1,d9.2)";
    static __thread char fmt_110[] = "(i6,i5,d11.3,2d10.2,3d9.1,d10.2)";
    static __thread char fmt_70[] = "(/\002    it   nf\002,6x,\002f\002,7x,\002reld\
f\002,3x,\002preldf\002,3x,\002reldx\002,2x,\002model  stppar\002,2x,\002d*s\
tep\002,2x,\002npreldf\002)";
    static __thread char fmt_80[] = "(/\002    it   nf\002,7x,\002f\002,8x,\002reld\
f\002,4x,\002preldf\002,4x,\002reldx\002,3x,\002stppar\002,3x,\002d*step\002\
,3x,\002npreldf\002)";
    static __thread char fmt_140[] = "(/\002 ***** x-convergence *****\002)";
    static __thread char fmt_160[] = "(/\002 ***** relative function convergence **\
***\002)";
    static __thread char fmt_180[] = "(/\002 ***** x- and relative function convergen\
ce *****\002)";
    static __thread char fmt_200[] = "(/\002 ***** absolute function convergence **\
***\002)";
    static __thread char fmt_220[] = "(/\002 ***** singular convergence *****\002)";
    static __thread char fmt_240[] = "(/\002 ***** false convergence *****\002)";
    static __thread char fmt_260[] = "(/\002 ***** function evaluation limit *****\
\002)";
    static __thread char fmt_280[] = "(/\002 ***** iteration limit *****\002)";
    static __thread char fmt_300[] = "(/\002 ***** stopx *****\002)";
    static __thread char fmt_320[] = "(/\002 ***** initial f(x) cannot be computed **\
***\002)";
    static __thread char fmt_340[] = "(/\002 ***** bad parameters to assess *****\002)"
	    ;
    static __thread char fmt_360[] = "(/\002 ***** gradient could not be computed ***\
**\002)";
    static __thread char fmt_380[] = "(/\002 ***** iv(1) =\002,i5,\002 *****\002)";
    static __thread char fmt_400[] = "(/\002     i     initial x(i)\002,8x,\002d(i\
)\002//(1x,i5,d17.6,d14.3))";
    static __thread char fmt_410[] = "(/\002     0    1\002,d10.3)";
    static __thread char fmt_420[] = "(/\002     0    1\002,d11.3)";
    static __thread char fmt_450[] = "(/\002 function\002,d17.6,\002   reldx\002,d17.\
3/\002 func. evals\002,i8,9x,\002grad. evals\002,i8/\002 preldf\002,d16.3,6x,\
\002npreldf\002,d15.3)";
    static __thread char fmt_460[] = "(/1x,i4,\002 extra func. evals for covariance a\
nd diagnostics.\002)";
    static __thread char fmt_470[] = "(1x,i4,\002 extra grad. evals for covariance an\
d diagnostics.\002)";
    static __thread char fmt_490[] = "(/\002     i      final x(i)\002,8x,\002d(i)\
\002,10x,\002g(i)\002/)";
    static __thread char fmt_510[] = "(1x,i5,d16.6,2d14.3)";
    static __thread char fmt_530[] = "(/\002 inconsistent dimensions\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static __thread integer i__, m, nf, ng, ol, pu, iv1, alg;
    static __thread doublereal oldf, reldf, nreldf, preldf;

    /* Fortran I/O blocks */
    static __thread cilist io___340 = { 0, 0, 0, fmt_30, 0 };
    static __thread cilist io___341 = { 0, 0, 0, fmt_40, 0 };
    static __thread cilist io___343 = { 0, 0, 0, fmt_100, 0 };
    static __thread cilist io___344 = { 0, 0, 0, fmt_110, 0 };
    static __thread cilist io___345 = { 0, 0, 0, fmt_70, 0 };
    static __thread cilist io___346 = { 0, 0, 0, fmt_80, 0 };
    static __thread cilist io___348 = { 0, 0, 0, fmt_100, 0 };
    static __thread cilist io___349 = { 0, 0, 0, fmt_110, 0 };
    static __thread cilist io___350 = { 0, 0, 0, fmt_140, 0 };
    static __thread cilist io___351 = { 0, 0, 0, fmt_160, 0 };
    static __thread cilist io___352 = { 0, 0, 0, fmt_180, 0 };
    static __thread cilist io___353 = { 0, 0, 0, fmt_200, 0 };
    static __thread cilist io___354 = { 0, 0, 0, fmt_220, 0 };
    static __thread cilist io___355 = { 0, 0, 0, fmt_240, 0 };
    static __thread cilist io___356 = { 0, 0, 0, fmt_260, 0 };
    static __thread cilist io___357 = { 0, 0, 0, fmt_280, 0 };
    static __thread cilist io___358 = { 0, 0, 0, fmt_300, 0 };
    static __thread cilist io___359 = { 0, 0, 0, fmt_320, 0 };
    static __thread cilist io___360 = { 0, 0, 0, fmt_340, 0 };
    static __thread cilist io___361 = { 0, 0, 0, fmt_360, 0 };
    static __thread cilist io___362 = { 0, 0, 0, fmt_380, 0 };
    static __thread cilist io___363 = { 0, 0, 0, fmt_400, 0 };
    static __thread cilist io___365 = { 0, 0, 0, fmt_30, 0 };
    static __thread cilist io___366 = { 0, 0, 0, fmt_40, 0 };
    static __thread cilist io___367 = { 0, 0, 0, fmt_70, 0 };
    static __thread cilist io___368 = { 0, 0, 0, fmt_80, 0 };
    static __thread cilist io___369 = { 0, 0, 0, fmt_410, 0 };
    static __thread cilist io___370 = { 0, 0, 0, fmt_420, 0 };
    static __thread cilist io___372 = { 0, 0, 0, fmt_450, 0 };
    static __thread cilist io___373 = { 0, 0, 0, fmt_460, 0 };
    static __thread cilist io___374 = { 0, 0, 0, fmt_470, 0 };
    static __thread cilist io___375 = { 0, 0, 0, fmt_490, 0 };
    static __thread cilist io___376 = { 0, 0, 0, fmt_510, 0 };
    static __thread cilist io___377 = { 0, 0, 0, fmt_530, 0 };



/*  ***  print iteration summary for ***sol (version 2.3)  *** */

/*  ***  parameter declarations  *** */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  local variables  *** */

/* /6 */
/* /7 */
/*     character*4 model1(6), model2(6) */
/* / */

/*  ***  intrinsic functions  *** */
/* /+ */
/* / */
/*  ***  no external functions or subroutines  *** */

/*  ***  subscripts for iv and v  *** */


/*  ***  iv subscript values  *** */

/* /6 */
    /* Parameter adjustments */
    --iv;
    --v;
    --x;
    --g;
    --d__;

    /* Function Body */
/* /7 */
/*     parameter (algsav=51, needhd=36, nfcall=6, nfcov=52, ngcall=30, */
/*    1           ngcov=53, niter=31, outlev=19, prntit=39, prunit=21, */
/*    2           solprt=22, statpr=23, sused=64, x0prt=24) */
/* / */

/*  ***  v subscript values  *** */

/* /6 */
/* /7 */
/*     parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6, preduc=7, */
/*    1           reldx=17, stppar=5) */
/* / */

/* /6 */
/* /7 */
/*     parameter (zero=0.d+0) */
/* / */
/* /6 */
/* /7 */
/*     data model1/'    ','    ','    ','    ','  g ','  s '/, */
/*    1     model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/ */
/* / */

/* -------------------------------  body  -------------------------------- */

    pu = iv[(0 + (0 + (prunit << 2))) / 4];
    if (pu == 0) {
	goto L999;
    }
    iv1 = iv[1];
    if (iv1 > 62) {
	iv1 += -51;
    }
    ol = iv[outlev];
    alg = iv[algsav];
    if (iv1 < 2 || iv1 > 15) {
	goto L370;
    }
    if (iv1 >= 12) {
	goto L120;
    }
    if (iv1 == 2 && iv[niter] == 0) {
	goto L390;
    }
    if (ol == 0) {
	goto L120;
    }
    if (iv1 >= 10 && iv[prntit] == 0) {
	goto L120;
    }
    if (iv1 > 2) {
	goto L10;
    }
    ++iv[prntit];
    if (iv[prntit] < abs(ol)) {
	goto L999;
    }
L10:
    nf = iv[nfcall] - (i__1 = iv[nfcov], abs(i__1));
    iv[prntit] = 0;
    reldf = zero;
    preldf = zero;
/* Computing MAX */
    d__3 = (d__1 = v[f0], abs(d__1)), d__4 = (d__2 = v[f], abs(d__2));
    oldf = max(d__3,d__4);
    if (oldf <= zero) {
	goto L20;
    }
    reldf = v[fdif] / oldf;
    preldf = v[preduc] / oldf;
L20:
    if (ol > 0) {
	goto L60;
    }

/*        ***  print short summary line  *** */

    if (iv[needhd] == 1 && alg == 1) {
	io___340.ciunit = pu;
	s_wsfe(&io___340);
	e_wsfe();
    }
    if (iv[needhd] == 1 && alg == 2) {
	io___341.ciunit = pu;
	s_wsfe(&io___341);
	e_wsfe();
    }
    iv[needhd] = 0;
    if (alg == 2) {
	goto L50;
    }
    m = iv[sused];
    io___343.ciunit = pu;
    s_wsfe(&io___343);
    do_fio(&c__1, (char *)&iv[niter], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[f], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[reldx], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&model1[m - 1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&model2[m - 1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&v[stppar], (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L120;

L50:
    io___344.ciunit = pu;
    s_wsfe(&io___344);
    do_fio(&c__1, (char *)&iv[niter], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[f], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[reldx], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[stppar], (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L120;

/*     ***  print long summary line  *** */

L60:
    if (iv[needhd] == 1 && alg == 1) {
	io___345.ciunit = pu;
	s_wsfe(&io___345);
	e_wsfe();
    }
    if (iv[needhd] == 1 && alg == 2) {
	io___346.ciunit = pu;
	s_wsfe(&io___346);
	e_wsfe();
    }
    iv[needhd] = 0;
    nreldf = zero;
    if (oldf > zero) {
	nreldf = v[nreduc] / oldf;
    }
    if (alg == 2) {
	goto L90;
    }
    m = iv[sused];
    io___348.ciunit = pu;
    s_wsfe(&io___348);
    do_fio(&c__1, (char *)&iv[niter], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[f], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[reldx], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&model1[m - 1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&model2[m - 1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&v[stppar], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[dstnrm], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nreldf, (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L120;

L90:
    io___349.ciunit = pu;
    s_wsfe(&io___349);
    do_fio(&c__1, (char *)&iv[niter], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[f], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[reldx], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[stppar], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[dstnrm], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nreldf, (ftnlen)sizeof(doublereal));
    e_wsfe();

L120:
    if (iv[statpr] < 0) {
	goto L430;
    }
    switch (iv1) {
	case 1:  goto L999;
	case 2:  goto L999;
	case 3:  goto L130;
	case 4:  goto L150;
	case 5:  goto L170;
	case 6:  goto L190;
	case 7:  goto L210;
	case 8:  goto L230;
	case 9:  goto L250;
	case 10:  goto L270;
	case 11:  goto L290;
	case 12:  goto L310;
	case 13:  goto L330;
	case 14:  goto L350;
	case 15:  goto L520;
    }

L130:
    io___350.ciunit = pu;
    s_wsfe(&io___350);
    e_wsfe();
    goto L430;

L150:
    io___351.ciunit = pu;
    s_wsfe(&io___351);
    e_wsfe();
    goto L430;

L170:
    io___352.ciunit = pu;
    s_wsfe(&io___352);
    e_wsfe();
    goto L430;

L190:
    io___353.ciunit = pu;
    s_wsfe(&io___353);
    e_wsfe();
    goto L430;

L210:
    io___354.ciunit = pu;
    s_wsfe(&io___354);
    e_wsfe();
    goto L430;

L230:
    io___355.ciunit = pu;
    s_wsfe(&io___355);
    e_wsfe();
    goto L430;

L250:
    io___356.ciunit = pu;
    s_wsfe(&io___356);
    e_wsfe();
    goto L430;

L270:
    io___357.ciunit = pu;
    s_wsfe(&io___357);
    e_wsfe();
    goto L430;

L290:
    io___358.ciunit = pu;
    s_wsfe(&io___358);
    e_wsfe();
    goto L430;

L310:
    io___359.ciunit = pu;
    s_wsfe(&io___359);
    e_wsfe();

    goto L390;

L330:
    io___360.ciunit = pu;
    s_wsfe(&io___360);
    e_wsfe();
    goto L999;

L350:
    io___361.ciunit = pu;
    s_wsfe(&io___361);
    e_wsfe();
    if (iv[niter] > 0) {
	goto L480;
    }
    goto L390;

L370:
    io___362.ciunit = pu;
    s_wsfe(&io___362);
    do_fio(&c__1, (char *)&iv[1], (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;

/*  ***  initial call on itsum  *** */

L390:
    if (iv[x0prt] != 0) {
	io___363.ciunit = pu;
	s_wsfe(&io___363);
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&d__[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
/*     *** the following are to avoid undefined variables when the */
/*     *** function evaluation limit is 1... */
    v[dstnrm] = zero;
    v[fdif] = zero;
    v[nreduc] = zero;
    v[preduc] = zero;
    v[reldx] = zero;
    if (iv1 >= 12) {
	goto L999;
    }
    iv[needhd] = 0;
    iv[prntit] = 0;
    if (ol == 0) {
	goto L999;
    }
    if (ol < 0 && alg == 1) {
	io___365.ciunit = pu;
	s_wsfe(&io___365);
	e_wsfe();
    }
    if (ol < 0 && alg == 2) {
	io___366.ciunit = pu;
	s_wsfe(&io___366);
	e_wsfe();
    }
    if (ol > 0 && alg == 1) {
	io___367.ciunit = pu;
	s_wsfe(&io___367);
	e_wsfe();
    }
    if (ol > 0 && alg == 2) {
	io___368.ciunit = pu;
	s_wsfe(&io___368);
	e_wsfe();
    }
    if (alg == 1) {
	io___369.ciunit = pu;
	s_wsfe(&io___369);
	do_fio(&c__1, (char *)&v[f], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (alg == 2) {
	io___370.ciunit = pu;
	s_wsfe(&io___370);
	do_fio(&c__1, (char *)&v[f], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/* 365  format(/11h     0    1,e11.3) */
    goto L999;

/*  ***  print various information requested on solution  *** */

L430:
    iv[needhd] = 1;
    if (iv[statpr] == 0) {
	goto L480;
    }
/* Computing MAX */
    d__3 = (d__1 = v[f0], abs(d__1)), d__4 = (d__2 = v[f], abs(d__2));
    oldf = max(d__3,d__4);
    preldf = zero;
    nreldf = zero;
    if (oldf <= zero) {
	goto L440;
    }
    preldf = v[preduc] / oldf;
    nreldf = v[nreduc] / oldf;
L440:
    nf = iv[nfcall] - iv[nfcov];
    ng = iv[ngcall] - iv[ngcov];
    io___372.ciunit = pu;
    s_wsfe(&io___372);
    do_fio(&c__1, (char *)&v[f], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[reldx], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ng, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nreldf, (ftnlen)sizeof(doublereal));
    e_wsfe();

    if (iv[nfcov] > 0) {
	io___373.ciunit = pu;
	s_wsfe(&io___373);
	do_fio(&c__1, (char *)&iv[nfcov], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (iv[ngcov] > 0) {
	io___374.ciunit = pu;
	s_wsfe(&io___374);
	do_fio(&c__1, (char *)&iv[ngcov], (ftnlen)sizeof(integer));
	e_wsfe();
    }

L480:
    if (iv[solprt] == 0) {
	goto L999;
    }
    iv[needhd] = 1;
    io___375.ciunit = pu;
    s_wsfe(&io___375);
    e_wsfe();
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___376.ciunit = pu;
	s_wsfe(&io___376);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&d__[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L500: */
    }
    goto L999;

L520:
    io___377.ciunit = pu;
    s_wsfe(&io___377);
    e_wsfe();
L999:
    return 0;
/*  ***  last card of itsum follows  *** */
} /* itsum_ */

#undef model2
#undef model1


/* Subroutine */ int litvmu_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* Initialized data */

    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static __thread integer i__, j, i0, ii, ij;
    static __thread doublereal xi;
    static __thread integer im1, np1;


/*  ***  solve  (l**t)*x = y,  where  l  is an  n x n  lower triangular */
/*  ***  matrix stored compactly by rows.  x and y may occupy the same */
/*  ***  storage.  *** */

/* /6 */
    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
/* /7 */
/*     parameter (zero=0.d+0) */
/* / */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = y[i__];
    }
    np1 = *n + 1;
    i0 = *n * (*n + 1) / 2;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = np1 - ii;
	xi = x[i__] / l[i0];
	x[i__] = xi;
	if (i__ <= 1) {
	    goto L999;
	}
	i0 -= i__;
	if (xi == zero) {
	    goto L30;
	}
	im1 = i__ - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    ij = i0 + j;
	    x[j] -= xi * l[ij];
/* L20: */
	}
L30:
	;
    }
L999:
    return 0;
/*  ***  last card of litvmu follows  *** */
} /* litvmu_ */

/* Subroutine */ int livmul_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* Initialized data */

    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static __thread integer i__, j, k;
    static __thread doublereal t;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);


/*  ***  solve  l*x = y, where  l  is an  n x n  lower triangular */
/*  ***  matrix stored compactly by rows.  x and y may occupy the same */
/*  ***  storage.  *** */

/* /6 */
    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
/* /7 */
/*     parameter (zero=0.d+0) */
/* / */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (y[k] != zero) {
	    goto L20;
	}
	x[k] = zero;
/* L10: */
    }
    goto L999;
L20:
    j = k * (k + 1) / 2;
    x[k] = y[k] / l[j];
    if (k >= *n) {
	goto L999;
    }
    ++k;
    i__1 = *n;
    for (i__ = k; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	t = dotprd_(&i__2, &l[j + 1], &x[1]);
	j += i__;
	x[i__] = (y[i__] - t) / l[j];
/* L30: */
    }
L999:
    return 0;
/*  ***  last card of livmul follows  *** */
} /* livmul_ */

/* Subroutine */ int parck_(integer *alg, doublereal *d__, integer *iv, 
	integer *liv, integer *lv, integer *n, doublereal *v)
{
    /* Initialized data */

    static __thread integer algsav = 51;
    static __thread integer dinit = 38;
    static __thread integer dtype = 16;
    static __thread integer dtype0 = 54;
    static __thread integer epslon = 19;
    static __thread integer inits = 25;
    static __thread integer ivneed = 3;
    static __thread integer lastiv = 44;
    static __thread integer lastv = 45;
    static __thread integer lmat = 42;
    static __thread integer nextiv = 46;
    static __thread integer nextv = 47;
    static __thread integer nvdflt = 50;
    static __thread integer oldn = 38;
    static __thread integer parprt = 20;
    static __thread integer parsav = 49;
    static __thread integer perm = 58;
    static __thread integer prunit = 21;
    static __thread integer vneed = 4;
    static __thread doublereal big = 0.;
    static __thread doublereal machep = -1.;
    static __thread doublereal tiny = 1.;
    static __thread doublereal zero = 0.;
    static __thread struct {
	char e_1[272];
	real e_2;
	} equiv_457 = { {'e', 'p', 's', 'l', 'o', 'n', '.', '.', 'p', 'h', 
		'm', 'n', 'f', 'c', '.', '.', 'p', 'h', 'm', 'x', 'f', 'c', 
		'.', '.', 'd', 'e', 'c', 'f', 'a', 'c', '.', '.', 'i', 'n', 
		'c', 'f', 'a', 'c', '.', '.', 'r', 'd', 'f', 'c', 'm', 'n', 
		'.', '.', 'r', 'd', 'f', 'c', 'm', 'x', '.', '.', 't', 'u', 
		'n', 'e', 'r', '1', '.', '.', 't', 'u', 'n', 'e', 'r', '2', 
		'.', '.', 't', 'u', 'n', 'e', 'r', '3', '.', '.', 't', 'u', 
		'n', 'e', 'r', '4', '.', '.', 't', 'u', 'n', 'e', 'r', '5', 
		'.', '.', 'a', 'f', 'c', 't', 'o', 'l', '.', '.', 'r', 'f', 
		'c', 't', 'o', 'l', '.', '.', 'x', 'c', 't', 'o', 'l', '.', 
		'.', '.', 'x', 'f', 't', 'o', 'l', '.', '.', '.', 'l', 'm', 
		'a', 'x', '0', '.', '.', '.', 'l', 'm', 'a', 'x', 's', '.', 
		'.', '.', 's', 'c', 't', 'o', 'l', '.', '.', '.', 'd', 'i', 
		'n', 'i', 't', '.', '.', '.', 'd', 't', 'i', 'n', 'i', 't', 
		'.', '.', 'd', '0', 'i', 'n', 'i', 't', '.', '.', 'd', 'f', 
		'a', 'c', '.', '.', '.', '.', 'd', 'l', 't', 'f', 'd', 'c', 
		'.', '.', 'd', 'l', 't', 'f', 'd', 'j', '.', '.', 'd', 'e', 
		'l', 't', 'a', '0', '.', '.', 'f', 'u', 'z', 'z', '.', '.', 
		'.', '.', 'r', 'l', 'i', 'm', 'i', 't', '.', '.', 'c', 'o', 
		's', 'm', 'i', 'n', '.', '.', 'h', 'u', 'b', 'e', 'r', 'c', 
		'.', '.', 'r', 's', 'p', 't', 'o', 'l', '.', '.', 's', 'i', 
		'g', 'm', 'i', 'n', '.', '.', 'e', 't', 'a', '0', '.', '.', 
		'.', '.', 'b', 'i', 'a', 's', '.', '.', '.', '.'}, (float)0. }
		;

#define vn ((real *)&equiv_457)

    static __thread doublereal vm[34] = { .001,-.99,.001,.01,1.2,.01,1.2,0.,0.,.001,
	    -1.,0.0,0.,0.0,0.,0.,0.0,0.0,0.,-10.,0.,0.,0.,0.0,0.0,0.0,1.01,
	    1e10,0.0,0.,0.,0.,0.0,0. };
    static __thread doublereal vx[34] = { .9,-.001,10.,.8,100.,.8,100.,.5,.5,1.,1.,0.0,
	    0.0,.1,1.,1.,0.0,0.0,1.,0.0,0.0,0.0,1.,1.,1.,1.,1e10,0.0,1.,0.0,
	    1.,1.,1.,1. };
    static __thread struct {
	char e_1[8];
	integer e_2;
	} equiv_458 = { {'p', ' ', ' ', ' ', 'n', ' ', ' ', ' '}, 0 };

#define varnm ((integer *)&equiv_458)

    static __thread struct {
	char e_1[8];
	integer e_2;
	} equiv_459 = { {'s', ' ', ' ', ' ', 'h', ' ', ' ', ' '}, 0 };

#define sh ((integer *)&equiv_459)

    static __thread struct {
	char e_1[12];
	real e_2;
	} equiv_460 = { {'-', '-', '-', 'c', 'h', 'a', 'n', 'g', 'e', 'd', 
		' ', 'v'}, (float)0. };

#define cngd ((real *)&equiv_460)

    static __thread struct {
	char e_1[12];
	real e_2;
	} equiv_461 = { {'n', 'o', 'n', 'd', 'e', 'f', 'a', 'u', 'l', 't', 
		' ', 'v'}, (float)0. };

#define dflt ((real *)&equiv_461)

    static __thread integer ijmp = 33;
    static __thread integer jlim[2] = { 0,24 };
    static __thread integer ndflt[2] = { 32,25 };
    static __thread integer miniv[2] = { 80,59 };

    /* Format strings */
    static __thread char fmt_20[] = "(/\002 the first parameter to deflt should be\
\002,i3,\002 rather than\002,i3)";
    static __thread char fmt_40[] = "(/\002 /// bad\002,a1,\002 =\002,i5)";
    static __thread char fmt_70[] = "(/\002 /// \002,1a1,\002 changed from \002,i5\
,\002 to \002,i5)";
    static __thread char fmt_90[] = "(/\002 ///  iv(1) =\002,i5,\002 should be betwee\
n 0 and 14.\002)";
    static __thread char fmt_130[] = "(/\002 ///  \002,2a4,\002.. v(\002,i2,\002) \
=\002,d11.3,\002 should\002,\002 be between\002,d11.3,\002 and\002,d11.3)";
    static __thread char fmt_160[] = "(/\002 iv(nvdflt) =\002,i5,\002 rather than \
\002,i5)";
    static __thread char fmt_180[] = "(/\002 ///  d(\002,i3,\002) =\002,d11.3,\002 sh\
ould be positive\002)";
    static __thread char fmt_220[] = "(/\002 nondefault values....\002/\002 init\002,\
a1,\002..... iv(25) =\002,i3)";
    static __thread char fmt_260[] = "(/\002 \002,3a4,\002alues....\002/)";
    static __thread char fmt_240[] = "(\002 dtype..... iv(16) =\002,i3)";
    static __thread char fmt_270[] = "(1x,2a4,\002.. v(\002,i2,\002) =\002,d15.7)";
    static __thread char fmt_310[] = "(/\002 /// liv =\002,i5,\002 must be at leas\
t\002,i5)";
    static __thread char fmt_330[] = "(/\002 /// lv =\002,i5,\002 must be at least\
\002,i5)";
    static __thread char fmt_350[] = "(/\002 /// alg =\002,i5,\002 must be 1 or 2\002)"
	    ;

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static __thread integer i__, j, k, l, m, ii;
    static __thread doublereal vk;
    static __thread integer pu, iv1, miv1, miv2;
    extern /* Subroutine */ int deflt_(integer *, integer *, integer *, 
	    integer *, doublereal *);
    static __thread real which[3];
    extern /* Subroutine */ int vdflt_(integer *, integer *, doublereal *), 
	    vcopy_(integer *, doublereal *, doublereal *);
    static __thread integer parsv1, ndfalt;
    extern doublereal rmdcon_(integer *);

    /* Fortran I/O blocks */
    static __thread cilist io___432 = { 0, 0, 0, fmt_20, 0 };
    static __thread cilist io___433 = { 0, 0, 0, fmt_40, 0 };
    static __thread cilist io___436 = { 0, 0, 0, fmt_70, 0 };
    static __thread cilist io___437 = { 0, 0, 0, fmt_90, 0 };
    static __thread cilist io___444 = { 0, 0, 0, fmt_130, 0 };
    static __thread cilist io___445 = { 0, 0, 0, fmt_160, 0 };
    static __thread cilist io___446 = { 0, 0, 0, fmt_180, 0 };
    static __thread cilist io___447 = { 0, 0, 0, fmt_220, 0 };
    static __thread cilist io___448 = { 0, 0, 0, fmt_260, 0 };
    static __thread cilist io___449 = { 0, 0, 0, fmt_240, 0 };
    static __thread cilist io___451 = { 0, 0, 0, fmt_260, 0 };
    static __thread cilist io___452 = { 0, 0, 0, fmt_270, 0 };
    static __thread cilist io___454 = { 0, 0, 0, fmt_310, 0 };
    static __thread cilist io___455 = { 0, 0, 0, fmt_330, 0 };
    static __thread cilist io___456 = { 0, 0, 0, fmt_350, 0 };



/*  ***  check ***sol (version 2.3) parameters, print changed values  *** */

/*  ***  alg = 1 for regression, alg = 2 for general unconstrained opt. */


/* rmdcon -- returns machine-dependent constants. */
/* vcopy  -- copies one vector to another. */
/* vdflt  -- supplies default parameter values to v alone. */
/* /+ */
/* / */

/*  ***  local variables  *** */

/* /6 */
/* /7 */
/*     character*1 varnm(2), sh(2) */
/*     character*4 cngd(3), dflt(3), vn(2,34), which(3) */
/* / */

/*  ***  iv and v subscripts  *** */



/* /6 */
    /* Parameter adjustments */
    --iv;
    --v;
    --d__;

    /* Function Body */
/* /7 */
/*     parameter (algsav=51, dinit=38, dtype=16, dtype0=54, epslon=19, */
/*    1           inits=25, ivneed=3, lastiv=44, lastv=45, lmat=42, */
/*    2           nextiv=46, nextv=47, nvdflt=50, oldn=38, parprt=20, */
/*    3           parsav=49, perm=58, prunit=21, vneed=4) */
/*     save big, machep, tiny */
/* / */

/* /6 */
/* /7 */
/*     data vn(1,1),vn(2,1)/'epsl','on..'/ */
/*     data vn(1,2),vn(2,2)/'phmn','fc..'/ */
/*     data vn(1,3),vn(2,3)/'phmx','fc..'/ */
/*     data vn(1,4),vn(2,4)/'decf','ac..'/ */
/*     data vn(1,5),vn(2,5)/'incf','ac..'/ */
/*     data vn(1,6),vn(2,6)/'rdfc','mn..'/ */
/*     data vn(1,7),vn(2,7)/'rdfc','mx..'/ */
/*     data vn(1,8),vn(2,8)/'tune','r1..'/ */
/*     data vn(1,9),vn(2,9)/'tune','r2..'/ */
/*     data vn(1,10),vn(2,10)/'tune','r3..'/ */
/*     data vn(1,11),vn(2,11)/'tune','r4..'/ */
/*     data vn(1,12),vn(2,12)/'tune','r5..'/ */
/*     data vn(1,13),vn(2,13)/'afct','ol..'/ */
/*     data vn(1,14),vn(2,14)/'rfct','ol..'/ */
/*     data vn(1,15),vn(2,15)/'xcto','l...'/ */
/*     data vn(1,16),vn(2,16)/'xfto','l...'/ */
/*     data vn(1,17),vn(2,17)/'lmax','0...'/ */
/*     data vn(1,18),vn(2,18)/'lmax','s...'/ */
/*     data vn(1,19),vn(2,19)/'scto','l...'/ */
/*     data vn(1,20),vn(2,20)/'dini','t...'/ */
/*     data vn(1,21),vn(2,21)/'dtin','it..'/ */
/*     data vn(1,22),vn(2,22)/'d0in','it..'/ */
/*     data vn(1,23),vn(2,23)/'dfac','....'/ */
/*     data vn(1,24),vn(2,24)/'dltf','dc..'/ */
/*     data vn(1,25),vn(2,25)/'dltf','dj..'/ */
/*     data vn(1,26),vn(2,26)/'delt','a0..'/ */
/*     data vn(1,27),vn(2,27)/'fuzz','....'/ */
/*     data vn(1,28),vn(2,28)/'rlim','it..'/ */
/*     data vn(1,29),vn(2,29)/'cosm','in..'/ */
/*     data vn(1,30),vn(2,30)/'hube','rc..'/ */
/*     data vn(1,31),vn(2,31)/'rspt','ol..'/ */
/*     data vn(1,32),vn(2,32)/'sigm','in..'/ */
/*     data vn(1,33),vn(2,33)/'eta0','....'/ */
/*     data vn(1,34),vn(2,34)/'bias','....'/ */
/* / */


/* /6 */
/* /7 */
/*     data varnm(1)/'p'/, varnm(2)/'n'/, sh(1)/'s'/, sh(2)/'h'/ */
/*     data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/, */
/*    1     dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/ */
/* / */

/* ...............................  body  ................................ */

    pu = 0;
    if (prunit <= *liv) {
	pu = iv[prunit];
    }
    if (*alg < 1 || *alg > 2) {
	goto L340;
    }
    if (iv[1] == 0) {
	deflt_(alg, &iv[1], liv, lv, &v[1]);
    }
    iv1 = iv[1];
    if (iv1 != 13 && iv1 != 12) {
	goto L10;
    }
    miv1 = miniv[*alg - 1];
    if (perm <= *liv) {
/* Computing MAX */
	i__1 = miv1, i__2 = iv[perm] - 1;
	miv1 = max(i__1,i__2);
    }
    if (ivneed <= *liv) {
/* Computing MAX */
	i__1 = iv[ivneed];
	miv2 = miv1 + max(i__1,0);
    }
    if (lastiv <= *liv) {
	iv[lastiv] = miv2;
    }
    if (*liv < miv1) {
	goto L300;
    }
    iv[ivneed] = 0;
/* Computing MAX */
    i__1 = iv[vneed];
    iv[lastv] = max(i__1,0) + iv[lmat] - 1;
    iv[vneed] = 0;
    if (*liv < miv2) {
	goto L300;
    }
    if (*lv < iv[lastv]) {
	goto L320;
    }
L10:
    if (*alg == iv[algsav]) {
	goto L30;
    }
    if (pu != 0) {
	io___432.ciunit = pu;
	s_wsfe(&io___432);
	do_fio(&c__1, (char *)&(*alg), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iv[algsav], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    iv[1] = 82;
    goto L999;
L30:
    if (iv1 < 12 || iv1 > 14) {
	goto L60;
    }
    if (*n >= 1) {
	goto L50;
    }
    iv[1] = 81;
    if (pu == 0) {
	goto L999;
    }
    io___433.ciunit = pu;
    s_wsfe(&io___433);
    do_fio(&c__1, (char *)&varnm[*alg - 1], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;
L50:
    if (iv1 != 14) {
	iv[nextiv] = iv[perm];
    }
    if (iv1 != 14) {
	iv[nextv] = iv[lmat];
    }
    if (iv1 == 13) {
	goto L999;
    }
    k = iv[parsav] - epslon;
    i__1 = *lv - k;
    vdflt_(alg, &i__1, &v[k + 1]);
    iv[dtype0] = 2 - *alg;
    iv[oldn] = *n;
    which[0] = dflt[0];
    which[1] = dflt[1];
    which[2] = dflt[2];
    goto L110;
L60:
    if (*n == iv[oldn]) {
	goto L80;
    }
    iv[1] = 17;
    if (pu == 0) {
	goto L999;
    }
    io___436.ciunit = pu;
    s_wsfe(&io___436);
    do_fio(&c__1, (char *)&varnm[*alg - 1], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&iv[oldn], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;

L80:
    if (iv1 <= 11 && iv1 >= 1) {
	goto L100;
    }
    iv[1] = 80;
    if (pu != 0) {
	io___437.ciunit = pu;
	s_wsfe(&io___437);
	do_fio(&c__1, (char *)&iv1, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L999;

L100:
    which[0] = cngd[0];
    which[1] = cngd[1];
    which[2] = cngd[2];

L110:
    if (iv1 == 14) {
	iv1 = 12;
    }
    if (big > tiny) {
	goto L120;
    }
    tiny = rmdcon_(&c__1);
    machep = rmdcon_(&c__3);
    big = rmdcon_(&c__6);
    vm[11] = machep;
    vx[11] = big;
    vx[12] = big;
    vm[13] = machep;
    vm[16] = tiny;
    vx[16] = big;
    vm[17] = tiny;
    vx[17] = big;
    vx[19] = big;
    vx[20] = big;
    vx[21] = big;
    vm[23] = machep;
    vm[24] = machep;
    vm[25] = machep;
    vx[27] = rmdcon_(&c__5);
    vm[28] = machep;
    vx[29] = big;
    vm[32] = machep;
L120:
    m = 0;
    i__ = 1;
    j = jlim[*alg - 1];
    k = epslon;
    ndfalt = ndflt[*alg - 1];
    i__1 = ndfalt;
    for (l = 1; l <= i__1; ++l) {
	vk = v[k];
	if (vk >= vm[i__ - 1] && vk <= vx[i__ - 1]) {
	    goto L140;
	}
	m = k;
	if (pu != 0) {
	    io___444.ciunit = pu;
	    s_wsfe(&io___444);
	    do_fio(&c__1, (char *)&vn[(i__ << 1) - 2], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&vn[(i__ << 1) - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&vk, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&vm[i__ - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&vx[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
L140:
	++k;
	++i__;
	if (i__ == j) {
	    i__ = ijmp;
	}
/* L150: */
    }

    if (iv[nvdflt] == ndfalt) {
	goto L170;
    }
    iv[1] = 51;
    if (pu == 0) {
	goto L999;
    }
    io___445.ciunit = pu;
    s_wsfe(&io___445);
    do_fio(&c__1, (char *)&iv[nvdflt], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ndfalt, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;
L170:
    if ((iv[dtype] > 0 || v[dinit] > zero) && iv1 == 12) {
	goto L200;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] > zero) {
	    goto L190;
	}
	m = 18;
	if (pu != 0) {
	    io___446.ciunit = pu;
	    s_wsfe(&io___446);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&d__[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
L190:
	;
    }
L200:
    if (m == 0) {
	goto L210;
    }
    iv[1] = m;
    goto L999;

L210:
    if (pu == 0 || iv[parprt] == 0) {
	goto L999;
    }
    if (iv1 != 12 || iv[inits] == *alg - 1) {
	goto L230;
    }
    m = 1;
    io___447.ciunit = pu;
    s_wsfe(&io___447);
    do_fio(&c__1, (char *)&sh[*alg - 1], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&iv[inits], (ftnlen)sizeof(integer));
    e_wsfe();
L230:
    if (iv[dtype] == iv[dtype0]) {
	goto L250;
    }
    if (m == 0) {
	io___448.ciunit = pu;
	s_wsfe(&io___448);
	do_fio(&c__3, (char *)&which[0], (ftnlen)sizeof(real));
	e_wsfe();
    }
    m = 1;
    io___449.ciunit = pu;
    s_wsfe(&io___449);
    do_fio(&c__1, (char *)&iv[dtype], (ftnlen)sizeof(integer));
    e_wsfe();
L250:
    i__ = 1;
    j = jlim[*alg - 1];
    k = epslon;
    l = iv[parsav];
    ndfalt = ndflt[*alg - 1];
    i__1 = ndfalt;
    for (ii = 1; ii <= i__1; ++ii) {
	if (v[k] == v[l]) {
	    goto L280;
	}
	if (m == 0) {
	    io___451.ciunit = pu;
	    s_wsfe(&io___451);
	    do_fio(&c__3, (char *)&which[0], (ftnlen)sizeof(real));
	    e_wsfe();
	}
	m = 1;
	io___452.ciunit = pu;
	s_wsfe(&io___452);
	do_fio(&c__1, (char *)&vn[(i__ << 1) - 2], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vn[(i__ << 1) - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&v[k], (ftnlen)sizeof(doublereal));
	e_wsfe();
L280:
	++k;
	++l;
	++i__;
	if (i__ == j) {
	    i__ = ijmp;
	}
/* L290: */
    }

    iv[dtype0] = iv[dtype];
    parsv1 = iv[parsav];
    vcopy_(&iv[nvdflt], &v[parsv1], &v[epslon]);
    goto L999;

L300:
    iv[1] = 15;
    if (pu == 0) {
	goto L999;
    }
    io___454.ciunit = pu;
    s_wsfe(&io___454);
    do_fio(&c__1, (char *)&(*liv), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&miv2, (ftnlen)sizeof(integer));
    e_wsfe();
    if (*liv < miv1) {
	goto L999;
    }
    if (*lv < iv[lastv]) {
	goto L320;
    }
    goto L999;

L320:
    iv[1] = 16;
    if (pu == 0) {
	goto L999;
    }
    io___455.ciunit = pu;
    s_wsfe(&io___455);
    do_fio(&c__1, (char *)&(*lv), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&iv[lastv], (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;

L340:
    iv[1] = 67;
    if (pu == 0) {
	goto L999;
    }
    io___456.ciunit = pu;
    s_wsfe(&io___456);
    do_fio(&c__1, (char *)&(*alg), (ftnlen)sizeof(integer));
    e_wsfe();

L999:
    return 0;
/*  ***  last card of parck follows  *** */
} /* parck_ */

#undef dflt
#undef cngd
#undef sh
#undef varnm
#undef vn


doublereal reldst_(integer *p, doublereal *d__, doublereal *x, doublereal *x0)
{
    /* Initialized data */

    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static __thread integer i__;
    static __thread doublereal t, emax, xmax;


/*  ***  compute and return relative difference between x and x0  *** */
/*  ***  nl2sol version 2.2  *** */

/* /+ */
/* / */
/* /6 */
    /* Parameter adjustments */
    --x0;
    --x;
    --d__;

    /* Function Body */
/* /7 */
/*     parameter (zero=0.d+0) */
/* / */

    emax = zero;
    xmax = zero;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = (d__1 = d__[i__] * (x[i__] - x0[i__]), abs(d__1));
	if (emax < t) {
	    emax = t;
	}
	t = d__[i__] * ((d__1 = x[i__], abs(d__1)) + (d__2 = x0[i__], abs(
		d__2)));
	if (xmax < t) {
	    xmax = t;
	}
/* L10: */
    }
    ret_val = zero;
    if (xmax > zero) {
	ret_val = emax / xmax;
    }
/* L999: */
    return ret_val;
/*  ***  last card of reldst follows  *** */
} /* reldst_ */

logical stopx_(integer *idummy)
{
    /* System generated locals */
    logical ret_val;

/*     *****parameters... */

/*     .................................................................. */

/*     *****purpose... */
/*     this function may serve as the stopx (asynchronous interruption) */
/*     function for the nl2sol (nonlinear least-squares) package at */
/*     those installations which do not wish to implement a */
/*     dynamic stopx. */

/*     *****algorithm notes... */
/*     at installations where the nl2sol system is used */
/*     interactively, this dummy stopx should be replaced by a */
/*     function that returns .true. if and only if the interrupt */
/*     (break) key has been pressed since the last call on stopx. */

/*     .................................................................. */

    ret_val = FALSE_;
    return ret_val;
} /* stopx_ */

/* Subroutine */ int vaxpy_(integer *p, doublereal *w, doublereal *a, 
	doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static __thread integer i__;


/*  ***  set w = a*x + y  --  w, x, y = p-vectors, a = scalar  *** */



    /* Parameter adjustments */
    --y;
    --x;
    --w;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	w[i__] = *a * x[i__] + y[i__];
    }
    return 0;
} /* vaxpy_ */

/* Subroutine */ int vcopy_(integer *p, doublereal *y, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static __thread integer i__;


/*  ***  set y = x, where x and y are p-vectors  *** */



    /* Parameter adjustments */
    --x;
    --y;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = x[i__];
    }
    return 0;
} /* vcopy_ */

/* Subroutine */ int vdflt_(integer *alg, integer *lv, doublereal *v)
{
    /* Initialized data */

    static __thread doublereal one = 1.;
    static __thread doublereal three = 3.;
    static __thread integer afctol = 31;
    static __thread integer bias = 43;
    static __thread integer cosmin = 47;
    static __thread integer decfac = 22;
    static __thread integer delta0 = 44;
    static __thread integer dfac = 41;
    static __thread integer dinit = 38;
    static __thread integer dltfdc = 42;
    static __thread integer dltfdj = 43;
    static __thread integer dtinit = 39;
    static __thread integer d0init = 40;
    static __thread integer epslon = 19;
    static __thread integer eta0 = 42;
    static __thread integer fuzz = 45;
    static __thread integer huberc = 48;
    static __thread integer incfac = 23;
    static __thread integer lmax0 = 35;
    static __thread integer lmaxs = 36;
    static __thread integer phmnfc = 20;
    static __thread integer phmxfc = 21;
    static __thread integer rdfcmn = 24;
    static __thread integer rdfcmx = 25;
    static __thread integer rfctol = 32;
    static __thread integer rlimit = 46;
    static __thread integer rsptol = 49;
    static __thread integer sctol = 37;
    static __thread integer sigmin = 50;
    static __thread integer tuner1 = 26;
    static __thread integer tuner2 = 27;
    static __thread integer tuner3 = 28;
    static __thread integer tuner4 = 29;
    static __thread integer tuner5 = 30;
    static __thread integer xctol = 33;
    static __thread integer xftol = 34;

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static __thread doublereal machep;
    extern doublereal rmdcon_(integer *);
    static __thread doublereal mepcrt, sqteps;


/*  ***  supply ***sol (version 2.3) default values to v  *** */

/*  ***  alg = 1 means regression constants. */
/*  ***  alg = 2 means general unconstrained optimization constants. */

/* /+ */
/* / */
/* rmdcon... returns machine-dependent constants */


/*  ***  subscripts for v  *** */


/* /6 */
    /* Parameter adjustments */
    --v;

    /* Function Body */
/* /7 */
/*     parameter (one=1.d+0, three=3.d+0) */
/* / */

/*  ***  v subscript values  *** */

/* /6 */
/* /7 */
/*     parameter (afctol=31, bias=43, cosmin=47, decfac=22, delta0=44, */
/*    1           dfac=41, dinit=38, dltfdc=42, dltfdj=43, dtinit=39, */
/*    2           d0init=40, epslon=19, eta0=42, fuzz=45, huberc=48, */
/*    3           incfac=23, lmax0=35, lmaxs=36, phmnfc=20, phmxfc=21, */
/*    4           rdfcmn=24, rdfcmx=25, rfctol=32, rlimit=46, rsptol=49, */
/*    5           sctol=37, sigmin=50, tuner1=26, tuner2=27, tuner3=28, */
/*    6           tuner4=29, tuner5=30, xctol=33, xftol=34) */
/* / */

/* -------------------------------  body  -------------------------------- */

    machep = rmdcon_(&c__3);
    v[afctol] = 1e-20;
    if (machep > 1e-10) {
/* Computing 2nd power */
	d__1 = machep;
	v[afctol] = d__1 * d__1;
    }
    v[decfac] = .5;
    sqteps = rmdcon_(&c__4);
    v[dfac] = .6;
    v[delta0] = sqteps;
    v[dtinit] = 1e-6;
    d__1 = one / three;
    mepcrt = pow_dd(&machep, &d__1);
    v[d0init] = 1.;
    v[epslon] = .1;
    v[incfac] = 2.;
    v[lmax0] = 1.;
    v[lmaxs] = 1.;
    v[phmnfc] = -.1;
    v[phmxfc] = .1;
    v[rdfcmn] = .1;
    v[rdfcmx] = 4.;
/* Computing MAX */
/* Computing 2nd power */
    d__3 = mepcrt;
    d__1 = 1e-10, d__2 = d__3 * d__3;
    v[rfctol] = max(d__1,d__2);
    v[sctol] = v[rfctol];
    v[tuner1] = .1;
    v[tuner2] = 1e-4;
    v[tuner3] = .75;
    v[tuner4] = .5;
    v[tuner5] = .75;
    v[xctol] = sqteps;
    v[xftol] = machep * 100.;

    if (*alg >= 2) {
	goto L10;
    }

/*  ***  regression  values */

/* Computing MAX */
    d__1 = 1e-6, d__2 = machep * 100.;
    v[cosmin] = max(d__1,d__2);
    v[dinit] = 0.;
    v[dltfdc] = mepcrt;
    v[dltfdj] = sqteps;
    v[fuzz] = 1.5;
    v[huberc] = .7;
    v[rlimit] = rmdcon_(&c__5);
    v[rsptol] = .001;
    v[sigmin] = 1e-4;
    goto L999;

/*  ***  general optimization values */

L10:
    v[bias] = .8;
    v[dinit] = -1.;
    v[eta0] = machep * 1e3;

L999:
    return 0;
/*  ***  last card of vdflt follows  *** */
} /* vdflt_ */

/* Subroutine */ int vscopy_(integer *p, doublereal *y, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static __thread integer i__;


/*  ***  set p-vector y to scalar s  *** */



    /* Parameter adjustments */
    --y;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = *s;
    }
    return 0;
} /* vscopy_ */

doublereal v2norm_(integer *p, doublereal *x)
{
    /* Initialized data */

    static __thread doublereal one = 1.;
    static __thread doublereal zero = 0.;
    static __thread doublereal sqteta = 0.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static __thread integer i__, j;
    static __thread doublereal r__, t, xi, scale;
    extern doublereal rmdcon_(integer *);


/*  ***  return the 2-norm of the p-vector x, taking  *** */
/*  ***  care to avoid the most likely underflows.    *** */


/* /+ */
/* / */

/* /6 */
    /* Parameter adjustments */
    --x;

    /* Function Body */
/* /7 */
/*     parameter (one=1.d+0, zero=0.d+0) */
/*     save sqteta */
/* / */

    if (*p > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L999;
L10:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] != zero) {
	    goto L30;
	}
/* L20: */
    }
    ret_val = zero;
    goto L999;

L30:
    scale = (d__1 = x[i__], abs(d__1));
    if (i__ < *p) {
	goto L40;
    }
    ret_val = scale;
    goto L999;
L40:
    t = one;
    if (sqteta == zero) {
	sqteta = rmdcon_(&c__2);
    }

/*     ***  sqteta is (slightly larger than) the square root of the */
/*     ***  smallest positive floating point number on the machine. */
/*     ***  the tests involving sqteta are done to prevent underflows. */

    j = i__ + 1;
    i__1 = *p;
    for (i__ = j; i__ <= i__1; ++i__) {
	xi = (d__1 = x[i__], abs(d__1));
	if (xi > scale) {
	    goto L50;
	}
	r__ = xi / scale;
	if (r__ > sqteta) {
	    t += r__ * r__;
	}
	goto L60;
L50:
	r__ = scale / xi;
	if (r__ <= sqteta) {
	    r__ = zero;
	}
	t = one + t * r__ * r__;
	scale = xi;
L60:
	;
    }

    ret_val = scale * sqrt(t);
L999:
    return ret_val;
/*  ***  last card of v2norm follows  *** */
} /* v2norm_ */

/* Subroutine */ int humsl_(integer *n, doublereal *d__, doublereal *x, S_fp 
	calcf, S_fp calcgh, integer *iv, integer *liv, integer *lv, 
	doublereal *v, integer *uiparm, doublereal *urparm, U_fp ufparm)
{
    /* Initialized data */

    static __thread integer nextv = 47;
    static __thread integer nfcall = 6;
    static __thread integer nfgcal = 7;
    static __thread integer g = 28;
    static __thread integer h__ = 56;
    static __thread integer toobig = 2;
    static __thread integer vneed = 4;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static __thread doublereal f;
    static __thread integer g1, h1, lh, nf, iv1;
    extern /* Subroutine */ int deflt_(integer *, integer *, integer *, 
	    integer *, doublereal *), humit_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *);


/*  ***  minimize general unconstrained objective function using   *** */
/*  ***  (analytic) gradient and hessian provided by the caller.   *** */

/*     dimension v(78 + n*(n+12)), uiparm(*), urparm(*) */

/* ------------------------------  discussion  --------------------------- */

/*        this routine is like sumsl, except that the subroutine para- */
/*     meter calcg of sumsl (which computes the gradient of the objec- */
/*     tive function) is replaced by the subroutine parameter calcgh, */
/*     which computes both the gradient and (lower triangle of the) */
/*     hessian of the objective function.  the calling sequence is... */
/*             call calcgh(n, x, nf, g, h, uiparm, urparm, ufparm) */
/*     parameters n, x, nf, g, uiparm, urparm, and ufparm are the same */
/*     as for sumsl, while h is an array of length n*(n+1)/2 in which */
/*     calcgh must store the lower triangle of the hessian at x.  start- */
/*     ing at h(1), calcgh must store the hessian entries in the order */
/*     (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), ... */
/*        the value printed (by itsum) in the column labelled stppar */
/*     is the levenberg-marquardt used in computing the current step. */
/*     zero means a full newton step.  if the special case described in */
/*     ref. 1 is detected, then stppar is negated.  the value printed */
/*     in the column labelled npreldf is zero if the current hessian */
/*     is not positive definite. */
/*        it sometimes proves worthwhile to let d be determined from the */
/*     diagonal of the hessian matrix by setting iv(dtype) = 1 and */
/*     v(dinit) = 0.  the following iv and v components are relevant... */

/* iv(dtol)..... iv(59) gives the starting subscript in v of the dtol */
/*             array used when d is updated.  (iv(dtol) can be */
/*             initialized by calling humsl with iv(1) = 13.) */
/* iv(dtype).... iv(16) tells how the scale vector d should be chosen. */
/*             iv(dtype) .le. 0 means that d should not be updated, and */
/*             iv(dtype) .ge. 1 means that d should be updated as */
/*             described below with v(dfac).  default = 0. */
/* v(dfac)..... v(41) and the dtol and d0 arrays (see v(dtinit) and */
/*             v(d0init)) are used in updating the scale vector d when */
/*             iv(dtype) .gt. 0.  (d is initialized according to */
/*             v(dinit), described in sumsl.)  let */
/*                  d1(i) = max(sqrt(abs(h(i,i))), v(dfac)*d(i)), */
/*             where h(i,i) is the i-th diagonal element of the current */
/*             hessian.  if iv(dtype) = 1, then d(i) is set to d1(i) */
/*             unless d1(i) .lt. dtol(i), in which case d(i) is set to */
/*                  max(d0(i), dtol(i)). */
/*             if iv(dtype) .ge. 2, then d is updated during the first */
/*             iteration as for iv(dtype) = 1 (after any initialization */
/*             due to v(dinit)) and is left unchanged thereafter. */
/*             default = 0.6. */
/* v(dtinit)... v(39), if positive, is the value to which all components */
/*             of the dtol array (see v(dfac)) are initialized.  if */
/*             v(dtinit) = 0, then it is assumed that the caller has */
/*             stored dtol in v starting at v(iv(dtol)). */
/*             default = 10**-6. */
/* v(d0init)... v(40), if positive, is the value to which all components */
/*             of the d0 vector (see v(dfac)) are initialized.  if */
/*             v(dfac) = 0, then it is assumed that the caller has */
/*             stored d0 in v starting at v(iv(dtol)+n).  default = 1.0. */

/*  ***  reference  *** */

/* 1. gay, d.m. (1981), computing optimal locally constrained steps, */
/*         siam j. sci. statist. comput. 2, pp. 186-197. */
/* . */
/*  ***  general  *** */

/*     coded by david m. gay (winter 1980).  revised sept. 1982. */
/*     this subroutine was written in connection with research supported */
/*     in part by the national science foundation under grants */
/*     mcs-7600324 and mcs-7906671. */

/* ----------------------------  declarations  --------------------------- */


/* deflt... provides default input values for iv and v. */
/* humit... reverse-communication routine that does humsl algorithm. */


/*  ***  subscripts for iv   *** */


/* /6 */
    /* Parameter adjustments */
    --x;
    --d__;
    --iv;
    --v;
    --uiparm;
    --urparm;

    /* Function Body */
/* /7 */
/*     parameter (nextv=47, nfcall=6, nfgcal=7, g=28, h=56, toobig=2, */
/*    1           vneed=4) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

    lh = *n * (*n + 1) / 2;
    if (iv[1] == 0) {
	deflt_(&c__2, &iv[1], liv, lv, &v[1]);
    }
    if (iv[1] == 12 || iv[1] == 13) {
	iv[vneed] += *n * (*n + 3) / 2;
    }
    iv1 = iv[1];
    if (iv1 == 14) {
	goto L10;
    }
    if (iv1 > 2 && iv1 < 12) {
	goto L10;
    }
    g1 = 1;
    h1 = 1;
    if (iv1 == 12) {
	iv[1] = 13;
    }
    goto L20;

L10:
    g1 = iv[g];
    h1 = iv[h__];

L20:
    humit_(&d__[1], &f, &v[g1], &v[h1], &iv[1], &lh, liv, lv, n, &v[1], &x[1])
	    ;
    if ((i__1 = iv[1] - 2) < 0) {
	goto L30;
    } else if (i__1 == 0) {
	goto L40;
    } else {
	goto L50;
    }

L30:
    nf = iv[nfcall];
    (*calcf)(n, &x[1], &nf, &f, &uiparm[1], &urparm[1], (U_fp)ufparm);
    if (nf <= 0) {
	iv[toobig] = 1;
    }
    goto L20;

L40:
    (*calcgh)(n, &x[1], &iv[nfgcal], &v[g1], &v[h1], &uiparm[1], &urparm[1], (
	    U_fp)ufparm);
    goto L20;

L50:
    if (iv[1] != 14) {
	goto L999;
    }

/*  ***  storage allocation */

    iv[g] = iv[nextv];
    iv[h__] = iv[g] + *n;
    iv[nextv] = iv[h__] + *n * (*n + 1) / 2;
    if (iv1 != 13) {
	goto L10;
    }

L999:
    return 0;
/*  ***  last card of humsl follows  *** */
} /* humsl_ */

/* Subroutine */ int humit_(doublereal *d__, doublereal *fx, doublereal *g, 
	doublereal *h__, integer *iv, integer *lh, integer *liv, integer *lv, 
	integer *n, doublereal *v, doublereal *x)
{
    /* Initialized data */

    static __thread integer cnvcod = 55;
    static __thread integer dg = 37;
    static __thread integer dtol = 59;
    static __thread integer dtype = 16;
    static __thread integer irc = 29;
    static __thread integer kagqt = 33;
    static __thread integer lmat = 42;
    static __thread integer mode = 35;
    static __thread integer model = 5;
    static __thread integer mxfcal = 17;
    static __thread integer mxiter = 18;
    static __thread integer nextv = 47;
    static __thread integer nfcall = 6;
    static __thread integer nfgcal = 7;
    static __thread integer ngcall = 30;
    static __thread integer niter = 31;
    static __thread integer radinc = 8;
    static __thread integer restor = 9;
    static __thread integer step = 40;
    static __thread integer stglim = 11;
    static __thread integer stlstg = 41;
    static __thread integer toobig = 2;
    static __thread integer vneed = 4;
    static __thread integer w = 34;
    static __thread integer xirc = 13;
    static __thread integer x0 = 43;
    static __thread integer dgnorm = 1;
    static __thread integer dinit = 38;
    static __thread integer dstnrm = 2;
    static __thread integer dtinit = 39;
    static __thread integer d0init = 40;
    static __thread integer f = 10;
    static __thread integer f0 = 13;
    static __thread integer fdif = 11;
    static __thread integer gtstep = 4;
    static __thread integer incfac = 23;
    static __thread integer lmax0 = 35;
    static __thread integer lmaxs = 36;
    static __thread integer preduc = 7;
    static __thread integer radfac = 16;
    static __thread integer radius = 8;
    static __thread integer rad0 = 9;
    static __thread integer reldx = 17;
    static __thread integer stppar = 5;
    static __thread integer tuner4 = 29;
    static __thread integer tuner5 = 30;
    static __thread doublereal one = 1.;
    static __thread doublereal onep2 = 1.2;
    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static __thread integer i__, j, k, l;
    static __thread doublereal t;
    static __thread integer w1, x01, dg1, nn1o2, temp1, step1;
    extern /* Subroutine */ int deflt_(integer *, integer *, integer *, 
	    integer *, doublereal *), parck_(integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *), dupdu_(
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *);
    static __thread integer dummy;
    extern /* Subroutine */ int assst_(integer *, integer *, integer *, 
	    doublereal *), vcopy_(integer *, doublereal *, doublereal *), 
	    itsum_(doublereal *, doublereal *, integer *, integer *, integer *
	    , integer *, doublereal *, doublereal *), gqtst_(doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *), vaxpy_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern logical stopx_(integer *);
    extern doublereal v2norm_(integer *, doublereal *), dotprd_(integer *, 
	    doublereal *, doublereal *), reldst_(integer *, doublereal *, 
	    doublereal *, doublereal *);
    static __thread integer lstgst;
    extern /* Subroutine */ int slvmul_(integer *, doublereal *, doublereal *,
	     doublereal *), vscopy_(integer *, doublereal *, doublereal *);


/*  ***  carry out humsl (unconstrained minimization) iterations, using */
/*  ***  hessian matrix provided by the caller. */

/*  ***  parameter declarations  *** */


/* --------------------------  parameter usage  -------------------------- */

/* d.... scale vector. */
/* fx... function value. */
/* g.... gradient vector. */
/* h.... lower triangle of the hessian, stored rowwise. */
/* iv... integer value array. */
/* lh... length of h = p*(p+1)/2. */
/* liv.. length of iv (at least 60). */
/* lv... length of v (at least 78 + n*(n+21)/2). */
/* n.... number of variables (components in x and g). */
/* v.... floating-point value array. */
/* x.... parameter vector. */

/*  ***  discussion  *** */

/*        parameters iv, n, v, and x are the same as the corresponding */
/*     ones to humsl (which see), except that v can be shorter (since */
/*     the part of v that humsl uses for storing g and h is not needed). */
/*     moreover, compared with humsl, iv(1) may have the two additional */
/*     output values 1 and 2, which are explained below, as is the use */
/*     of iv(toobig) and iv(nfgcal).  the value iv(g), which is an */
/*     output value from humsl, is not referenced by humit or the */
/*     subroutines it calls. */

/* iv(1) = 1 means the caller should set fx to f(x), the function value */
/*             at x, and call humit again, having changed none of the */
/*             other parameters.  an exception occurs if f(x) cannot be */
/*             computed (e.g. if overflow would occur), which may happen */
/*             because of an oversized step.  in this case the caller */
/*             should set iv(toobig) = iv(2) to 1, which will cause */
/*             humit to ignore fx and try a smaller step.  the para- */
/*             meter nf that humsl passes to calcf (for possible use by */
/*             calcgh) is a copy of iv(nfcall) = iv(6). */
/* iv(1) = 2 means the caller should set g to g(x), the gradient of f at */
/*             x, and h to the lower triangle of h(x), the hessian of f */
/*             at x, and call humit again, having changed none of the */
/*             other parameters except perhaps the scale vector d. */
/*                  the parameter nf that humsl passes to calcg is */
/*             iv(nfgcal) = iv(7).  if g(x) and h(x) cannot be evaluated, */
/*             then the caller may set iv(nfgcal) to 0, in which case */
/*             humit will return with iv(1) = 65. */
/*                  note -- humit overwrites h with the lower triangle */
/*             of  diag(d)**-1 * h(x) * diag(d)**-1. */
/* . */
/*  ***  general  *** */

/*     coded by david m. gay (winter 1980).  revised sept. 1982. */
/*     this subroutine was written in connection with research supported */
/*     in part by the national science foundation under grants */
/*     mcs-7600324 and mcs-7906671. */

/*        (see sumsl and humsl for references.) */

/* +++++++++++++++++++++++++++  declarations  ++++++++++++++++++++++++++++ */

/*  ***  local variables  *** */


/*     ***  constants  *** */


/*  ***  no intrinsic functions  *** */

/*  ***  external functions and subroutines  *** */


/* assst.... assesses candidate step. */
/* deflt.... provides default iv and v input values. */
/* dotprd... returns inner product of two vectors. */
/* dupdu.... updates scale vector d. */
/* gqtst.... computes optimally locally constrained step. */
/* itsum.... prints iteration summary and info on initial and final x. */
/* parck.... checks validity of input iv and v values. */
/* reldst... computes v(reldx) = relative step size. */
/* slvmul... multiplies symmetric matrix times vector, given the lower */
/*             triangle of the matrix. */
/* stopx.... returns .true. if the break key has been pressed. */
/* vaxpy.... computes scalar times one vector plus another. */
/* vcopy.... copies one vector to another. */
/* vscopy... sets all elements of a vector to a scalar. */
/* v2norm... returns the 2-norm of a vector. */

/*  ***  subscripts for iv and v  *** */


/*  ***  iv subscript values  *** */

/* /6 */
    /* Parameter adjustments */
    --h__;
    --iv;
    --v;
    --x;
    --g;
    --d__;

    /* Function Body */
/* /7 */
/*     parameter (cnvcod=55, dg=37, dtol=59, dtype=16, irc=29, kagqt=33, */
/*    1           lmat=42, mode=35, model=5, mxfcal=17, mxiter=18, */
/*    2           nextv=47, nfcall=6, nfgcal=7, ngcall=30, niter=31, */
/*    3           radinc=8, restor=9, step=40, stglim=11, stlstg=41, */
/*    4           toobig=2, vneed=4, w=34, xirc=13, x0=43) */
/* / */

/*  ***  v subscript values  *** */

/* /6 */
/* /7 */
/*     parameter (dgnorm=1, dinit=38, dstnrm=2, dtinit=39, d0init=40, */
/*    1           f=10, f0=13, fdif=11, gtstep=4, incfac=23, lmax0=35, */
/*    2           lmaxs=36, preduc=7, radfac=16, radius=8, rad0=9, */
/*    3           reldx=17, stppar=5, tuner4=29, tuner5=30) */
/* / */

/* /6 */
/* /7 */
/*     parameter (one=1.d+0, onep2=1.2d+0, zero=0.d+0) */
/* / */

/* +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++ */

    i__ = iv[1];
    if (i__ == 1) {
	goto L30;
    }
    if (i__ == 2) {
	goto L40;
    }

/*  ***  check validity of iv and v input values  *** */

    if (iv[1] == 0) {
	deflt_(&c__2, &iv[1], liv, lv, &v[1]);
    }
    if (iv[1] == 12 || iv[1] == 13) {
	iv[vneed] = iv[vneed] + *n * (*n + 21) / 2 + 7;
    }
    parck_(&c__2, &d__[1], &iv[1], liv, lv, n, &v[1]);
    i__ = iv[1] - 2;
    if (i__ > 12) {
	goto L999;
    }
    nn1o2 = *n * (*n + 1) / 2;
    if (*lh >= nn1o2) {
	switch (i__) {
	    case 1:  goto L210;
	    case 2:  goto L210;
	    case 3:  goto L210;
	    case 4:  goto L210;
	    case 5:  goto L210;
	    case 6:  goto L210;
	    case 7:  goto L160;
	    case 8:  goto L120;
	    case 9:  goto L160;
	    case 10:  goto L10;
	    case 11:  goto L10;
	    case 12:  goto L20;
	}
    }
    iv[1] = 66;
    goto L350;

/*  ***  storage allocation  *** */

L10:
    iv[dtol] = iv[lmat] + nn1o2;
    iv[x0] = iv[dtol] + (*n << 1);
    iv[step] = iv[x0] + *n;
    iv[stlstg] = iv[step] + *n;
    iv[dg] = iv[stlstg] + *n;
    iv[w] = iv[dg] + *n;
    iv[nextv] = iv[w] + (*n << 2) + 7;
    if (iv[1] != 13) {
	goto L20;
    }
    iv[1] = 14;
    goto L999;

/*  ***  initialization  *** */

L20:
    iv[niter] = 0;
    iv[nfcall] = 1;
    iv[ngcall] = 1;
    iv[nfgcal] = 1;
    iv[mode] = -1;
    iv[model] = 1;
    iv[stglim] = 1;
    iv[toobig] = 0;
    iv[cnvcod] = 0;
    iv[radinc] = 0;
    v[rad0] = zero;
    v[stppar] = zero;
    if (v[dinit] >= zero) {
	vscopy_(n, &d__[1], &v[dinit]);
    }
    k = iv[dtol];
    if (v[dtinit] > zero) {
	vscopy_(n, &v[k], &v[dtinit]);
    }
    k += *n;
    if (v[d0init] > zero) {
	vscopy_(n, &v[k], &v[d0init]);
    }
    iv[1] = 1;
    goto L999;

L30:
    v[f] = *fx;
    if (iv[mode] >= 0) {
	goto L210;
    }
    iv[1] = 2;
    if (iv[toobig] == 0) {
	goto L999;
    }
    iv[1] = 63;
    goto L350;

/*  ***  make sure gradient could be computed  *** */

L40:
    if (iv[nfgcal] != 0) {
	goto L50;
    }
    iv[1] = 65;
    goto L350;

/*  ***  update the scale vector d  *** */

L50:
    dg1 = iv[dg];
    if (iv[dtype] <= 0) {
	goto L70;
    }
    k = dg1;
    j = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j += i__;
	v[k] = h__[j];
	++k;
/* L60: */
    }
    dupdu_(&d__[1], &v[dg1], &iv[1], liv, lv, n, &v[1]);

/*  ***  compute scaled gradient and its norm  *** */

L70:
    dg1 = iv[dg];
    k = dg1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = g[i__] / d__[i__];
	++k;
/* L80: */
    }
    v[dgnorm] = v2norm_(n, &v[dg1]);

/*  ***  compute scaled hessian  *** */

    k = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = one / d__[i__];
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    h__[k] = t * h__[k] / d__[j];
	    ++k;
/* L90: */
	}
/* L100: */
    }

    if (iv[cnvcod] != 0) {
	goto L340;
    }
    if (iv[mode] == 0) {
	goto L300;
    }

/*  ***  allow first step to have scaled 2-norm at most v(lmax0)  *** */

    v[radius] = v[lmax0];

    iv[mode] = 0;


/* -----------------------------  main loop  ----------------------------- */


/*  ***  print iteration summary, check iteration limit  *** */

L110:
    itsum_(&d__[1], &g[1], &iv[1], liv, lv, n, &v[1], &x[1]);
L120:
    k = iv[niter];
    if (k < iv[mxiter]) {
	goto L130;
    }
    iv[1] = 10;
    goto L350;

L130:
    iv[niter] = k + 1;

/*  ***  initialize for start of next iteration  *** */

    dg1 = iv[dg];
    x01 = iv[x0];
    v[f0] = v[f];
    iv[irc] = 4;
    iv[kagqt] = -1;

/*     ***  copy x to x0  *** */

    vcopy_(n, &v[x01], &x[1]);

/*  ***  update radius  *** */

    if (k == 0) {
	goto L150;
    }
    step1 = iv[step];
    k = step1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = d__[i__] * v[k];
	++k;
/* L140: */
    }
    v[radius] = v[radfac] * v2norm_(n, &v[step1]);

/*  ***  check stopx and function evaluation limit  *** */

L150:
    if (! stopx_(&dummy)) {
	goto L170;
    }
    iv[1] = 11;
    goto L180;

/*     ***  come here when restarting after func. eval. limit or stopx. */

L160:
    if (v[f] >= v[f0]) {
	goto L170;
    }
    v[radfac] = one;
    k = iv[niter];
    goto L130;

L170:
    if (iv[nfcall] < iv[mxfcal]) {
	goto L190;
    }
    iv[1] = 9;
L180:
    if (v[f] >= v[f0]) {
	goto L350;
    }

/*        ***  in case of stopx or function evaluation limit with */
/*        ***  improved v(f), evaluate the gradient at x. */

    iv[cnvcod] = iv[1];
    goto L290;

/* . . . . . . . . . . . . .  compute candidate step  . . . . . . . . . . */

L190:
    step1 = iv[step];
    dg1 = iv[dg];
    l = iv[lmat];
    w1 = iv[w];
    gqtst_(&d__[1], &v[dg1], &h__[1], &iv[kagqt], &v[l], n, &v[step1], &v[1], 
	    &v[w1]);
    if (iv[irc] == 6) {
	goto L210;
    }

/*  ***  check whether evaluating f(x0 + step) looks worthwhile  *** */

    if (v[dstnrm] <= zero) {
	goto L210;
    }
    if (iv[irc] != 5) {
	goto L200;
    }
    if (v[radfac] <= one) {
	goto L200;
    }
    if (v[preduc] <= onep2 * v[fdif]) {
	goto L210;
    }

/*  ***  compute f(x0 + step)  *** */

L200:
    x01 = iv[x0];
    step1 = iv[step];
    vaxpy_(n, &x[1], &one, &v[step1], &v[x01]);
    ++iv[nfcall];
    iv[1] = 1;
    iv[toobig] = 0;
    goto L999;

/* . . . . . . . . . . . . .  assess candidate step  . . . . . . . . . . . */

L210:
    x01 = iv[x0];
    v[reldx] = reldst_(n, &d__[1], &x[1], &v[x01]);
    assst_(&iv[1], liv, lv, &v[1]);
    step1 = iv[step];
    lstgst = iv[stlstg];
    if (iv[restor] == 1) {
	vcopy_(n, &x[1], &v[x01]);
    }
    if (iv[restor] == 2) {
	vcopy_(n, &v[lstgst], &v[step1]);
    }
    if (iv[restor] != 3) {
	goto L220;
    }
    vcopy_(n, &v[step1], &v[lstgst]);
    vaxpy_(n, &x[1], &one, &v[step1], &v[x01]);
    v[reldx] = reldst_(n, &d__[1], &x[1], &v[x01]);

L220:
    k = iv[irc];
    switch (k) {
	case 1:  goto L230;
	case 2:  goto L260;
	case 3:  goto L260;
	case 4:  goto L260;
	case 5:  goto L230;
	case 6:  goto L240;
	case 7:  goto L250;
	case 8:  goto L250;
	case 9:  goto L250;
	case 10:  goto L250;
	case 11:  goto L250;
	case 12:  goto L250;
	case 13:  goto L330;
	case 14:  goto L300;
    }

/*     ***  recompute step with new radius  *** */

L230:
    v[radius] = v[radfac] * v[dstnrm];
    goto L150;

/*  ***  compute step of length v(lmaxs) for singular convergence test. */

L240:
    v[radius] = v[lmaxs];
    goto L190;

/*  ***  convergence or false convergence  *** */

L250:
    iv[cnvcod] = k - 4;
    if (v[f] >= v[f0]) {
	goto L340;
    }
    if (iv[xirc] == 14) {
	goto L340;
    }
    iv[xirc] = 14;

/* . . . . . . . . . . . .  process acceptable step  . . . . . . . . . . . */

L260:
    if (iv[irc] != 3) {
	goto L290;
    }
    temp1 = lstgst;

/*     ***  prepare for gradient tests  *** */
/*     ***  set  temp1 = hessian * step + g(x0) */
/*     ***             = diag(d) * (h * step + g(x0)) */

/*        use x0 vector as temporary. */
    k = x01;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = d__[i__] * v[step1];
	++k;
	++step1;
/* L270: */
    }
    slvmul_(n, &v[temp1], &h__[1], &v[x01]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[temp1] = d__[i__] * v[temp1] + g[i__];
	++temp1;
/* L280: */
    }

/*  ***  compute gradient and hessian  *** */

L290:
    ++iv[ngcall];
    iv[1] = 2;
    goto L999;

L300:
    iv[1] = 2;
    if (iv[irc] != 3) {
	goto L110;
    }

/*  ***  set v(radfac) by gradient tests  *** */

    temp1 = iv[stlstg];
    step1 = iv[step];

/*     ***  set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))  *** */

    k = temp1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = (v[k] - g[i__]) / d__[i__];
	++k;
/* L310: */
    }

/*     ***  do gradient tests  *** */

    if (v2norm_(n, &v[temp1]) <= v[dgnorm] * v[tuner4]) {
	goto L320;
    }
    if (dotprd_(n, &g[1], &v[step1]) >= v[gtstep] * v[tuner5]) {
	goto L110;
    }
L320:
    v[radfac] = v[incfac];
    goto L110;

/* . . . . . . . . . . . . . .  misc. details  . . . . . . . . . . . . . . */

/*  ***  bad parameters to assess  *** */

L330:
    iv[1] = 64;
    goto L350;

/*  ***  print summary of final iteration and other requested items  *** */

L340:
    iv[1] = iv[cnvcod];
    iv[cnvcod] = 0;
L350:
    itsum_(&d__[1], &g[1], &iv[1], liv, lv, n, &v[1], &x[1]);

L999:
    return 0;

/*  ***  last card of humit follows  *** */
} /* humit_ */

/* Subroutine */ int dupdu_(doublereal *d__, doublereal *hdiag, integer *iv, 
	integer *liv, integer *lv, integer *n, doublereal *v)
{
    /* Initialized data */

    static __thread integer dfac = 41;
    static __thread integer dtol = 59;
    static __thread integer dtype = 16;
    static __thread integer niter = 31;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static __thread integer i__;
    static __thread doublereal t;
    static __thread integer d0i;
    static __thread doublereal vdfac;
    static __thread integer dtoli;


/*  ***  update scale vector d for humsl  *** */

/*  ***  parameter declarations  *** */


/*  ***  local variables  *** */


/*  ***  intrinsic functions  *** */
/* /+ */
/* / */
/*  ***  subscripts for iv and v  *** */

/* /6 */
    /* Parameter adjustments */
    --iv;
    --v;
    --hdiag;
    --d__;

    /* Function Body */
/* /7 */
/*     parameter (dfac=41, dtol=59, dtype=16, niter=31) */
/* / */

/* -------------------------------  body  -------------------------------- */

    i__ = iv[(0 + (0 + (dtype << 2))) / 4];
    if (i__ == 1) {
	goto L10;
    }
    if (iv[niter] > 0) {
	goto L999;
    }

L10:
    dtoli = iv[dtol];
    d0i = dtoli + *n;
    vdfac = v[dfac];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = sqrt((d__1 = hdiag[i__], abs(d__1))), d__3 = vdfac * d__[i__];
	t = max(d__2,d__3);
	if (t < v[dtoli]) {
/* Computing MAX */
	    d__1 = v[dtoli], d__2 = v[d0i];
	    t = max(d__1,d__2);
	}
	d__[i__] = t;
	++dtoli;
	++d0i;
/* L20: */
    }

L999:
    return 0;
/*  ***  last card of dupdu follows  *** */
} /* dupdu_ */

/* Subroutine */ int gqtst_(doublereal *d__, doublereal *dig, doublereal *
	dihdi, integer *ka, doublereal *l, integer *p, doublereal *step, 
	doublereal *v, doublereal *w)
{
    /* Initialized data */

    static __thread integer dgnorm = 1;
    static __thread integer dstnrm = 2;
    static __thread integer dst0 = 3;
    static __thread integer epslon = 19;
    static __thread integer gtstep = 4;
    static __thread integer nreduc = 6;
    static __thread integer phmnfc = 20;
    static __thread integer phmxfc = 21;
    static __thread integer preduc = 7;
    static __thread integer radius = 8;
    static __thread integer rad0 = 9;
    static __thread integer stppar = 5;
    static __thread doublereal epsfac = 50.;
    static __thread doublereal four = 4.;
    static __thread doublereal half = .5;
    static __thread doublereal kappa = 2.;
    static __thread doublereal negone = -1.;
    static __thread doublereal one = 1.;
    static __thread doublereal p001 = .001;
    static __thread doublereal six = 6.;
    static __thread doublereal three = 3.;
    static __thread doublereal two = 2.;
    static __thread doublereal zero = 0.;
    static __thread doublereal big = 0.;
    static __thread doublereal dgxfac = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static __thread integer i__, j, k, q;
    static __thread doublereal t;
    static __thread integer x, k1, q0;
    static __thread doublereal t1, t2, lk, si, sk, uk, wi, sw;
    static __thread integer im1, lk0, uk0;
    static __thread doublereal aki, akk, rad;
    static __thread integer inc, irc;
    static __thread doublereal phi, eps, dst;
    static __thread integer diag, emin, emax;
    static __thread doublereal root;
    static __thread integer diag0;
    static __thread doublereal delta;
    static __thread integer kalim, kamin;
    static __thread doublereal radsq, gtsta;
    extern /* Subroutine */ int lsqrt_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    extern doublereal v2norm_(integer *, doublereal *);
    static __thread doublereal alphak, psifac;
    static __thread integer dggdmx;
    static __thread doublereal oldphi;
    extern doublereal rmdcon_(integer *);
    static __thread doublereal phimin, phimax;
    static __thread integer phipin;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);
    static __thread integer dstsav;
    extern /* Subroutine */ int livmul_(integer *, doublereal *, doublereal *,
	     doublereal *);
    extern doublereal lsvmin_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int litvmu_(integer *, doublereal *, doublereal *,
	     doublereal *);
    static __thread logical restrt;
    static __thread doublereal twopsi;


/*  *** compute goldfeld-quandt-trotter step by more-hebden technique *** */
/*  ***  (nl2sol version 2.2), modified a la more and sorensen  *** */

/*  ***  parameter declarations  *** */

/*     dimension dihdi(p*(p+1)/2), l(p*(p+1)/2), w(4*p+7) */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  purpose  *** */

/*        given the (compactly stored) lower triangle of a scaled */
/*     hessian (approximation) and a nonzero scaled gradient vector, */
/*     this subroutine computes a goldfeld-quandt-trotter step of */
/*     approximate length v(radius) by the more-hebden technique.  in */
/*     other words, step is computed to (approximately) minimize */
/*     psi(step) = (g**t)*step + 0.5*(step**t)*h*step  such that the */
/*     2-norm of d*step is at most (approximately) v(radius), where */
/*     g  is the gradient,  h  is the hessian, and  d  is a diagonal */
/*     scale matrix whose diagonal is stored in the parameter d. */
/*     (gqtst assumes  dig = d**-1 * g  and  dihdi = d**-1 * h * d**-1.) */

/*  ***  parameter description  *** */

/*     d (in)  = the scale vector, i.e. the diagonal of the scale */
/*              matrix  d  mentioned above under purpose. */
/*   dig (in)  = the scaled gradient vector, d**-1 * g.  if g = 0, then */
/*              step = 0  and  v(stppar) = 0  are returned. */
/* dihdi (in)  = lower triangle of the scaled hessian (approximation), */
/*              i.e., d**-1 * h * d**-1, stored compactly by rows., i.e., */
/*              in the order (1,1), (2,1), (2,2), (3,1), (3,2), etc. */
/*    ka (i/o) = the number of hebden iterations (so far) taken to deter- */
/*              mine step.  ka .lt. 0 on input means this is the first */
/*              attempt to determine step (for the present dig and dihdi) */
/*              -- ka is initialized to 0 in this case.  output with */
/*              ka = 0  (or v(stppar) = 0)  means  step = -(h**-1)*g. */
/*     l (i/o) = workspace of length p*(p+1)/2 for cholesky factors. */
/*     p (in)  = number of parameters -- the hessian is a  p x p  matrix. */
/*  step (i/o) = the step computed. */
/*     v (i/o) contains various constants and variables described below. */
/*     w (i/o) = workspace of length 4*p + 6. */

/*  ***  entries in v  *** */

/* v(dgnorm) (i/o) = 2-norm of (d**-1)*g. */
/* v(dstnrm) (output) = 2-norm of d*step. */
/* v(dst0)   (i/o) = 2-norm of d*(h**-1)*g (for pos. def. h only), or */
/*             overestimate of smallest eigenvalue of (d**-1)*h*(d**-1). */
/* v(epslon) (in)  = max. rel. error allowed for psi(step).  for the */
/*             step returned, psi(step) will exceed its optimal value */
/*             by less than -v(epslon)*psi(step).  suggested value = 0.1. */
/* v(gtstep) (out) = inner product between g and step. */
/* v(nreduc) (out) = psi(-(h**-1)*g) = psi(newton step)  (for pos. def. */
/*             h only -- v(nreduc) is set to zero otherwise). */
/* v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step */
/*             (more*s sigma).  the error v(dstnrm) - v(radius) must lie */
/*             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius). */
/* v(phmxfc) (in)  (see v(phmnfc).) */
/*             suggested values -- v(phmnfc) = -0.25, v(phmxfc) = 0.5. */
/* v(preduc) (out) = psi(step) = predicted obj. func. reduction for step. */
/* v(radius) (in)  = radius of current (scaled) trust region. */
/* v(rad0)   (i/o) = value of v(radius) from previous call. */
/* v(stppar) (i/o) is normally the marquardt parameter, i.e. the alpha */
/*             described below under algorithm notes.  if h + alpha*d**2 */
/*             (see algorithm notes) is (nearly) singular, however, */
/*             then v(stppar) = -alpha. */

/*  ***  usage notes  *** */

/*     if it is desired to recompute step using a different value of */
/*     v(radius), then this routine may be restarted by calling it */
/*     with all parameters unchanged except v(radius).  (this explains */
/*     why step and w are listed as i/o).  on an initial call (one with */
/*     ka .lt. 0), step and w need not be initialized and only compo- */
/*     nents v(epslon), v(stppar), v(phmnfc), v(phmxfc), v(radius), and */
/*     v(rad0) of v must be initialized. */

/*  ***  algorithm notes  *** */

/*        the desired g-q-t step (ref. 2, 3, 4, 6) satisfies */
/*     (h + alpha*d**2)*step = -g  for some nonnegative alpha such that */
/*     h + alpha*d**2 is positive semidefinite.  alpha and step are */
/*     computed by a scheme analogous to the one described in ref. 5. */
/*     estimates of the smallest and largest eigenvalues of the hessian */
/*     are obtained from the gerschgorin circle theorem enhanced by a */
/*     simple form of the scaling described in ref. 7.  cases in which */
/*     h + alpha*d**2 is nearly (or exactly) singular are handled by */
/*     the technique discussed in ref. 2.  in these cases, a step of */
/*     (exact) length v(radius) is returned for which psi(step) exceeds */
/*     its optimal value by less than -v(epslon)*psi(step).  the test */
/*     suggested in ref. 6 for detecting the special case is performed */
/*     once two matrix factorizations have been done -- doing so sooner */
/*     seems to degrade the performance of optimization routines that */
/*     call this routine. */

/*  ***  functions and subroutines called  *** */

/* dotprd - returns inner product of two vectors. */
/* litvmu - applies inverse-transpose of compact lower triang. matrix. */
/* livmul - applies inverse of compact lower triang. matrix. */
/* lsqrt  - finds cholesky factor (of compactly stored lower triang.). */
/* lsvmin - returns approx. to min. sing. value of lower triang. matrix. */
/* rmdcon - returns machine-dependent constants. */
/* v2norm - returns 2-norm of a vector. */

/*  ***  references  *** */

/* 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive */
/*             nonlinear least-squares algorithm, acm trans. math. */
/*             software, vol. 7, no. 3. */
/* 2.  gay, d.m. (1981), computing optimal locally constrained steps, */
/*             siam j. sci. statist. computing, vol. 2, no. 2, pp. */
/*             186-197. */
/* 3.  goldfeld, s.m., quandt, r.e., and trotter, h.f. (1966), */
/*             maximization by quadratic hill-climbing, econometrica 34, */
/*             pp. 541-551. */
/* 4.  hebden, m.d. (1973), an algorithm for minimization using exact */
/*             second derivatives, report t.p. 515, theoretical physics */
/*             div., a.e.r.e. harwell, oxon., england. */
/* 5.  more, j.j. (1978), the levenberg-marquardt algorithm, implemen- */
/*             tation and theory, pp.105-116 of springer lecture notes */
/*             in mathematics no. 630, edited by g.a. watson, springer- */
/*             verlag, berlin and new york. */
/* 6.  more, j.j., and sorensen, d.c. (1981), computing a trust region */
/*             step, technical report anl-81-83, argonne national lab. */
/* 7.  varga, r.s. (1965), minimal gerschgorin sets, pacific j. math. 15, */
/*             pp. 719-729. */

/*  ***  general  *** */

/*     coded by david m. gay. */
/*     this subroutine was written in connection with research */
/*     supported by the national science foundation under grants */
/*     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and */
/*     mcs-7906671. */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  local variables  *** */


/*     ***  constants  *** */

/*  ***  intrinsic functions  *** */
/* /+ */
/* / */
/*  ***  external functions and subroutines  *** */


/*  ***  subscripts for v  *** */

/* /6 */
    /* Parameter adjustments */
    --dihdi;
    --l;
    --step;
    --dig;
    --d__;
    --v;
    --w;

    /* Function Body */
/* /7 */
/*     parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19, gtstep=4, */
/*    1           nreduc=6, phmnfc=20, phmxfc=21, preduc=7, radius=8, */
/*    2           rad0=9, stppar=5) */
/* / */

/* /6 */
/* /7 */
/*     parameter (epsfac=50.0d+0, four=4.0d+0, half=0.5d+0, */
/*    1     kappa=2.0d+0, negone=-1.0d+0, one=1.0d+0, p001=1.0d-3, */
/*    2     six=6.0d+0, three=3.0d+0, two=2.0d+0, zero=0.0d+0) */
/*     save dgxfac */
/* / */

/*  ***  body  *** */

/*     ***  store largest abs. entry in (d**-1)*h*(d**-1) at w(dggdmx). */
    dggdmx = *p + 1;
/*     ***  store gerschgorin over- and underestimates of the largest */
/*     ***  and smallest eigenvalues of (d**-1)*h*(d**-1) at w(emax) */
/*     ***  and w(emin) respectively. */
    emax = dggdmx + 1;
    emin = emax + 1;
/*     ***  for use in recomputing step, the final values of lk, uk, dst, */
/*     ***  and the inverse derivative of more*s phi at 0 (for pos. def. */
/*     ***  h) are stored in w(lk0), w(uk0), w(dstsav), and w(phipin) */
/*     ***  respectively. */
    lk0 = emin + 1;
    phipin = lk0 + 1;
    uk0 = phipin + 1;
    dstsav = uk0 + 1;
/*     ***  store diag of (d**-1)*h*(d**-1) in w(diag),...,w(diag0+p). */
    diag0 = dstsav;
    diag = diag0 + 1;
/*     ***  store -d*step in w(q),...,w(q0+p). */
    q0 = diag0 + *p;
    q = q0 + 1;
/*     ***  allocate storage for scratch vector x  *** */
    x = q + *p;
    rad = v[radius];
/* Computing 2nd power */
    d__1 = rad;
    radsq = d__1 * d__1;
/*     ***  phitol = max. error allowed in dst = v(dstnrm) = 2-norm of */
/*     ***  d*step. */
    phimax = v[phmxfc] * rad;
    phimin = v[phmnfc] * rad;
/* Computing 2nd power */
    d__1 = rad;
    psifac = two * v[epslon] / (three * (four * (v[phmnfc] + one) * (kappa + 
	    one) + kappa + two) * (d__1 * d__1));
/*     ***  oldphi is used to detect limits of numerical accuracy.  if */
/*     ***  we recompute step and it does not change, then we accept it. */
    oldphi = zero;
    eps = v[epslon];
    irc = 0;
    restrt = FALSE_;
    kalim = *ka + 50;

/*  ***  start or restart, depending on ka  *** */

    if (*ka >= 0) {
	goto L290;
    }

/*  ***  fresh start  *** */

    k = 0;
    uk = negone;
    *ka = 0;
    kalim = 50;
    v[dgnorm] = v2norm_(p, &dig[1]);
    v[nreduc] = zero;
    v[dst0] = zero;
    kamin = 3;
    if (v[dgnorm] == zero) {
	kamin = 0;
    }

/*     ***  store diag(dihdi) in w(diag0+1),...,w(diag0+p)  *** */

    j = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j += i__;
	k1 = diag0 + i__;
	w[k1] = dihdi[j];
/* L10: */
    }

/*     ***  determine w(dggdmx), the largest element of dihdi  *** */

    t1 = zero;
    j = *p * (*p + 1) / 2;
    i__1 = j;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = (d__1 = dihdi[i__], abs(d__1));
	if (t1 < t) {
	    t1 = t;
	}
/* L20: */
    }
    w[dggdmx] = t1;

/*  ***  try alpha = 0  *** */

L30:
    lsqrt_(&c__1, p, &l[1], &dihdi[1], &irc);
    if (irc == 0) {
	goto L50;
    }
/*        ***  indef. h -- underestimate smallest eigenvalue, use this */
/*        ***  estimate to initialize lower bound lk on alpha. */
    j = irc * (irc + 1) / 2;
    t = l[j];
    l[j] = one;
    i__1 = irc;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	w[i__] = zero;
    }
    w[irc] = one;
    litvmu_(&irc, &w[1], &l[1], &w[1]);
    t1 = v2norm_(&irc, &w[1]);
    lk = -t / t1 / t1;
    v[dst0] = -lk;
    if (restrt) {
	goto L210;
    }
    goto L70;

/*     ***  positive definite h -- compute unmodified newton step.  *** */
L50:
    lk = zero;
    t = lsvmin_(p, &l[1], &w[q], &w[q]);
    if (t >= one) {
	goto L60;
    }
    if (big <= zero) {
	big = rmdcon_(&c__6);
    }
    if (v[dgnorm] >= t * t * big) {
	goto L70;
    }
L60:
    livmul_(p, &w[q], &l[1], &dig[1]);
    gtsta = dotprd_(p, &w[q], &w[q]);
    v[nreduc] = half * gtsta;
    litvmu_(p, &w[q], &l[1], &w[q]);
    dst = v2norm_(p, &w[q]);
    v[dst0] = dst;
    phi = dst - rad;
    if (phi <= phimax) {
	goto L260;
    }
    if (restrt) {
	goto L210;
    }

/*  ***  prepare to compute gerschgorin estimates of largest (and */
/*  ***  smallest) eigenvalues.  *** */

L70:
    k = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wi = zero;
	if (i__ == 1) {
	    goto L90;
	}
	im1 = i__ - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    ++k;
	    t = (d__1 = dihdi[k], abs(d__1));
	    wi += t;
	    w[j] += t;
/* L80: */
	}
L90:
	w[i__] = wi;
	++k;
/* L100: */
    }

/*  ***  (under-)estimate smallest eigenvalue of (d**-1)*h*(d**-1)  *** */

    k = 1;
    t1 = w[diag] - w[1];
    if (*p <= 1) {
	goto L120;
    }
    i__1 = *p;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = diag0 + i__;
	t = w[j] - w[i__];
	if (t >= t1) {
	    goto L110;
	}
	t1 = t;
	k = i__;
L110:
	;
    }

L120:
    sk = w[k];
    j = diag0 + k;
    akk = w[j];
    k1 = k * (k - 1) / 2 + 1;
    inc = 1;
    t = zero;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == k) {
	    goto L130;
	}
	aki = (d__1 = dihdi[k1], abs(d__1));
	si = w[i__];
	j = diag0 + i__;
	t1 = half * (akk - w[j] + si - aki);
	t1 += sqrt(t1 * t1 + sk * aki);
	if (t < t1) {
	    t = t1;
	}
	if (i__ < k) {
	    goto L140;
	}
L130:
	inc = i__;
L140:
	k1 += inc;
/* L150: */
    }

    w[emin] = akk - t;
    uk = v[dgnorm] / rad - w[emin];
    if (v[dgnorm] == zero) {
	uk = uk + p001 + p001 * uk;
    }
    if (uk <= zero) {
	uk = p001;
    }

/*  ***  compute gerschgorin (over-)estimate of largest eigenvalue  *** */

    k = 1;
    t1 = w[diag] + w[1];
    if (*p <= 1) {
	goto L170;
    }
    i__1 = *p;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = diag0 + i__;
	t = w[j] + w[i__];
	if (t <= t1) {
	    goto L160;
	}
	t1 = t;
	k = i__;
L160:
	;
    }

L170:
    sk = w[k];
    j = diag0 + k;
    akk = w[j];
    k1 = k * (k - 1) / 2 + 1;
    inc = 1;
    t = zero;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == k) {
	    goto L180;
	}
	aki = (d__1 = dihdi[k1], abs(d__1));
	si = w[i__];
	j = diag0 + i__;
	t1 = half * (w[j] + si - aki - akk);
	t1 += sqrt(t1 * t1 + sk * aki);
	if (t < t1) {
	    t = t1;
	}
	if (i__ < k) {
	    goto L190;
	}
L180:
	inc = i__;
L190:
	k1 += inc;
/* L200: */
    }

    w[emax] = akk + t;
/* Computing MAX */
    d__1 = lk, d__2 = v[dgnorm] / rad - w[emax];
    lk = max(d__1,d__2);

/*     ***  alphak = current value of alpha (see alg. notes above).  we */
/*     ***  use more*s scheme for initializing it. */
    alphak = (d__1 = v[stppar], abs(d__1)) * v[rad0] / rad;

    if (irc != 0) {
	goto L210;
    }

/*  ***  compute l0 for positive definite h  *** */

    livmul_(p, &w[1], &l[1], &w[q]);
    t = v2norm_(p, &w[1]);
    w[phipin] = dst / t / t;
/* Computing MAX */
    d__1 = lk, d__2 = phi * w[phipin];
    lk = max(d__1,d__2);

/*  ***  safeguard alphak and add alphak*i to (d**-1)*h*(d**-1)  *** */

L210:
    ++(*ka);
    if (-v[dst0] >= alphak || alphak < lk || alphak >= uk) {
/* Computing MAX */
	d__1 = p001, d__2 = sqrt(lk / uk);
	alphak = uk * max(d__1,d__2);
    }
    if (alphak <= zero) {
	alphak = half * uk;
    }
    if (alphak <= zero) {
	alphak = uk;
    }
    k = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k += i__;
	j = diag0 + i__;
	dihdi[k] = w[j] + alphak;
/* L220: */
    }

/*  ***  try computing cholesky decomposition  *** */

    lsqrt_(&c__1, p, &l[1], &dihdi[1], &irc);
    if (irc == 0) {
	goto L240;
    }

/*  ***  (d**-1)*h*(d**-1) + alphak*i  is indefinite -- overestimate */
/*  ***  smallest eigenvalue for use in updating lk  *** */

    j = irc * (irc + 1) / 2;
    t = l[j];
    l[j] = one;
    i__1 = irc;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L230: */
	w[i__] = zero;
    }
    w[irc] = one;
    litvmu_(&irc, &w[1], &l[1], &w[1]);
    t1 = v2norm_(&irc, &w[1]);
    lk = alphak - t / t1 / t1;
    v[dst0] = -lk;
    goto L210;

/*  ***  alphak makes (d**-1)*h*(d**-1) positive definite. */
/*  ***  compute q = -d*step, check for convergence.  *** */

L240:
    livmul_(p, &w[q], &l[1], &dig[1]);
    gtsta = dotprd_(p, &w[q], &w[q]);
    litvmu_(p, &w[q], &l[1], &w[q]);
    dst = v2norm_(p, &w[q]);
    phi = dst - rad;
    if (phi <= phimax && phi >= phimin) {
	goto L270;
    }
    if (phi == oldphi) {
	goto L270;
    }
    oldphi = phi;
    if (phi < zero) {
	goto L330;
    }

/*  ***  unacceptable alphak -- update lk, uk, alphak  *** */

L250:
    if (*ka >= kalim) {
	goto L270;
    }
/*     ***  the following dmin1 is necessary because of restarts  *** */
    if (phi < zero) {
	uk = min(uk,alphak);
    }
/*     *** kamin = 0 only iff the gradient vanishes  *** */
    if (kamin == 0) {
	goto L210;
    }
    livmul_(p, &w[1], &l[1], &w[q]);
    t1 = v2norm_(p, &w[1]);
    alphak += phi / t1 * (dst / t1) * (dst / rad);
    lk = max(lk,alphak);
    goto L210;

/*  ***  acceptable step on first try  *** */

L260:
    alphak = zero;

/*  ***  successful step in general.  compute step = -(d**-1)*q  *** */

L270:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = q0 + i__;
	step[i__] = -w[j] / d__[i__];
/* L280: */
    }
    v[gtstep] = -gtsta;
    v[preduc] = half * (abs(alphak) * dst * dst + gtsta);
    goto L410;


/*  ***  restart with new radius  *** */

L290:
    if (v[dst0] <= zero || v[dst0] - rad > phimax) {
	goto L310;
    }

/*     ***  prepare to return newton step  *** */

    restrt = TRUE_;
    ++(*ka);
    k = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k += i__;
	j = diag0 + i__;
	dihdi[k] = w[j];
/* L300: */
    }
    uk = negone;
    goto L30;

L310:
    kamin = *ka + 3;
    if (v[dgnorm] == zero) {
	kamin = 0;
    }
    if (*ka == 0) {
	goto L50;
    }

    dst = w[dstsav];
    alphak = (d__1 = v[stppar], abs(d__1));
    phi = dst - rad;
    t = v[dgnorm] / rad;
    uk = t - w[emin];
    if (v[dgnorm] == zero) {
	uk = uk + p001 + p001 * uk;
    }
    if (uk <= zero) {
	uk = p001;
    }
    if (rad > v[rad0]) {
	goto L320;
    }

/*        ***  smaller radius  *** */
    lk = zero;
    if (alphak > zero) {
	lk = w[lk0];
    }
/* Computing MAX */
    d__1 = lk, d__2 = t - w[emax];
    lk = max(d__1,d__2);
    if (v[dst0] > zero) {
/* Computing MAX */
	d__1 = lk, d__2 = (v[dst0] - rad) * w[phipin];
	lk = max(d__1,d__2);
    }
    goto L250;

/*     ***  bigger radius  *** */
L320:
    if (alphak > zero) {
/* Computing MIN */
	d__1 = uk, d__2 = w[uk0];
	uk = min(d__1,d__2);
    }
/* Computing MAX */
    d__1 = zero, d__2 = -v[dst0], d__1 = max(d__1,d__2), d__2 = t - w[emax];
    lk = max(d__1,d__2);
    if (v[dst0] > zero) {
/* Computing MAX */
	d__1 = lk, d__2 = (v[dst0] - rad) * w[phipin];
	lk = max(d__1,d__2);
    }
    goto L250;

/*  ***  decide whether to check for special case... in practice (from */
/*  ***  the standpoint of the calling optimization code) it seems best */
/*  ***  not to check until a few iterations have failed -- hence the */
/*  ***  test on kamin below. */

L330:
/* Computing MIN */
    d__1 = zero, d__2 = v[dst0];
    delta = alphak + min(d__1,d__2);
    twopsi = alphak * dst * dst + gtsta;
    if (*ka >= kamin) {
	goto L340;
    }
/*     *** if the test in ref. 2 is satisfied, fall through to handle */
/*     *** the special case (as soon as the more-sorensen test detects */
/*     *** it). */
    if (delta >= psifac * twopsi) {
	goto L370;
    }

/*  ***  check for the special case of  h + alpha*d**2  (nearly) */
/*  ***  singular.  use one step of inverse power method with start */
/*  ***  from lsvmin to obtain approximate eigenvector corresponding */
/*  ***  to smallest eigenvalue of (d**-1)*h*(d**-1).  lsvmin returns */
/*  ***  x and w with  l*w = x. */

L340:
    t = lsvmin_(p, &l[1], &w[x], &w[1]);

/*     ***  normalize w  *** */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L350: */
	w[i__] = t * w[i__];
    }
/*     ***  complete current inv. power iter. -- replace w by (l**-t)*w. */
    litvmu_(p, &w[1], &l[1], &w[1]);
    t2 = one / v2norm_(p, &w[1]);
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L360: */
	w[i__] = t2 * w[i__];
    }
    t = t2 * t;

/*  ***  now w is the desired approximate (unit) eigenvector and */
/*  ***  t*x = ((d**-1)*h*(d**-1) + alphak*i)*w. */

    sw = dotprd_(p, &w[q], &w[1]);
    t1 = (rad + dst) * (rad - dst);
    root = sqrt(sw * sw + t1);
    if (sw < zero) {
	root = -root;
    }
    si = t1 / (sw + root);

/*  ***  the actual test for the special case... */

/* Computing 2nd power */
    d__1 = t2 * si;
/* Computing 2nd power */
    d__2 = dst;
    if (d__1 * d__1 <= eps * (d__2 * d__2 + alphak * radsq)) {
	goto L380;
    }

/*  ***  update upper bound on smallest eigenvalue (when not positive) */
/*  ***  (as recommended by more and sorensen) and continue... */

    if (v[dst0] <= zero) {
/* Computing MIN */
/* Computing 2nd power */
	d__3 = t2;
	d__1 = v[dst0], d__2 = d__3 * d__3 - alphak;
	v[dst0] = min(d__1,d__2);
    }
/* Computing MAX */
    d__1 = lk, d__2 = -v[dst0];
    lk = max(d__1,d__2);

/*  ***  check whether we can hope to detect the special case in */
/*  ***  the available arithmetic.  accept step as it is if not. */

/*     ***  if not yet available, obtain machine dependent value dgxfac. */
L370:
    if (dgxfac == zero) {
	dgxfac = epsfac * rmdcon_(&c__3);
    }

    if (delta > dgxfac * w[dggdmx]) {
	goto L250;
    }
    goto L270;

/*  ***  special case detected... negate alphak to indicate special case */

L380:
    alphak = -alphak;
    v[preduc] = half * twopsi;

/*  ***  accept current step if adding si*w would lead to a */
/*  ***  further relative reduction in psi of less than v(epslon)/3. */

    t1 = zero;
    t = si * (alphak * sw - half * si * (alphak + t * dotprd_(p, &w[x], &w[1])
	    ));
    if (t < eps * twopsi / six) {
	goto L390;
    }
    v[preduc] += t;
    dst = rad;
    t1 = -si;
L390:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = q0 + i__;
	w[j] = t1 * w[i__] - w[j];
	step[i__] = w[j] / d__[i__];
/* L400: */
    }
    v[gtstep] = dotprd_(p, &dig[1], &w[q]);

/*  ***  save values for use in a possible restart  *** */

L410:
    v[dstnrm] = dst;
    v[stppar] = alphak;
    w[lk0] = lk;
    w[uk0] = uk;
    v[rad0] = rad;
    w[dstsav] = dst;

/*     ***  restore diagonal of dihdi  *** */

    j = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j += i__;
	k = diag0 + i__;
	dihdi[j] = w[k];
/* L420: */
    }

/* L999: */
    return 0;

/*  ***  last card of gqtst follows  *** */
} /* gqtst_ */

/* Subroutine */ int lsqrt_(integer *n1, integer *n, doublereal *l, 
	doublereal *a, integer *irc)
{
    /* Initialized data */

    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static __thread integer i__, j, k;
    static __thread doublereal t;
    static __thread integer i0, j0, ij, ik, jk;
    static __thread doublereal td;
    static __thread integer im1, jm1;


/*  ***  compute rows n1 through n of the cholesky factor  l  of */
/*  ***  a = l*(l**t),  where  l  and the lower triangle of  a  are both */
/*  ***  stored compactly by rows (and may occupy the same storage). */
/*  ***  irc = 0 means all went well.  irc = j means the leading */
/*  ***  principal  j x j  submatrix of  a  is not positive definite -- */
/*  ***  and  l(j*(j+1)/2)  contains the (nonpos.) reduced j-th diagonal. */

/*  ***  parameters  *** */

/*     dimension l(n*(n+1)/2), a(n*(n+1)/2) */

/*  ***  local variables  *** */


/*  ***  intrinsic functions  *** */
/* /+ */
/* / */
/* /6 */
    /* Parameter adjustments */
    --a;
    --l;

    /* Function Body */
/* /7 */
/*     parameter (zero=0.d+0) */
/* / */

/*  ***  body  *** */

    i0 = *n1 * (*n1 - 1) / 2;
    i__1 = *n;
    for (i__ = *n1; i__ <= i__1; ++i__) {
	td = zero;
	if (i__ == 1) {
	    goto L40;
	}
	j0 = 0;
	im1 = i__ - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    t = zero;
	    if (j == 1) {
		goto L20;
	    }
	    jm1 = j - 1;
	    i__3 = jm1;
	    for (k = 1; k <= i__3; ++k) {
		ik = i0 + k;
		jk = j0 + k;
		t += l[ik] * l[jk];
/* L10: */
	    }
L20:
	    ij = i0 + j;
	    j0 += j;
	    t = (a[ij] - t) / l[j0];
	    l[ij] = t;
	    td += t * t;
/* L30: */
	}
L40:
	i0 += i__;
	t = a[i0] - td;
	if (t <= zero) {
	    goto L60;
	}
	l[i0] = sqrt(t);
/* L50: */
    }

    *irc = 0;
    goto L999;

L60:
    l[i0] = t;
    *irc = i__;

L999:
    return 0;

/*  ***  last card of lsqrt  *** */
} /* lsqrt_ */

doublereal lsvmin_(integer *p, doublereal *l, doublereal *x, doublereal *y)
{
    /* Initialized data */

    static __thread doublereal half = .5;
    static __thread doublereal one = 1.;
    static __thread doublereal r9973 = 9973.;
    static __thread doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Local variables */
    static __thread doublereal b;
    static __thread integer i__, j;
    static __thread doublereal t;
    static __thread integer j0, ii, ji, jj, ix, jm1, pm1, jjj;
    static __thread doublereal splus;
    extern /* Subroutine */ int vaxpy_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static __thread doublereal xplus;
    extern doublereal v2norm_(integer *, doublereal *), dotprd_(integer *, 
	    doublereal *, doublereal *);
    static __thread doublereal sminus, xminus;


/*  ***  estimate smallest sing. value of packed lower triang. matrix l */

/*  ***  parameter declarations  *** */

/*     dimension l(p*(p+1)/2) */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  purpose  *** */

/*     this function returns a good over-estimate of the smallest */
/*     singular value of the packed lower triangular matrix l. */

/*  ***  parameter description  *** */

/*  p (in)  = the order of l.  l is a  p x p  lower triangular matrix. */
/*  l (in)  = array holding the elements of  l  in row order, i.e. */
/*             l(1,1), l(2,1), l(2,2), l(3,1), l(3,2), l(3,3), etc. */
/*  x (out) if lsvmin returns a positive value, then x is a normalized */
/*             approximate left singular vector corresponding to the */
/*             smallest singular value.  this approximation may be very */
/*             crude.  if lsvmin returns zero, then some components of x */
/*             are zero and the rest retain their input values. */
/*  y (out) if lsvmin returns a positive value, then y = (l**-1)*x is an */
/*             unnormalized approximate right singular vector correspond- */
/*             ing to the smallest singular value.  this approximation */
/*             may be crude.  if lsvmin returns zero, then y retains its */
/*             input value.  the caller may pass the same vector for x */
/*             and y (nonstandard fortran usage), in which case y over- */
/*             writes x (for nonzero lsvmin returns). */

/*  ***  algorithm notes  *** */

/*     the algorithm is based on (1), with the additional provision that */
/*     lsvmin = 0 is returned if the smallest diagonal element of l */
/*     (in magnitude) is not more than the unit roundoff times the */
/*     largest.  the algorithm uses a random number generator proposed */
/*     in (4), which passes the spectral test with flying colors -- see */
/*     (2) and (3). */

/*  ***  subroutines and functions called  *** */

/*        v2norm - function, returns the 2-norm of a vector. */

/*  ***  references  *** */

/*     (1) cline, a., moler, c., stewart, g., and wilkinson, j.h.(1977), */
/*         an estimate for the condition number of a matrix, report */
/*         tm-310, applied math. div., argonne national laboratory. */

/*     (2) hoaglin, d.c. (1976), theoretical properties of congruential */
/*         random-number generators --  an empirical view, */
/*         memorandum ns-340, dept. of statistics, harvard univ. */

/*     (3) knuth, d.e. (1969), the art of computer programming, vol. 2 */
/*         (seminumerical algorithms), addison-wesley, reading, mass. */

/*     (4) smith, c.s. (1971), multiplicative pseudo-random number */
/*         generators with prime modulus, j. assoc. comput. mach. 18, */
/*         pp. 586-593. */

/*  ***  history  *** */

/*     designed and coded by david m. gay (winter 1977/summer 1978). */

/*  ***  general  *** */

/*     this subroutine was written in connection with research */
/*     supported by the national science foundation under grants */
/*     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989. */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  local variables  *** */


/*  ***  constants  *** */


/*  ***  intrinsic functions  *** */
/* /+ */
/* / */
/*  ***  external functions and subroutines  *** */


/* /6 */
    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
/* /7 */
/*     parameter (half=0.5d+0, one=1.d+0, r9973=9973.d+0, zero=0.d+0) */
/* / */

/*  ***  body  *** */

    ix = 2;
    pm1 = *p - 1;

/*  ***  first check whether to return lsvmin = 0 and initialize x  *** */

    ii = 0;
    j0 = *p * pm1 / 2;
    jj = j0 + *p;
    if (l[jj] == zero) {
	goto L110;
    }
    ix = ix * 3432 % 9973;
    b = half * (one + (real) ix / r9973);
    xplus = b / l[jj];
    x[*p] = xplus;
    if (*p <= 1) {
	goto L60;
    }
    i__1 = pm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii += i__;
	if (l[ii] == zero) {
	    goto L110;
	}
	ji = j0 + i__;
	x[i__] = xplus * l[ji];
/* L10: */
    }

/*  ***  solve (l**t)*x = b, where the components of b have randomly */
/*  ***  chosen magnitudes in (.5,1) with signs chosen to make x large. */

/*     do j = p-1 to 1 by -1... */
    i__1 = pm1;
    for (jjj = 1; jjj <= i__1; ++jjj) {
	j = *p - jjj;
/*       ***  determine x(j) in this iteration. note for i = 1,2,...,j */
/*       ***  that x(i) holds the current partial sum for row i. */
	ix = ix * 3432 % 9973;
	b = half * (one + (real) ix / r9973);
	xplus = b - x[j];
	xminus = -b - x[j];
	splus = abs(xplus);
	sminus = abs(xminus);
	jm1 = j - 1;
	j0 = j * jm1 / 2;
	jj = j0 + j;
	xplus /= l[jj];
	xminus /= l[jj];
	if (jm1 == 0) {
	    goto L30;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ji = j0 + i__;
	    splus += (d__1 = x[i__] + l[ji] * xplus, abs(d__1));
	    sminus += (d__1 = x[i__] + l[ji] * xminus, abs(d__1));
/* L20: */
	}
L30:
	if (sminus > splus) {
	    xplus = xminus;
	}
	x[j] = xplus;
/*       ***  update partial sums  *** */
	if (jm1 > 0) {
	    vaxpy_(&jm1, &x[1], &xplus, &l[j0 + 1], &x[1]);
	}
/* L50: */
    }

/*  ***  normalize x  *** */

L60:
    t = one / v2norm_(p, &x[1]);
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L70: */
	x[i__] = t * x[i__];
    }

/*  ***  solve l*y = x and return lsvmin = 1/twonorm(y)  *** */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j - 1;
	j0 = j * jm1 / 2;
	jj = j0 + j;
	t = zero;
	if (jm1 > 0) {
	    t = dotprd_(&jm1, &l[j0 + 1], &y[1]);
	}
	y[j] = (x[j] - t) / l[jj];
/* L100: */
    }

    ret_val = one / v2norm_(p, &y[1]);
    goto L999;

L110:
    ret_val = zero;
L999:
    return ret_val;
/*  ***  last card of lsvmin follows  *** */
} /* lsvmin_ */

/* Subroutine */ int slvmul_(integer *p, doublereal *y, doublereal *s, 
	doublereal *x)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static __thread integer i__, j, k;
    static __thread doublereal xi;
    static __thread integer im1;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);


/*  ***  set  y = s * x,  s = p x p symmetric matrix.  *** */
/*  ***  lower triangle of  s  stored rowwise.         *** */

/*  ***  parameter declarations  *** */

/*     dimension s(p*(p+1)/2) */

/*  ***  local variables  *** */


/*  ***  no intrinsic functions  *** */

/*  ***  external function  *** */


/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --x;
    --y;
    --s;

    /* Function Body */
    j = 1;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = dotprd_(&i__, &s[j], &x[1]);
	j += i__;
/* L10: */
    }

    if (*p <= 1) {
	goto L999;
    }
    j = 1;
    i__1 = *p;
    for (i__ = 2; i__ <= i__1; ++i__) {
	xi = x[i__];
	im1 = i__ - 1;
	++j;
	i__2 = im1;
	for (k = 1; k <= i__2; ++k) {
	    y[k] += s[j] * xi;
	    ++j;
/* L30: */
	}
/* L40: */
    }

L999:
    return 0;
/*  ***  last card of slvmul follows  *** */
} /* slvmul_ */

doublereal d1mach_(integer *i__)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static __thread cilist io___712 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___713 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___714 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___715 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___716 = { 0, 6, 0, "(a,i12)", 0 };
    static __thread cilist io___717 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___718 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___719 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___720 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___721 = { 0, 6, 0, "(a,i12)", 0 };


/* *********************************************************************72 */

/* c D1MACH returns double precision real machine-dependent constants. */

/*  Discussion: */

/*    D1MACH can be used to obtain machine-dependent parameters */
/*    for the local machine environment.  It is a function */
/*    with one input argument, and can be called as follows: */

/*      D = D1MACH ( I ) */

/*    where I=1,...,5.  The output value of D above is */
/*    determined by the input value of I:. */

/*    D1MACH ( 1) = B**(EMIN-1), the smallest positive magnitude. */
/*    D1MACH ( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*    D1MACH ( 3) = B**(-T), the smallest relative spacing. */
/*    D1MACH ( 4) = B**(1-T), the largest relative spacing. */
/*    D1MACH ( 5) = LOG10(B) */

/*  Licensing: */

/*    This code is distributed under the GNU LGPL license. */

/*  Modified: */

/*    25 April 2007 */

/*  Author: */

/*    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer. */
/*    This FORTRAN77 version by John Burkardt. */

/*  Reference: */

/*    Phyllis Fox, Andrew Hall, Norman Schryer, */
/*    Algorithm 528: */
/*    Framework for a Portable Library, */
/*    ACM Transactions on Mathematical Software, */
/*    Volume 4, Number 2, June 1978, page 176-188. */

/*  Parameters: */

/*    Input, integer I, the index of the desired constant. */

/*    Output, double precision D1MACH, the value of the constant. */

    if (*i__ < 1) {
	s_wsfe(&io___712);
	do_fio(&c__1, " ", (ftnlen)1);
	e_wsfe();
	s_wsfe(&io___713);
	do_fio(&c__1, "D1MACH - Fatal error!", (ftnlen)21);
	e_wsfe();
	s_wsfe(&io___714);
	do_fio(&c__1, "  The input argument I is out of bounds.", (ftnlen)40);
	e_wsfe();
	s_wsfe(&io___715);
	do_fio(&c__1, "  Legal values satisfy 1 <= I <= 5.", (ftnlen)35);
	e_wsfe();
	s_wsfe(&io___716);
	do_fio(&c__1, "  I = ", (ftnlen)6);
	do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	e_wsfe();
	ret_val = 0.;
	s_stop("", (ftnlen)0);
    } else if (*i__ == 1) {
	ret_val = 4.450147717014403e-308;
    } else if (*i__ == 2) {
	ret_val = 8.988465674311579e307;
    } else if (*i__ == 3) {
	ret_val = 1.110223024625157e-16;
    } else if (*i__ == 4) {
	ret_val = 2.220446049250313e-16;
    } else if (*i__ == 5) {
	ret_val = .301029995663981;
    } else if (5 < *i__) {
	s_wsfe(&io___717);
	do_fio(&c__1, " ", (ftnlen)1);
	e_wsfe();
	s_wsfe(&io___718);
	do_fio(&c__1, "D1MACH - Fatal error!", (ftnlen)21);
	e_wsfe();
	s_wsfe(&io___719);
	do_fio(&c__1, "  The input argument I is out of bounds.", (ftnlen)40);
	e_wsfe();
	s_wsfe(&io___720);
	do_fio(&c__1, "  Legal values satisfy 1 <= I <= 5.", (ftnlen)35);
	e_wsfe();
	s_wsfe(&io___721);
	do_fio(&c__1, "  I = ", (ftnlen)6);
	do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	e_wsfe();
	ret_val = 0.;
	s_stop("", (ftnlen)0);
    }
    return ret_val;
} /* d1mach_ */

integer i1mach_(integer *i__)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static __thread cilist io___722 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___723 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___724 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___725 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___726 = { 0, 6, 0, "(a,i12)", 0 };
    static __thread cilist io___727 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___728 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___729 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___730 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___731 = { 0, 6, 0, "(a,i12)", 0 };


/* *********************************************************************72 */

/* c I1MACH returns integer machine dependent constants. */

/*  Discussion: */

/*    Input/output unit numbers. */

/*      I1MACH(1) = the standard input unit. */
/*      I1MACH(2) = the standard output unit. */
/*      I1MACH(3) = the standard punch unit. */
/*      I1MACH(4) = the standard error message unit. */

/*    Words. */

/*      I1MACH(5) = the number of bits per integer storage unit. */
/*      I1MACH(6) = the number of characters per integer storage unit. */

/*    Integers. */

/*    Assume integers are represented in the S digit base A form: */

/*      Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0)) */

/*    where 0 <= X(1:S-1) < A. */

/*      I1MACH(7) = A, the base. */
/*      I1MACH(8) = S, the number of base A digits. */
/*      I1MACH(9) = A**S-1, the largest integer. */

/*    Floating point numbers */

/*    Assume floating point numbers are represented in the T digit */
/*    base B form: */

/*      Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) ) */

/*    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX. */

/*      I1MACH(10) = B, the base. */

/*    Single precision */

/*      I1MACH(11) = T, the number of base B digits. */
/*      I1MACH(12) = EMIN, the smallest exponent E. */
/*      I1MACH(13) = EMAX, the largest exponent E. */

/*    Double precision */

/*      I1MACH(14) = T, the number of base B digits. */
/*      I1MACH(15) = EMIN, the smallest exponent E. */
/*      I1MACH(16) = EMAX, the largest exponent E. */

/*  Licensing: */

/*    This code is distributed under the GNU LGPL license. */

/*  Modified: */

/*    25 April 2007 */

/*  Author: */

/*    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer. */
/*    This FORTRAN77 version by John Burkardt. */

/*  Reference: */

/*    Phyllis Fox, Andrew Hall, Norman Schryer, */
/*    Algorithm 528, */
/*    Framework for a Portable Library, */
/*    ACM Transactions on Mathematical Software, */
/*    Volume 4, Number 2, June 1978, page 176-188. */

/*  Parameters: */

/*    Input, integer I, chooses the parameter to be returned. */
/*    1 <= I <= 16. */

/*    Output, integer I1MACH, the value of the chosen parameter. */

    if (*i__ < 1) {
	s_wsfe(&io___722);
	do_fio(&c__1, " ", (ftnlen)1);
	e_wsfe();
	s_wsfe(&io___723);
	do_fio(&c__1, "I1MACH - Fatal error!", (ftnlen)21);
	e_wsfe();
	s_wsfe(&io___724);
	do_fio(&c__1, "  The input argument I is out of bounds.", (ftnlen)40);
	e_wsfe();
	s_wsfe(&io___725);
	do_fio(&c__1, "  Legal values satisfy 1 <= I <= 16.", (ftnlen)36);
	e_wsfe();
	s_wsfe(&io___726);
	do_fio(&c__1, "  I = ", (ftnlen)6);
	do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	e_wsfe();
	ret_val = 0;
	s_stop("", (ftnlen)0);
    } else if (*i__ == 1) {
	ret_val = 5;
    } else if (*i__ == 2) {
	ret_val = 6;
    } else if (*i__ == 3) {
	ret_val = 7;
    } else if (*i__ == 4) {
	ret_val = 6;
    } else if (*i__ == 5) {
	ret_val = 32;
    } else if (*i__ == 6) {
	ret_val = 4;
    } else if (*i__ == 7) {
	ret_val = 2;
    } else if (*i__ == 8) {
	ret_val = 31;
    } else if (*i__ == 9) {
	ret_val = 2147483647;
    } else if (*i__ == 10) {
	ret_val = 2;
    } else if (*i__ == 11) {
	ret_val = 24;
    } else if (*i__ == 12) {
	ret_val = -125;
    } else if (*i__ == 13) {
	ret_val = 128;
    } else if (*i__ == 14) {
	ret_val = 53;
    } else if (*i__ == 15) {
	ret_val = -1021;
    } else if (*i__ == 16) {
	ret_val = 1024;
    } else if (16 < *i__) {
	s_wsfe(&io___727);
	do_fio(&c__1, " ", (ftnlen)1);
	e_wsfe();
	s_wsfe(&io___728);
	do_fio(&c__1, "I1MACH - Fatal error!", (ftnlen)21);
	e_wsfe();
	s_wsfe(&io___729);
	do_fio(&c__1, "  The input argument I is out of bounds.", (ftnlen)40);
	e_wsfe();
	s_wsfe(&io___730);
	do_fio(&c__1, "  Legal values satisfy 1 <= I <= 16.", (ftnlen)36);
	e_wsfe();
	s_wsfe(&io___731);
	do_fio(&c__1, "  I = ", (ftnlen)6);
	do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	e_wsfe();
	ret_val = 0;
	s_stop("", (ftnlen)0);
    }
    return ret_val;
} /* i1mach_ */

doublereal r1mach_(integer *i__)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static __thread cilist io___732 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___733 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___734 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___735 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___736 = { 0, 6, 0, "(a,i12)", 0 };
    static __thread cilist io___737 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___738 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___739 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___740 = { 0, 6, 0, "(a)", 0 };
    static __thread cilist io___741 = { 0, 6, 0, "(a,i12)", 0 };


/* *********************************************************************72 */

/* c R1MACH returns single precision real machine constants. */

/*  Discussion: */

/*    Assume that single precision real numbers are stored with a mantissa */
/*    of T digits in base B, with an exponent whose value must lie */
/*    between EMIN and EMAX.  Then for values of I between 1 and 5, */
/*    R1MACH will return the following values: */

/*      R1MACH(1) = B**(EMIN-1), the smallest positive magnitude. */
/*      R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude. */
/*      R1MACH(3) = B**(-T), the smallest relative spacing. */
/*      R1MACH(4) = B**(1-T), the largest relative spacing. */
/*      R1MACH(5) = log10(B) */

/*  Licensing: */

/*    This code is distributed under the GNU LGPL license. */

/*  Modified: */

/*    25 April 2007 */

/*  Author: */

/*    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer. */
/*    This FORTRAN77 version by John Burkardt. */

/*  Reference: */

/*    Phyllis Fox, Andrew Hall, Norman Schryer, */
/*    Algorithm 528, */
/*    Framework for a Portable Library, */
/*    ACM Transactions on Mathematical Software, */
/*    Volume 4, Number 2, June 1978, page 176-188. */

/*  Parameters: */

/*    Input, integer I, chooses the parameter to be returned. */
/*    1 <= I <= 5. */

/*    Output, real R1MACH, the value of the chosen parameter. */

    if (*i__ < 1) {
	s_wsfe(&io___732);
	do_fio(&c__1, " ", (ftnlen)1);
	e_wsfe();
	s_wsfe(&io___733);
	do_fio(&c__1, "R1MACH - Fatal error!", (ftnlen)21);
	e_wsfe();
	s_wsfe(&io___734);
	do_fio(&c__1, "  The input argument I is out of bounds.", (ftnlen)40);
	e_wsfe();
	s_wsfe(&io___735);
	do_fio(&c__1, "  Legal values satisfy 1 <= I <= 5.", (ftnlen)35);
	e_wsfe();
	s_wsfe(&io___736);
	do_fio(&c__1, "  I = ", (ftnlen)6);
	do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	e_wsfe();
	ret_val = (float)0.;
	s_stop("", (ftnlen)0);
    } else if (*i__ == 1) {
	ret_val = (float)1.1754944e-38;
    } else if (*i__ == 2) {
	ret_val = (float)3.4028235e38;
    } else if (*i__ == 3) {
	ret_val = (float)5.9604645e-8;
    } else if (*i__ == 4) {
	ret_val = (float)1.1920929e-7;
    } else if (*i__ == 5) {
	ret_val = (float).30103;
    } else if (5 < *i__) {
	s_wsfe(&io___737);
	do_fio(&c__1, " ", (ftnlen)1);
	e_wsfe();
	s_wsfe(&io___738);
	do_fio(&c__1, "R1MACH - Fatal error!", (ftnlen)21);
	e_wsfe();
	s_wsfe(&io___739);
	do_fio(&c__1, "  The input argument I is out of bounds.", (ftnlen)40);
	e_wsfe();
	s_wsfe(&io___740);
	do_fio(&c__1, "  Legal values satisfy 1 <= I <= 5.", (ftnlen)35);
	e_wsfe();
	s_wsfe(&io___741);
	do_fio(&c__1, "  I = ", (ftnlen)6);
	do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	e_wsfe();
	ret_val = (float)0.;
	s_stop("", (ftnlen)0);
    }
    return ret_val;
} /* r1mach_ */

#ifdef __cplusplus
	}
#endif
