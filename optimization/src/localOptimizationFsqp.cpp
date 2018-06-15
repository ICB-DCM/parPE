#include "localOptimizationFsqp.h"

#include <cmath>
#include <iostream>
#include <assert.h>
#include <cstring>
#include <stdexcept>
#include <optimizationResultWriter.h>
#include <misc.h>
#include <logging.h>

extern "C" {
#include <f2c.h>
#undef max
}

// callback functions types for FFSQP, see below
using objType    = void (*) (integer &nparam, integer &j, doublereal *x, doublereal &fj);
using constrType = void (*) (integer &nparam, integer &j, doublereal *x, doublereal &gj);
using gradobType = void (*) (integer &nparam, integer &j, doublereal *x, doublereal *gradfj, doublereal *dummy);
using gradcnType = void (*) (integer &nparam, integer &j, doublereal *x, doublereal *gradgj, doublereal *dummy);


extern "C" {
// These two functions are needed for linking by ql0001, so far they were never called though
integer lnblnk_(char *a, ftnlen) {
    throw std::runtime_error("FSQP: lnblnk_ was called but is not implemented");

    return strlen(a);
}
int basout_(integer *, integer *, char *, ftnlen) {
    // TODO: should print?
    throw std::runtime_error("FSQP: basout_ was called but is not implemented");

    return 0;
}

/**
 * @brief The ffsqp routine. See OptimizerFsqp::optimize for documentation.
 * @param nparam
 * @param nf
 * @param nineqn
 * @param nineq
 * @param neqn
 * @param neq
 * @param mode
 * @param iprint
 * @param miter
 * @param inform__
 * @param bigbnd
 * @param eps
 * @param epseqn
 * @param udelta
 * @param bl
 * @param bu
 * @param x
 * @param f
 * @param g
 * @param iw
 * @param iwsize
 * @param w
 * @param nwsize
 * @param obj
 * @param constr
 * @param gradob
 * @param gradcn
 * @return
 */
int ffsqp_(integer &nparam, integer &nf, integer &nineqn,
           integer &nineq, integer &neqn, integer &neq, integer &mode, integer &
           iprint, integer &miter, integer &inform__, doublereal &bigbnd,
           doublereal &eps, doublereal &epseqn, doublereal &udelta, doublereal *
           bl, doublereal *bu, doublereal *x, doublereal *f, doublereal *g,
           integer *iw, integer &iwsize, doublereal *w, integer &nwsize, objType
           obj, constrType constr, gradobType gradob, gradcnType gradcn);

//#include "ffsqp.c"
}


namespace parpe {

extern "C" {
// FFSQP callback functions
void obj(integer &nparam, integer &j, doublereal *x, doublereal &fj);
void constr (integer &nparam, integer &j, doublereal *x, doublereal &gj);
void gradob (integer &nparam, integer &j, doublereal *x, doublereal *gradfj, doublereal *dummy);
void gradcn (integer &nparam, integer &j, doublereal *x, doublereal *gradgj, doublereal *dummy);
}

// FFSQP print level options
enum class iprintPrintLevel {
    nothing,
    firstLast,
    eachIteration,
    eachIterationDetailed
};

// FFSQP exit status
enum class informExitStatus {
    converged,
    infeasibleStartingPointByLinearContraints,
    infeasibleStartingPointByNonLinearContraints,
    maxIterationsReached,
    lineSearchStepTooSmall,
    failureQPSolverD0,
    failureQPSolverD1,
    failureInconsistentInput,
    equivalentIterates,
    bigbndExceeded
};

// make sure float sizes match; otherwise data must be copied to new containers first
static_assert(sizeof(double) == sizeof(doublereal), "");

/** Mutex for managing access to FFSQP routines which are not thread-safe */
typedef std::recursive_mutex mutexFsqpType;
static mutexFsqpType mutexFsqp {};

std::unique_lock<mutexFsqpType> fsqpGetLock()
{
    return std::unique_lock<mutexFsqpType>(mutexFsqp);
}

InverseUniqueLock<mutexFsqpType> fsqpReleaseLock()
{
    return InverseUniqueLock<mutexFsqpType>(&mutexFsqp);
}


/**
 * @brief Wrapper for a FFSQP optimization problem.
 *
 * An instance of this class is slipped into the FFSQP callback functions.
 *
 * NOTE: FFSQP is not thread-safe, therefore we lock a mutex when passing control to
 * FFSQP routines. Not sure if strictly necessary, but in FFSQP code consider
 * changing `static` to `static __thread` for non-`fmt_` variables.
 */
class FsqpProblem {
public:
    FsqpProblem(OptimizationProblem *problem)
        : problem(problem),
          reporter(problem->getReporter()),
          nparam(problem->costFun->numParameters()),
          miter(problem->getOptimizationOptions().maxOptimizerIterations),
          bl(std::vector<doublereal>(nparam)),
          bu(std::vector<doublereal>(nparam)),
          x(std::vector<doublereal>(nparam)),
          f(std::vector<doublereal>(std::max(1L, nf))),
          g(std::vector<doublereal>(std::max(1L, nineq + neq))),
          iwsize(6 * nparam + 8 * std::max(1L, nineq + neq) + 7 * std::max(1L, nf) + 30),
          iw(std::vector<integer>(iwsize))

    {
        nwsize = 4 * nparam * nparam
                + 5 * std::max(1L, nineq + neq) * nparam
                + 3 * std::max(1L, nf) * nparam
                + 26 * (nparam + std::max(1L, nf))
                + 45 * std::max(1L, nineq + neq) + 100;

        // Reserve one more to hide a pointer to this instance behind the last parameter
        // so we can use it in objective and constraint functions
        w.resize(nwsize + 1);
        // make sure it fits
        auto thisthis = this;
        assert(sizeof(doublereal) >= sizeof(&thisthis));
        memcpy(w.data(), &thisthis, 1 * sizeof(&thisthis));

        logmessage(LOGLVL_DEBUG, "w0 %p", &w.data()[0]);

        problem->fillInitialParameters(x);
        problem->fillParametersMin(bl);
        problem->fillParametersMax(bu);

    }

    std::tuple<int, double, std::vector<double> > optimize()
    {
        if(reporter)
            reporter->starting(x);

        // lock while ffsqp has control
        auto lock = fsqpGetLock();

        ffsqp_(nparam, nf, nineqn, nineq, neq, neqn, mode, iprint, miter, inform,
               bigbnd, eps, epseqn, udelta, bl.data(), bu.data(),
               x.data(), f.data(), g.data(), iw.data(), iwsize, &w[1], nwsize,
                parpe::obj, parpe::constr, parpe::gradob, parpe::gradcn);

        if(reporter)
            reporter->finished(f[0], x, inform);

        std::cout<<"Final cost "<<f[0]<<std::endl;

        return std::tuple<int, double, std::vector<double> >(
                    inform,
                    f[0],
                x);
    }

    void obj(integer &nparam, integer &j, doublereal *x, doublereal &fj) {
        gradientDummy.resize(nparam);
        reporter->evaluate(gsl::span<double const>(x, nparam), fj, gradientDummy);

        std::cout<<"np:"<<nparam<<" j:"<<j<<" x:"<<x[0]<<" fj:"<<fj<<std::endl;
    }

    // Once we want contraints: void constr (integer &nparam, integer &j, doublereal *x, doublereal &gj) {    }

    void gradob (integer &nparam, integer &j, doublereal *x, doublereal *gradfj, doublereal *dummy) {
        static_assert(sizeof(double) == sizeof(doublereal), "");

        double fvalDummy = NAN;
        reporter->evaluate(gsl::span<double const>(x, nparam), fvalDummy, gsl::span<double>(gradfj, nparam));
        reporter->iterationFinished(gsl::span<double const>(x, nparam), fvalDummy, gsl::span<double>(gradfj, nparam));
        std::cout<<"np:"<<nparam<<" j:"<<j<<" x:"<<x[0]<<" gradfj:"<<gradfj[0]<<std::endl;
    }

    // Once we want contraints: void gradcn (integer &nparam, integer &j, doublereal *x, doublereal *gradgj, doublereal *dummy) {}

    OptimizationProblem *problem = nullptr;
    std::unique_ptr<OptimizationReporter> reporter;

    /* Do gradient evaluation always during `obj` and save results,
     * so that we can fail in obj if current parameters are infeasible.
     * Otherwise ffsqp will terminate "successfully".
     * Actual caching is done in OptimizationReporter, but need to pass an array
     * so that the gradient is evaluated.
     */
    std::vector<doublereal> gradientDummy;

    // Number of optimization variables
    integer nparam = 0;
    // Number of objective functions
    integer nf = 1;
    // Number of nonlinear inequality constraints
    integer nineqn = 0;
    // Total number of inequality constraints
    integer nineq = 0;
    // Number of nonlinear equality constraints
    integer neqn = 0;
    // Total number of equality constraints
    integer neq = 0;
    // TODO enum
    integer mode = 0 + 10 * 0 + 100 * 1;
    // Print level
    integer iprint = static_cast<integer>(iprintPrintLevel::eachIteration);
    // Maximum allowed iterations
    integer miter = 0;
    // Status (see enum inform above)
    integer inform = 0;

    // Infinite bounds
    doublereal bigbnd = INFINITY;
    // Final norm requirement for the Newton direction
    doublereal eps = 1e-14;
    // Maximum allowed nonlinear equality constraint violation at optimum
    doublereal epseqn;
    // perturbation size for finite differences
    doublereal udelta;

    // Lower bounds for x
    std::vector<doublereal> bl;
    // Upper bounds for x
    std::vector<doublereal> bu;
    // Optimization parameters
    std::vector<doublereal> x;

    // Final cost
    std::vector<doublereal> f;
    // Final values of constraints
    std::vector<doublereal> g;

    // Number of elements in iw
    integer iwsize = 0;
    // Workspace
    std::vector<integer> iw;
    // Number of elements in w
    integer nwsize = 0;
    // Workspace
    std::vector<doublereal> w;
};

std::tuple<int, double, std::vector<double> > OptimizerFsqp::optimize(parpe::OptimizationProblem *problem)
{
    FsqpProblem p(problem);
    return p.optimize();
}

/**
 * @brief Get the position of ff in w (see FFSQP source)
 * @param nparam
 * @param j
 * @return nwff
 */
constexpr int getNwff(int nparam, int j) {
    return 1 + nparam*nparam /* nwhes1*/ + (nparam+1)*(nparam+1) /* nwff */ + j;
}

/**
 * @brief Get the position of nwgrf in w (see FFSQP source)
 * @param nparam
 * @param j
 * @return nwgrf
 */
constexpr int getNwgrf(int nparam, int j) {
    return getNwff(nparam, j) + (1 + 1) /*nwx TODO: nobj*/
            + 3 * (nparam + 1) /* nwdi, nwd, nwgm */
            + 1 /* nwgrg */ +  (0 * nparam + 1) /* nwgrf */;
}

/**
 * @brief Recover out data that was hidden in front of FFSQP's w. Find from provided fj.
 * @param fj
 * @param nparam
 * @param j
 * @return
 */
FsqpProblem *getProblemFromFj(doublereal &fj, integer nparam, integer j) {
    // need to go relative to fj because location of x changes; for position see "nwff"

    parpe::FsqpProblem *problem = nullptr;
    int nwff = getNwff(nparam, j);

    logmessage(LOGLVL_DEBUG, "w0 obj: %p", &fj - nwff + 1);

    memcpy(&problem, &fj - nwff + 1, sizeof(problem));

    RELEASE_ASSERT(problem, "Error retrieving FsqpProblem.");

    return problem;
}

/**
 * @brief Recover out data that was hidden in front of FFSQP's w. Find from provided gradfj.
 * @param gradfj
 * @param nparam
 * @param j
 * @return
 */
FsqpProblem *getProblemFromGradFj(doublereal *gradfj, integer nparam, integer j) {
    parpe::FsqpProblem *problem = nullptr;
    int nwgrf = getNwgrf(nparam, j);

    logmessage(LOGLVL_DEBUG, "w0 gradobj: %p", gradfj - nwgrf + 1);

    // NOTE: Will have to change that once we want to include constraints
    memcpy(&problem, gradfj - nwgrf + 1, sizeof(problem));

    RELEASE_ASSERT(problem, "Error retrieving FsqpProblem.");

    return problem;
}


/**
 * @brief Objective function to be passed to FFSQP.
 * @param nparam Length of x
 * @param j Objective index
 * @param x Parameters
 * @param fj (out) Function value for objective j
 */
void obj(integer &nparam, integer &j, doublereal *x, doublereal &fj) {
    RELEASE_ASSERT(j == 1, "Error: j > 1. Only a single objective is currently supported.");

    auto unlockFsqp = fsqpReleaseLock();

    auto problem = getProblemFromFj(fj, nparam, j);
    problem->obj(nparam, j, x, fj);
}

/**
 * @brief Constraints function to be passed to FFSQP.
 * @param nparam Length of x
 * @param j Constraint index
 * @param x Parameters
 * @param gj Value of constraint j
 */
void constr (integer &nparam, integer &j, doublereal *x, doublereal &gj) {
    //auto unlockFsqp = fsqpReleaseLock();

    // no constraints currently supported
}

/**
 * @brief Objective function gradient to be passed to FFSQP.
 * @param nparam Length of x
 * @param j Objective index
 * @param x Parameters
 * @param gradfj (out) Gradient value for objective j at x
 * @param dummy Passed to gradob for forward difference calculation
 */
void gradob (integer &nparam, integer &j, doublereal *x, doublereal *gradfj, doublereal *dummy) {
    RELEASE_ASSERT(j == 1, "Error: j > 1. Only a single objective is currently supported.");

    auto unlockFsqp = fsqpReleaseLock();

    auto problem = getProblemFromGradFj(gradfj, nparam, j);
    problem->gradob(nparam, j, x, gradfj, dummy);
}

/**
 * @brief Constraint gradient to be passed to FFSQP.
 * @param nparam Length of x
 * @param j Objective index
 * @param x Parameters
 * @param gradgj (out) Gradient value for constraint j at x
 * @param dummy Passed to gradob for forward difference calculation
 */
void gradcn (integer &nparam, integer &j, doublereal *x, doublereal *gradgj, doublereal *dummy) {
    //auto unlockFsqp = fsqpReleaseLock();

    // no constraints currently supported
}


} // namespace parpe

