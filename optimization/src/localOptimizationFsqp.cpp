#include "localOptimizationFsqp.h"

#include <cmath>
#include <iostream>
#include <assert.h>
#include <cstring>
#include <stdexcept>
#include <optimizationResultWriter.h>


extern "C" {
#include <f2c.h>
#undef max

// callback functions types for FFSQP, see below
using objType    = void (*) (integer &nparam, integer &j, doublereal *x, doublereal &fj);
using constrType = void (*) (integer &nparam, integer &j, doublereal *x, doublereal &gj);
using gradobType = void (*) (integer &nparam, integer &j, doublereal *x, doublereal *gradfj, doublereal *dummy);
using gradcnType = void (*) (integer &nparam, integer &j, doublereal *x, doublereal *gradgj, doublereal *dummy);


// FFSQP callback functions
void obj(integer &nparam, integer &j, doublereal *x, doublereal &fj);
void constr (integer &nparam, integer &j, doublereal *x, doublereal &gj);
void gradob (integer &nparam, integer &j, doublereal *x, doublereal *gradfj, doublereal *dummy);
void gradcn (integer &nparam, integer &j, doublereal *x, doublereal *gradgj, doublereal *dummy);

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

enum class iprint { nothing, firstLast, eachIteration, eachIterationDetailed };
enum class inform {
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

namespace parpe {

class FsqpProblem {
public:
    FsqpProblem(OptimizationProblem *problem)
        : problem(problem),
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

        // std::cerr<<"w0 "<<&w.data()[0]<<std::endl;

        problem->fillInitialParameters(x.data());
        problem->fillParametersMin(bl.data());
        problem->fillParametersMax(bu.data());

    }

    std::tuple<int, double, std::vector<double> > optimize()
    {
        if(reporter)
            reporter->starting(nparam, x.data());

        ffsqp_(nparam, nf, nineqn, nineq, neq, neqn, mode, iprint, miter, inform,
               bigbnd, eps, epseqn, udelta, bl.data(), bu.data(),
               x.data(), f.data(), g.data(), iw.data(), iwsize, &w[1], nwsize,
                obj, constr, gradob, gradcn);

        if(reporter)
            reporter->finished(f[0], x.data(), inform);

        std::cout<<"Final cost "<<f[0]<<std::endl;

        return std::tuple<int, double, std::vector<double> >(
                    inform,
                    f[0],
                x);

    }

    OptimizationProblem *problem = nullptr;
    std::unique_ptr<OptimizationReporter> reporter;

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
    integer iprint = static_cast<integer>(iprint::eachIteration);
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

} // namespace parpe


/**
 * @brief Objective function to be passed to FFSQP.
 * @param nparam Length of x
 * @param j Objective index
 * @param x Parameters
 * @param fj (out) Function value for objective j
 */
void obj(integer &nparam, integer &j, doublereal *x, doublereal &fj) {
    // Recover user-data:
    parpe::FsqpProblem *problem = nullptr;
    // need to go relative to fj because location of x changes; for position see "nwff"
    int nwff = 1 + nparam*nparam + (nparam+1)*(nparam+1) + j;
    // std::cerr<<"obj "<<&fj - nwff + 1<<std::endl;
    memcpy(&problem, &fj - nwff + 1, sizeof(problem));

    if(problem->reporter && problem->reporter->beforeCostFunctionCall(nparam, x) != 0)
        return;

    problem->problem->costFun->evaluate(x, fj, nullptr);

    if(problem->reporter && problem->reporter->afterCostFunctionCall(nparam, x, fj, nullptr))
        return;

    std::cout<<"np:"<<nparam<<" j:"<<j<<" x:"<<x[0]<<" fj:"<<fj<<std::endl;
}

/**
 * @brief Constraints function to be passed to FFSQP.
 * @param nparam Length of x
 * @param j Constraint index
 * @param x Parameters
 * @param gj Value of constraint j
 */
void constr (integer &nparam, integer &j, doublereal *x, doublereal &gj) {
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
    parpe::FsqpProblem *problem = nullptr;
    // See nwgrf
    int nwff = 1 + nparam*nparam + (nparam+1)*(nparam+1) + j;
    // NOTE: Will have to change that once we want to include constraints
    int nwgrf = nwff + 1 + 1 + 3 * (nparam + 1) + 1 + 1 * nparam + 1 - nparam; // Where does the last "- nparam" come from?
    //std::cerr<<"gradfj "<<gradfj-nwgrf+1<<std::endl;
    memcpy(&problem, gradfj-nwgrf+1, sizeof(problem));

    if(problem->reporter && problem->reporter->beforeCostFunctionCall(nparam, x) != 0)
        return;

    static_assert(sizeof(double) == sizeof(doublereal), "");
    doublereal fj;
    problem->problem->costFun->evaluate(x, fj, gradfj);

    if(problem->reporter && problem->reporter->afterCostFunctionCall(nparam, x, fj, gradfj))
        return;

    if(problem->reporter && problem->reporter->iterationFinished(nparam, x, fj, gradfj))
        return;


    std::cout<<"np:"<<nparam<<" j:"<<j<<" x:"<<x[0]<<" gradfj:"<<gradfj<<std::endl;
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
    // no constraints currently supported
}

namespace parpe {

std::tuple<int, double, std::vector<double> > OptimizerFsqp::optimize(parpe::OptimizationProblem *problem)
{
    FsqpProblem p(problem);
    return p.optimize();
}

}
