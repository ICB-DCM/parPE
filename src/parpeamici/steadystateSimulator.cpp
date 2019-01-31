#include <parpeamici/steadystateSimulator.h>

#include <parpecommon/logging.h>
#include <parpecommon/parpeException.h>

#include <cassert>
#include <cmath>
#include <cstring>

namespace parpe {

#define XDOT_REL_TOLERANCE 1e-6

std::unique_ptr<amici::ReturnData> SteadystateSimulator::getSteadystateSolution(amici::Model &model,
                                                         amici::Solver &solver,
                                                         amici::ExpData *edata,
                                                         int *status,
                                                         int *iterationDone) {
    if (model.nt() > 1)
        throw(ParPEException("SteadystateSimulator::getSteadystateSolution "
                             "works only with nt == 1"));

    std::unique_ptr<amici::ReturnData> rdata;
    bool inSteadyState = FALSE;
    int iterations = 0;

    while (!inSteadyState) {
        ++iterations;

        rdata = amici::runAmiciSimulation(solver, edata, model);

        if (rdata->status < 0) {
            error("Failed to integrate."); // TODO add dataset info,
                                           // case/control, celline
            return rdata;
        }

        inSteadyState = reachedSteadyState(rdata->xdot.data(), rdata->x.data(), rdata->nt,
                                           rdata->nx, XDOT_REL_TOLERANCE);

        if (inSteadyState) {
            break;
        }

        if (iterations >= 100) {
            logmessage(LOGLVL_WARNING, "getSteadystateSolutionForExperiment: "
                                       "no steady after %d iterations... "
                                       "aborting...",
                       iterations);
            *status = -1;
            break;
        }

        if (iterations % 10 == 0) {
            logmessage(LOGLVL_DEBUG, "getSteadystateSolutionForExperiment: no "
                                     "steady state after %d iterations... "
                                     "trying on...",
                       iterations);
        }

        // use previous solution as initial conditions
        model.setInitialStates(rdata->x);
    }
    // logmessage(LOGLVL_DEBUG, "getSteadystateSolutionForExperiment:
    // steadystate after %d iterations", iterations);

    *iterationDone = iterations;

    return rdata;
}


bool SteadystateSimulator::reachedSteadyState(const double *xdot,
                                              const double *x,
                                              int numTimepoints, int numStates,
                                              double tolerance) {
    assert(numTimepoints == 1 &&
           "SteadystateSimulator currently only supports one timepoint!");
    for (int state = 0; state < numStates; ++state) {
        double sensitivity = fabs(xdot[state]) / (fabs(x[state]) + tolerance);
        if (sensitivity > tolerance) {
            // logmessage(LOGLVL_DEBUG, "No steady state: %d: x %e xdot %e
            // relxdot %e s %e\n", state,  (x[state]), (xdot[state]) ,
            // (xdot[state]) / (x[state]), sensitivity);
            return FALSE;
        }
    }
    return TRUE;
}

} // namespace parpe
