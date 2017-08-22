#include "steadystateSimulator.h"
#include <amici_interface_cpp.h>
#include <cassert>
#include <cmath>
#include <cstring>
#include <logging.h>
#include <mpi.h>
#include <parpeException.h>

#define XDOT_REL_TOLERANCE 1e-6

ReturnData *SteadystateSimulator::getSteadystateSolution(Model *model,
                                                         UserData *udata,
                                                         ExpData *edata,
                                                         int *status,
                                                         int *iterationDone) {
    if (udata->nt > 1)
        throw(ParPEException("SteadystateSimulator::getSteadystateSolution "
                             "works only with nt == 1"));

    ReturnData *rdata = NULL;
    bool inSteadyState = FALSE;
    int iterations = 0;

    while (!inSteadyState) {
        ++iterations;

        rdata = getSimulationResults(model, udata, edata);

        if (*rdata->status < 0) {
            error("Failed to integrate."); // TODO add dataset info,
                                           // case/control, celline
            return rdata;
        }

        inSteadyState = reachedSteadyState(rdata->xdot, rdata->x, rdata->nt,
                                           rdata->nx, XDOT_REL_TOLERANCE);

        if (inSteadyState) {
            break;
        } else if (iterations >= 100) {
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
        updateInitialConditions(udata->x0data, rdata->x, rdata->nx);

        delete rdata;
    }
    // logmessage(LOGLVL_DEBUG, "getSteadystateSolutionForExperiment:
    // steadystate after %d iterations", iterations);

    *iterationDone = iterations;

    return rdata;
}

void SteadystateSimulator::updateInitialConditions(double destination[],
                                                   const double src[],
                                                   int count) {
    memcpy(destination, src, count * sizeof(double));
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
