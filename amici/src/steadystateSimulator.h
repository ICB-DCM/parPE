#ifndef STEADYSTATESIMULATOR_H
#define STEADYSTATESIMULATOR_H

#include <udata.h>
#include <edata.h>
#include <rdata.h>


class SteadystateSimulator
{
public:

    /**
     * @brief SteadystateSimulator::getSteadystateSolution
     * @param udata Note: udata->x0 will be overwritten.
     * @param edata
     * @param status
     * @param iterationDone
     * @return
     */

    static ReturnData *getSteadystateSolution(UserData *udata, ExpData *edata, int *status, int *iterationDone);


    /**
     * @brief getSteadystateSolution Simulate the model until steady state
     * @param udata: Model data and options
     * @param edata: ExpData TODO: could remove in case of no events?
     * @param status: return status, 0 if successful
     * @param iterationDone: Iterations of simulation to t=10^9 until steady-state is reached
     * @return
     */
    static bool reachedSteadyState(const double *xdot, const double *x, int numTimepoints, int numStates, double tolerance);

protected:
    static void updateInitialConditions(double destination[], const double src[], int count);

};

#endif // STEADYSTATESIMULATOR_H
