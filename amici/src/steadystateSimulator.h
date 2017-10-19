#ifndef STEADYSTATESIMULATOR_H
#define STEADYSTATESIMULATOR_H

#include <edata.h>
#include <rdata.h>
#include <udata.h>

class Model;

namespace parpe {


/**
 * @brief The SteadystateSimulator class runs an AMICI simulation until a
 * steady-state is reached.
 */
class SteadystateSimulator {
  public:
    /**
     * @brief SteadystateSimulator::getSteadystateSolution
     * @param udata Note: Content of udata->x0 will be overwritten.
     * @param edata
     * @param status
     * @param iterationDone
     * @return
     */

    static ReturnData *getSteadystateSolution(Model *model, UserData *udata,
                                              ExpData *edata, int *status,
                                              int *iterationDone);

    /**
     * @brief getSteadystateSolution Simulate the model until steady state
     * @param udata: Model data and options
     * @param edata: ExpData TODO: could remove in case of no events?
     * @param status: return status, 0 if successful
     * @param iterationDone: Iterations of simulation to t=10^9 until
     * steady-state is reached
     * @return
     */
    static bool reachedSteadyState(const double *xdot, const double *x,
                                   int numTimepoints, int numStates,
                                   double tolerance);

  protected:
    /**
     * @brief Replace old UserData::x0 by result of last simulation
     * @param destination
     * @param src
     * @param count
     */
    static void updateInitialConditions(double destination[],
                                        const double src[], int count);
};

} // namespace parpe

#endif // STEADYSTATESIMULATOR_H
