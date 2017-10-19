#include "SteadyStateMultiConditionProblem.h"
#include "optimizationOptions.h"
#include "wrapfunctions.h"

#include <LoadBalancerWorker.h>
#include <hdf5Misc.h>
#include <logging.h>
#include <multiConditionProblemResultWriter.h>
#include <optimizationApplication.h>
#include <misc.h>
#include <iostream>
#include <unistd.h>

/** @brief This example demonstrates the use of the loadbalancer / queue for
 * parallel ODE simulation.
 * The example is based on the `steadystate` example included in AMICI.
 *
 * This model has 5 parameters, 3 states, 4 condition specific fixed parameters.
 *
 * To run, e.g.: mpiexec -np 4
 * ../parPE-build/amici/examples/steadystate/example_steadystate_multi -o
 * steadystate_`date +%F` amici/examples/steadystate/data.h5
 */

class SteadystateApplication : public parpe::OptimizationApplication {
  public:
    using OptimizationApplication::OptimizationApplication;

    virtual void initProblem(std::string inFileArgument,
                             std::string outFileArgument) override {
        model = std::unique_ptr<Model>(getModel());
        dataProvider = std::make_unique<SteadyStateMultiConditionDataProvider>(model.get(), inFileArgument);

        problem =
            new SteadyStateMultiConditionProblem(dataProvider.get(), &loadBalancer);

        parpe::JobIdentifier id;
        resultWriter =
            new parpe::MultiConditionProblemResultWriter(outFileArgument, true, id);
        problem->resultWriter = resultWriter;
    }

    virtual int runSingleMpiProcess() override;

    virtual ~SteadystateApplication() {
        delete resultWriter;
        delete problem;
    }

    std::unique_ptr<SteadyStateMultiConditionDataProvider> dataProvider;
    std::unique_ptr<Model> model;
};

class SteadystateLocalOptimizationApplication : public SteadystateApplication {
  public:
    SteadystateLocalOptimizationApplication() : SteadystateApplication() {}

    virtual int runMaster() override {
        // Single optimization
        return getLocalOptimum(problem);
    }
};

class SteadystateMultiStartOptimizationApplication
    : public SteadystateApplication {
  public:
    using SteadystateApplication::SteadystateApplication;

    virtual int runMaster() override {

        int status = 0;

        // Multistart optimization
        std::unique_ptr<parpe::OptimizationOptions> options(parpe::OptimizationOptions::fromHDF5(
                                                                 dataProvider->getHdf5FileId()));
        // if numStarts > 1: need to use multiple MPI
        // workers, otherwise simulation crashes due
        // to CVODES threading issues

        parpe::MultiConditionProblemMultiStartOptimization multiStartOptimization(
            options->numStarts, options->retryOptimization);
        multiStartOptimization.options = options.get();
        multiStartOptimization.resultWriter = problem->resultWriter;
        multiStartOptimization.dp = dataProvider.get();
        multiStartOptimization.loadBalancer = &loadBalancer;

        parpe::logmessage(parpe::LOGLVL_DEBUG, multiStartOptimization.options->toString());

        multiStartOptimization.run();

        return status;
    }
};

int main(int argc, char **argv) {
    int status = 0;

    // SteadystateLocalOptimizationApplication app(argc, argv);
    SteadystateMultiStartOptimizationApplication app(argc, argv);
    status = app.run();

    return status;
}

int SteadystateApplication::runSingleMpiProcess() {
    return getLocalOptimum(problem);
}
