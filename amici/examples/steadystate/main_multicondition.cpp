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

        problem = std::unique_ptr<parpe::MultiConditionProblem>(
                    new SteadyStateMultiConditionProblem(dataProvider.get(), &loadBalancer));

        parpe::JobIdentifier id;
        resultWriter = std::make_unique<parpe::MultiConditionProblemResultWriter>(outFileArgument, true, id);

        problem->resultWriter = std::make_unique<parpe::MultiConditionProblemResultWriter>(*resultWriter);
        problem->resultWriter->setJobId(id);

    }

    virtual int runSingleMpiProcess() override {
        parpe::MultiConditionProblemMultiStartOptimization ms(
            1,
            problem->getOptimizationOptions().retryOptimization);
        ms.options = problem->getOptimizationOptions();
        ms.resultWriter = problem->resultWriter.get();
        ms.dp = problem->getDataProvider();
        ms.loadBalancer = &loadBalancer;
        ms.run();
        return 0;
        //return getLocalOptimum(problem);
    }

    std::unique_ptr<SteadyStateMultiConditionDataProvider> dataProvider;
    std::unique_ptr<Model> model;
};

class SteadystateLocalOptimizationApplication : public SteadystateApplication {
  public:

    using SteadystateApplication::SteadystateApplication;

    virtual int runMaster() override {
        // Single optimization
        return getLocalOptimum(problem.get());
    }
};

int main(int argc, char **argv) {
    int status = 0;

//     SteadystateLocalOptimizationApplication app(argc, argv);
    SteadystateApplication app(argc, argv);
    status = app.run();

    return status;
}


