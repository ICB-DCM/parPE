#include "steadyStateMultiConditionDataprovider.h"
#include "wrapfunctions.h"

#ifdef PARPE_ENABLE_MPI
#include <parpeloadbalancer/loadBalancerWorker.h>
#endif

#include <parpecommon/parpeConfig.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpecommon/hdf5Misc.h>
#include <parpecommon/logging.h>
#include <parpeamici/optimizationApplication.h>
#include <parpecommon/misc.h>

#include <iostream>
#include <unistd.h>

/** @file
 *
 * This example demonstrates the use of the loadbalancer / queue for
 * parallel ODE simulation.
 * The example is based on the `steadystate` example included in AMICI.
 *
 * This model has 5 parameters, 3 states, 4 condition specific fixed parameters.
 *
 * To run, e.g.: mpiexec -np 4
 * ../parPE-build/amici/examples/steadystate/example_steadystate_multi -o
 * steadystate_`date +%F` amici/examples/steadystate/data.h5
 */

/**
 * @brief The SteadystateApplication class subclasses parpe::OptimizationApplication
 * which provides a frame for a standalone program to solve a multi-start local optimization
 * problem.
 */

class SteadystateApplication : public parpe::OptimizationApplication {
public:
    using OptimizationApplication::OptimizationApplication;

    ~SteadystateApplication() override = default;

    void initProblem(std::string const& inFileArgument,
                     std::string const& outFileArgument) override
    {

        // The same file should only be opened/created once, an then only be reopened
        h5File = parpe::hdf5CreateFile(outFileArgument, true);
        logParPEVersion(h5File);

        dataProvider = std::make_unique<SteadyStateMultiConditionDataProvider>(
                    amici::generic_model::getModel(), inFileArgument);

        // read options from file
        auto optimizationOptions = parpe::OptimizationOptions::fromHDF5(
                    dataProvider->getHdf5File());

        // Create one instance for the problem, one for the application for clear ownership
        auto multiCondProb = new parpe::MultiConditionProblem(
                    dataProvider.get(),
                    &loadBalancer,
                    std::make_unique<parpe::Logger>(),
                    // TODO remove this resultwriter
                    std::make_unique<parpe::OptimizationResultWriter>(
                        h5File,
                        std::string("/multistarts/"))
                    );

        // If hierarchical optimization was requested, wrap the original problem
        if(optimizationOptions->hierarchicalOptimization) {
            problem.reset(new parpe::HierarchicalOptimizationProblemWrapper(
                              std::unique_ptr<parpe::MultiConditionProblem>(multiCondProb),
                              dataProvider.get())
                          );
        } else {
            problem.reset(multiCondProb);
        }

        problem->setOptimizationOptions(*optimizationOptions);

        // On master, copy input data to result file
        if(parpe::getMpiRank() < 1)
            dataProvider->copyInputData(h5File);

        // TODO: we can set the correct start?
        auto ms = new parpe::MultiConditionProblemMultiStartOptimizationProblem(
                    dataProvider.get(),
                    problem->getOptimizationOptions(),
                    multiCondProb->getResultWriter(),
                    &loadBalancer,
                    std::make_unique<parpe::Logger>()
                    );
        multiStartOptimizationProblem.reset(ms);
    }

    std::unique_ptr<SteadyStateMultiConditionDataProvider> dataProvider;
};

/**
 * @brief The SteadystateLocalOptimizationApplication class overrides the multi-start optimization
 * in the base class and performs only a single optimization run. This is mostly for debugging.
 */
class SteadystateLocalOptimizationApplication : public SteadystateApplication {
public:

    using SteadystateApplication::SteadystateApplication;

    void runMaster() override {
        // Single optimization
        getLocalOptimum(problem.get());
    }
};


int main(int argc, char **argv) {
    int status = EXIT_SUCCESS;

    // SteadystateLocalOptimizationApplication app(argc, argv);
    SteadystateApplication app;
    status = app.run(argc, argv);

    return status;
}


