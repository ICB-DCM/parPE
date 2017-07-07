#include "wrapfunctions.h"
#include "SteadyStateMultiConditionProblem.h"
#include "optimizationOptions.h"

#include <logging.h>
#include <loadBalancerWorker.h>
#include <hdf5Misc.h>
#include <multiConditionProblemResultWriter.h>
#include <optimizationApplication.h>

#include <stdio.h>
#include <unistd.h>

/*
 * This example demonstrates the use of the loadbalancer / queue for parallel ODE simulation.
 *
 * To run, e.g.: mpiexec -np 4 ../parPE-build/amici/examples/steadystate/example_steadystate_multi -o steadystate_`date +%F` amici/examples/steadystate/data.h5
 */

class SteadystateApplication : public OptimizationApplication {
public:
    SteadystateApplication() : OptimizationApplication() {}

    virtual void initProblem(const char *inFileArgument, const char *outFileArgument) {
        dataProvider = new SteadyStateMultiConditionDataProvider(inFileArgument);
        problem = new SteadyStateMultiConditionProblem(dataProvider);

        JobIdentifier id = {0};
        resultWriter = new MultiConditionProblemResultWriter(outFileArgument, true, id);
        problem->resultWriter = resultWriter;
    }

    virtual void destroyProblem() {
        delete resultWriter;
        delete problem;
        delete dataProvider;
    }

    virtual void runWorker() {
        if(getMpiCommSize() > 1)
            loadBalancerWorkerRun(handleWorkPackage, problem);
    }

    SteadyStateMultiConditionDataProvider *dataProvider;
};

class SteadystateLocalOptimizationApplication : public SteadystateApplication {
public:
    SteadystateLocalOptimizationApplication() : SteadystateApplication() {}

    virtual int runMaster() {
            // Single optimization
            return getLocalOptimum(problem);
    }
};

class SteadystateMultiStartOptimizationApplication : public SteadystateApplication {
public:
    SteadystateMultiStartOptimizationApplication() : SteadystateApplication() {}

    virtual int runMaster() {

        int status = 0;

        // Multistart optimization
        MultiConditionProblemGeneratorForMultiStart generator;
        generator.options = OptimizationOptions::fromHDF5(dataProvider->fileId); // if numStarts > 1: need to use multiple MPI workers, otherwise simulation crashes due to CVODES threading issues
        generator.resultWriter = reinterpret_cast<MultiConditionProblemResultWriter *>(problem->resultWriter);
        generator.dp = dataProvider;

        std::cout<<generator.options->toString();

        runParallelMultiStartOptimization(&generator, generator.options->numStarts, generator.options->retryOptimization);

        return status;
    }

};


int main(int argc, char **argv)
{
    int status = 0;

    //SteadystateLocalOptimizationApplication app(argc, argv);
    SteadystateMultiStartOptimizationApplication app;
    app.init(argc, argv);
    status = app.run();

    return status;
}
