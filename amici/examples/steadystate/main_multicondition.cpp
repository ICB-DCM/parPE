#include <stdio.h>
#include "steadystateProblemParallel.h"
#include <logging.h>
#include <loadBalancerWorker.h>
#include "hdf5Misc.h"
#include <unistd.h>
#include "SteadyStateMultiConditionProblem.h"
#include "wrapfunctions.h"
#include "simulationWorkerAmici.h"
#include "multiConditionProblemResultWriter.h"
#include "optimizationApplication.h"
/*
 * This example demonstrates the use of the loadbalancer / queue for parallel ODE simulation.
 *
 * To run, e.g.: mpiexec -np 4 ../parPE-build/amici/examples/steadystate/example_steadystate_multi -o steadystate_`date +%F` amici/examples/steadystate/data.h5
 */

class SteadystateApplication : public OptimizationApplication {
public:
    SteadystateApplication(int argc, char **argv) : OptimizationApplication(argc, argv) {}

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
        // TODO : move to base class; need wrapper; ParallelProblem interface?
        if(getMpiCommSize() > 1)
            loadBalancerWorkerRun(handleWorkPackage, problem);
    }

    SteadyStateMultiConditionDataProvider *dataProvider;
};

class SteadystateLocalOptimizationApplication : public SteadystateApplication {
public:
    SteadystateLocalOptimizationApplication(int argc, char **argv) : SteadystateApplication(argc, argv) {}

    virtual int runMaster() {
            // Single optimization
            return getLocalOptimum(problem);
    }
};

class SteadystateMultiStartOptimizationApplication : public SteadystateApplication {
public:
    SteadystateMultiStartOptimizationApplication(int argc, char **argv) : SteadystateApplication(argc, argv) {}

    virtual int runMaster() {
        return getLocalOptimum(problem);

        int status = 0;
        // Multistart optimization
        OptimizationOptions options;
        options.maxOptimizerIterations = 1;
        options.numStarts = 1; // if numStarts > 1: need to use multiple MPI workers, otherwise simulation crashes due to CVODES threading issues

        MultiConditionProblemGeneratorForMultiStart generator;
        generator.options = &options;
        generator.resultWriter = reinterpret_cast<MultiConditionProblemResultWriter *>(problem->resultWriter);
        generator.dp = dataProvider;

        runParallelMultiStartOptimization(&generator, options.numStarts, options.retryOptimization);

        return status;
    }

};


int main(int argc, char **argv)
{
    int status = 0;

    SteadystateLocalOptimizationApplication app(argc, argv);
    status = app.run();

    return status;
}
