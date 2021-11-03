#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/multiConditionProblem.h>
#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpeoptimization/optimizationOptions.h>
#include <amici/model.h>
#include <amici/hdf5.h>
#include <iostream>
#include <memory>

// to avoid including model-specific header files
namespace amici::generic_model {
std::unique_ptr<amici::Model> getModel();
}
using namespace parpe;

int main(int argc, char **argv) {
#ifndef NDEBUG
    // Set stdout to unbuffered when debugging
    setbuf(stdout, NULL);
#endif

    std::string inFileArgument = "";
    std::string outFileArgument = "";

    // Check command line arguments
    if (argc == 2) {
        inFileArgument = argv[1];
    } else if (argc == 3) {
        inFileArgument = argv[1];
        outFileArgument = argv[2];
    } else {
        std::stringstream ss;
        ss << "Error: USAGE: "<< argv[0]
           << " HDF5_INPUT_FILE [HDF5_OUTPUT_FILE]\n";
        fprintf(stderr, "%s", ss.str().c_str());
        return 1;
    }


    parpe::logmessage(parpe::loglevel::info,
                      "Reading options and data from '%s'.",
                      inFileArgument.c_str());

    // setup data and problem
    MultiConditionDataProviderHDF5 dataProvider(
        amici::generic_model::getModel(), inFileArgument);
    auto options = OptimizationOptions::fromHDF5(dataProvider.getHdf5File());

    std::unique_ptr<OptimizationResultWriter> rw;
    if(!outFileArgument.empty()) {
        rw = std::make_unique<OptimizationResultWriter>(
            outFileArgument, true, "/");
    }

    MultiConditionProblem problem { &dataProvider, nullptr, nullptr,
                                  std::move(rw)};

    // Read nominal parameters
    auto optimizationParams = amici::hdf5::getDoubleDataset1D(
                dataProvider.getHdf5File(), "/parameters/nominalValues");

    double fval = NAN;
    std::vector<double> gradient(optimizationParams.size(), NAN);
    problem.cost_fun_->evaluate(optimizationParams, fval, gradient);

    std::cout<<gradient<<std::endl;
    std::for_each(gradient.begin(), gradient.end(), [](double &d){ d = std::fabs(d); });
    auto maxAbsGradient = *std::max_element(gradient.begin(), gradient.end());
    std::for_each(gradient.begin(), gradient.end(), [fval](double &d){ d /= std::fabs(fval); });
    auto maxRelGradient = *std::max_element(gradient.begin(), gradient.end());
    parpe::logmessage(loglevel::info, "Max(abs(grad)) = " + std::to_string(maxAbsGradient));
    parpe::logmessage(loglevel::info, "Max(abs(grad)/fval) = " + std::to_string(maxRelGradient));

    parpe::logmessage(loglevel::info, "Likelihood: " + std::to_string(fval));
}
