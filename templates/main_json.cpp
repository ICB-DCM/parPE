/*
 * Executable for evaluating the objective function for parameters specified
 * in a JSON file. The results will be written to a JSON file.
 *
 * The input is a list of objects with the following fields:
 * - data_file: the name of the file containing the conditions
 * - parameters: object with parameterId => parameterValue
 * - gradient: boolean indicating whether the gradient should be computed
 *
 * The output is a list of objects with the following fields:
 * - data_file: the name of the file containing the conditions
 * - parameters: object with parameterId => parameterValue
 * - fval: the objective function value
 * - gradient: the gradient (object with parameterId => gradient entry),
 *     or `false` if not computed
 */

#include <parpecommon/parpeConfig.h>

#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/multiConditionProblem.h>
#include <parpeamici/standaloneSimulator.h>
#include <parpecommon/misc.h>
#include <parpecommon/parpeConfig.h>

#include <cstdio> // remove
#include <iostream>
#include <fstream>
#include <string_view>
#include <nlohmann/json.hpp>
#ifdef PARPE_ENABLE_MPI
#include <mpi.h>
#endif

using json = nlohmann::json;

// fields in the input/output json
constexpr std::string_view DATA_FILE = "data_file";
constexpr std::string_view PARAMETERS = "parameters";
constexpr std::string_view FVAL = "fval";
constexpr std::string_view GRADIENT = "gradient";
constexpr std::string_view STATUS = "status";
constexpr std::string_view ERROR = "error";


namespace amici::generic_model {
std::unique_ptr<amici::Model> getModel();
}


void printUsage() {
    std::cerr<<"Error: wrong number of arguments.\n";
    std::cerr<<"Usage: ... infile outfile\n";
}

int main(int argc, char **argv) {
    int status = EXIT_SUCCESS;

    if(argc != 3) {
        printUsage();
        return EXIT_FAILURE;
    }

    std::string inFileArgument = argv[1];
    std::string outFileArgument = argv[2];

    std::ifstream istream(inFileArgument);
    json jobs = json::parse(istream);

    json results = json::array();
    std::unique_ptr<parpe::MultiConditionDataProvider> dataProvider;
    std::string h5file_name;

    for (auto& job : jobs) {
        std::cout<<job<<std::endl;
        auto cur_result = job;
        try {
            if (job[DATA_FILE] != h5file_name) {
                h5file_name = job[DATA_FILE];
                dataProvider = std::make_unique<parpe::MultiConditionDataProviderHDF5>(
                    amici::generic_model::getModel(), h5file_name);
            }
            parpe::MultiConditionProblem problem(
                dataProvider.get(), nullptr,
                std::make_unique<parpe::Logger>(),
                nullptr);

                   // parameter values
            auto objective = dynamic_cast<parpe::SummedGradientFunctionGradientFunctionAdapter<int>*>(problem.cost_fun_.get());

            auto parameter_ids = objective->getParameterIds();
            std::vector<double> parameters;
            parameters.reserve(parameter_ids.size());

            for (auto& parameter_id: parameter_ids) {
                parameters.push_back(job[PARAMETERS][parameter_id]);
            }

            double fval;
            std::vector<double> gradient(parameter_ids.size());
            auto fv = objective->evaluate(parameters, fval, job[GRADIENT]?gsl::make_span(gradient):gsl::span<double>(nullptr, 0));
            cur_result[FVAL] = fval;
            if (job[GRADIENT]) {
                cur_result[GRADIENT] = json::object();
                for (size_t i = 0; i < parameter_ids.size(); ++i) {
                    cur_result[GRADIENT][parameter_ids[i]] = gradient[i];
                }
            }
            cur_result[STATUS] = fv;
        } catch(std::exception& e) {
            std::cerr<<e.what()<<std::endl;
            cur_result[ERROR] = e.what();

        }
        results.push_back(cur_result);
        // not efficient, but let's save it after every simulation for now
        std::ofstream ostream(outFileArgument);
        ostream << std::setw(4) << results << std::endl;
    }

    return status;
}
