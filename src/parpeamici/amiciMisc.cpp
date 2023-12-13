#include <parpeamici/amiciMisc.h>
#include <amici/amici.h>
#include <amici/model.h>
#include <amici/edata.h>
#include <amici/logging.h>
#include <parpecommon/logging.h>

#include <memory>

namespace parpe {

using amici::ReturnData;

std::unique_ptr<ReturnData> run_amici_simulation(
    amici::Solver &solver,
    const amici::ExpData *edata,
    amici::Model &model,
    bool rethrow,
    Logger *logger) {

    auto rdata = amici::runAmiciSimulation(solver, edata, model, rethrow);

    if (logger != nullptr) {
        // TODO: subclass amici::Logger to print messages without delay

        // for now, print collected messages after simulation
        for(auto const& log_item: rdata->messages) {
            auto lvl = loglevel::debug;
            switch (log_item.severity) {
            case amici::LogSeverity::debug:
                lvl = loglevel::debug;
                break;
            case amici::LogSeverity::warning:
                lvl = loglevel::warning;
                break;
            case amici::LogSeverity::error:
                lvl = loglevel::error;
                break;
            }
            if(!log_item.identifier.empty()) {
                logger->logmessage(lvl, "[" + log_item.identifier + "] " + log_item.message);
            } else {
                logger->logmessage(loglevel::warning, log_item.message);
            }

        }
    }
    return rdata;
}

} // namespace parpe

