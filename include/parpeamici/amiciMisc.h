#include "amici/rdata.h"
#include "parpecommon/logging.h"
#include <amici/defines.h>
#include <amici/misc.h>
#include <memory>

namespace parpe {

using amici::getUnscaledParameter;

using amici::getScaledParameter;

std::unique_ptr<amici::ReturnData>
run_amici_simulation(amici::Solver& solver,
                     amici::ExpData const* edata,
                     amici::Model& model,
                     bool rethrow = false,
                     Logger* logger = nullptr);
}
