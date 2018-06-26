#include "parpeException.h"

namespace parpe {

ParPEException::ParPEException(const char *message) : message(message) {}

ParPEException::ParPEException(std::string message) : message(std::move(message)) {}

const char *ParPEException::what() const noexcept { return message.c_str(); }

} // namespace parpe
