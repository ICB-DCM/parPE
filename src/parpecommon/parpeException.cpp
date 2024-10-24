#include <parpecommon/parpeException.h>

namespace parpe {

ParPEException::ParPEException(char const* message)
    : message(message) {}

ParPEException::ParPEException(std::string message)
    : message(std::move(message)) {}

char const* ParPEException::what() const noexcept { return message.c_str(); }

} // namespace parpe
