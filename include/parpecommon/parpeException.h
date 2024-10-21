#ifndef PARPEEXCEPTION_H
#define PARPEEXCEPTION_H

#include <exception>
#include <string>

namespace parpe {

class ParPEException : public std::exception {
  public:
    explicit ParPEException(char const* message);

    explicit ParPEException(std::string message);

    ~ParPEException() throw() override = default;

    char const* what() const noexcept override;

  private:
    std::string message;
};

} // namespace parpe
#endif // PARPEEXCEPTION_H
