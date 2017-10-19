#ifndef PARPEEXCEPTION_H
#define PARPEEXCEPTION_H

#include <exception>
#include <string>

namespace parPE {

class ParPEException : public std::exception {
  public:
    ParPEException(const char *message);

    ParPEException(const std::string &message);

    virtual ~ParPEException() throw() {}

    virtual const char *what() const noexcept override;

  protected:
    std::string message;
};

} // namespace parPE
#endif // PARPEEXCEPTION_H
