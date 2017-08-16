#ifndef PARPEEXCEPTION_H
#define PARPEEXCEPTION_H

#include <exception>
#include <string>

class ParPEException : public std::exception {
  public:
    ParPEException(const char *message);

    ParPEException(const std::string &message);

    virtual ~ParPEException() throw() {}

    virtual const char *what() const noexcept;

  protected:
    std::string message;
};

#endif // PARPEEXCEPTION_H
