#ifndef __EXCEPTIONMESSENGER__
#define __EXCEPTIONMESSENGER__

#include <exception>
#include "TObject.h"


class ExceptionMessenger : public std::exception {
public:
  explicit ExceptionMessenger(const std::string &msg)
      : std::exception(), message(msg){};
  ~ExceptionMessenger() throw(){};
  const char *what() const throw() { return message.c_str(); };

private:
  std::string message;
  ClassDef(ExceptionMessenger, 1)
};


#endif
