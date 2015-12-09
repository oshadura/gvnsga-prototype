#ifndef __EXCEPTIONMESSENGER__
#define __EXCEPTIONMESSENGER__

#include <exception>

class ExceptionMessenger : public std::exception
{
public:
	explicit ExceptionMessenger(const std::string& msg): std::exception(), message(msg){};
	~ExceptionMessenger() throw () {};
	const char* ReturnBody() const throw(){
		return message.c_str();
	};
private:
	std::string message;
};

#endif