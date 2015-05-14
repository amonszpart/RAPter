// Class borrowed from Nicolas Mellado<nmellad0@gmail.com>

#ifndef GF2_EXCEPTION_H
#define GF2_EXCEPTION_H

#include <string>

#define GENERATE_CLASS_EXCEPTION(ClassName)                                    \
class Exception : public Utilities::Exception {                                \
public :                                                                       \
    Exception (const std::string & msg) :                                      \
    Utilities::Exception (ClassName, msg) {}                                   \
};                                                                             \


namespace Utilities{
    //! Base class for Exceptions
    class Exception {
    private:
        std::string _className;
        std::string _message;
    protected:
        Exception (const std::string className, const std::string & msg) :
                _className(className),
                _message (msg) {}
    public :
         virtual  ~Exception () {}
        const std::string getMessage () const{
            return "[" + _className + "] : " + _message;
        }
    };
}

#endif // GF2_EXCEPTION_H
