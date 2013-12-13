/*
 * Version.h
 *
 *  Created on: April 17, 2012
 *      Author: aidan
 */

#ifndef __gpgmp_Version_h__
#define __gpgmp_Version_h__


#include <string>
#include <sstream>


namespace gpgmp {
    
    
class Version {
    
public:
    static const int major    = 1;
    static const int minor    = 0;
    static const int patch    = 0;
    static const int revision; // set by Version.cpp (which is generated at build)

    inline static std::string Version::toString() {
        std::ostringstream ss;
        ss << major << '.' << minor << '.' << patch << '.' << revision;
        return ss.str();
    }
};

} // namespace gpgmp


#endif // !__gpgmp_Version_h__
