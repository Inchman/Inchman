/*
 * Base.h
 *
 *  Created on: 08/09/2010
 *      Author: aidan
 */

#ifndef __gpgmp_Base_h__
#define __gpgmp_Base_h__


#include <Globals.h>

#include <string>


namespace gpgmp {
    
/**
 * Base class for all Inchman classes.
 */
class Base {
    
public:
    /**
     * Constructs the object with the given id.
     * @param id Id of the new object.
     */
    Base(const std::string &id);
    
    /**
     * Destructor.
     */
    virtual ~Base();
    
    /**
     * @return Gives the ID
     */
    inline const std::string &id() const;
    
    
private:
    std::string m_id; ///< Object id
};
    
    
inline const std::string & Base::id() const { return m_id; }


} // namespace gpgmp


#endif // !__gpgmp_Base_h__
