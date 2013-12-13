/*
 *  InchmanDiffusionModel.h
 *  inchman
 *
 *  Created by Aidan Lane on Feb 1, 2011.
 *  Copyright 2011-2012. All rights reserved.
 *
 */

#ifndef __InchmanDiffusionModel_h__
#define __InchmanDiffusionModel_h__


#include <boost/python.hpp>

#include <DiffusionModel.h>
#include <sbml/SBMLTypes.h>

namespace LIBSBML_CPP_NAMESPACE{
class Parameter;
class XMLNode;
class Compartment;
class Reaction;
class Species;
}

namespace gpgmp {

class InchmanDiffusionModel : public DiffusionModel {

    typedef DiffusionModel Super;
    typedef InchmanDiffusionModel Self;
    
public:
    InchmanDiffusionModel();
    
    static InchmanDiffusionModel *
    fromFile(const std::string &filename,
             const std::map<std::string, double> &paramOverrides);

    bool load(const std::string &filename,
              const std::map<std::string, double> &paramOverrides);
    
    
protected:
    virtual void createOutputFile();
    
    using Super::setParameter;
    void setParameter(LIBSBML_CPP_NAMESPACE::Parameter *sbmlParameter);
    
    using Super::addCompartment;
    gpgmp::Compartment *addCompartment(LIBSBML_CPP_NAMESPACE::Compartment *sbmlCompartment);
    
    using Super::addSpecies;
    gpgmp::Species *addSpecies(LIBSBML_CPP_NAMESPACE::Species *sbmlSpecies);
    
    using Super::addReaction;
    gpgmp::Reaction *addReaction(LIBSBML_CPP_NAMESPACE::Reaction *sbmlReaction);

    
private:
    bool tryEvalChildContentsString(LIBSBML_CPP_NAMESPACE::XMLNode *node, const std::string &childName, std::string *out);
    bool tryEvalNodeContentsString(LIBSBML_CPP_NAMESPACE::XMLNode *node, std::string *out);
    
    template <class T>
    bool tryEvalChildContentsMathML(LIBSBML_CPP_NAMESPACE::XMLNode *node, const std::string &childName, T *out);
    
    template <class T>
    bool tryEvalNodeContentsMathML(LIBSBML_CPP_NAMESPACE::XMLNode *node, T *out);
    
    template <class C, class R, class T>
    bool tryEvalChildContentsMathML(LIBSBML_CPP_NAMESPACE::XMLNode *node, const std::string &childName, C *outObj, R (C::*outSetter)(T));
    
    template <class C, class R, class T>
    bool tryEvalNodeContentsMathML(LIBSBML_CPP_NAMESPACE::XMLNode *node, C *outObj, R (C::*outSetter)(T));
    
    template <class C, class R, class T>
    bool tryGetNodeContents(LIBSBML_CPP_NAMESPACE::XMLNode *node, const std::string &childName, C *outObj, R (C::*outSetter)(T));
    
    template <class C, class R, class T>
    bool tryGetNodeContents(LIBSBML_CPP_NAMESPACE::XMLNode *node, C *outObj, R (C::*outSetter)(T));

    std::string c_inchmanFile;
    
    boost::python::object c_pyMainModule;
    boost::python::object c_pyMainNamespace;
};

    
} // namespace gpgmp

#endif // !__InchmanDiffusionModel_h__w
