/*
 *  InchmanDiffusionModel.cpp
 *  inchman
 *
 *  Created by Aidan Lane on Feb 1, 2011.
 *  Copyright 2011-2012. All rights reserved.
 *
 */

#include "InchmanDiffusionModel.h"

#include <Compartment.h>
#include <Reaction.h>
#include <Species.h>

#include <sbml/SBMLTypes.h>
#include <sbml/math/ASTNode.h>
#include <sbml/xml/XMLNode.h>
#include <sbml/common/libsbml-namespace.h>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <functional>

#include <hdf5_hl.h>

#include <boost/log/trivial.hpp>

using namespace std;
using namespace gpgmp;

namespace gpgmp {

LIBSBML_CPP_NAMESPACE::XMLNode *findChild(LIBSBML_CPP_NAMESPACE::XMLNode *node, const string &name);
LIBSBML_CPP_NAMESPACE::XMLNode *findAnnotation(LIBSBML_CPP_NAMESPACE::SBase *sb, const string &name);
LIBSBML_CPP_NAMESPACE::XMLNode *findAnnotationChild(LIBSBML_CPP_NAMESPACE::SBase *sb, const string &annotationName, const string &childName);


template <class T>
bool tryReadChildContentsValue(LIBSBML_CPP_NAMESPACE::XMLNode *node, const string &childName, T &out) {
    try {
        LIBSBML_CPP_NAMESPACE::XMLNode *c = findChild(node, childName);
        if (c && c->getNumChildren() > 0)
            out = boost::lexical_cast<T>(c->getChild(0).getCharacters());
        return true;
    } catch (bad_cast &e) {
        // leave "out" as it was
        return false;
    }
}

InchmanDiffusionModel::InchmanDiffusionModel()
/*   Can't null-ify these here, as the Python instance needs to exist first
      - there not just pointers!
:    c_pyMainModule(0),
     c_pyMainNamespace(0)
*/
{}
    
InchmanDiffusionModel *
InchmanDiffusionModel::fromFile(const string &filename,
                                const map<string, double> &paramOverrides) {
    InchmanDiffusionModel *model = new InchmanDiffusionModel();
    if (!model->load(filename, paramOverrides)) {
        delete model;
        model = 0;
    }
    return model;
}

bool InchmanDiffusionModel::load(const string &filename,
                                 const map<string, double> &paramOverrides)
{
    BOOST_LOG_TRIVIAL(info) << "Loading file: " << filename;
    
    loadTextFile(c_inchmanFile, filename);
    
    // Read the annotated SBML file
    LIBSBML_CPP_NAMESPACE::SBMLDocument *document = LIBSBML_CPP_NAMESPACE::SBMLReader().readSBMLFromString(c_inchmanFile);
    if (document->getNumErrors() > 0) {
        document->printErrors(cerr);
        delete document;
        return false;
    }
    
    LIBSBML_CPP_NAMESPACE::Model *sbmlModel = document->getModel();

    // Note: DiffusionModel manages the python instance
    c_pyMainModule = boost::python::import("__main__");
    c_pyMainNamespace = c_pyMainModule.attr("__dict__");
    
    // TODO: this is an experimental approach to allowing "true" and "false" to be used in MathML expressions (e.g. compartment isAllowed)
    c_pyMainNamespace["false"] = 0.0;
    c_pyMainNamespace["true"]  = 1.0;
        
    // Parameters - DO THIS FIRST, so that all MathML expressions can refer to them
    BOOST_LOG_TRIVIAL(debug) << "Parameters from model:";
    if (sbmlModel->getNumParameters() > 0) {
        for (unsigned i=0; i < sbmlModel->getNumParameters(); i++)
            setParameter(sbmlModel->getParameter(i));
    } else {
        BOOST_LOG_TRIVIAL(debug)  << "  (none)";
    }
    
    // Parameter overrides
    BOOST_LOG_TRIVIAL(debug) << "Parameters from config:";
    if (paramOverrides.size() > 0) {
        typedef pair<string, double> ParamPair;
        BOOST_FOREACH(const ParamPair &p, paramOverrides) {
            setParameter(p.first, p.second);
        }
    } else {
        BOOST_LOG_TRIVIAL(debug) << "  (none)";
    }
    
    // SBML
    LIBSBML_CPP_NAMESPACE::XMLNode *sbml_a = findAnnotation(document, "sbml");
    if (sbml_a) {
        
        tryGetNodeContents<Super, void, const std::string &>(sbml_a, "solver", (Super *)this, &Super::setSolver);
        tryEvalChildContentsMathML(sbml_a, "runs",  (Super *)this, &Super::setNumRuns);
        tryEvalChildContentsMathML(sbml_a, "steps", (Super *)this, &Super::setMaxSteps);
        tryEvalChildContentsMathML(sbml_a, "time",  (Super *)this, &Super::setMaxTime);
        tryEvalChildContentsMathML(sbml_a, "outputInterval", (Super *)this, &Super::setOutputInterval);
        
        // Get init script
        LIBSBML_CPP_NAMESPACE::XMLNode *init_n = findChild(sbml_a, "init");
        if (init_n) {
            // TODO: ensure that the type is "python"
            LIBSBML_CPP_NAMESPACE::XMLNode *script_n = findChild(init_n, "script");
            if (script_n
                && script_n->getNumChildren() == 1
                && script_n->getChild(0).isText()) {
                setInitScriptContents(script_n->getChild(0).getCharacters());
            }
        }
        
        // Get events script
        LIBSBML_CPP_NAMESPACE::XMLNode *events_n = findChild(sbml_a, "events");
        if (events_n) {
            // TODO: ensure that the type is "python"
            LIBSBML_CPP_NAMESPACE::XMLNode *script_n = findChild(events_n, "script");
            if (script_n
                && script_n->getNumChildren() == 1
                && script_n->getChild(0).isText()) {
                setEventsScriptContents(script_n->getChild(0).getCharacters());
            }
            
            // Get the time stamps
            LIBSBML_CPP_NAMESPACE::XMLNode *listOfTimeStamps_n = findChild(events_n, "listOfTimeStamps");
            if (listOfTimeStamps_n) {
                for (unsigned int i=0; i < listOfTimeStamps_n->getNumChildren(); i++) {
                    LIBSBML_CPP_NAMESPACE::XMLNode &timeStampNode = listOfTimeStamps_n->getChild(i);
                    tryEvalNodeContentsMathML(&timeStampNode, (Super *)this, &Super::addEventTime);
                }
            }
        }
        
        // Get Drift Diffusivity Method and associated attributes
        LIBSBML_CPP_NAMESPACE::XMLNode *driftDiffusivity_n = findChild(sbml_a, "driftDiffusivity");
        if (driftDiffusivity_n)
        {
            // Method
            LIBSBML_CPP_NAMESPACE::XMLNode *method_n = findChild(driftDiffusivity_n, "method");
            if (method_n
                && method_n->getNumChildren() == 1
                && method_n->getChild(0).isText()) {
                setComputeDriftDiffusivityMethod(method_n->getChild(0).getCharacters());
            }
            
            // Non-Linear
            string nonLinearStr = driftDiffusivity_n->getAttrValue("nonLinear");
            setNonlinearDiffusivity(boost::iequals(nonLinearStr, "true")
                                    || boost::iequals(nonLinearStr, "1"));
            
            // Compute Moments
            string computeMomentsStr = driftDiffusivity_n->getAttrValue("computeMoments");
            setComputeMoments(boost::iequals(computeMomentsStr, "true")
                              || boost::iequals(computeMomentsStr, "1"));
        }

        // get setField Method
        LIBSBML_CPP_NAMESPACE::XMLNode *setField_n = findChild(sbml_a, "setField");
        if (setField_n) {
            // Method
            LIBSBML_CPP_NAMESPACE::XMLNode *method_n = findChild(setField_n, "method");
            if (method_n
                && method_n->getNumChildren() == 1
                && method_n->getChild(0).isText()) {
                setFieldMethod(method_n->getChild(0).getCharacters());
            }
        }

        // get new individuals
        LIBSBML_CPP_NAMESPACE::XMLNode *newIndividuals_n = findChild(sbml_a, "newIndividualsMethod");
        if (newIndividuals_n) {
            // Method
            LIBSBML_CPP_NAMESPACE::XMLNode *method_n = findChild(newIndividuals_n, "method");
            if (method_n
                && method_n->getNumChildren() == 1
                && method_n->getChild(0).isText()) {
                setNewIndividualsMethod(method_n->getChild(0).getCharacters());
            }
        }

    }
    
    // Model
    LIBSBML_CPP_NAMESPACE::XMLNode *model_a = findAnnotation(sbmlModel, "model");
    if (model_a) {
        // TODO: move to using MathML
        tryEvalChildContentsMathML(model_a, "gridWidth",      (Super *)this, &Self::setGridModelWidth);
        tryEvalChildContentsMathML(model_a, "gridHeight",     (Super *)this, &Self::setGridModelHeight);
        tryEvalChildContentsMathML(model_a, "physicalWidth",  (Super *)this, &Self::setPhysicalLength);
        //tryEvalChildContentsMathML(model_a, "physicalHeight", ?);
    }
    
    // Compartments
    for (unsigned i=0; i < sbmlModel->getNumCompartments(); i++)
        addCompartment(sbmlModel->getCompartment(i));
    
    // Species - add AFTER compartments, due to references to them
    for (unsigned i=0; i < sbmlModel->getNumSpecies(); i++)
        addSpecies(sbmlModel->getSpecies(i));
    
    // Reactions - add AFTER compartments and species, due to references to them
    for (unsigned i=0; i < sbmlModel->getNumReactions(); i++)
        addReaction(sbmlModel->getReaction(i));
    
    delete document;
    
    return true;
}


void InchmanDiffusionModel::createOutputFile()
{
    Super::createOutputFile();
    
    if (outputFormat() != gpgmp::OutputHdf5)
        return;
    
    
    herr_t status;
    hid_t hdf5File = H5Fopen(outputFilename().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    
    /*
     * Write the Job's Parameters
     */
    hid_t parentOfParams = hdf5File;
    hid_t paramsRoot = H5Gcreate2(parentOfParams, "Parameters",
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    typedef std::map<std::string, Real> Parameters_t;
    BOOST_FOREACH(const Parameters_t::value_type &i, parameters()) { // output the absolute initial/input param values, before any scripts are run
        H5LTset_attribute_float(parentOfParams, "Parameters", i.first.c_str(), &i.second, 1);
    }
    status = H5Gclose(paramsRoot);
    if (parentOfParams != hdf5File) {
        status = H5Gclose(parentOfParams);
    }
    
    
    /*
     * Write the Model's InchmanFile
     */
    hid_t modelRoot = H5Gopen2(hdf5File, "Model", H5P_DEFAULT);
    status = H5LTmake_dataset_string(modelRoot, "InchmanFile", c_inchmanFile.c_str());
    
    
    H5Fclose(hdf5File);
}


void InchmanDiffusionModel::setParameter(LIBSBML_CPP_NAMESPACE::Parameter *sbmlParameter) {
    // default value
    std::string domain("single");

    // check if it's a field parameter
    LIBSBML_CPP_NAMESPACE::XMLNode *a = findAnnotation(sbmlParameter, "parameter");
    if (a)
        domain = a->getAttrValue("domain");

    if (domain.compare("field") == 0) {
       setFieldParameter(sbmlParameter->getId());

    } else if (domain.compare("single") == 0) {

        // and set it
        setParameter(sbmlParameter->getId(), sbmlParameter->getValue());
    } else {
        BOOST_LOG_TRIVIAL(error) << "Parameter "<<sbmlParameter->getId()<< " is of type domain - not supported yet!";
    }
}


Compartment * InchmanDiffusionModel::addCompartment(LIBSBML_CPP_NAMESPACE::Compartment *sbmlCompartment)
{
    Compartment *imCompartment = addCompartment(sbmlCompartment->getId());
    LIBSBML_CPP_NAMESPACE::XMLNode *a = findAnnotation(sbmlCompartment, "compartment");
    if (a) {
        tryEvalChildContentsMathML(a, "x",      imCompartment, &Compartment::setX0);
        tryEvalChildContentsMathML(a, "y",      imCompartment, &Compartment::setY0);
        tryEvalChildContentsMathML(a, "width",  imCompartment, &Compartment::setWidth);
        tryEvalChildContentsMathML(a, "height", imCompartment, &Compartment::setHeight);
    }
    BOOST_LOG_TRIVIAL(debug) <<"Creating compartment "<<*imCompartment;
    return imCompartment;
}


Species * InchmanDiffusionModel::addSpecies(LIBSBML_CPP_NAMESPACE::Species *sbmlSpecies)
{
    Species *imSpecies = addSpecies(sbmlSpecies->getId());
    
    // Get diffusionConstant AND Set initial amounts
    assert(sbmlSpecies->getModel()->getNumCompartments() == 0
           || (compartments().size() == sbmlSpecies->getModel()->getNumCompartments())); // MUST add compartments before species
    LIBSBML_CPP_NAMESPACE::XMLNode *a = findAnnotation(sbmlSpecies, "species");
    if (a)
    {
        tryEvalChildContentsMathML(a, "diffusionConstant", imSpecies, &Species::setDiffusionConstant);
        
        LIBSBML_CPP_NAMESPACE::XMLNode *list = findChild(a, "listOfCompartmentParameters");
        if (list) {
            for (unsigned int i=0; i < list->getNumChildren(); i++) {
                LIBSBML_CPP_NAMESPACE::XMLNode &child = list->getChild(i);
                if (child.getName() == "compartmentParameters") {
                    Compartment *c = compartment(child.getAttrValue("compartment"));
                    if (c) {
                        int initialAmount = 0;
                        tryEvalChildContentsMathML(&child, "initialAmount", &initialAmount);
                        c->setInitialAmount(imSpecies, initialAmount);
                    }
                    else {
                        // TODO: handle this case - warn user
                    }
                }
            }
        }

        // is the individual flag set?
        LIBSBML_CPP_NAMESPACE::XMLNode *individual = findChild(a, "individual");
        if (individual) {
            imSpecies->setHasIndividualProperties(true);
            setHasContinuousParameterSpecies(true);
        }

    }
    
    return imSpecies;
}


Reaction * InchmanDiffusionModel::addReaction(LIBSBML_CPP_NAMESPACE::Reaction *sbmlReaction)
{
    Reaction *imReaction = addReaction(sbmlReaction->getId());    
    
    /*
     Set Kinetic Law
     */
    if (sbmlReaction->isSetKineticLaw())
        imReaction->setCustomLaw(sbmlReaction->getKineticLaw()->getFormula());
    
    
    /*
     Add reactants
     */
    unsigned int numReactants = sbmlReaction->getNumReactants();
    for (unsigned int i=0; i < numReactants; ++i) {
        Species *s = species(sbmlReaction->getReactant(i)->getSpecies());
        if (s) {
            // first get stoichiometry
            int st = sbmlReaction->getReactant(i)->getStoichiometry();
            if (st < 0) BOOST_LOG_TRIVIAL(warning) << "Stoichiometry for reactant " << s->id() <<" in reaction " << imReaction->id() <<" is negative.";
            //else BOOST_LOG_TRIVIAL(warning) << "Stoichiometry for reactant  " << s->id() <<" in reaction " << imReaction->id() <<" is "<< st <<".";

            // and increase overall stoichiometry
            imReaction->increaseReactantStoichiometry(s, st);
        }
        else {
            // TODO: handle this case - warn user
        }
    }
    
    /*
     Add products
     */
    unsigned int numProducts = sbmlReaction->getNumProducts();
    for (unsigned int i=0; i < numProducts; ++i) {
        Species *s = species(sbmlReaction->getProduct(i)->getSpecies());
        if (s) {
            // first get stoichiometry
            int st = sbmlReaction->getProduct(i)->getStoichiometry();
            if (st < 0) BOOST_LOG_TRIVIAL(warning) << "Stoichiometry for product " << s->id() <<" in reaction " << imReaction->id() <<" is negative.";
            //else BOOST_LOG_TRIVIAL(warning) << "Stoichiometry for product  " << s->id() <<" in reaction " << imReaction->id() <<" is "<< st <<".";

            // and increase overall stoichiometry
            imReaction->increaseProductStoichiometry(s, st);
        }
        else {
            // TODO: handle this case - warn user
        }
    }
    
    
    /*
     Add compartment localizations
     */
    assert(sbmlReaction->getModel()->getNumCompartments() == 0
           || (compartments().size() == sbmlReaction->getModel()->getNumCompartments())); // MUST add compartments before reactions
    LIBSBML_CPP_NAMESPACE::XMLNode *list = findAnnotationChild(sbmlReaction, "reaction", "listOfCompartmentParameters");
    if (list) {
        unsigned int numChildren = list->getNumChildren();
        for (unsigned int i=0; i < numChildren; i++) {
            LIBSBML_CPP_NAMESPACE::XMLNode &child = list->getChild(i);
            if (child.getName() == "compartmentParameters")
            {
                // TODO: verify this code!
                double isAllowed = 0.0;
                tryEvalChildContentsMathML(&child, "isAllowed", &isAllowed);
                if (isAllowed != 0.0) {
                    Compartment *c = compartment(child.getAttrValue("compartment"));
                    if (c) {
                        imReaction->addCompartment(c); // IMPORTANT: call this only AFTER the reaction's been added to the model
                    }
                    else {
                        // TODO: handle this case - warn user
                    }
                }
            }
        }
    }
    
    return imReaction;
}


template<class T>
struct PrimitiveSetter {
    PrimitiveSetter(T *ptr) : ptr(ptr) {}
    void set(T value) { *ptr = value; }
    T *ptr;
};
    
bool InchmanDiffusionModel::tryEvalChildContentsString(LIBSBML_CPP_NAMESPACE::XMLNode *node, const std::string &childName, std::string *out) {
    return tryGetNodeContents<std::string, std::string &, const std::string &>(node, childName, out, &std::string::assign);
}
bool InchmanDiffusionModel::tryEvalNodeContentsString(LIBSBML_CPP_NAMESPACE::XMLNode *node, std::string *out) {
    return tryGetNodeContents<std::string, std::string &, const std::string &>(node, out, &std::string::assign);
}

template <class T>
bool InchmanDiffusionModel::tryEvalChildContentsMathML(LIBSBML_CPP_NAMESPACE::XMLNode *node, const string &childName, T *out) {
    PrimitiveSetter<T> setter(out);
    return tryEvalChildContentsMathML(node, childName, &setter, &PrimitiveSetter<T>::set);
}

template <class T>
bool InchmanDiffusionModel::tryEvalNodeContentsMathML(LIBSBML_CPP_NAMESPACE::XMLNode *node, T *out) {
    PrimitiveSetter<T> setter(out);
    return tryEvalNodeContentsMathML(node, &setter, &PrimitiveSetter<T>::set);
}

template <class C, class R, class T>
bool InchmanDiffusionModel::tryEvalChildContentsMathML(LIBSBML_CPP_NAMESPACE::XMLNode *node, const string &childName, C *outObj, R (C::*outSetter)(T)) {
    LIBSBML_CPP_NAMESPACE::XMLNode *child = findChild(node, childName);
    return tryEvalNodeContentsMathML(child, outObj, outSetter);
}

template <class C, class R, class T>
bool InchmanDiffusionModel::tryEvalNodeContentsMathML(LIBSBML_CPP_NAMESPACE::XMLNode *node, C *outObj, R (C::*outSetter)(T))
{
    bool success = false;
    if (node
        && node->getNumChildren() > 0)
    {
        string xml = LIBSBML_CPP_NAMESPACE::XMLNode::convertXMLNodeToString(&node->getChild(0)); // use the child/contents of the node
        LIBSBML_CPP_NAMESPACE::ASTNode *math = LIBSBML_CPP_NAMESPACE::readMathMLFromString(xml.c_str());
        if (math) {
            char *formula = LIBSBML_CPP_NAMESPACE::SBML_formulaToString(math);
            if (formula) {
                try {
                    boost::python::object result = boost::python::eval(formula, c_pyMainNamespace);
                    std::mem_fun(outSetter)(outObj, static_cast<T>(boost::python::extract<double>(result))); // MUST extract as double
                    success = true;
                } catch (boost::python::error_already_set const &) {
#if 1
                    cerr << "ERROR: failed to evaluate formula: \"" << formula << "\"" << endl;
                    exit(1);
#endif
                    success = false; // just to be sure...
                }
                safe_free(formula);
            }
        }
        delete math;
    }
    return success;
}

template <class C, class R, class T>
bool InchmanDiffusionModel::tryGetNodeContents(LIBSBML_CPP_NAMESPACE::XMLNode *node, const string &childName, C *outObj, R (C::*outSetter)(T)) {
    LIBSBML_CPP_NAMESPACE::XMLNode *child = findChild(node, childName);
    return tryGetNodeContents(child, outObj, outSetter);
}

template <class C, class R, class T>
bool InchmanDiffusionModel::tryGetNodeContents(LIBSBML_CPP_NAMESPACE::XMLNode *node, C *outObj, R (C::*outSetter)(T))
{
    bool success = false;
    if (node
        && node->getNumChildren() > 0)
    {
        string xml = LIBSBML_CPP_NAMESPACE::XMLNode::convertXMLNodeToString(&node->getChild(0)); // use the child/contents of the node
        std::mem_fun(outSetter)(outObj, xml);
        success = true;
    }
    return success;
}


LIBSBML_CPP_NAMESPACE::XMLNode *findChild(LIBSBML_CPP_NAMESPACE::XMLNode *node, const string &name)
{
    if (node) {
        for (unsigned int i=0; i < node->getNumChildren(); i++) {
            LIBSBML_CPP_NAMESPACE::XMLNode &child = node->getChild(i);
            if (child.getName() == name)
                return &child;
        }
    }
    return 0;
}

LIBSBML_CPP_NAMESPACE::XMLNode *findAnnotation(LIBSBML_CPP_NAMESPACE::SBase *sb, const string &name)
{
    if (sb
        && sb->isSetAnnotation())
    {
        LIBSBML_CPP_NAMESPACE::XMLNode *node = sb->getAnnotation();
        if (node)
            return findChild(node, name);
    }
    return 0;
}

LIBSBML_CPP_NAMESPACE::XMLNode *findAnnotationChild(LIBSBML_CPP_NAMESPACE::SBase *sb, const string &annotationName, const string &childName) {
    return findChild(findAnnotation(sb, annotationName), childName);
}


} // namespace gpgmp
