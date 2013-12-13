/*
 * Species.h
 *
 *  Created on: 08/09/2010
 *      Author: aidan
 */

#include "Species.h"
#include "IndividualSolver.h"

namespace gpgmp {

Species::Species(const std::string &id, Real diffusionConstant)
    :   Base(id),
      m_diffusionConstant(diffusionConstant),
      m_hasIndividualProperties(false)
{
    // create the individual properties vector
    m_rx = new std::vector<Real>(gpgmp::IndividualSolver::c_maxNumIndividuals);
    m_ry = new std::vector<Real>(gpgmp::IndividualSolver::c_maxNumIndividuals);
    m_dx = new std::vector<Real>(gpgmp::IndividualSolver::c_maxNumIndividuals);
    m_dy = new std::vector<Real>(gpgmp::IndividualSolver::c_maxNumIndividuals);

    // and add them to the dictionary
    m_individualProperties["rx"] = *m_rx;
    m_individualProperties["ry"] = *m_ry;
    m_individualProperties["diffx"] = *m_dx;
    m_individualProperties["diffy"] = *m_dy;

    delete m_rx;
    delete m_ry;
    delete m_dx;
    delete m_dy;

}

Species::~Species()
{
}


void Species::setDiffusionConstant(Real diffusionConstant) {
    m_diffusionConstant = diffusionConstant;
}


} // namespace gpgmp
