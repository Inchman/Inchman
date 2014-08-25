/*
* SimpleDiffusion.cpp
*
*  Created on: 05/03/2010
*      Author: matthias
*/

#include "DiffusionProblem.h"

#include "Compartment.h"

#include <iostream>


namespace gpgmp {

/**
    Defines two compartments, "World" and "Source". Registers
    the species "A" and sets the initial amount for compartment "Source".

    If a non-vanishing drift field or an anisotropic diffusivity
    is specified, we set up an inhomogeneous
    simulation. The drift-diffusion method is set to gpgmp::DT_HOMOGENEOUS and
    the kernel parameters are set accordingly. Both, diffusivity and drift can be anisotropic
    but are homogeneous over the domain.

    \see setComputeDriftDiffusivityMethod
    \see setParameter
    \see AnnihilationProblem
    */
DiffusionProblem::DiffusionProblem(Real length, int dx, int dy, Real diffusionConstantX, int numMolecules,
                                   Real diffusionConstantY, Real rx, Real ry)
:   DiffusionModel(length, dx, dy)
{
    // add one species
    Species *speciesA = addSpecies("A", diffusionConstantX);
    
    // set boundary mask to be empty/clear
    clearBoundaryMasks();

    // add World compartment
    addCompartment("World", 0, 0, dx-1, dy-1);
    
    // add Source compartment
    int x0 = dx/2;
    int y0 = dy/2;
    Compartment *source = addCompartment("Source", x0, y0, x0, y0);
    source->setInitialAmount(speciesA, numMolecules, gpgmp::HomogeneousDistribution);

    if (diffusionConstantY != 0. || rx != 0. || ry != 0.) {
        // add user parameters for kernel
        setParameter("diffX", diffusionConstantX);
        setParameter("diffY", diffusionConstantY);
        setParameter("driftX", rx);
        setParameter("driftY", ry);

        // set homogeneous diffusivity and drift
        setComputeDriftDiffusivityMethod(gpgmp::DT_HOMOGENEOUS);
    }
}


} // namespace gpgmp
