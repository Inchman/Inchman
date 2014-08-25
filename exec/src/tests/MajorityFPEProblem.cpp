#include "MajorityFPEProblem.h"
#include <sstream>
#include "Compartment.h"

namespace gpgmp {

	/**
	*/
	MajorityFPEProblem::MajorityFPEProblem(Real length, int dx, int dy,
		int numMolecules)
		: DiffusionModel(length, dx, dy)
	{
		// one species
		Species *speciesA = addSpecies("A", 1.);

		// set boundary mask
		const cl_int4 boundaryMask = { { 1, 1, 1, 1 } };
		const cl_int4 sourceMask = { { 0, 0, 0, 0 } };
		setBoundaryMasks(boundaryMask, sourceMask);

		// add source compartment
		Compartment *compartmentSource = new Compartment("TrueSource", dx / 2, 0, dx / 2, dy);
		compartmentSource->setInitialAmount(speciesA, numMolecules * dy, gpgmp::HomogeneousDistribution);
		addCompartment(compartmentSource);

		// add user parameters for kernel
		setFieldParameter("diffFieldX");
		setFieldParameter("driftFieldX");

		// set diffusivity/drift method
		std::ostringstream ss;
		ss <<
			"// sets the diffusivity and drift to the constants given by the host parameters\n"
			"All->DiffusivityX = diffFieldX;\n"
			"All->DiffusivityY = 0.;\n"
			"\n"
			"All->DriftX = driftFieldX;\n"
			"All->DriftY = 0.;\n";

		// set non-linear diffusivity and drift
		setComputeDriftDiffusivityMethod(ss.str());

		// set init script
		loadInitScript("C:\\Users\\Matthias\\Documents\\Projects\\data\\majority\\flame\\virtual\\avoidance_large\\results\\radius_0.01\\initFPESolver.py");

		setNonlinearDiffusivity(false);
		setComputeMoments(false);
	} // constructor

} // namespace gpgmp
