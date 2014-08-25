#ifndef MAJORITYFPEPROBLEM_H
#define MAJORITYFPEPROBLEM_H
#include "DiffusionModel.h"

namespace gpgmp {
	/**
	**/

	class MajorityFPEProblem : public DiffusionModel
	{
	public:
		/**
		Sets up the model.
		*/
		MajorityFPEProblem(Real length, int dx, int dy,
			int numMolecules);
	};

}// namespace gpgmp
#endif // NONLINEARPROBLEM_H
