#ifndef MAJORITYVOTEPROBLEM_H
#define MAJORITYVOTEPROBLEM_H

#include "DiffusionModel.h"
#include "Globals.h"

#include <hdf5.h>

namespace gpgmp {

class MajorityVoteProblem : public DiffusionModel
{
public:
    //MajorityVoteProblem(Real length, int dx, int dy, uint numMolecules, uint nReds, uint numEncounters, Real avoidanceRadius); // this one for individuals

    MajorityVoteProblem(Real length, int dx, int dy, Real diffX, Real diffY, Real mux, Real muy, uint numMolecules, uint nReds, uint numEncounters, Real reactionProb, const std::string &initScriptsPath);

    MajorityVoteProblem(Real length, int dx, int dy, uint numMolecules, uint nReds, uint numEncounters, Real reactionProb, Real recognitionError, const std::string &initScriptsPath);

    MajorityVoteProblem(Real length, int dx, int dy, uint numMolecules, uint nReds, Real reactionProb, const std::string &initScriptsPath);

    void writeSingleTimeHDF5(hid_t currentDataGroup, const char *buffer, size_t bufferStep, hid_t bufferType);

private:
    unsigned int getSpeciesIndexFrom2D(int i, int j, int n) {/*std::cout <<i<<", "<<j<<", "<<n<<" yields "<<-i*(i-2*n-1)/2+j<<".\n"<<std::flush;*/ return -i*(i-2*n-1)/2+j;}
    unsigned int getChangeoverLimit(uint pos, uint n) {return (4*n + 3*n*n - 7)/8;}
    bool isMarginal(uint pos, uint n);
    unsigned int getNextRedPosition(uint pos, uint n);
private:
    bool m_individual;
    bool m_isNonDiffusive;
};

} // namespace gpgmp
#endif // MAJORITYVOTEPROBLEM_H
