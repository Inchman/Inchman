#ifndef INDIVIDUALSOLVER_H
#define INDIVIDUALSOLVER_H

// these are to combine keys and indices in single unsigned long
// we assume that unsigned long has 64 bits
#define KEYS_MULT  4294967296
#define KEY_INDEX_PAIR(key, index) (key + index *KEYS_MULT)
#define GET_KEY(kvp) (kvp & 0xffffffff)
#define GET_INDEX(kvp) (kvp >> 32)

#include "Globals.h"
#include "DiffusionModel.h"
#include "Kernel.h"

#include <map>
#include <vector>

namespace gpgmp {

// forward declaration of Solver class
class Solver;

/**
 * This class serves to implement individuals. Each individual species
 * is represented by a list of long unsigned integers, where the first
 * 32 bits give the index of the grid cell (as laid out in the state
 * array) or 0xffffffff if the individual is not active.
 * The last 32 bits denote the individuals id. For each species, we
 * reserve a maximum number of individuals. The list will always be sorted
 * according to their cell id.
 *
 * This class is not derived from Solver, instead it is used in addition
 * to all solvers.
 */
class IndividualSolver
{
public:
    /**
     * Constructor.
     *
     * @param model The reference to the underlying diffusion model
     * @param solver Reference to the asssociated solver
     * @param context Open-CL context
     * @param queue Open-CL queue
     * @param contextDevice Open-CL device context
     */
    IndividualSolver(DiffusionModel &model,
                     const Solver &solver,
                     cl::Context *context,
                     cl::CommandQueue *queue,
                     cl::Device &contextDevice);
    ~IndividualSolver();

    /**
     * Computes the diffusion constants. This is done by updating the state array
     * according to the position of each
     * individual and then just computing the diffusivity/drift using the standard method.
     *
     * @param speciesIndex Index of the species
     * @param d_diffusionConstantsX Device array for diffusivity in x-dir
     * @param d_diffusionConstantsY Device array for diffusivity in y-dir
     * @param d_rx Device array for drift in x-dir
     * @param d_ry Device array for drift in y-dir
     * @param d_state The state array that will be filled in
     * @param simTime Simulation time (needed for custom diffusivity/drift methods)
     * @param d_errors Device array to allow error handling
     * @param d_sumMoments Device array to compute averages
     * @param d_sumStates Device array to compute totals
     */
    void computeDiffusionConstants(int speciesIndex, cl::Buffer *d_diffusionConstantsX,
                                   cl::Buffer *d_diffusionConstantsY,
                                   cl::Buffer *d_rx, cl::Buffer *d_ry,
                                   cl::Buffer *d_state,
                                   Real simTime, cl::Buffer *d_errors,
                                   cl::Buffer *d_sumMoments,
                                   cl::Buffer *d_sumStates);
    /**
     * Performs the diffusion sweep in x direction. Works like the standard algorithm
     * only that it changes the cell indices of the individuals.
     * @param seed Random seed
     * @param speciesIndex Index of the species
     * @param dt Length of time step
     * @param d_state State array
     * @param d_errors Device array for error handling
     */
    void diffuseX(int seed, int speciesIndex, float dt, cl::Buffer *d_state, cl::Buffer *d_errors);
    /**
     * Performs the diffusion sweep in y direction. Works like the standard algorithm
     * only that it changes the cell indices of the individuals.
     * @param seed Random seed
     * @param speciesIndex Index of the species
     * @param dt Length of time step
     * @param d_state State array
     * @param d_errors Device array for error handling
     */
    void diffuseY(int seed, int speciesIndex, float dt, cl::Buffer *d_state, cl::Buffer *d_errors);

    // reactions
    /**
     * After the Gillespie step, this method is called to maintain order in the
     * indivduals list. New particles could have been added or particles could have
     * been destroyed in reactions. In order to avoid resorting of the list, we
     * insert/destroy individuals after the Gillespie step using stream compactification.
     *
     * @param speciesIndex Species index
     * @param seed RNG seet
     * @param d_state state array
     */
    void maintainPopulation(int speciesIndex, int seed, cl::Buffer *d_state);

    /**
     * Fetches the position list for a particular species.
     *
     * @param speciesIndex Index of the species
     * @return Combined position/id list for this species
     */
    std::vector<cl_ulong> getPositionList(int speciesIndex);
    /**
     * Fetches the custom property list for all species.
     *
     * @param speciesIndex Species Index
     * @return The property list
     */
    std::map<int, std::vector<Real> > getPropertiesList(int speciesIndex);

    // output
    /** @name Output functions for debugging purposes.
     *  This set of functions can be used for debugging to output the
     *  various internal states to text files.
     */
    ///@{
    void writePositions(int speciesIndex, std::string filename, int index=0, int num=0);
    void writeAll(int speciesIndex);
    void writeDiffusionConstants(std::string filename,
                                     cl::Buffer *d_diffusionConstantsX, cl::Buffer *d_diffusionConstantsY,
                                     cl::Buffer *d_rx, cl::Buffer *d_ry);
    void writeProperties(std::string filename, int speciesIndex);
    void writeCellIndices(int speciesIndex, std::string filename);
    void writeReductionOffsets(std::string filename, int index=0);
    void checkPositionListConsistency(int speciesIndex, cl::Buffer *d_state);
    ///@}

private:    
    void writeLostParticles(std::string filename);
    void writeSegmentedScanBuffers(std::string filename);
    void writeState(std::string filename, cl::Buffer *d_state);

    /** @name Parallel radix sort functions.  The radix sort algorithm
     * is an adapted version of the one found in the NVIDIA toolkit.
     *
     */
    ///@{

    /**
     * Sorts the species using parallel radix sort.
     * @param species The species index to sort.
     * @param nIndividuals Number of active individuals.
     */
    void sortSpecies(int species, uint nIndividuals);
    void radixSortStep(uint startbit, int species, uint nkeys, uint totalBlocks);
    void scanExclusive(cl::Buffer *d_Dst, cl::Buffer *d_Src, uint batchSize, uint arrayLength);
    unsigned int factorRadix2(unsigned int &log2L, unsigned int L);
    unsigned int iSnapUp(unsigned int dividend, unsigned int divisor);
    ///@}

    /**
     * Computes the start (left) and end (right) indices of the particles in each cell.
     * This method assumes that the individual list is sorted, of course.
     *
     * @param speciesIndex Species index
     * @param nIndividuals Number of active individauls
     */
    void computeCellIndices(int speciesIndex, uint nIndividuals);

    /**
     * Performs a segmented scan on the data array. The array can be 1D or 2D.
     *
     * @param d_data Data array
     * @param d_flags Flags array
     * @param arrayLengthX X dimension
     * @param arrayLengthY Y dimension
     */
    void segmentedScan(cl::Buffer *d_data, cl::Buffer *d_flags, uint arrayLengthX, uint arrayLengthY);
    /**
     * Unit test for the segmented scan method.
     */
    void testSegmentedScan();

    /** @name Helper functions for handling device buffers.
     */
    ///@{
    // helper functions copied from diffusion model
    // todo: would be nice to have a proper buffer class..
    template <class T> cl::Buffer *oclCreateBuffer(size_t count, cl_mem_flags flags = CL_MEM_READ_WRITE) const;
    template <class ContainerT> cl::Buffer *oclCreateBufferFrom(ContainerT &from, cl_mem_flags flags = CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR) const;
    template <class T> cl::Buffer *oclCreateBufferFrom(T *from, size_t count, cl_mem_flags flags) const;
    template <class ContainerT> void oclCopy(const cl::Buffer &from, ContainerT &to);
    template <class T> void oclCopy(const cl::Buffer &from, T *to, size_t count);
    template <class ContainerT> void oclCopy(const ContainerT &from, cl::Buffer &to);
    template <class T> void oclCopy(const T *from, cl::Buffer &to, size_t count);
    void oclCopyDeviceBuffers(cl::Buffer &from, cl::Buffer &to, size_t count);
    ///@}
    ///

    /**
     * Builds the kernels from the template.
     */
    void buildKernels();

public:
    static const unsigned int c_maxNumIndividuals = 262144; ///< Maximum number of individuals per species

private:
    typedef cl_ulong type_key; ///>  type for key-value pairs


    DiffusionModel &m_diffusionModel; ///< Reference to the associated diffusion model
    const Solver &m_solver; ///< Reference to the associated solver class

    cl::Context *m_context; ///< Pointer to the GPU context
    cl::CommandQueue *m_commandQueue; ///< Pointer to the current command queue
    cl::Device &m_contextDevice; ///< Pointer to associated device

    // these buffers holds all individual properties and positions/keys, indexed per species
    std::map<int, cl::Buffer *> individualProperties; ///< Handle to device buffer for individual properties, indexed per species.
    std::map<int, cl::Buffer *> individualPositions; ///< Handle to device buffer for individual positions, indexed per species.
    std::map<int, uint> m_nIndividuals; ///< Number of individuals indexed per species.

    // these buffers are for the cell indices. We need one for each species.
    std::map<int,cl::Buffer *> d_cellIndicesLeft; ///< Buffer to index the start of the sorted list for each cell.
    std::map<int,cl::Buffer *> d_cellIndicesRight; ///< Buffer to index the end of the sorted list for each cell.

    // various helper buffers .. will be re-used
    cl::Buffer *d_tempKeys; ///< Temporary buffer to hold keys. Will be re-used for each species
    cl::Buffer *d_mCounters; ///< Holds counters for radix sort.
    cl::Buffer *d_mCountersSum; ///< Holds sum of counters for radix sort.
    cl::Buffer *d_mBlockOffsets; ///< Holds block offsets for radix sort counters.
    cl::Buffer *d_buffer; ///< a buffer (todo: find out what it's for :))
    cl::Buffer *d_reductionOffsets; ///< Used for stream reduction after removing particles
    cl::Buffer *d_reductionOffsetsBlocks; ///< Used for stream reduction after removing particles
    cl::Buffer *d_reductionOffsetsAdding; ///< Used for stream reduction after adding particles
    cl::Buffer *d_reductionOffsetsBlocksAdding; ///< Used for stream reduction after adding particles
    cl::Buffer *d_total; ///< Single value holding total count for species
    cl::Buffer *d_totalAdded; ///< Single value holding total added count
    cl::Buffer *d_lostParticles; ///< keeps track of lost individual particles
    cl::Buffer *d_particlesBuffer; ///< temporary buffer to hold individuals
    cl::Buffer *d_domainReductionOffsets; ///< To reduce over number of new particles
    cl::Buffer *d_newParticles; ///< Contains new particles per cell

    // buffers for segmented scan
    cl::Buffer *d_ssFlags; ///< Flags for the segmented scan
    cl::Buffer *d_ssDataBlocks; ///< Segmented-scan block results
    cl::Buffer *d_ssFlagsBlocks; ///< Segmented-scan block OR results
    cl::Buffer *d_ssFirstBlockFlag; ///< Segmented-scan first flag of each block
    cl::Buffer *d_ssPartialFlags; ///< Segmented-scan partial OR tree


    // kernels
    Kernel m_kernel_radixSortBlocksKeysOnly; ///< Kernel to sort keys
    Kernel m_kernel_findRadixOffsets;
    Kernel m_kernel_scanNaive;
    Kernel m_kernel_reorderDataKeysOnly;
    Kernel m_kernel_scanExclusiveLocal1;
    Kernel m_kernel_scanExclusiveLocal2;
    Kernel m_kernel_uniformUpdate;
    Kernel m_kernel_computeCellIndices;
    Kernel m_kernel_computeDiffusionConstants;
    Kernel m_kernel_individualDiffuseX;
    Kernel m_kernel_individualDiffuseY;
    Kernel m_kernel_padKeys;
    Kernel m_kernel_maintainPopulation;
    Kernel m_kernel_scanIndividuals;
    Kernel m_kernel_scanIndividualsBlocks;
    Kernel m_kernel_scatterIndividuals;
    Kernel m_kernel_writeLostParticles;
    Kernel m_kernel_segmentedScanBlock;
    Kernel m_kernel_segmentedScanUpSweep;
    Kernel m_kernel_segmentedScanDownSweep;
    Kernel m_kernel_segmentedScanBlock_single;
    Kernel m_kernel_computeReductionOffsets;
    Kernel m_kernel_clearReductionOffsetsAdding;
    Kernel m_kernel_updateNewParticles;
    Kernel m_kernel_copyTempToKeys;

    // static constants
    static const uint c_keyBits = 32; ///< Number of significant bits for sorting. todo: can adapt depending on grid size ..

    // Topology for radix sort and scan
    static const unsigned int c_warp_size = 32; ///< Warp size
    static const unsigned int c_bitStep = 4; ///< Number of bits per radix step
    static const int c_ctaSize = 128; ///< Workgroup size for radix
    uint m_radixNumBlocks, m_radixNumBlocks2; ///< Block topology for radix sort
    size_t m_wgSizeScan; ///<  Work group size for scan.

    static const int c_ssWgSize = 16; ///< Workgroup size for segmented scan
};

template <class ContainerT>
void IndividualSolver::oclCopy(const cl::Buffer &from, ContainerT &to) {
    assert(!to.empty());
    oclCopy(from, &to.front(), to.size());
}
template <class T>
void IndividualSolver::oclCopy(const cl::Buffer &from, T *to, size_t count) {
    m_commandQueue->enqueueReadBuffer(from, CL_TRUE, 0, count * sizeof(T), (void *)to);
    m_commandQueue->flush();
}

template <class ContainerT>
void IndividualSolver::oclCopy(const ContainerT &from, cl::Buffer &to) {
    assert(!from.empty());
    oclCopy(&from.front(), to, from.size());
}
template <class T>
void IndividualSolver::oclCopy(const T* from, cl::Buffer &to, size_t count) {
    m_commandQueue->enqueueWriteBuffer(to, CL_TRUE, 0, count * sizeof(T), (const void *)from);
    m_commandQueue->flush();
}

template <class T>
cl::Buffer *IndividualSolver::oclCreateBuffer(size_t count, cl_mem_flags flags) const {
    return oclCreateBufferFrom((T *)0, count, flags);
}

template <class ContainerT>
cl::Buffer *IndividualSolver::oclCreateBufferFrom(ContainerT &from, cl_mem_flags flags) const {
    return oclCreateBufferFrom(&from.front(), from.size(), flags);
}

template <class T>
cl::Buffer *IndividualSolver::oclCreateBufferFrom(T *from, size_t count, cl_mem_flags flags) const
{
    return new cl::Buffer(*m_context, flags, count * sizeof(T), (void *)from);
}
}
#endif // INDIVIDUALSOLVER_H
