#include "InhomogeneousSolver.h"
#include "Species.h"

#include <sstream>

namespace gpgmp {
InhomogeneousSolver::InhomogeneousSolver(DiffusionModel &diffusionModel,
                                         std::vector<int> *runState,
                                         std::vector<Real> *deterministicState, SolverType solverType,
                                         bool deterministic, bool hasIndividualSpecies, bool computeMoments,
                                         bool nonlinearDiffusivity)
    : Solver(diffusionModel, runState, deterministicState, solverType, deterministic, hasIndividualSpecies),
      m_computeMoments(computeMoments), m_nonlinearDiffusivity(nonlinearDiffusivity),
      m_smallestDtX(diffusionModel.numSpecies()), m_smallestDtY(diffusionModel.numSpecies()),
      m_nextDiffusionTimeX(diffusionModel.numSpecies(),0.), m_nextDiffusionTimeY(diffusionModel.numSpecies(),0.),
      m_diffusionDtX(diffusionModel.numSpecies()), m_diffusionDtY(diffusionModel.numSpecies()),
      m_meanX(diffusionModel.numSpecies()), m_sumState(diffusionModel.numSpecies())
{
    // initialize lost particles buffer
    if (!deterministic)
        d_lostParticles = oclCreateBuffer<cl_int4>(diffusionModel.gridArea());

    // Now build program and create kernels
    std::ostringstream ss;

    // Generate CL code and build the program
    // we need that also for the deterministic solver.. it's got the computeDiffusionConstants
    // and computeMoments kernels.
    oclGenerateHeader(ss);
    oclReadAndReplaceTemplateFile(ss, "StochasticInhomogeneousDriftDiffusion.cl");

    m_stochasticProgram = oclBuildProgram(ss);

    // build deterministic diffusion program
    if (deterministic) {
        std::ostringstream sss;

        oclGenerateHeader(sss);
        oclReadAndReplaceTemplateFile(sss, "Deterministic.cl");
        m_deterministicProgram = oclBuildProgram(sss);
    }

    // Create kernels
    std::cout <<"Creating kernels..\n"<<std::flush;

    // we need these kernels for deterministic and stochastic
    // we use the stochastic time-step and directional splitting for deterministic
    // this is probably overkill..
    computeDiffusionConstants_forSpecies.resize(m_diffusionModel.numSpecies());
    for (size_t i=0; i<m_diffusionModel.numSpecies(); ++i) {
        std::string name = "inhomogeneous_computeDiffusionConstants_forSpecies_" + m_diffusionModel.species()[i]->id();
        computeDiffusionConstants_forSpecies[i] = Kernel(m_stochasticProgram, name.c_str(),
                                                         cl::EnqueueArgs(*m_commandQueue, m_globalRange));
        computeDiffusionConstants_forSpecies[i].setLocalSize(m_localRange);
    }

    computeMomentsX = Kernel(m_stochasticProgram,
                             "inhomogeneous_computeMomentsX",
                             cl::EnqueueArgs(*m_commandQueue, m_globalRange));
    computeMomentsX.setLocalSize(m_localRange);
    computeMomentsX.instance().setArg(4, 2*LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(Real),0);

    reduceMomentsBlocks = Kernel(m_stochasticProgram,
                                 "inhomogeneous_reduceMomentsBlocks",
                                 cl::EnqueueArgs(*m_commandQueue, m_globalRangeReduce));
    reduceMomentsBlocks.setLocalSize(m_globalRangeReduce);
    reduceMomentsBlocks.instance().setArg(6, 2*m_diffusionModel.gridWidth()*m_diffusionModel.gridHeight()/(LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK)*sizeof(Real),0);

    computeIndividualTimestepX = Kernel(m_stochasticProgram,
                                        "inhomogeneous_computeIndividualTimestepX",
                                        cl::EnqueueArgs(*m_commandQueue, m_globalRange));
    computeIndividualTimestepX.setLocalSize(m_localRange);
    computeIndividualTimestepX.instance().setArg(10, LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(Real),0);

    computeIndividualTimestepY = Kernel(m_stochasticProgram,
                                        "inhomogeneous_computeIndividualTimestepY",
                                        cl::EnqueueArgs(*m_commandQueue, m_globalRange));
    computeIndividualTimestepY.setLocalSize(m_localRange);
    computeIndividualTimestepY.instance().setArg(10, LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(Real),0);

    reduceBlocks = Kernel(m_stochasticProgram,
                          "inhomogeneous_reduceBlocks",
                          cl::EnqueueArgs(*m_commandQueue, m_globalRangeReduce));
    reduceBlocks.setLocalSize(m_globalRangeReduce);
    reduceBlocks.instance().setArg(3, m_localSizeReduce, 0);


    // these are only needed for the stochastic solver
    if (!deterministic) {
        diffuseX = Kernel(m_stochasticProgram,
                          "inhomogeneous_diffuseX",
                          cl::EnqueueArgs(*m_commandQueue, m_globalRange));
        diffuseX.setLocalSize(m_localRange);

        diffuseY = Kernel(m_stochasticProgram,
                          "inhomogeneous_diffuseY",
                          cl::EnqueueArgs(*m_commandQueue, m_globalRange));
        diffuseY.setLocalSize(m_localRange);

        updateStateX = Kernel(m_stochasticProgram,
                              "inhomogeneous_updateStateX",
                              cl::EnqueueArgs(*m_commandQueue, m_globalRange));
        updateStateX.setLocalSize(m_localRange);

        updateStateY = Kernel(m_stochasticProgram,
                              "inhomogeneous_updateStateY",
                              cl::EnqueueArgs(*m_commandQueue, m_globalRange));
        updateStateY.setLocalSize(m_localRange);

    }// create stochastic kernels
    else if (deterministic)
        // all deterministic kernels
    {
        // inhomogeneous diffusion kernel
        inhomogeneous_deterministic_diffuseX = Kernel(m_deterministicProgram,
                                                      "inhomogeneous_deterministic_diffuseX",
                                                      cl::EnqueueArgs(*m_commandQueue, m_globalRange));
        inhomogeneous_deterministic_diffuseX.setLocalSize(m_localRange);
        inhomogeneous_deterministic_diffuseY = Kernel(m_deterministicProgram,
                                                      "inhomogeneous_deterministic_diffuseY",
                                                      cl::EnqueueArgs(*m_commandQueue, m_globalRange));
        inhomogeneous_deterministic_diffuseY.setLocalSize(m_localRange);
    } // deterministic kernels

    std::cout <<"Finished creating kernels. Allocating buffers..\n"<<std::flush;

    // TODO: We might want to bundle these? less kernel parameters to pass!
    // device array for diffusion constants and drift vector
    d_diffusionConstantsX =
            oclCreateBuffer<Real>(m_diffusionModel.gridArea() * m_diffusionModel.numSpecies());
    d_diffusionConstantsY =
            oclCreateBuffer<Real>(m_diffusionModel.gridArea() * m_diffusionModel.numSpecies());
    d_rx                  =
            oclCreateBuffer<Real>(m_diffusionModel.gridArea() * m_diffusionModel.numSpecies());
    d_ry                  =
            oclCreateBuffer<Real>(m_diffusionModel.gridArea() * m_diffusionModel.numSpecies());

    // device array for individual (cell) time steps
    std::cout <<"Allocating "<< sizeof(Real)* m_diffusionModel.gridArea() * m_diffusionModel.numSpecies() << " bytes for minTimeBlock.\n";
    d_timeSteps = oclCreateBuffer<Real>(m_diffusionModel.gridArea() * m_diffusionModel.numSpecies());

    // device array for minimum time step per block
    std::cout <<"Allocating "<< sizeof(Real)*m_globalRangeReduce[0]<< " bytes for minTimeBlock.\n";
    d_minTimeBlock =
            oclCreateBuffer<Real>(m_globalRangeReduce[0]);

    // device array for global minimum time step
    d_globalMinTimestep = oclCreateBuffer<Real>(1);

    // device array for mean X per block
    // todo: can probably put in one - less memory copying
    d_sumMomentsBlocks =
            oclCreateBuffer<Real>(m_globalRangeReduce[0]);
    d_sumStatesBlocks  =
            oclCreateBuffer<Real>(m_globalRangeReduce[0]);

    // device array for global minimum time step for each species
    d_sumMoments = oclCreateBufferFrom(m_meanX);
    d_sumStates  = oclCreateBufferFrom(m_sumState);

    // device array for transition probabilites
    d_transitionProbabilities =
            oclCreateBuffer<cl_float4>(m_diffusionModel.gridArea() * m_diffusionModel.numSpecies());

    std::cout <<"Finished allocating device buffers.\n"<<std::flush;

    // compute moments for each species
    if (m_computeMoments) {
        // compute them
        for (cl_int speciesIndex=0; speciesIndex<m_diffusionModel.numSpecies(); speciesIndex++) {
            this->computeMoments(speciesIndex);
        }
        // and copy back onto device
        oclCopy(*d_sumMoments, m_meanX);
        oclCopy(*d_sumStates, m_sumState);

        // then print our the details :)
        for (size_t speciesIndex=0; speciesIndex<m_diffusionModel.numSpecies(); speciesIndex++) {
            std::cout <<"Finding total of "<<m_sumState[speciesIndex]<<" for species "
                     << speciesIndex <<" and a x moment of "<<m_meanX[speciesIndex]<<"."
                     << " That gives a <x>="<<m_meanX[speciesIndex]/m_sumState[speciesIndex]<<".\n";
            std::cout << std::flush;
        }
        std::cout <<"done with computing moments.\n" << std::flush;
    } // compute moments

    std::cout <<"Computing first time steps.\n"<<std::flush;

    // compute diffusivity/drift field for each species and update time step
    for (cl_int speciesIndex=0; speciesIndex<m_diffusionModel.numSpecies(); speciesIndex++)
    {
        // diffusivity
        computeDiffusivityDrift(speciesIndex, 0.);

        // timestep X
        double dtx = computeTimestepX(speciesIndex, 0.);
        std::cout <<"Finding first time step (X) "<< dtx
                 <<" for species " << speciesIndex << ".\n",
        std::cout << std::flush;

        // timestep Y
        double dty = computeTimestepY(speciesIndex, 0.);
        std::cout <<"Finding first time step (Y) "<< dty
                 <<" for species " << speciesIndex << ".\n",
        std::cout << std::flush;

        // and also save into relative and absolute arrays
        m_diffusionDtX[speciesIndex] = dtx;
        m_diffusionDtY[speciesIndex] = dty;
        m_nextDiffusionTimeX[speciesIndex] = dtx;
        m_nextDiffusionTimeY[speciesIndex] = dty;
    }
}

InhomogeneousSolver::~InhomogeneousSolver()
{

    delete d_lostParticles;

    delete d_diffusionConstantsX;
    delete d_diffusionConstantsY;
    delete d_rx;
    delete d_ry;

    delete d_timeSteps;
    delete d_minTimeBlock;
    delete d_globalMinTimestep;
    delete d_sumMomentsBlocks;
    delete d_sumStatesBlocks;
    delete d_sumMoments;
    delete d_sumStates;
    delete d_transitionProbabilities;

    // todo: iterate through indivdual properties buffers and delete them
    // todo: release kernels etc!
}

std::vector<Real> InhomogeneousSolver::getNextDiffusionTime()
{
    // return vector
    std::vector<Real> ret(m_diffusionModel.numSpecies());

    for (int i=0; i<m_diffusionModel.numSpecies(); i++)
    {
        if (m_nextDiffusionTimeX[i] < m_nextDiffusionTimeY[i])
            ret[i] = m_nextDiffusionTimeX[i];
        else
            ret[i] = m_nextDiffusionTimeY[i];
    }

    return ret;
}

void InhomogeneousSolver::diffuse(int speciesIndex, Real currentTime, int &seed)
{
    // diffuse X or Y or both?
    if (m_nextDiffusionTimeX[speciesIndex] <= currentTime)
    {
        //printf("Diffusing species %d in x-direction at time=%g. (nextDiffusionTimeX=%g.)\n",
        //       speciesIndex, currentTime,  m_nextDiffusionTimeX[speciesIndex]);

        // diffuse if necessary
        if (!m_deterministic) {
            // check if it is an individual species
            if (m_individualSpecies[speciesIndex]) {
                //std::cout <<"Diffusing individual species "<<speciesIndex<<" in X direction.\n"<<std::flush;
                // individual diffusion
                m_individualSolver -> diffuseX(seed, speciesIndex, m_diffusionDtX[speciesIndex],  d_state, d_errors);
                seed++;
                // and we need to recompute the diffusion constants .. but will do that below
            } else {
                // use normal diffusion kernel
                diffuseX(*d_state, *d_transitionProbabilities,
                         *d_lostParticles, seed, speciesIndex,
                         m_diffusionDtX[speciesIndex], *d_errors,
                         *d_rx, *d_ry,
                         *d_diffusionConstantsX, *d_diffusionConstantsY);
                seed++;
                // and update the state
                updateStateX(*d_state, *d_lostParticles, speciesIndex);
            }
        } else if (m_deterministic) {
            inhomogeneous_deterministic_diffuseX(*d_state, m_diffusionDtX[speciesIndex],
                                                 *d_diffusionConstantsX,
                                                 *d_rx, speciesIndex);
        } // diffuse X

        // now compute the next time step for x

        // for some problems we need to compute the moments first
        if (m_computeMoments)
        {
            // compute mean x
            computeMomentsX(*d_state, *d_sumMomentsBlocks, *d_sumStatesBlocks, speciesIndex/*, *d_errors*/);

            // reduce over blocks
            reduceMomentsBlocks(*d_sumMomentsBlocks, *d_sumStatesBlocks,
                                *d_sumMoments, *d_sumStates, speciesIndex, *d_errors);

        } // compute new diffusion constants

        // we might need to compute the new diffusion constants first
        // todo: these three had simTime instead of currentTime previously.. hmmm..
        if (m_nonlinearDiffusivity || m_individualSpecies[speciesIndex]) {
            if (m_individualSpecies[speciesIndex]) {
                m_individualSolver->computeDiffusionConstants(speciesIndex,
                                                            d_diffusionConstantsX,
                                                            d_diffusionConstantsY, d_rx, d_ry,
                                                            d_state, currentTime, d_errors,
                                                            d_sumMoments, d_sumStates);
            } else {
                // compute the according diffusion constant
                if (m_hasFields)
                    computeDiffusionConstants_forSpecies[speciesIndex](*d_diffusionConstantsX,
                                                            *d_diffusionConstantsY,
                                                            *d_rx, *d_ry, *d_state, *d_fields, currentTime,
                                                            *d_errors,*d_sumMoments, *d_sumStates);

                else
                    computeDiffusionConstants_forSpecies[speciesIndex](*d_diffusionConstantsX,
                                                        *d_diffusionConstantsY,
                                                        *d_rx, *d_ry, *d_state, currentTime,
                                                        *d_errors,*d_sumMoments, *d_sumStates);
            } // compute diffusion constants

            // and compute new timestep (X)
            computeIndividualTimestepX(*d_timeSteps,
                                       *d_minTimeBlock,
                                       *d_diffusionConstantsX, *d_diffusionConstantsY,
                                       *d_rx, *d_ry,
                                       *d_transitionProbabilities, speciesIndex, currentTime,
                                       *d_errors);

            // reduce over blocks
            reduceBlocks(*d_errors, *d_minTimeBlock, *d_globalMinTimestep);

            // and copy it back into time steps array
            oclCopy(*d_globalMinTimestep, &m_smallestDtX[speciesIndex], 1);
        } // re-compute diffusivity

        // set it as next diffusion time
        m_nextDiffusionTimeX[speciesIndex] += m_smallestDtX[speciesIndex];
        m_diffusionDtX[speciesIndex] = m_smallestDtX[speciesIndex];

        //std::cout <<"New nextDiffusiontimeY is "<<m_nextDiffusionTimeY[speciesIndex]<<"\n"<<std::flush;
        //std::cout <<"New diffusionDtY is "<<m_diffusionDtY[speciesIndex]<<"\n"<<std::flush;

    } // diffusion sweep X

    // do Y sweep if needed
    if (m_nextDiffusionTimeY[speciesIndex] <= currentTime)
    {
        //printf("Diffusing species %d in y-direction at time=%g. (nextDiffusionTimeY=%g.)\n",
        //       speciesIndex, currentTime, m_nextDiffusionTimeY[speciesIndex]);

        if (!m_deterministic) {
            // check if it is an individual species
            if (m_individualSpecies[speciesIndex]) {
                //std::cout <<"Diffusing individual species "<<speciesIndex<<" in Y direction.\n"<<std::flush;
                // indidivual diffusion
                m_individualSolver -> diffuseY(seed, speciesIndex, m_diffusionDtY[speciesIndex],d_state, d_errors);
                seed++;
                // and we need to recompute the diffusion constants .. but will do that below
            } else {
                // diffuse Y if necessary
                diffuseY(*d_state, *d_transitionProbabilities,
                         *d_lostParticles, seed, speciesIndex, m_diffusionDtY[speciesIndex],
                         *d_errors, *d_rx, *d_ry,
                         *d_diffusionConstantsX, *d_diffusionConstantsY);
                seed++;
                // and update the state
                updateStateY(*d_state, *d_lostParticles, speciesIndex);
            } // diffuse stochastically (Y)
        } else if (m_deterministic) {
            inhomogeneous_deterministic_diffuseY(*d_state, m_diffusionDtY[speciesIndex],
                                                 *d_diffusionConstantsY,
                                                 *d_ry, speciesIndex);
        } // diffuse Y

        // todo: implement compute moments Y

        // compute the new diffusion constants here when needed
        // todo: these three had simTime instead of currentTime previously.. hmmm..
        if (m_nonlinearDiffusivity || m_individualSpecies[speciesIndex]) {
            if (m_individualSpecies[speciesIndex]) {
                m_individualSolver->computeDiffusionConstants(speciesIndex,
                                                            d_diffusionConstantsX,
                                                            d_diffusionConstantsY, d_rx, d_ry,
                                                            d_state, currentTime, d_errors,
                                                            d_sumMoments, d_sumStates);
            } else {
                // compute the according diffusion constant
                if (m_hasFields)
                    computeDiffusionConstants_forSpecies[speciesIndex](*d_diffusionConstantsX,
                                                            *d_diffusionConstantsY,
                                                            *d_rx, *d_ry, *d_state, *d_fields, currentTime,
                                                            *d_errors, *d_sumMoments, *d_sumStates);
                else
                    computeDiffusionConstants_forSpecies[speciesIndex](*d_diffusionConstantsX,
                                                                            *d_diffusionConstantsY,
                                                                            *d_rx, *d_ry, *d_state, currentTime,
                                                                            *d_errors, *d_sumMoments, *d_sumStates);

            } // compute new diffusion constants

            // and compute the next time step for y
            computeIndividualTimestepY(*d_timeSteps,
                                       *d_minTimeBlock,
                                       *d_diffusionConstantsX, *d_diffusionConstantsY,
                                       *d_rx, *d_ry, *d_transitionProbabilities,
                                       speciesIndex, currentTime, *d_errors);

            // reduce over blocks
            reduceBlocks(*d_errors, *d_minTimeBlock, *d_globalMinTimestep);

            // and copy it back into time steps array
            oclCopy(*d_globalMinTimestep, &m_smallestDtY[speciesIndex], 1);
        }

        // set it as next diffusion time
        m_nextDiffusionTimeY[speciesIndex] += m_smallestDtY[speciesIndex];
        m_diffusionDtY[speciesIndex] = m_smallestDtY[speciesIndex];

        //std::cout <<"New nextDiffusiontimeY is "<<m_nextDiffusionTimeY[speciesIndex]<<"\n"<<std::flush;
        //std::cout <<"New diffusionDtY is "<<m_diffusionDtY[speciesIndex]<<"\n"<<std::flush;
    } // do Y sweep
}

void InhomogeneousSolver::computeMoments(cl_int speciesIndex)
{
    // Compute Mean X
    computeMomentsX(*d_state, *d_sumMomentsBlocks, *d_sumStatesBlocks, speciesIndex/*, *d_errors*/);
    // reduce over blocks
    reduceMomentsBlocks(*d_sumMomentsBlocks, *d_sumStatesBlocks, *d_sumMoments, *d_sumStates, speciesIndex, *d_errors);
}

void InhomogeneousSolver::computeDiffusivityDrift(cl_int speciesIndex, Real simTime)
{
    // check if species is individual-based
    if (m_individualSpecies[speciesIndex]) {
        std::cout <<"Computing (individual) diffusion constants for species " << speciesIndex <<".\n"<<std::flush;
        m_individualSolver->computeDiffusionConstants(speciesIndex,
                                                      d_diffusionConstantsX,
                                                      d_diffusionConstantsY, d_rx, d_ry,
                                                      d_state, simTime, d_errors,
                                                      d_sumMoments, d_sumStates);
        /*
        std::string fn = "diffusionConstants_" + oss.str() + ".dat";
        individualSolver->writeDiffusionConstants(fn,
                                d_diffusionConstantsX, d_diffusionConstantsY,
                                d_rx, d_ry);
        */
    } else {
        std::cout <<"Computing (non-individual) diffusion constants for species " << speciesIndex <<".\n"<<std::flush;
        if (m_hasFields)
            computeDiffusionConstants_forSpecies[speciesIndex](*d_diffusionConstantsX,
                                                               *d_diffusionConstantsY, *d_rx, *d_ry,
                                                               *d_state, *d_fields, simTime, *d_errors,
                                                               *d_sumMoments, *d_sumStates);
            else
            computeDiffusionConstants_forSpecies[speciesIndex](*d_diffusionConstantsX,
                                                               *d_diffusionConstantsY, *d_rx, *d_ry,
                                                               *d_state, simTime, *d_errors,
                                                               *d_sumMoments, *d_sumStates);
    }
}

Real InhomogeneousSolver::computeTimestepX(cl_int speciesIndex, Real simTime)
{
    // computes the individual time step - x
    computeIndividualTimestepX(*d_timeSteps, *d_minTimeBlock,
                               *d_diffusionConstantsX, *d_diffusionConstantsY,
                               *d_rx, *d_ry,
                               *d_transitionProbabilities, speciesIndex, simTime, *d_errors);

    // reduce over blocks
    reduceBlocks(*d_errors, *d_minTimeBlock, *d_globalMinTimestep);

    // and copy it back into time steps array
    oclCopy(*d_globalMinTimestep, &m_smallestDtX[speciesIndex], 1);

    return m_smallestDtX[speciesIndex];
}

Real InhomogeneousSolver::computeTimestepY(cl_int speciesIndex, Real simTime)
{
    // computes the individual time step - y
    computeIndividualTimestepY(*d_timeSteps, *d_minTimeBlock,
                               *d_diffusionConstantsX, *d_diffusionConstantsY,
                               *d_rx, *d_ry,
                               *d_transitionProbabilities, speciesIndex, simTime, *d_errors);

    // reduce over blocks
    reduceBlocks(*d_errors, *d_minTimeBlock, *d_globalMinTimestep);

    // and copy it back into time steps array
    oclCopy(*d_globalMinTimestep, &m_smallestDtY[speciesIndex], 1);

    return m_smallestDtY[speciesIndex];
}


} // namespace gpgmp
