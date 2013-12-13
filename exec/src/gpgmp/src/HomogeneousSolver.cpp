#include "HomogeneousSolver.h"
#include "Species.h"

#include <sstream>
#include <boost/log/trivial.hpp>

namespace gpgmp {
HomogeneousSolver::HomogeneousSolver(DiffusionModel &diffusionModel, std::vector<int> *runState, std::vector<Real> *deterministicState, SolverType solverType, bool deterministic, bool hasIndividualSpecies)
    : Solver(diffusionModel, runState, deterministicState, solverType, deterministic, hasIndividualSpecies),
      m_dt(diffusionModel.numSpecies())
{
    // initialize lost particles buffer
    if (!deterministic)
        d_lostParticles = oclCreateBuffer<cl_int4>(diffusionModel.gridArea());

    // Build programs
    std::ostringstream ss;
    oclGenerateHeader(ss);

    if (!deterministic) {
        oclReadAndReplaceTemplateFile(ss, "StochasticHomogeneousDiffusion.cl");
        m_diffusionProgram = oclBuildProgram(ss);
    } else {
        oclReadAndReplaceTemplateFile(ss, "Deterministic.cl");
        m_diffusionProgram = oclBuildProgram(ss);
    }

    // Create kernels
    if (!deterministic) {
        // diffusion kernel
        stochastic_diffuse = Kernel(m_diffusionProgram,
                                    "homogeneous_stochastic_diffuse",
                                    cl::EnqueueArgs(*m_commandQueue, m_globalRange));

        //optimizeKernelLocalSize(stochastic_diffuse);
        stochastic_diffuse.setLocalSize(m_localRange);

        // state update kernel
        //optimizeKernelLocalSizeAndMemory<cl_int4>(stochastic_updateState, 1, 2); // local mem area = (width+2)*(height+2)
        stochastic_updateState= Kernel(m_diffusionProgram,
                                       "homogeneous_stochastic_updateState",
                                       cl::EnqueueArgs(*m_commandQueue, m_globalRange));
        stochastic_updateState.setLocalSize(m_localRange);

        // for the source boundaries (if any)
        cl::NDRange globalRange1D(m_diffusionModel.gridWidth());
        stochastic_source_boundary= Kernel(m_diffusionProgram,
                                           "homogeneous_stochastic_source_boundary",
                                           cl::EnqueueArgs(*m_commandQueue, globalRange1D));
    } else {
        // Homogeneous diffusion kernel
        homogeneous_deterministic_diffuse = Kernel(m_diffusionProgram,
                                                   "homogeneous_deterministic_diffuse",
                                                   cl::EnqueueArgs(*m_commandQueue, m_globalRange));
        homogeneous_deterministic_diffuse.setLocalSize(m_localRange);
        homogeneous_deterministic_diffuse.setLocalMemory(sizeof(Real)*m_diffusionModel.species().size());

        // Reaction kernel for RK4
        if (m_diffusionModel.reactions().size()>0) {
            if (m_solverType==gpgmp::deterministic_homogeneous_RK4) {
                deterministic_performReaction_rk4 = Kernel(m_diffusionProgram,
                                                           "deterministic_performReaction_rk4",
                                                           cl::EnqueueArgs(*m_commandQueue, m_globalRange));
                deterministic_performReaction_rk4.setLocalSize(m_localRange);
                deterministic_performReaction_rk4.setLocalMemory(LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(Real)*5*m_diffusionModel.species().size());
            }
            else if (m_solverType==gpgmp::deterministic_homogeneous_aqss) {
                deterministic_performReaction_aqss = Kernel(m_diffusionProgram,
                                                            "deterministic_performReaction_aqss",
                                                            cl::EnqueueArgs(*m_commandQueue, m_globalRange));
                deterministic_performReaction_aqss.setLocalSize(m_localRange);
                deterministic_performReaction_aqss.setLocalMemory(LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(Real)*(3*m_diffusionModel.species().size()+1));
            }
        } // reaction kernels
    } // deterministic kernels

    // compute time steps (since this is the homogeneous solver, we only need to do that once)

    // physical cell size
    Real h = m_diffusionModel.physicalLength() / Real(m_diffusionModel.gridWidth());

    // compute the time steps
    if (!deterministic)
    {
        // our time step is chosen such that on average a fraction p0 of the particles
        // will be redistributed per time step.
        for (size_t i=0; i<m_diffusionModel.numSpecies(); i++) {
            if (m_diffusionModel.species()[i]->diffusionConstant()>0) {
                m_dt[i] = m_diffusionModel.p0()/((Real)2.*GMP_DIMENSIONALITY)*h*h/m_diffusionModel.species()[i]->diffusionConstant();
                BOOST_LOG_TRIVIAL(debug) <<"Species " << i <<" has a time step of "<< m_dt[i];
            } else {
                m_dt[i] = m_diffusionModel.maxTime();
                BOOST_LOG_TRIVIAL(debug) <<"Species " << i <<" is not diffusing.";
            }
        }
    }
    else
    {
        Real dtdiff;

        for (size_t i=0; i<m_diffusionModel.numSpecies(); i++) {
            if (m_diffusionModel.species()[i]->diffusionConstant()>0) {                
                m_dt[i]=(Real)1e-2f*h*h/(4.f*m_diffusionModel.species()[i]->diffusionConstant());
               BOOST_LOG_TRIVIAL(debug) << "Species " << i << " has a time step of " << m_dt[i] << ".";
            } else {
                m_dt[i] = m_diffusionModel.maxTime();
                BOOST_LOG_TRIVIAL(debug) << "Species " << i << " is not diffusing.";
            }
            // find smallest diffusion time step
            if (i==0) {
                dtdiff = m_dt[0];
            } else {
                if (m_dt[i] < dtdiff) dtdiff = m_dt[i];
            }
        }
        BOOST_LOG_TRIVIAL(debug) << "Diffusion time step is " << dtdiff <<".";

        // and we put it in for all species ..
        // todo: why do we need a common one for all species?
        m_dt.assign(m_diffusionModel.numSpecies(), dtdiff);
    }

    // and set it as next diffusion time as well
    m_nextDiffusionTime = m_dt;
}

HomogeneousSolver::~HomogeneousSolver()
{
    delete d_lostParticles;
}

void HomogeneousSolver::diffuse(int speciesIndex, Real currentTime, int &seed)
{
    // diffuse species
    if (!m_deterministic) {
        stochastic_diffuse(*d_state, *d_lostParticles, seed , speciesIndex);
        stochastic_updateState(*d_state, *d_lostParticles, speciesIndex);

        // if source boundaries present, we need to update again
        int boundary;
        if (m_diffusionModel.sourceMask().s[0]>0) {
            boundary = 0;
            stochastic_source_boundary(*d_state, seed, speciesIndex, boundary, m_diffusionModel.sourceMask().s[1]);
        }
        if (m_diffusionModel.sourceMask().s[1]>0) {
            boundary = 1;
            stochastic_source_boundary(*d_state, seed, speciesIndex, boundary, m_diffusionModel.sourceMask().s[1]);
        }
        if (m_diffusionModel.sourceMask().s[2]>0) {
            boundary = 2;
            stochastic_source_boundary(*d_state, seed, speciesIndex, boundary, m_diffusionModel.sourceMask().s[2]);
        }
        if (m_diffusionModel.sourceMask().s[3]>0) {
            boundary=3;
            stochastic_source_boundary(*d_state, seed, speciesIndex, boundary, m_diffusionModel.sourceMask().s[3]);
        }
    } else {
        // todo: implement source boundaries for deterministic here
        homogeneous_deterministic_diffuse(*d_state, m_dt[speciesIndex]);
    }

    // and update the next diffusion time
    m_nextDiffusionTime[speciesIndex] += m_dt[speciesIndex];
}


} // namespace gpgmp
