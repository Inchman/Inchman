#include "helpers.cl"

// Would be nice to have templates..
// kernel to compute the moments
__kernel void inhomogeneous_computeMomentsX(__global GPGMP_STATE_TYPE *d_state,
                                            __global Real *d_sumMomentsBlocks,
                                            __global Real *d_sumStatesBlocks,
                                            int speciesIndex/*, __global Real *d_errors*/,
                                            __local Real *s_local)
{
    // we compute the moment with respect to the cell coordinate
    Real xc = (Real)get_global_id(0);

    // compute moment and #particles in cell
    const size_t local_id = get_local_id_2d;
    s_local[local_id] = xc * d_state[get_global_id_for_species_2d(speciesIndex)];
    s_local[local_id+get_local_area_2d] = d_state[get_global_id_for_species_2d(speciesIndex)];

    // add up all moments for this block - need another kernel to reduce
    // over all blocks
    // Up-sweep phase (see Harris et al. in GPU Gems)
    int offset = 1;

    for (size_t d = get_local_area_2d>>1; d > 0; d >>= 1)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (local_id < d) {
            int ai = offset*(2*local_id+1)-1;
            int bi = offset*(2*local_id+2)-1;
            s_local[bi] += s_local[ai];
            s_local[bi+get_local_area_2d] += s_local[ai+get_local_area_2d];
        }
        offset *= 2;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // now the sum for this block/group should be in the last element
    if (local_id == get_local_area_2d-1) {
        d_sumMomentsBlocks[get_block_id_2d] = s_local[local_id];
        d_sumStatesBlocks[get_block_id_2d] = s_local[local_id+get_local_area_2d];
    }
} // inhomogeneous_compueMomentsX

__kernel void inhomogeneous_reduceMomentsBlocks(__global Real *d_sumMomentsBlocks, __global Real *d_sumStatesBlocks,
                                                __global Real *d_sumMoments, __global Real *d_sumStates,
                                                int speciesIndex, __global Real *d_errors, __local Real *s_local)
{
    // shared data worker array
    // TODO:  local ID is global here.  change..
    //size_t local_id = get_local_id(0) + get_local_id(1)*get_local_size(0);
    size_t local_id = get_local_id(0);

    // copy to shared array
    s_local[local_id] = d_sumMomentsBlocks[local_id];
    s_local[local_id+get_local_size(0)] = d_sumStatesBlocks[local_id];

    // Up-sweep phase (see Harris et al. in GPU Gems)
    int offset = 1;
//    for (size_t d = get_global_area_2d>>1; d > 0; d >>= 1)
    for (size_t d = get_local_size(0)>>1; d > 0; d >>= 1)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (local_id < d)
        {
            int ai = offset*(2*local_id+1)-1;
            int bi = offset*(2*local_id+2)-1;
            s_local[bi] += s_local[ai];
            s_local[bi+get_local_size(0)] += s_local[ai+get_local_size(0)];
        }
        offset *= 2;
    }

    // and write back to device
    //if (local_id == get_global_area_2d-1) {
    if (local_id == get_local_size(0)-1) {
//        d_sumMoments[speciesIndex] = s_local[get_global_area_2d-1];
//        d_sumStates[speciesIndex] = s_local[get_global_area_2d-1+get_global_area_2d];
        d_sumMoments[speciesIndex] = s_local[get_local_size(0)-1];
        d_sumStates[speciesIndex] = s_local[get_local_size(0)-1+get_local_size(0)];
    }
}//inhomogeneousReduceMomentsBlock



// <<<! inhomogeneous_computeDiffusionConstants_forSpecies_ALL !>>>

__kernel void inhomogeneous_computeIndividualTimestepX(
        __global Real *d_timeSteps,
        __global Real *d_minTimeBlock,
        __global Real *d_diffusionConstantsX,
        __global Real *d_diffusionConstantsY,
        __global Real *d_rx, __global Real *d_ry,
        __global float4 *d_transitionProbabilities,
        int speciesIndex, Real simTime,
        __global Real *d_errors, __local Real *s_local)
{
    // get absolute position on the device
//    const size_t device_x = get_global_id(0);
//    const size_t device_y = get_global_id(1);
    const size_t index = get_global_id_for_species_2d(speciesIndex);

    const size_t local_id = get_local_id_2d;

    // diffusion constant and drift field in this cell
    Real diffHere, rx;

    const Real deltax = PhysicalCellWidth;
    const Real deltay = PhysicalCellHeight;
#pragma unused(deltax)
#pragma unused(deltay)

    // set error code to good - todo: Make optional
    // d_errors[0] = 0;

    diffHere = d_diffusionConstantsX[index];
    rx = d_rx[index];

    // compute time step for this cell
    // TODO Here we assume that |r/h| < D/h^2 .. not necessarily true!
    // But I don't know how to handle the other case reasonably..
    Real timeStep;

#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS_FPE

    timeStep = (Real)StayProbability/6./
            (diffHere/(deltax*deltax)+fabs(rx/(2.*deltax)));

#endif
#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS

    if (rx==0) {
        timeStep = deltax*deltax/(2.*diffHere)*(1.-2./3.);
    } else {
      if (diffHere == 0) {
	timeStep = deltax/fabs(rx);
      } else {
        Real epsx = deltax*rx/(2.*diffHere);
        Real temp = (1.+(-1./(tanh(epsx)*epsx)+1./(sinh(epsx)*sinh(epsx))));
        // TODO: ok this is a hack.
        // if the bracketed expression is <=0 we assume eps~0 so take
        // pure diffusion time step. There must be a better
        // way to write the expression for the time step to
        // avoid that. For now we go with that.
        if (temp <= 0)
            temp = 1./3.;

        timeStep = deltax/rx * tanh(epsx)*temp;
      }
    } // zero field

#endif  
#ifdef GPGMP_SOLVER_DETERMINISTIC_INHOMOGENEOUS_RK4
    // TODO: We can probably use a different time step here..
    if (rx==0) {
        timeStep = deltax*deltax/(2.*diffHere)*(1.-2./3.);
    } else {
        Real epsx = deltax*rx/(2.*diffHere);
        Real temp = (1.+(-1./(tanh(epsx)*epsx)+1./(sinh(epsx)*sinh(epsx))));
        // TODO: ok this is a hack.
        // if the bracketed expression is <=0 we assume eps~0 so take
        // pure diffusion time step. There must be a better
        // way to write the expression for the time step to
        // avoid that. For now we go with that.
        if (temp <= 0)
            temp = 1./3.;

        timeStep = 0.1*deltax/rx * tanh(epsx)*temp;
    } // zero field
#endif
#ifdef GPGMP_SOLVER_DETERMINISTIC_INHOMOGENEOUS_AQSS
    // TODO: We can probably use a different time step here..
    if (rx==0) {
        timeStep = deltax*deltax/(2.*diffHere)*(1.-2./3.);
    } else {
        Real epsx = deltax*rx/(2.*diffHere);
        Real temp = (1.+(-1./(tanh(epsx)*epsx)+1./(sinh(epsx)*sinh(epsx))));
        // TODO: ok this is a hack.
        // if the bracketed expression is <=0 we assume eps~0 so take
        // pure diffusion time step. There must be a better
        // way to write the expression for the time step to
        // avoid that. For now we go with that.
        if (temp <= 0)
            temp = 1./3.;

        timeStep = deltax/rx * tanh(epsx)*temp;
    } // zero field

#endif

    // and save the smaller one to the worker array
    s_local[local_id] = timeStep;

    // compute transition rates .. do that only for discrete FPE algorithm
    // for the CA we do it later (since we need the smallest time step)

#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS_FPE

    // to go to the right (i+1)
    d_transitionProbabilities[index].x = diffHere/(deltax*deltax) + rx/(2.*deltax);
    // to go to the left (i-1)
    d_transitionProbabilities[index].y = diffHere/(deltax*deltax) - rx/(2.*deltax);

#endif

    // Find minimum time step of this block (will be reduced over all
    // blocks in the reduceBlocks kernel)

    // Up-sweep phase (see Harris et al. in GPU Gems)
    int offset = 1;

    for (size_t d = get_local_area_2d>>1; d > 0; d >>= 1)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (local_id < d) {
            int ai = offset*(2*local_id+1)-1;
            int bi = offset*(2*local_id+2)-1;
            s_local[bi] = min(s_local[ai], s_local[bi]);
        }
        offset *= 2;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // now the minimum for this block/group should be in the last element
    if (local_id==get_local_area_2d-1) {
        d_minTimeBlock[get_block_id_2d] = s_local[local_id];
    }

} // computeIndividualTimestepX

__kernel void inhomogeneous_computeIndividualTimestepY(
        __global Real *d_timeSteps, __global Real *d_minTimeBlock,
        __global Real *d_diffusionConstantsX,
        __global Real *d_diffusionConstantsY,
        __global Real *d_rx, __global Real *d_ry,
        __global float4 *d_transitionProbabilities,
        int speciesIndex, Real simTime,
        __global Real *d_errors, __local Real *s_local)
{
    // get absolute position on the device
//    const size_t device_x = get_global_id(0);
//    const size_t device_y = get_global_id(1);
    const size_t index = get_global_id_for_species_2d(speciesIndex);

    const size_t local_id = get_local_id_2d;

    // diffusion constant and drift field in this cell
    // hopefully compiler optimizes it
    // TODO: should only need deltay here .. but check
    const Real deltax = PhysicalCellWidth;
    const Real deltay = PhysicalCellHeight;
#pragma unused(deltax)
#pragma unused(deltay)

    // set error code to good - todo: Make optional
    // d_errors[0] = 0;

    // this one for y sweeps
    Real diffHere, rx;

    diffHere = d_diffusionConstantsY[index];
    rx = d_ry[index];

    // compute time step for this cell
    // TODO Here we assume that |r/h| < D/h^2 .. not necessarily true!
    // But I don't know how to handle the other case reasonably..
    Real timeStep;

#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS_FPE

    timeStep = (Real)StayProbability/6./
            (diffHere/(deltax*deltax)+fabs(rx/(2.*deltax)));

#endif
#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS

    if (rx==0) {
        timeStep = deltax*deltax/(2.*diffHere)*(1.-2./3.);
    } else {
      if (diffHere == 0) {
	timeStep = deltax/fabs(rx);
      } else {
        Real epsx = deltax*rx/(2.*diffHere);
        Real temp = (1.+(-1./(tanh(epsx)*epsx)+1./(sinh(epsx)*sinh(epsx))));
        // TODO: ok this is a hack.
        // if the bracketed expression is <=0 we assume eps~0 so take
        // pure diffusion time step. There must be a better
        // way to write the expression for the time step to
        // avoid that. For now we go with that.
        if (temp <= 0)
            temp = 1./3.;

        timeStep = deltax/rx * tanh(epsx)*temp;
      }
    } // zero field

#endif  
#ifdef GPGMP_SOLVER_DETERMINISTIC_INHOMOGENEOUS_RK4
    // todo: use a better time step for the deterministic solvers!
    if (rx==0) {
        timeStep = deltax*deltax/(2.*diffHere)*(1.-2./3.);
    } else {
        Real epsx = deltax*rx/(2.*diffHere);
        Real temp = (1.+(-1./(tanh(epsx)*epsx)+1./(sinh(epsx)*sinh(epsx))));
        // TODO: ok this is a hack.
        // if the bracketed expression is <=0 we assume eps~0 so take
        // pure diffusion time step. There must be a better
        // way to write the expression for the time step to
        // avoid that. For now we go with that.
        if (temp <= 0)
            temp = 1./3.;

        timeStep = 0.1*deltax/rx * tanh(epsx)*temp;
    } // zero field
#endif
#ifdef GPGMP_SOLVER_DETERMINISTIC_INHOMOGENEOUS_AQSS
    // todo: use a better time step for the deterministic solvers!
    if (rx==0) {
        timeStep = deltax*deltax/(2.*diffHere)*(1.-2./3.);
    } else {
        Real epsx = deltax*rx/(2.*diffHere);
        Real temp = (1.+(-1./(tanh(epsx)*epsx)+1./(sinh(epsx)*sinh(epsx))));
        // TODO: ok this is a hack.
        // if the bracketed expression is <=0 we assume eps~0 so take
        // pure diffusion time step. There must be a better
        // way to write the expression for the time step to
        // avoid that. For now we go with that.
        if (temp <= 0)
            temp = 1./3.;

        timeStep = deltax/rx * tanh(epsx)*temp;
    } // zero field
#endif

    // and save the smaller one to the worker array
    s_local[local_id] = timeStep;

    // compute transition rates .. do that only for discrete FPE algorithm
    // for the CA we do it later (since we need the smallest time step)

#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS_FPE

    // to go up (j+1)
    d_transitionProbabilities[index].z = diffHere/(deltay*deltay) + rx/(2.*deltay);
    // to go down (j-1)
    d_transitionProbabilities[index].w = diffHere/(deltay*deltay) - rx/(2.*deltay);

#endif

    // Find minimum time step of this block (will be reduced over all
    // blocks in the reduceBlocks kernel)

    // Up-sweep phase (see Harris et al. in GPU Gems)
    int offset = 1;

    for (size_t d = get_local_area_2d>>1; d > 0; d >>= 1)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (local_id < d) {
            int ai = offset*(2*local_id+1)-1;
            int bi = offset*(2*local_id+2)-1;
            s_local[bi] = min(s_local[ai], s_local[bi]);
        }
        offset *= 2;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // now the minimum for this block/group should be in the last element
    if (local_id==get_local_area_2d-1) {
        d_minTimeBlock[get_block_id_2d] = s_local[local_id];
    }

    if (index==0) {
      d_errors[0]=0.5;
      d_errors[1]=d_minTimeBlock[get_block_id_2d];
    }
} // computeIndividualTimestepY


__kernel void inhomogeneous_reduceBlocks(__global Real *d_errors,
        __global Real *d_minTimeBlock,
        __global Real *d_minTimeStep,
        __local Real *s_local)
{    

    // shared data worker array
    // we assume that we only have a 1D geometry (global_size=local_size)
    // that covers all blocks!
    // todo: Change that to 2D but in a way that works.
    //__local Real s_local[inhomogeneous_reduceBlocks_local_area];

    //size_t local_id = get_local_id(0) + get_local_id(1)*get_local_size(0);
    size_t local_id = get_local_id(0);

    // copy to shared array
    s_local[local_id] = d_minTimeBlock[local_id];

    // Up-sweep phase (see Harris et al. in GPU Gems)
    int offset = 1;
//    for (size_t d = get_global_area_2d>>1; d > 0; d >>= 1)
    for (size_t d = get_local_size(0)>>1; d > 0; d >>= 1)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (local_id < d) {
            int ai = offset*(2*local_id+1)-1;
            int bi = offset*(2*local_id+2)-1;
            s_local[bi] = min(s_local[ai], s_local[bi]);
        }
        offset *= 2;
    }

    // and write back to device
//    if (local_id == get_global_area_2d-1) {
    if (local_id == get_local_size(0)-1) {
//        *d_minTimeStep = s_local[get_global_area_2d-1];
        *d_minTimeStep = s_local[get_local_size(0)-1];
        //		printf("Thread %d finds smallest time step %g.\n",
        //				local_id, *d_minTimeStep);
    }

} // reduceBlocks

// main kernel to implement diffusion
__kernel void inhomogeneous_diffuseX(__global int *d_state,
        __global float4 *d_transitionProbabilities,
        __global int4 *d_lostParticles, int seed,
        int speciesIndex, Real dt, __global Real *d_errors,
        __global Real *d_rx, __global Real *d_ry,
        __global Real *d_diffusionConstantsX,
        __global Real *d_diffusionConstantsY)
{
    // TODO: we can probably optimize a bit better now that we do separate
    // x/y sweeps. But we leave that for now.

    const Real deltax = PhysicalCellWidth;
    const Real deltay = PhysicalCellHeight;
#pragma unused(deltax)
#pragma unused(deltay)

    // set good error code
    // d_errors[0] = 0;

    // get absolute position on the device
    const size_t device_x = get_global_id(0);
    const size_t device_y = get_global_id(1);   
    const size_t global_id = get_global_id_for_species_2d(speciesIndex);

    // Initialize counters and keys for RNG
    threefry2x32_key_t k = {{global_id, 0xdecafbad+seed}};
    threefry2x32_ctr_t c = {{0, 0xf00dcafe}};
    union {
        threefry2x32_ctr_t c;
        int2 r;
    } u;

    // initialize RNG counter - note: we can't seed here!
    c.v[0]=0;

    // fill in my species number
    int my_state = d_state[global_id];

    // this is a mask - decides which directions we can go from here
    // we use an int object, where 0 means -x, 1 means +x, 2 means
    // -y and 3 means +y
    size_t directionMask[4];

    // decide if we're at an outer boundary where we can't go everywhere
    // boundaryMask for that direction needs to be:
    // 0 (reflective boundary or source) TODO check out the source boundary rule
    // 1 (absorbing boundary, i.e. particle sink)
    // TODO maybe add extra kernel for every BC to avoid conditionals
    directionMask[0] = (device_x == 0) ? BoundaryMaskX1 : 1;
    directionMask[1] = (device_x == GridModelWidth-1) ? BoundaryMaskX2 : 1;
    directionMask[2] = (device_y == 0) ? BoundaryMaskY1 : 1;
    directionMask[3] = (device_y == GridModelHeight-1) ? BoundaryMaskY2 : 1;    

    // this counts the number of particles that wandered off into each direction
    // TODO: Can we use a four vector here? Would be better..
    int lostParticles[4];
    lostParticles[0] = 0;
    lostParticles[1] = 0;
    lostParticles[2] = 0;
    lostParticles[3] = 0;

    // values for the direction picking.

    // we need only two since we do two independent sweeps!
    // p1 is for right/up, p2 for left/down
    Real p1, p2, ps;

#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS_FPE    

    // fetch transition probabilites
    float4 transProbs = d_transitionProbabilities[get_global_id_for_species_2d(speciesIndex)];

    // go right (plus x)
    p1 = transProbs.x*dt;
    // go left (minus x)
    p2 = p1 + transProbs.y*dt;

#endif
#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS

    // fetch drift field and diffusion constant
    Real diffHere, rx, dx;
    diffHere =	d_diffusionConstantsX[get_global_id_for_species_2d(speciesIndex)];
    rx = d_rx[get_global_id_for_species_2d(speciesIndex)];
    dx = deltax;

    // this one for y sweep
    //    diffHere =	d_diffusionConstantsY[get_global_id_for_species_2d(speciesIndex)];
    //    rx = d_ry[get_global_id_for_species_2d(speciesIndex)];
    //    dx = deltay;

    // compute probabilities
    ps = 1. - (dt * (2.*diffHere + rx*rx*dt))/(dx*dx);
    p1 = (1.-ps)*(1. + (dx*rx)/(2.*diffHere + rx*rx*dt))/2.;
    p2 = p1 + (1.-ps)*(1. - (dx*rx)/(2.*diffHere + rx*rx*dt))/2.;
#endif       

    // main loop. For every particle decide if it stays in the cell or moves.
    for (int i=0; i<my_state; i++)
    {
        // draw two random numbers
        // We only need to draw every second round since we only need one number
        // shouldn't give any divergence I think.. but check!
        float rand;
        if (! (i%2)) {
            c.v[0]++;
            u.c = threefry2x32(c, k);
            rand = u.r.x*M_RAN_INVM32+0.5;
        } else {
            rand = (u.r.y*M_RAN_INVM32+0.5);
        }

        // which direction?
        size_t dir=0;
        dir += (rand<=p1);

        // this one for y
        //if (!sweepX) dir+=2; // for y sweeps

        // only lose it if it doesn't stay
        lostParticles[dir] += directionMask[dir] * (rand<=p2);
    }

    // done. We now save the lost particles for each thread on the device memory
    d_lostParticles[get_global_id_2d].x = lostParticles[0];
    d_lostParticles[get_global_id_2d].y = lostParticles[1];
    d_lostParticles[get_global_id_2d].z = lostParticles[2];
    d_lostParticles[get_global_id_2d].w = lostParticles[3];

    // and need to substract the ones we lost
    d_state[get_global_id_for_species_2d(speciesIndex)] -= (lostParticles[0] + lostParticles[1] + lostParticles[2] + lostParticles[3]);
} // diffuseX

__kernel void inhomogeneous_diffuseY(__global int *d_state,
        __global float4 *d_transitionProbabilities,
        __global int4 *d_lostParticles, int seed,
        int speciesIndex, Real dt, __global Real *d_errors,
        __global Real *d_rx, __global Real *d_ry,
        __global Real *d_diffusionConstantsX,
        __global Real *d_diffusionConstantsY)
{
    // TODO: we can probably optimize a bit better now that we do separate
    // x/y sweeps. But we leave that for now.

    const Real deltax = PhysicalCellWidth;
    const Real deltay = PhysicalCellHeight;
#pragma unused(deltax)
#pragma unused(deltay)

    // set good error code
    // d_errors[0] = 0;

    // get absolute position on the device
    const size_t device_x = get_global_id(0);
    const size_t device_y = get_global_id(1);   
    const size_t global_id = get_global_id_for_species_2d(speciesIndex);

    // Initialize counters and keys for RNG
    threefry2x32_key_t k = {{global_id, 0xdecafbad+seed}};
    threefry2x32_ctr_t c = {{0, 0xf00dcafe}};
    union {
        threefry2x32_ctr_t c;
        int2 r;
    } u;

    // initialize RNG counter - note: we can't seed here!
    c.v[0]=0;

    // fill in my species number
    int my_state = d_state[global_id];

    // this is a mask - decides which directions we can go from here
    // we use an int object, where 0 means -x, 1 means +x, 2 means
    // -y and 3 means +y
    size_t directionMask[4];

    // decide if we're at an outer boundary where we can't go everywhere
    // boundaryMask for that direction needs to be:
    // 0 (reflective boundary or source) TODO check out the source boundary rule
    // 1 (absorbing boundary, i.e. particle sink)
    // TODO maybe add extra kernel for every BC to avoid conditionals
    directionMask[0] = (device_x == 0) ? BoundaryMaskX1 : 1;
    directionMask[1] = (device_x == GridModelWidth-1) ? BoundaryMaskX2 : 1;
    directionMask[2] = (device_y == 0) ? BoundaryMaskY1 : 1;
    directionMask[3] = (device_y == GridModelHeight-1) ? BoundaryMaskY2 : 1;    

    // this counts the number of particles that wandered off into each direction
    // TODO: Can we use a four vector here? Would be better..
    int lostParticles[4];
    lostParticles[0] = 0;
    lostParticles[1] = 0;
    lostParticles[2] = 0;
    lostParticles[3] = 0;

    // values for the direction picking.

    // we need only two since we do two independent sweeps!
    // p1 is for right/up, p2 for left/down
    Real p1, p2, ps;

#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS_FPE    

    // fetch transition probabilites
    float4 transProbs = d_transitionProbabilities[get_global_id_for_species_2d(speciesIndex)];

    // go up (plus y)
    p1 = transProbs.z*dt;

    // go down (minus y)
    p2 = p1 + transProbs.w*dt;
#endif

#ifdef GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS

    // fetch drift field and diffusion constant
    Real diffHere, rx, dx;

    diffHere =	d_diffusionConstantsY[get_global_id_for_species_2d(speciesIndex)];
    rx = d_ry[get_global_id_for_species_2d(speciesIndex)];
    dx = deltay;

    // compute probabilities
    ps = 1. - (dt * (2.*diffHere + rx*rx*dt))/(dx*dx);
    p1 = (1.-ps)*(1. + (dx*rx)/(2.*diffHere + rx*rx*dt))/2.;
    p2 = p1 + (1.-ps)*(1. - (dx*rx)/(2.*diffHere + rx*rx*dt))/2.;
#endif       

    // main loop. For every particle decide if it stays in the cell or moves.
    for (int i=0; i<my_state; i++)
    {
        // draw two random numbers
        // We only need to draw every second round since we only need one number
        // shouldn't give any divergence I think.. but check!
        float rand;
        if (! (i%2)) {
            c.v[0]++;
            u.c = threefry2x32(c, k);
            rand = u.r.x*M_RAN_INVM32+0.5;
        } else {
            rand = (u.r.y*M_RAN_INVM32+0.5);
        }

        // which direction?
        size_t dir=2;
        dir += (rand<=p1);

        // only lose it if it doesn't stay
        lostParticles[dir] += directionMask[dir] * (rand<=p2);
    }

    // done. We now save the lost particles for each thread on the device memory
    d_lostParticles[get_global_id_2d].x = lostParticles[0];
    d_lostParticles[get_global_id_2d].y = lostParticles[1];
    d_lostParticles[get_global_id_2d].z = lostParticles[2];
    d_lostParticles[get_global_id_2d].w = lostParticles[3];

    // and need to substract the ones we lost
    d_state[get_global_id_for_species_2d(speciesIndex)] -=
            (lostParticles[0] + lostParticles[1] + lostParticles[2] + lostParticles[3]);
} // diffuseY

__kernel void inhomogeneous_updateStateX(__global int *d_state,
                                         __global int4 *d_lostParticles, int speciesIndex)
{
    // get absolute position on the device
    const size_t device_x = get_global_id(0);
    const size_t device_y = get_global_id(1);

    // Now we need to update the global memory for each thread
    int state = d_state[get_global_id_2d + speciesIndex*get_global_area_2d];

    if (device_x > 0)
        state += d_lostParticles[(device_x-1)+device_y*get_global_size(0)].y;
    if (device_x < get_global_size(0)-1)
        state += d_lostParticles[(device_x+1)+device_y*get_global_size(0)].x;

    // and write it
    d_state[get_global_id_2d + speciesIndex*get_global_area_2d] = state;
}

__kernel void inhomogeneous_updateStateY(__global int *d_state,
        __global int4 *d_lostParticles, int speciesIndex)
{
    // get absolute position on the device
    const size_t device_x = get_global_id(0);
    const size_t device_y = get_global_id(1);

    // Now we need to update the global memory for each thread
    int state = d_state[get_global_id_2d + speciesIndex*get_global_area_2d];

    if (device_y > 0)
        state += d_lostParticles[device_x+(device_y-1)*get_global_size(0)].w;
    if (device_y < get_global_size(1)-1)
        state += d_lostParticles[device_x+(device_y+1)*get_global_size(0)].z;

    // and write it
    d_state[get_global_id_2d + speciesIndex*get_global_area_2d] = state;
}
