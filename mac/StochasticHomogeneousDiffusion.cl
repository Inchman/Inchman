// Function Prototypes - Apple's OpenCL compiler will complain if they're not specified (Mac OS X 10.7.3)
__kernel void homogeneous_stochastic_diffuse(__global int *d_state, __global int4 *d_lostParticles,  int seed, int speciesIndex);
__kernel void homogeneous_stochastic_updateState(__global int *d_state, __global int4 *d_lostParticles, int speciesIndex);
__kernel void homogeneous_stochastic_source_boundary(__global int *d_state, int seed, int speciesIndex, int boundary, int source);


// Diffusion kernel
__kernel void homogeneous_stochastic_diffuse(__global int *d_state, __global int4 *d_lostParticles,
                                             int seed, int speciesIndex)
{

    // get absolute position on the device
    const size_t device_x = get_global_id(0);
    const size_t device_y = get_global_id(1);
    const int global_id = get_global_id_2d;
    const int state_offset = speciesIndex*get_global_area_2d;

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
    int my_state = d_state[global_id + state_offset];

    // this is a mask - decides which directions we can go from here
    // we use an int object, where 0 means -x, 1 means +x, 2 means
    // -y and 3 means +y
    size_t directionMask[4];//={1,1,1,1};
    
    // decide if we're at an outer boundary where we can't go everywhere
    // BoundaryMask* for that direction needs to be:
    // 0 (reflective boundary or source) TODO check out the source boundary rule
    // 1 (absorbing boundary, i.e. particle sink)
    // TODO Select operator should be quick but if thread divergence is high we can use bit shift..
    directionMask[0] = (device_x == 0) ? BoundaryMaskX1 : 1;
    directionMask[1] = (device_x == GridModelWidth-1) ? BoundaryMaskX2 : 1;   
    directionMask[2] = (device_y == 0) ? BoundaryMaskY1 : 1;
    directionMask[3] = (device_y == GridModelHeight-1) ? BoundaryMaskY2 : 1;
    
    // this counts the number of particles that wandered off into each direction
    // need to initialize individually to be compiler-friendly
    // TODO: Can we use a four vector here? Would be better..
    int lostParticles[4];
    lostParticles[0] = 0;
    lostParticles[1] = 0;
    lostParticles[2] = 0;
    lostParticles[3] = 0;
    
    // range values for picking directions
    // GPGMP_DIMENSIONALITY is 3 (for 2D slab)
    const Real pi = (1.-StayProbability)/(2.*GPGMP_DIMENSIONALITY);
    const Real p1 = StayProbability + pi;
    const Real p2 = p1 + pi;
    const Real p3 = p2 + pi;
    const Real p4 = p3 + pi;
    // We dont need p5 and p6 yet! (only in 3d)
    //const Real p5 = p4 + pi;
    //const Real p6 = p5 + pi;
    
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

      // pick direction - try to avoid conditionals
      size_t dir=3;

      // we deviate from the Rodriguez implementation now..
      dir -= (rand<=p3);
      dir -= (rand<=p2);
      dir -= (rand<=p1);
      
      // only lose the particle if it doesn't stay..
      lostParticles[dir] += directionMask[dir] * (rand > StayProbability && rand <= p4);
    }
    
    // done. We now save the lost particles for each thread on the device memory
    d_lostParticles[global_id].x = lostParticles[0];
    d_lostParticles[global_id].y = lostParticles[1];
    d_lostParticles[global_id].z = lostParticles[2];
    d_lostParticles[global_id].w = lostParticles[3];
    
    // and need to substract the ones we lost
    //if (my_state > 0) // this check is worth it, as it'll often be 0 for entire warps, saving many calculations. TODO: really??
    d_state[global_id + state_offset] -= (lostParticles[0] + lostParticles[1] + lostParticles[2] + lostParticles[3]);
}// diffusion kernel

__kernel void homogeneous_stochastic_updateState(__global int *d_state, __global int4 *d_lostParticles, int speciesIndex)
{
    // get absolute position on the device
    const size_t device_x = get_global_id(0);
    const size_t device_y = get_global_id(1);   
    
    // Now we need to update the global memory for each thread
    int state = d_state[get_global_id_2d + speciesIndex*get_global_area_2d];
    if (device_x > 0)
        state += d_lostParticles[(device_x-1)+device_y*GridModelWidth].y;
    if (device_x < GridModelWidth-1)
        state += d_lostParticles[(device_x+1)+device_y*GridModelWidth].x;
    if (device_y > 0)
        state += d_lostParticles[device_x+(device_y-1)*GridModelWidth].w;
    if (device_y < GridModelHeight-1)
        state += d_lostParticles[device_x+(device_y+1)*GridModelWidth].z;
    
    // and write it
    d_state[get_global_id_2d + speciesIndex*get_global_area_2d]=state;
} // update state kernel

// This kernel is used to implement source boundaries.
// It's a kernel on its own because we hardly ever need it ..
__kernel void homogeneous_stochastic_source_boundary(__global int *d_state, int seed,
                                                     int speciesIndex, int boundary, int source)
{
  // which boundary? 0: x=0, 1:x=xmax, 2:y=0, 3:y=ymax
  // absolute position on device - we only have 1D topology
  const size_t device = get_global_id(0);

    // Initialize counters and keys for RNG
  threefry2x32_key_t k = {{device, 0xdecafbad+seed}};
  threefry2x32_ctr_t c = {{0, 0xf00dcafe}};
  union {
    threefry2x32_ctr_t c;
    int2 r;
  } u;

  // initialize RNG counter
  c.v[0]=0;

  int global_id;
  switch (boundary) { // shouldn't give us any divergence
  case 0: // x=0
    global_id = ((device * GridModelWidth));
    break;
  case 1: // x=xmax
    global_id = ((GridModelWidth-1) + (device * GridModelWidth));
    break;
  case 2: // y=0
    global_id = (device);
    break;
  case 3: // y=ymax
    global_id = (get_global_id(0) + ((GridModelHeight-1) * GridModelWidth));
    break;
  }

  const int state_offset = speciesIndex*get_global_area_2d;

  const Real pi = (1.-StayProbability)/(2.*GPGMP_DIMENSIONALITY);
  const Real p1 = StayProbability + pi;

  for (int i=0; i<source; i++) {
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

    // only lose the particle if it doesn't stay..
    d_state[global_id + state_offset] += (rand > StayProbability && rand <= p1);
  }
} // homogeneous_stochastic_source_boundary

