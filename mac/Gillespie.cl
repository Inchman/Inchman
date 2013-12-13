// main kernel to perform a series of Gillespie steps
// this is a template for a dynamically generated Gillespie kernel
#if (GPGMP_NUM_REACTIONS > 0)
__kernel void performGillespie(__global int *_globalState,
                               int _seed, float _deltat,

#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
                               // we only need this if there are localized reactions in the model
                               image3d_t _reactionMask,
#endif
#ifdef GPGMP_DEBUG_KERNEL
			       __global float *_error,
#endif
                               __local int *_localState)
{
    
    // IMPORTANT: All varaibles are prefixed with "_", so as to mark them as "internal",
    //            as the user's variable names will be used directly as given.

    // get local and global ids
    int _localId = get_local_id(0) + get_local_id(1) * get_local_size(0); // LOCAL x + LOCAL y * LOCAL width
    int GridX = get_global_id(0); // GLOBAL x
    int GridY = get_global_id(1); // GLOBAL y
    int _globalId = GridX + GridY * get_global_size(0); // GLOBAL x + GLOBAL y * GLOBAL width

    // offset for species
    const int _speciesOffsetGlobal = get_global_area_2d;
    const int _speciesOffsetLocal = get_local_area_2d;

    // and copy for each species
    // poor coalescence for some reason..
    // we need to put species count in shared mem not registers..
    // otherwise we're gonna run out of registers!

// <<<! defineUserParameters !>>>
    
// <<<! copyStateToLocal !>>>
// <<<! defineSpeciesFromLocal !>>>

    // Initialize counters and keys .. note the funny names for init!
    threefry2x32_key_t _rng_k = {{_globalId, 0xdecafbad + _seed}};
    threefry2x32_ctr_t _rng_c = {{0, 0xf00dcafe}};

    // Get a random number
    union {
      threefry2x32_ctr_t c;
      int2 r;
    } _rng_u;

    // initialize RNG counter - note: we can't seed here!
    _rng_c.v[0]=0;

    // total propensity
    float _h[GPGMP_NUM_REACTIONS+1]; // local propensities
    float _time = 0.;

    // main loop
    while (_time < _deltat) {
      // compute reaction hazards and h0

// <<<! computeReactionHazards !>>>

      // get out if no reactions left to do
      if (_h[0] == 0) break;

      // draw two random numbers
      // increment counter .. we have to do that every time we draw
      _rng_c.v[0]++;
      _rng_u.c = threefry2x32(_rng_c, _rng_k);

      // this give us random numbers between 0 and 1
      float _r1 = _rng_u.r.x*M_RAN_INVM32+0.5;
      float _r2 = (_rng_u.r.y*M_RAN_INVM32+0.5) * _h[0]; // needs to be multiplied with total propensity

      // compute new time step
      float _dt = log(1.f/_r1) / _h[0];
      _time += _dt;

      // only perform reaction if time step is smaller than diffusion time step
      if (_time < _deltat) {

// <<<! performReactions !>>>

      }//if (time < deltat)

#ifdef GPGMP_DEBUG_KERNEL
      // check if all species are positive
      _error[get_global_id_2d * GPGMP_NUM_ERRORS] = 0;
      for (uint ii=0; ii<NumSpecies; ii++) {
	if ((_localState[_localId + ii * _speciesOffsetLocal]) <0) {
	  // set global error value for this thread
	  _error[get_global_id_2d * GPGMP_NUM_ERRORS] = 1;
	  _error[get_global_id_2d * GPGMP_NUM_ERRORS + 1] = ii;
	  // the index of the reaction should be in _error[2] now .. 
	  return;
	}
      }
#endif
    }// main loop

    // Copy state back to device

// <<<! undefSpeciesFromLocal !>>>
// <<<! copyLocalToState !>>>

}
#endif // (GPGMP_NUM_REACTIONS > 0)

__kernel void updateFields( Real dt,
                            Real PhysicalSimTime,
                            __global Real *_field,
                            __global GPGMP_STATE_TYPE *_state)
{
// <<<! updateFields !>>>
}
