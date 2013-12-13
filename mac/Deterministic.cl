// Function Prototypes - Apple's OpenCL compiler will complain if they're not specified (Mac OS X 10.7.3)
inline void deterministic_rhs_ode_rk4(__local Real *_localState, __local Real *_deriv,
                                      Real _offset, __local Real *_rk
#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
                                      // we only need this if there are localized reactions in the model
                                      , image3d_t _reactionMask
#endif
                                      );
inline void deterministic_rhs_ode_aqss(__local Real *_localState, __local Real *_q, __local Real *_d
#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
                                       , image3d_t _reactionMask
#endif
                                       );


// kernel to compute Laplacian and diffuse (homogeneous, no drift)
__kernel void homogeneous_deterministic_diffuse(__global Real *d_state, Real dtdiff,
                                                __local Real *_physicalDiffusionConstants)
{
    const Real dx2 = PhysicalCellWidth * PhysicalCellWidth;

    // we put the physical diffusion constants into shared memory
    // TODO: that might not be that great - all threads trying to write to same memory?

// <<<! deterministicStoreDiffusionConstants !>>>


    // get absolute position on the device
    size_t dx = get_global_id(0);
    size_t dy = get_global_id(1);

    // compute laplacian using registers
    for (int i=0; i<GPGMP_NUM_SPECIES; i++) {

        // get neighbours - we use mirror BCs

        Real sleft = (dx > 0
                      ? d_state[(dx-1)+dy*GridModelWidth+i*get_global_area_2d]
                      : d_state[(dx+1)+dy*GridModelWidth+i*get_global_area_2d]);

        Real sright = (dx < GridModelWidth-1
                       ? d_state[(dx+1)+dy*GridModelWidth+i*get_global_area_2d]
                       : d_state[(dx-1)+dy*GridModelWidth+i*get_global_area_2d]);

        Real stop = (dy > 0
                     ? d_state[dx+(dy-1)*GridModelWidth+i*get_global_area_2d]
                     : d_state[dx+(dy+1)*GridModelWidth+i*get_global_area_2d]);

        Real sbottom = (dy < GridModelWidth-1
                        ? d_state[dx+(dy+1)*GridModelWidth+i*get_global_area_2d]
                        : d_state[dx+(dy-1)*GridModelWidth+i*get_global_area_2d]);

        Real shere = d_state[dx+dy*GridModelWidth+i*get_global_area_2d];

        // compute laplacian
        Real slaplace = stop + sbottom + sleft + sright - 4.*shere;

        // diffusion constant
        Real diff = _physicalDiffusionConstants[i];

        // and store the state array
        d_state[dx+dy*GridModelWidth+i*get_global_area_2d] = shere + dtdiff*diff*slaplace/dx2;
    }

} // deterministic diffuse kernel (homogeneous)

// Kernel to compute Laplacian and diffuse in x direction (inhomogeneous with drift)
// TODO: The accuracy of the deterministic drift-diffusion solver is pretty bad..
// TODO: we'd probably need a better algorithm for the drift-diffusion part!
// i denotes the species index to diffuse
__kernel void inhomogeneous_deterministic_diffuseX(__global Real *d_state, Real dtdiff,
                                                   __global Real *d_diffx,
                                                   __global Real *d_rx,
                                                   int i)
{
    const Real dx2 = PhysicalCellWidth * PhysicalCellWidth;

    // get absolute position on the device
    size_t dx = get_global_id(0);
    size_t dy = get_global_id(1);

    // get neighbours - we use mirror BCs (zero-gradient, outflow)
    // todo: allow reflective BCs
    int indexleft   = (dx > 0 ? (dx-1)+dy*GridModelWidth+i*get_global_area_2d : (dx+1)+dy*GridModelWidth+i*get_global_area_2d);
    int indexright  = (dx < GridModelWidth-1 ? (dx+1)+dy*GridModelWidth+i*get_global_area_2d : (dx-1)+dy*GridModelWidth+i*get_global_area_2d);

    // compute laplacian (centered difference, second order)
    Real sleft     = d_state[indexleft]*d_diffx[indexleft];
    Real sright    = d_state[indexright]*d_diffx[indexright];
    Real sherex    = d_state[dx+dy*GridModelWidth+i*get_global_area_2d]*d_diffx[dx+dy*GridModelWidth+i*get_global_area_2d];
    Real slaplace  =  sleft + sright - 2.*sherex;

    // and drift field
    // we do centered difference here as well (second order)
    Real sdriftx  = d_rx[indexright]*d_state[indexright] - d_rx[indexleft]*d_state[indexleft];

    // and store the state array
    // remember that the drift field has negative sign (convention)
    d_state[dx+dy*GridModelWidth+i*get_global_area_2d] = d_state[dx+dy*GridModelWidth+i*get_global_area_2d]
                                                             + dtdiff*slaplace/dx2 - dtdiff*sdriftx/(2.*PhysicalCellWidth);
} // deterministic diffuse kernel (inhomogeneous)

// Kernel to compute Laplacian and diffuse in y direction (inhomogeneous with drift)
__kernel void inhomogeneous_deterministic_diffuseY(__global Real *d_state, Real dtdiff,
                                                   __global Real *d_diffy,
                                                   __global Real *d_ry,
                                                   int i)
// i denotes the species index to diffuse
{
    const Real dx2 = PhysicalCellWidth * PhysicalCellWidth;

    // get absolute position on the device
    size_t dx = get_global_id(0);
    size_t dy = get_global_id(1);

    // get neighbours - we use mirror BCs (zero-gradient, outflow)
    // todo: allow reflective BCs
    int indextop    = (dy > 0 ? dx+(dy-1)*GridModelWidth+i*get_global_area_2d : dx+(dy+1)*GridModelWidth+i*get_global_area_2d);
    int indexbottom = (dy < GridModelWidth-1 ? dx+(dy+1)*GridModelWidth+i*get_global_area_2d : dx+(dy-1)*GridModelWidth+i*get_global_area_2d);

    // compute laplacian (centered difference, second order)
    Real stop      = d_state[indextop]*d_diffy[indextop];
    Real sbottom   = d_state[indexbottom]*d_diffy[indexbottom];
    Real sherey    = d_state[dx+dy*GridModelWidth+i*get_global_area_2d]*d_diffy[dx+dy*GridModelWidth+i*get_global_area_2d];
    Real slaplace = stop + sbottom - 2.*sherey;

    // and drift field
    // we do centered difference here as well (second order)
    // remember that +y means towards GridModelHeight (indexbottom)
    Real sdrifty  = d_ry[indexbottom]*d_state[indexbottom] - d_ry[indextop]*d_state[indextop];

    // and store the state array
    // remember that the drift field has negative sign (convention)
    d_state[dx+dy*GridModelWidth+i*get_global_area_2d] = d_state[dx+dy*GridModelWidth+i*get_global_area_2d]
                                                             + dtdiff*slaplace/dx2 - dtdiff*sdrifty/(2.*PhysicalCellWidth);
} // deterministic diffuse kernel (inhomogeneous)


// Deterministic RK4 kernels to compute reactions
#if (GPGMP_NUM_REACTIONS > 0)

// returns the rhs of the chemical reaction ode for species i
inline void deterministic_rhs_ode_rk4(__local Real *_localState, __local Real *_deriv,
                                      Real _offset, __local Real *_rk
#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
                                      // we only need this if there are localized reactions in the model
                                      , image3d_t _reactionMask
#endif
                                 )
{
    // clean production/destruction arrays and add offset
    for (int j=0; j<GPGMP_NUM_SPECIES; j++) {
        _deriv[j] = 0.;
        _localState[j] += _offset * _rk[j];
    }

    // todo: are these actually global_id(x) ?
//    int device_x = get_group_id(0)*get_local_size(0)+get_local_id(0);
//    int device_y = get_group_id(1)*get_local_size(1)+get_local_id(1);
    const size_t GridX = get_global_id(0);
    const size_t GridY = get_global_id(1);
    
// <<<! defineUserParameters !>>>
    
    // compute ODE RHS and store it in derivative array
    const int _localId = get_local_id_2d; // used by defineSpeciesFromLocal
    const int _speciesOffsetLocal = 0;    // used by defineSpeciesFromLocal
// <<<! defineSpeciesFromLocalRK4 !>>>
// <<<! computeODERHS !>>>
// <<<! undefSpeciesFromLocal !>>>

    // remove offset
    for (int j=0; j<GPGMP_NUM_SPECIES; j++) {
        _localState[j] -= _offset * _rk[j];
    }
} // homogeneous_deterministic_rhs_ode

// Reaction kernel for RK4 Solver
__kernel void deterministic_performReaction_rk4(__global Real *d_state, Real deltat,

#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
                                                image3d_t _reactionMask,
#endif
                                                __local Real *s_local)
{
    // shared memory will be used to store the current state for this thread

    // LOCAL thread id
    int local_id = get_local_id_2d;

    // state array - will occupy first get_local_area_2d*sizeof(Real) bytes
    // pointer arithmetic makes it easier to access local memory just per s_state[i] (with i being species)
    __local Real *s_state = &s_local[NumSpecies*local_id];
    // these will be the arrays for the RK4 derivatives
    __local Real *s_rkk1 = &s_local[1*NumSpecies*get_local_area_2d + NumSpecies*local_id];
    __local Real *s_rkk2 = &s_local[2*NumSpecies*get_local_area_2d + NumSpecies*local_id];
    __local Real *s_rkk3 = &s_local[3*NumSpecies*get_local_area_2d + NumSpecies*local_id];
    __local Real *s_rkk4 = &s_local[4*NumSpecies*get_local_area_2d + NumSpecies*local_id];

    // copy current state over into shared memory
    for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
        s_state[i] = d_state[get_global_id_for_species_2d(i)];
    }

    // compute timestep - todo: adaptive step size control is the way to go
    // todo: or at least allow user to change nmax
    int nmax = 20;
    Real deltatchem = deltat / (Real)nmax;

    // do reaction loop
    for (int i=0; i<nmax; i++)
    {
        // advance one RK4 step
#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
        deterministic_rhs_ode_rk4(s_state, s_rkk1, 0,             s_rkk1, _reactionMask);
        deterministic_rhs_ode_rk4(s_state, s_rkk2, deltatchem/2., s_rkk1, _reactionMask);
        deterministic_rhs_ode_rk4(s_state, s_rkk3, deltatchem/2., s_rkk2, _reactionMask);
        deterministic_rhs_ode_rk4(s_state, s_rkk4, deltatchem,    s_rkk3, _reactionMask);
#else
        deterministic_rhs_ode_rk4(s_state, s_rkk1, 0,             s_rkk1);
        deterministic_rhs_ode_rk4(s_state, s_rkk2, deltatchem/2., s_rkk1);
        deterministic_rhs_ode_rk4(s_state, s_rkk3, deltatchem/2., s_rkk2);
        deterministic_rhs_ode_rk4(s_state, s_rkk4, deltatchem,    s_rkk3);
#endif

        // and add it to current state (loop over species)
        for (int j=0; j<GPGMP_NUM_SPECIES; j++) {
            s_state[j] += deltatchem/6. * (s_rkk1[j]+2.*s_rkk2[j]+2.*s_rkk3[j]+s_rkk4[j]);
            if (s_state[j]<0.)
                s_state[j]=0.;
        }
    }

    // Finally, we need to copy the state back to the device
    for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
        d_state[get_global_id_for_species_2d(i)] = s_state[i];
    }
} // deterministic performReaction (RK4)
#endif // NUM_REACTIONS>0

#if (GPGMP_NUM_REACTIONS > 0)
// compute RHS of the ODE for aqss solver
inline void deterministic_rhs_ode_aqss(__local Real *_localState, __local Real *_q, __local Real *_d
#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
                                       , image3d_t _reactionMask
#endif
                                       )
{
    // clean production/destruction arrays
    for (int j=0; j<GPGMP_NUM_SPECIES; j++) {
        _d[j]=0.;
        _q[j]=0.;
    }

    // TODO: Do we still need this ugly hack? I don't think so!
    _localState[-1] = 1.; // should be like that already!

    // get device coordinates
    int GridX = get_global_id(0);
    int GridY = get_global_id(1);

// <<<! defineUserParameters !>>>
   
    const int _localId = get_local_id_2d; // used by defineSpeciesFromLocal
    const int _speciesOffsetLocal = 0;    // used by defineSpeciesFromLocal
// <<<! defineSpeciesFromLocalRK4 !>>>
// <<<! computeODERHSaqss !>>>
// <<<! undefSpeciesFromLocal !>>>

    /**
    // cycle through reactions to get propensity
    for (int i=0; i<GPGMP_NUM_REACTIONS; i++)
    {
        // compute propensity for the reaction
        Real h0 = computePropensity(i, s_state, _reactionMask);

        // now we've got to go through the stoichiometry array and add it
        // to production/destruction array, according to the sign
        for (int j=0; j<GPGMP_NUM_SPECIES; j++) {
            int stoi = c_stoichiometry[i*NumSpecies+j];
            if (stoi < 0)
                s_d[j] += (Real)abs(stoi) * h0;
            else
                s_q[j] += (Real)stoi * h0;

        }
    }**/
} // homogeneous_deterministic_rhs_ode_aqss
#endif // NUM_REACTIONS>0

// main kernel to perform a series of reaction steps
// this one uses an alpha-QSS solver for stiff reactions
#if (GPGMP_NUM_REACTIONS > 0)
__kernel void deterministic_performReaction_aqss(__global Real *d_state, Real dtg,
                                                 __global Real *d_y0, __global Real *d_ys, __global Real *d_y1,
                                                 __global Real *d_ym1, __global Real *d_ym2,
                                                 __global Real *d_scrarray, __global Real *d_qs, __global Real *d_rtaus,
#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
                                                 image3d_t _reactionMask,
#endif
                                                 __local Real *s_local)
{
    // shared memory will be used to store the current state for this thread
    //__local Real s_local[homogeneous_deterministic_performStiffReaction_local_area * (3*NumSpecies + 1)];

    // LOCAL thread id
    const int local_id = get_local_id_2d;

    // y - will occupy first get_local_area_2d*sizeof(Real) bytes
    // the element s_y[-1] will contain value 1
    __local Real *s_y = &s_local[(NumSpecies+1)*local_id+1];

    // q array
    __local Real *s_q = &s_local[(NumSpecies+1)*get_local_area_2d + NumSpecies*local_id];

    // d array
    __local Real *s_d = &s_local[(2*NumSpecies+1)*get_local_area_2d + NumSpecies*local_id];

    // parameters
    Real const epsmin = 1e-2;
    Real const sqreps = 0.5;
    int  const itermax = 2;
    Real const epscl = 100.;
    Real const epsmax = 10.;
    Real const tfd = 1.000008;
    Real const dtmin = 1e-15;
    Real const ymin = 1e-20;

    // -- this is where the algorithm starts
    //initialize control parameter
    //tn is current time
    Real tn = 0.;

    // copy state array for current thread in shared memory and
    // keep bigger than ymin
    for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
        s_y[i] = d_state[get_global_id_for_species_2d(i)];
        d_y0[get_global_id_for_species_2d(i)] = s_y[i];
        s_y[i] = max(s_y[i], ymin);
    }

    // we also need s_y[-1]=1 for the rate laws (as in the stochastic simulator)
    s_y[-1] = 1.;

    // get first derivative value
    // q is production, d is consumption of the species.
#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
    deterministic_rhs_ode_aqss(s_y, s_q, s_d, _reactionMask);
#else
    deterministic_rhs_ode_aqss(s_y, s_q, s_d);
#endif

    // get initial stepsize
    // Following Mott & Oran (2001), we pick:
    // - strongly increasing functions (q >> d) use time step proportional
    //   to time reached to get to equilibrium
    // - decreasing functions take characteristic time as step size
    Real scrtch = 1e-25;

    for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
        Real ascr = fabs(s_q[i]);
        Real scr2 = copysign((Real)1./s_y[i], (Real)0.1*epsmin*ascr-s_d[i]);
        Real scr1 = scr2*s_d[i];
        scrtch = fmax(fmax(scr1, -fabs(ascr-s_d[i])*scr2), scrtch);
    }

    Real dt = min(sqreps/scrtch, dtg);
    Real ts, dtcurrent, eps;

    // line 100
    while (dtg > tn*tfd) {
        ts = tn;
        dtcurrent = dt;

        // this is r in the memo
        for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
            d_ys[get_global_id_for_species_2d(i)] = s_y[i];
            d_qs[get_global_id_for_species_2d(i)] = s_q[i];
            d_rtaus[get_global_id_for_species_2d(i)] = dtcurrent*s_d[i]/s_y[i];
        }

        // line 101
        while (dtg > tn*tfd) {
            // find predictor
            for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
                Real rtaui = dtcurrent*s_d[i]/s_y[i];
                Real alpha =
                (180.+rtaui*(60.+rtaui*(11.+rtaui)))/(360.+rtaui*(60.+rtaui*(12.+rtaui)));
                d_scrarray[get_global_id_for_species_2d(i)] =
                (s_q[i]-s_d[i])/(1.+alpha*rtaui);
            }

            // iterate
            for (int iteration=0; iteration<itermax; iteration++) {

                // limit decreasing functions to their minimum
                for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
                    d_ym2[get_global_id_for_species_2d(i)] = d_ym1[get_global_id_for_species_2d(i)];
                    d_ym1[get_global_id_for_species_2d(i)] = s_y[i];
                    s_y[i] = max(d_ys[get_global_id_for_species_2d(i)]
                                 +dt*d_scrarray[get_global_id_for_species_2d(i)],ymin);
                }

                // first corrector step
                if (iteration == 0) {
                    tn = ts + dt;
                    for (int i=0; i<GPGMP_NUM_SPECIES; i++)
                        d_y1[get_global_id_for_species_2d(i)] = s_y[i];
                } // end first iteration

                // evaluate derivative for corrector
                // q and d are now qp, dp
#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
                deterministic_rhs_ode_aqss(s_y, s_q, s_d, _reactionMask);
#else
                deterministic_rhs_ode_aqss(s_y, s_q, s_d);
#endif

                eps = 1e-10;

                for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
                    // this is pbar * dt (13)
                    Real rtaub = 0.5*(d_rtaus[get_global_id_for_species_2d(i)]+dt*s_d[i]/s_y[i]);
                    Real alpha = (180.*+rtaub*(60.+rtaub*(11.+rtaub)))/(360.+rtaub*(60.+rtaub*(12.+rtaub)));

                    // this is qtilde (11)
                    Real qt = d_qs[get_global_id_for_species_2d(i)]*(1.-alpha)+s_q[i]*alpha;

                    // this is pbar
                    Real pb = rtaub/dt;

                    // again the new corrected value is
                    // yc = y0 + dt*scrarray
                    d_scrarray[get_global_id_for_species_2d(i)] = (qt - d_ys[get_global_id_for_species_2d(i)]*pb)/(1.+alpha*rtaub);
                }
            } // end iteration

            // compute new f, check for convergence and limit decreasing functions.
            for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
                Real scr2 = fmax(d_ys[get_global_id_for_species_2d(i)]+dt*d_scrarray[get_global_id_for_species_2d(i)], (Real)0.);

                // scr1 = |yc-yp|
                Real scr1 = fabs(scr2-d_y1[get_global_id_for_species_2d(i)]);
                // y = yc
                s_y[i] = max(scr2, ymin);
                d_ym2[get_global_id_for_species_2d(i)] = d_ym1[get_global_id_for_species_2d(i)];
                d_ym1[get_global_id_for_species_2d(i)] = s_y[i];

                // compute only for the species that haven't died yet
                if (0.25*(d_ys[get_global_id_for_species_2d(i)]+s_y[i]) > ymin) {
                    // scr1 = |yc-yp|/yc
                    scr1 = scr1/s_y[i];

                    eps = fmax((Real)0.5*(scr1 + fmin(fabs(s_q[i]-s_d[i]) / (s_q[i]+s_d[i]+(Real)1e-30), scr1)),
                               eps);
                }
            }

            eps = eps*epscl;

            // check if step size is too small
            if (dt <= dtmin+1e-16*tn) {
                dt = dtg - ts;
                dt = fmin(dtmin, fabs(dt));
            }

            // stability check
            Real stab = 0.01;
            if (itermax >= 3) {
                for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
                    stab =
                    fmax(stab,
                         fabs(s_y[i]-d_ym1[get_global_id_for_species_2d(i)])
                         /(fabs(d_ym1[get_global_id_for_species_2d(i)] - d_ym2[get_global_id_for_species_2d(i)])+(Real)1e-20*s_y[i]));
                }
            }

            // check if step is valid
            if (eps <= epsmax && stab <= 1.) {
                // step is valid. Return if max time is reached
                if (dtg <= tn*tfd)
                    break;

            } else {
                // invalid step
                tn = ts;
            }

            // perform step size modification
            Real rteps = sqrt(eps);

            Real dto = dt;
            dt = fmin(fmin(dt*((Real)1.0/rteps + (Real).005), tfd*(dtg-tn)), dto/(stab+(Real)0.001));

            // begin new step if previous step converged
            if (eps <=epsmax && stab <= 1) {
                break;
            }

            // this is an unsuccessful step
            // goto 101
            dto = dt/dto;
            for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
                d_rtaus[get_global_id_for_species_2d(i)] *= dto;
            }
        } // end 101


        // successful step, get source term for next step
#ifdef GPGMP_HAS_LOCALIZED_REACTIONS
        deterministic_rhs_ode_aqss(s_y, s_q, s_d, _reactionMask);
#else
        deterministic_rhs_ode_aqss(s_y, s_q, s_d);
#endif
        // goto 100
    }// end 100

    // Finally, we need to copy the state back to the device
    for (int i=0; i<GPGMP_NUM_SPECIES; i++) {
        d_state[get_global_id_for_species_2d(i)] = s_y[i];
    }
} // homogeneous_deterministic_performStiffReaction

#endif // (GPGMP_NUM_REACTIONS > 0)
