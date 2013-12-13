typedef struct {
    const GPGMP_STATE_TYPE State;
    __global int *SpeciesState;
    const Real MeanX;
    Real DiffusivityX;
    Real DiffusivityY;
    Real DriftX;
    Real DriftY;
} _species_t;

// Function Prototypes - Apple's OpenCL compiler will complain if they're not specified (Mac OS X 10.7.3)
inline Real getMeanX(__global Real *sumMoments, __global Real *sumStates, int ispecies);
inline Real getPhysicalX(void /* C wants the void- makes the prototype work properly*/);
inline Real getPhysicalY(void /* C wants the void- makes the prototype work properly*/);
inline Real getCentralDifferenceX(_species_t *species);
inline Real getCentralDifferenceY(_species_t *species);

inline Real getMeanX(__global Real *sumMoments, __global Real *sumStates, int ispecies) {
    Real tmp = sumMoments[ispecies] / sumStates[ispecies];
    return PhysicalModelWidth / (Real)GridModelWidth * (tmp-(Real)GridModelWidth/2);
}


// TODO: maybe remove ifdefs and use two different source files for CA and FPE solvers?

// Helper function to compute a physical coordinate from the thread ID
// TODO: Add functionality to establish a physical coordinate system (including units)
//       in the model which is then automatically reproduced here.
inline Real getPhysicalX(void) {
    //return (Real)PhysicalModelWidth / (Real) GridModelWidth * ((int)get_global_id(0)-GridModelWidth/2);
    return (Real)PhysicalModelWidth / (Real) GridModelWidth * ((Real)get_global_id(0)-(Real)GridModelWidth/2.);
}
inline Real getPhysicalY(void) {
    //return (Real)PhysicalModelHeight / (Real) GridModelHeight * ((int)get_global_id(1)-GridModelHeight/2);
    return (Real)PhysicalModelHeight / (Real) GridModelHeight * ((Real)get_global_id(1)-(Real)GridModelHeight/2.);
}

// Helper function to compute the central difference in X direction.
inline Real getCentralDifferenceX(_species_t *species)
{
    __global int *d_state = species->SpeciesState;

    const size_t dx = get_global_id(0);
    const size_t dy = get_global_id(1);

    int sleft, sright;

    if (dx > 0) {
        sleft = d_state[(dx-1)+dy*get_global_size(0)];
    } else {
        sleft = d_state[(dx+1)+dy*get_global_size(0)];
    }

    if (dx < get_global_size(0)-1) {
        sright = d_state[(dx+1)+dy*get_global_size(0)];
    } else {
        sright = d_state[(dx-1)+dy*get_global_size(0)];
    }

    return sright-sleft;
}

// Helper function to compute the central difference in Y direction
inline Real getCentralDifferenceY(_species_t *species)
{
    __global int *d_state = species->SpeciesState;

    const size_t dx = get_global_id(0);
    const size_t dy = get_global_id(1);

    int stop, sbottom;

    if (dy > 0) {
        stop = d_state[dx+(dy-1)*get_global_size(0) ];
    } else {
        stop = d_state[dx+(dy+1)*get_global_size(0) ];
    }

    if (dy < get_global_size(0)-1) {
        sbottom = d_state[dx+(dy+1)*get_global_size(0) ];
    } else {
        sbottom = d_state[dx+(dy-1)*get_global_size(0) ];
    }

    // this is due to the definition of the coordinate system for the velocity vector
    return sbottom-stop;
}
