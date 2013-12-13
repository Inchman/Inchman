#!/bin/python
import numpy
import random
import matplotlib.pyplot as plt
import scipy.spatial.distance

# for analytic solutions
import gpgmp.test.homogeneous_drift_diffusion
import gpgmp.test.multiplicative_noise
import gpgmp.constants
import gpgmp.helpers
from gpgmp.constants import na as avogadro

import visual

# quickly displays results of simple A+B->A reaction
def simple_ab():
    # plot reference
    #tot=numpy.genfromtxt("total.dat.reference")
    #time=numpy.mgrid[0:200000]*1e-6
    #plt.plot(time, tot-500, "g-")

    # read in total..
    totsuper=numpy.genfromtxt("total.dat")
    
    # take time points from output
    time = totsuper[:,0]

    # compute reaction rate from reaction radius
    k1 = 8e-3*4*numpy.pi*(1+1)
    volsub = 2**3

    # and overplot it
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500, "b-")

    # plot variance lines
    var=numpy.exp(-k1/volsub*500.*time)*500.*(1.-numpy.exp(-k1/volsub*500.*time))
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500+numpy.sqrt(var), "b--")
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500-numpy.sqrt(var), "b--")
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500+2.*numpy.sqrt(var), "b:")
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500-2.*numpy.sqrt(var), "b:")

    # plot superstel
    plt.plot(totsuper[:,0], totsuper[:,2]-500, "r-", linewidth=2)

    # overplot smoldyn solution
    8e-3*4*numpy.pi*(1+1)
    volsub = 2**3

    # and overplot it
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500, "b-")

    # plot variance lines
    var=numpy.exp(-k1/volsub*500.*time)*500.*(1.-numpy.exp(-k1/volsub*500.*time))
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500+numpy.sqrt(var), "b--")
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500-numpy.sqrt(var), "b--")
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500+2.*numpy.sqrt(var), "b:")
    plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500-2.*numpy.sqrt(var), "b:")

    # plot superstel
    plt.plot(totsuper[:,0], totsuper[:,2]-500, "r-", linewidth=2)

    # overplot smoldyn solution
    smtot = numpy.genfromtxt('smoldyn.dat')
    plt.plot(smtot[:,0], smtot[:,1], "g-")

# reads in the protective zones and displays them on the screen
def show_protective(nparticles):
    # create scene
    scene = visual.display(title='Protective zones', width=600, height=400, center=(0,0,0))
    scene.select()

    # read in positions, state and protective zones
    positions = (numpy.genfromtxt('position.dat'))[:nparticles]
    state = (numpy.genfromtxt('state.dat'))[:nparticles]
    zones = (numpy.genfromtxt('fpt.dat'))[:nparticles]

    # go through particles and display them
    for i in range(nparticles):
        # color the spheres according to state
        cball = visual.color.green

        if (state[i,0]==1):
            cball = visual.color.red
        else:
            cball = visual.color.blue

        if (state[i,3]==1):
            cball = visual.color.yellow
            
        if (state[i,3]==2):
            cball = visual.color.orange

        visual.sphere(pos=(positions[i,0], positions[i,1], positions[i,2]), radius=zones[i], color=cball)
        visual.label(pos=(positions[i,0], positions[i,1], positions[i,2]), text='{0:d}'.format(i))

def read_ascii(filename='species', nparticles=-1):
    data = numpy.genfromtxt(filename, unpack=True)

    return (data[0, :nparticles], data[1, :nparticles], data[2, :nparticles])

def read_grid(filename='grid'):
    data = numpy.genfromtxt(filename, unpack=True)
    d0 = numpy.reshape(data, ((64*64*64)), order='F')
    dn = numpy.reshape(data, (64,64,64), order='C')

    # should have row major now.. first index is x, second index is z, third is y? or sth like that
    # this one agrees with the binning in plot_image_diffusion
    return numpy.sum(dn, axis=1)

def plot_image_diffusion(filename = 'species', nparticles = -1):
    x, y, z = read_ascii(filename, nparticles)
    # compute edges of bins
    minx = -10.
    maxx = 10.
    ncells = 64.
    dx = (maxx-minx)/ncells
    edges=numpy.mgrid[0:ncells+1]*dx+minx

    # bin it
    H, xedges, yedges = numpy.histogram2d(x, y, bins=edges)

    # plot 
    plt.imshow(H, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

    return H, xedges, yedges

def plot_rw_diffusion(tmax, filename='species', sliceX=0, sliceY=0, sx=1., sy=1., mux=0., muy=0.):
    # read data
    x, y, z = read_ascii(filename)
    
    # compute # particles
    nparticles = numpy.shape(x)[0]

    # bin it
    H, xedges, yedges = numpy.histogram2d(x, y, bins=(100, 100))

    # plot 
    dx = (numpy.diff(xedges))[0]
    dy = (numpy.diff(yedges))[0]

    xc = xedges[:-1]+0.5*numpy.diff(xedges)
    yc = yedges[:-1]+0.5*numpy.diff(yedges)

    fig=plt.figure(1, figsize=(gpgmp.constants.width2c,gpgmp.constants.height1c), dpi=150)
    plt.clf()

    fig.subplots_adjust(hspace=0.1, bottom=0.15, right=0.95, left=0.1)
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, sharex=ax1, sharey=ax1)
    
    # get slice position
    ix, vx = gpgmp.helpers.find_nearest(xc, sliceX)
    iy, vy = gpgmp.helpers.find_nearest(yc, sliceY)

    # plot slices
    ax1.plot(xc, H[ix, :], "ro")
    ax2.plot(yc, H[:, iy], "ro")
    ana = gpgmp.test.homogeneous_drift_diffusion.biased_diffusion_2d(xc, yc, tmax, sx**2/2., sy**2/2., mux, muy)*dx*dy*nparticles

    ax1.plot(xc, ana[ix, :], "b-")
    ax2.plot(yc, ana[:, iy], "b-")

def find_closest(x, y, z):

    # create matrix
    points=numpy.transpose(numpy.array([x,y,z]))

    # compute distance and transform to square matrix
    dist = scipy.spatial.distance.pdist(points, 'euclidean')
    sfdist = (scipy.spatial.distance.squareform(dist))

    # omit self distance
    sfdist[numpy.where(sfdist==0)]=numpy.inf

    # get closest indices and corresponding distance
    closest = numpy.argmin(sfdist, axis=0)
    cdist = numpy.min(sfdist, axis=0)

    return closest, cdist

def get_rradius(k, dt):
    # compute reaction radius (in mum)
    return (3.*k*dt/(4.*numpy.pi*1e3*avogadro))**(1./3.)*1e6

def zo_reaction(x, y, z, state, dt, box):
    # perform zeroth order reactions

    # probability for a 0-th order reaction to happen during time step delta t
    # is given by Poisson process with mean, variance = lambda dt
    # [compare Floyd (2007), p 13]

    # 0 -> A
    k0 = 1.2
    nna = numpy.random.poisson(k0*dt)
    nxa = numpy.random.uniform(low=box[0], high=box[1], size=(nna))
    nya = numpy.random.uniform(low=box[2], high=box[3], size=(nna))
    nza = numpy.random.uniform(low=box[4], high=box[5], size=(nna))
    nsa = numpy.ones((nna))
    x = numpy.append(x, nxa)
    y = numpy.append(y, nya)
    z = numpy.append(z, nza)
    state = numpy.append(state, nsa)

    # 0 -> B
    k1 = 1.
    nnb = numpy.random.poisson(k1*dt)
    nxb = numpy.random.uniform(low=box[0], high=box[1], size=(nnb))
    nyb = numpy.random.uniform(low=box[2], high=box[3], size=(nnb))
    nzb = numpy.random.uniform(low=box[4], high=box[5], size=(nnb))
    nsb = numpy.ones((nnb))*2
    x = numpy.append(x, nxb)
    y = numpy.append(y, nyb)
    z = numpy.append(z, nzb)
    state = numpy.append(state, nsb)

    return x, y, z, state

def bimolecular(x, y, z, state, dt):

    # TODO: find_closest is not sufficient
    # particle might not be able to react with its closest
    # neighbour but can still react with a different one
    # solution: find_closest for certain type

    # A + A -> 0
    rr1 = get_rradius(287158083.91571045, dt)

    # A + B -> 0
    #rr2 = get_rradius(2871580839.1571045, dt)
    rr2 = get_rradius(287158083.9157105, dt)

    # get closest neighbor for each particle
    closest, dist = find_closest(x,y,z)

    # check if smaller than reaction radius and fits reactant stoich *AND* is reciprocal closest
    # and perform it right away
    nparts = (numpy.size(x))
    identity = numpy.mgrid[0:nparts]

    ri1 = numpy.where(numpy.logical_and(numpy.logical_and(dist <= rr1,
                                                          numpy.logical_or(numpy.logical_and(state==1, state[closest]==1),
                                                                           numpy.logical_and(state==1, state[closest]==1))),
                                        closest[closest[identity]]==identity))
    state[ri1] = 0

    ri2 = numpy.where(numpy.logical_and(numpy.logical_and(dist <= rr2,
                                                          numpy.logical_or(numpy.logical_and(state==1, state[closest]==2),
                                                                           numpy.logical_and(state==2, state[closest]==1))),
                                        closest[closest[identity]]==identity))

#    ri2 = numpy.where(numpy.logical_and(numpy.logical_and(dist <= rr2,
#                                                          numpy.logical_and(state==1, state[closest]==2)),
#                                        closest[closest[identity]]==identity))
    state[ri2] = 0


def explicit_Euler(x, dt, ddm):
    
    # compute diffusivity and drift
    a, b, bp = ddm(x)

    # get Wiener increments
    # dw ~ N(0, dt) (note that dt=sigma^2)
    #dw = random.gauss(0, numpy.sqrt(dt))
    dw = numpy.random.normal(loc=0., scale=numpy.sqrt(dt), size=(numpy.shape(x)[0]))

    # returns dx
    return a*dt+b*dw

def explicit_Milstein(x, dt, ddm):
    # compute diffusivity and drift
    a, b, bp = ddm(x)
    cn = b*bp/2.

    # get Wiener increment
    # dw ~ N(0, dt) (note that dt=sigma^2)
    #dw = random.gauss(0, numpy.sqrt(dt))
    dw = numpy.random.normal(loc=0., scale=numpy.sqrt(dt), size=(numpy.shape(x)[0]))

    # returns dx
    return a*dt+b*dw+cn*(dw**2-dt)

def run_multiple():
    nsteps=1000
    naf = numpy.zeros((10,nsteps), dtype=float)
    nbf = numpy.zeros((10,nsteps), dtype=float)

    for avr in range(10):
        t, na, nb = main(tmax=2., nparticles=20000, nsteps=nsteps, prob='simpleab')
        naf[avr, :]=na
        nbf[avr, :]=nb

        print avr

    return t, naf, nbf

def apply_boundary_conditions(x,y,z,dx,dy,dz,box):

        # check if outside box and put back in (reflective)
        lxil = numpy.where(x < box[0])
        x[lxil] = x[lxil] - 2.*dx[lxil]
        lxir = numpy.where(x > box[1])
        x[lxir] = x[lxir] - 2.*dx[lxir]
        lyil = numpy.where(y < box[2])
        y[lyil] = y[lyil] - 2.*dy[lyil]
        lyir = numpy.where(y > box[3])
        y[lyir] = y[lyir] - 2.*dy[lyir]
        lzil = numpy.where(z < box[4])
        z[lzil] = z[lzil] - 2.*dz[lzil]
        lzir = numpy.where(z > box[5])
        z[lzir] = z[lzir] - 2.*dz[lzir]

def init_particles(prob, box, nparticles):
    if prob=='simpleab':
        x = numpy.random.uniform(low=box[0], high=box[1], size=(nparticles))
        y = numpy.random.uniform(low=box[2], high=box[3], size=(nparticles))
        z = numpy.random.uniform(low=box[4], high=box[5], size=(nparticles))
        state=numpy.zeros((nparticles), dtype=int)
        state[:nparticles/2]=1
        state[nparticles/2:]=2

    return x,y,z,state

def main(tmax=10., nsteps=100, nparticles=1, alg='explicit_Euler', prob='diffusion'):
    
    # seed RNG
    random.seed()

    # compute time step
    dt = tmax/nsteps

    # computational volume
    # (xl, xr, yl, yr)
    box = (-5., 5, -5, 5)

    # set up problem
    if prob=='diffusion':
        x0, y0 = (0., 0.)
        mux, sx = (0, 0.5)
        muy, sy = (0, 2)
        funcx = lambda x: (mux, sx, 0)        
        funcy = lambda y: (muy, sy, 0)        
        sliceX, sliceY = (tmax*mux, tmax*muy)

    elif prob=='multiplicative':
        x0, y0 = (16., 16.)
        mux, sx = (0.1, 5)
        muy, sy = (0.2, 4)
        funcx = lambda x: (mux*x, sx*x, sx)
        funcy = lambda y: (muy*y, sy*y, sy)
        sliceX, sliceY = (8., 8.)

    elif prob=='simpleab':
        diffc = 0.00001
        mux, sx = (0, numpy.sqrt(2.*diffc))
        muy, sy = (0, numpy.sqrt(2.*diffc))
        muz, sz = (0, numpy.sqrt(2.*diffc))
        funcx = lambda x: (mux, sx, 0)        
        funcy = lambda y: (muy, sy, 0)        
        funcz = lambda y: (muz, sz, 0)        
        # gives sidelength of 100./128. mum
        # to agree with Vigelius et al. (2010)       
        box = (-0.390625, 0.390625, -0.390625, 0.390625, -0.390625, 0.390625)
        # output counts
        nao = numpy.zeros((nsteps), dtype=int)
        nbo = numpy.zeros((nsteps), dtype=int)
        
    # create array
    x     = numpy.zeros((nparticles), dtype=float)
    y     = numpy.zeros((nparticles), dtype=float)
    z     = numpy.zeros((nparticles), dtype=float)
    state = numpy.zeros((nparticles), dtype=int)


    x,y,z,state = init_particles(prob, box)

    for i in range(nsteps):
        # compute dx
        dx=explicit_Euler(x, dt, funcx )
        dy=explicit_Euler(y, dt, funcy)      
        dz=explicit_Euler(z, dt, funcz)      
       
        # compute new position
        x = x + dx
        y = y + dy
        z = z + dz

        # apply BCs
        apply_boundary_conditions(x,y,z,dx,dy,dz,box)

        if prob=='simpleab':
            # count molecules
            na = numpy.sum(state==1)
            nb = numpy.sum(state==2)
            ntot = numpy.size(x)

            # and store them
            nao[i] = na
            nbo[i] = nb

            if (i % (nsteps/10.)) == 0:
                print "t={0:g}: A={1:d}, B={2:d}, tot={3:d}".format(i*dt, na, nb, ntot)

            # perform bimolecular reactions
            if ntot>0:
                bimolecular(x, y, z, state, dt)

            # remove all zero states
            zerind = numpy.where(state==0)
            x=numpy.delete(x, zerind)
            y=numpy.delete(y, zerind)
            z=numpy.delete(z, zerind)
            state=numpy.delete(state, zerind)

            # perform creation reactions
            #x, y, z, state = zo_reaction(x, y, z, state, dt, box)


    if prob=='diffusion' or prob=='multiplicative':
        # create histogram
        nx, ny = (50, 50)
        H, xedges, yedges = numpy.histogram2d(x, y, bins=(nx, ny))

        # plot 
        dx = (numpy.diff(xedges))[0]
        dy = (numpy.diff(yedges))[0]

        xc = xedges[:-1]+0.5*numpy.diff(xedges)
        yc = yedges[:-1]+0.5*numpy.diff(yedges)

        fig=plt.figure(1, figsize=(gpgmp.constants.width2c,gpgmp.constants.height1c), dpi=150)
        plt.clf()

        fig.subplots_adjust(hspace=0.1, bottom=0.15, right=0.95, left=0.1)
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2, sharex=ax1, sharey=ax1)

        # get slice position
        ix, vx = gpgmp.helpers.find_nearest(xc, sliceX)
        iy, vy = gpgmp.helpers.find_nearest(yc, sliceY)

        # plot slices
        ax1.plot(xc, H[ix, :], "ro")
        ax2.plot(yc, H[:, iy], "ro")
        if prob == 'diffusion':
            ana = gpgmp.test.homogeneous_drift_diffusion.biased_diffusion_2d(xc, yc, tmax, sx**2/2., sy**2/2., mux, muy)*dx*dy*nparticles
        elif prob == 'multiplicative':
            ana = gpgmp.test.multiplicative_noise.analytic_2d(xc, 8, yc, 8, sx, sy, tmax, mux, muy)*dx*dy*nparticles

            ax1.plot(xc, ana[ix, :], "b-")
            ax2.plot(yc, ana[:, iy], "b-")

        return x,y, H

    elif prob=='simpleab':
        t = numpy.mgrid[0:nsteps]*dt

        return t, nao, nbo


