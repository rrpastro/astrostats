
#rishi tried to edit this code
# Changes by Gizis to SDR's code
# change to use import numpy as np, import matplotlib as mpl
# and matplotlib.pyplot as plt

import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt

'''
Figure out total mass contained by Centaur population
'''

Mpluto = 1.309e25 # grams
kmcm = 1.0e5
fourthirds = 4./3.

# Chariklo (largest Centaur) diameters in km:
cs1 = 236.8 # Fornasier et al. (2013)
cs2 = 257. # Stansberry et al. (2008)
cs3 = 273. # Altenhoff et al. (2001)
cs4 = 302. # Jewitt & Kalas (1998)
cs5 = 236. # Groussin et al. (2004)
Chariklo = kmcm * np.array([cs1, cs2, cs3, cs4, cs5])

# Orcus diameter (Fornasier et al. 2013, Brown et al. 2010)
Orcus = kmcm * np.array([917., 900.])
# Orcus density (g/cm^3, same sources)
Odens = np.array([1.4, 1.53])
# Salacia diameter (Fornasier et al. 2013, Stansberry et al. 2012)
Salacia = kmcm * np.array([854., 905.])
# Salacia density
Sdens = np.array([1.29, 1.16])

# Threshold for separating 'large' from 'small' KBOs: large
# are those that obey N(>R) propto R^(1-q), where q = 4.
lthresh = 30. * kmcm

# Minimum size to consider
minsize = 0.01 * kmcm

# Thresholds for further breaks in the power laws (km):
t2 = 10. * kmcm
t3 = 2. * kmcm
t4 = 0.1 * kmcm

# Exponents for the power-law segments (large -> small):
q1 = 4.0
q2 = 2.0
q3 = 5.8
q4 = 2.5
q5 = 3.7 # This is for < 0.1 km -- how much mass?

'''
Find the total mass of Centaurs given a maximum Centaur radius, a
density, and the power-law exponents and thresholds given above
(integrate Schlichting et al. (2013) broken power law)
'''
def integrate_power(Rmax,dens):

    # First find the total number of LARGE Centaurs (> 30 km)
    # given the maximum size
    nlarge = (Rmax/lthresh)**(q1-1.)

    # Now, setting R0 = lthresh to scale the power law, find N0
    # for the largest bodies (q=4)
    N0_large = nlarge * (q1-1.)
    print 'N0_large', N0_large

    # Find the total mass contained in large bodies
    quant = fourthirds * np.pi * dens
    mlarge = lthresh**(q1-1.) * N0_large * quant * \
             (np.log(Rmax)-np.log(lthresh))
    print 'Large body mass', mlarge/Mpluto

    '''
    Now go up the rest of the broken power law, finding masses
    contained in each segment
    '''
    # Second segment: q=2
    N0_med = N0_large * (q2-1.)/(q1-1.)
    print 'N0_med', N0_med
    # Here we do not have a log function appear in the integral,
    # so we can write it in general terms of q2
    mmed = N0_med * (lthresh**(q2-1.)) * quant * \
           (lthresh**(4.-q2) - t2**(4.-q2)) / (4.-q2)
    print 'Medium body mass', mmed/Mpluto
    # Now rescale the second segment of the broken power law to
    # R0=t2 instead of R0=lthresh
    N0_med_new = N0_med * (t2/lthresh)**(1.-q2)
    print 'N0_med_new', N0_med_new
    # Third segment: q=5.8. Connect to previous segment at t2
    N0_small = N0_med_new * (q3-1.) / (q2-1.)
    print 'N0_small', N0_small
    msmall = N0_small * (t2**(q3-1.)) * quant* \
             (t2**(4.-q3) - t3**(4.-q3)) / (4.-q3)
    print 'Small body mass', msmall/Mpluto
    # Rescale third power law segment to R0=t3
    N0_small_new = N0_small * (t3/t2)**(1.-q3)
    print 'N0_small_new', N0_small_new
    # Now attach the fourth segment of the power law
    N0_little = N0_small_new * (q4-1.) / (q3-1.)
    print 'N0_little', N0_little
    mlittle = N0_little * (t3**(q4-1.)) * quant* \
              (t3**(4.-q4) - t4**(4.-q4)) / (4.-q4)
    print 'Little body mass (0.1-2 km)', mlittle/Mpluto
    # Rescale fourth power law segment to R0=t4
    N0_little_new = N0_little * (t4/t3)**(1.-q4)
    # Attach the final power law segment
    N0_tiny = N0_little_new * (q5-1.) / (q4-1.)
    print 'N0_tiny', N0_tiny
    mtiny = N0_tiny * (t4**(q5-1.)) * quant* \
              (t4**(4.-q5) - minsize**(4.-q5)) / (4.-q5)
    print 'Little body mass (minsize-0.1 km)', mtiny/Mpluto

    # Total number of Centaurs
    Ntot = (N0_tiny / (q5-1.)) * (minsize/t4)**(1.-q5)
    print 'Total number of Centaurs above r (km)', minsize/kmcm, Ntot

    # Add up total mass
    mtot = mlarge + mmed + msmall + mlittle + mtiny
    print 'Total mass', mtot/Mpluto

    # Get arrays to add to plot of N>(r) and then one of M<(r)
    r = np.logspace(np.log10(minsize), np.log10(Rmax), num=500, endpoint=True)
    ngr = np.zeros_like(r)
    mlr = np.zeros_like(r)
    r1 = np.where((r >= lthresh))
    r2 = np.where((r < lthresh) & (r >= t2))
    r3 = np.where((r < t2) & (r >= t3))
    r4 = np.where((r < t3) & (r >= t4))
    r5 = np.where((r < t4) & (r >= minsize))
    ngr[(r1)] = N0_large/(q1-1.)*(r[(r1)]/lthresh)**(1.-q1)
    ngr[(r2)] = N0_med/(q2-1.)*(r[(r2)]/lthresh)**(1.-q2)
    ngr[(r3)] = N0_small/(q3-1.)*(r[(r3)]/t2)**(1.-q3)
    ngr[(r4)] = N0_little/(q4-1.)*(r[(r4)]/t3)**(1.-q4)
    ngr[(r5)] = N0_tiny/(q5-1.)*(r[(r5)]/t4)**(1.-q5)
    mlr[(r5)] = (t4**(q5-1.)) * N0_tiny * quant * \
                (r[(r5)]**(4.-q5) - minsize**(4.-q5)) / (4.-q5)
    mlr[(r4)] = max(mlr[(r5)]) + (t3**(q4-1.)) * N0_little * quant * \
                (r[(r4)]**(4.-q4) - t4**(4.-q4)) / (4.-q4)
    mlr[(r3)] = max(mlr[(r4)]) + (t2**(q3-1.)) * N0_small * quant * \
                (r[(r3)]**(4.-q3) - t3**(4.-q3)) / (4.-q3)
    mlr[(r2)] = max(mlr[(r3)]) + (lthresh**(q2-1.)) * N0_med * quant * \
                (r[(r2)]**(4.-q2) - t2**(4.-q2)) / (4.-q2)
    mlr[(r1)] = max(mlr[(r2)]) + (lthresh**(q1-1.)) * N0_large * quant * \
                (np.log(r[(r1)]) - np.log(lthresh))
    print 'Sanity check total mass', max(mlr)/Mpluto
    r = r / kmcm
    
    return mtot, r, ngr, mlr


'''
Compute the total mass of Centaurs for several
Rmax/density combinations and plot results
'''
def mtot(detfrac):

    # Number of large Centaurs simulated by Tiscareno & Malhotra
    nobs = 53.

    # Total number of large Centaurs above lthresh, computed from
    # detectability metrics above
    ntot = nobs / detfrac

    # Power-law normalization for very largest objects, using R0
    # = lthresh
    N0_large = (q1 - 1.) * ntot
    print N0_large

    # Find estimate of largest Centaur size from Schlichting et
    # al. (2013) size distribution such that N(>30 km) =
    # Nobs/detfrac
    maxsize = lthresh * ntot**(1./(q1-1.))
    print maxsize
    # Set a minimum size of R = 0.01 km
    minsize = 0.01 * kmcm

    # Try three more estimates of the largest Centaur size:
    # Let's put Chariklo at its largest observed size -- what is the
    # total mass in Centaurs? Then minimum Chariklo size, then
    # maximum Orcus size (similar also to Salacia)
    Crad = np.array([max(Orcus)/2., maxsize, max(Chariklo)/2., \
            min(Chariklo)/2.]) 

    # Density values: try max(Orcus), min(Salacia) -- little
    # geological alteration, plus 67P nucleus
    dens = np.array([min(Sdens), max(Odens), 0.5])

    # Start the figure
    fig_size=(7,8.5)
    plt.rcParams["figure.figsize"] = fig_size
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    f.subplots_adjust(left=0.14,right=0.95,bottom=0.09,top=0.95, \
                      hspace=0.1)
    # Loop through all the max radius and density possibilities
    cols = ['b', 'r', 'g', 'k']
    lins = ['-', '--', '-.']
    labels = ['Orcus', 'TM03', 'Chariklo 1', 'Chariklo 2']
    denslabels = ['Salacia', 'Orcus', '67P']
    for i in range(len(Crad)):
        for j in range(len(dens)):
            mtot, r, ngr, mlr = integrate_power(Crad[i],dens[j])
            # Average mass of Centaur = mtot / N>(minimum size)
            mave = mtot / ngr[0]
            print 'Average mass', mave
            # Number of Centaurs > 1 km in size
            g1km = np.where((r >= 1.0))
            print 'Size threshold, N>:', r[(g1km[(0)][0])], \
                   ngr[(g1km[(0)][0])]
            if (j == 0):
                ax1.loglog(r, ngr, color=cols[i], linewidth=1.7, \
                           label=labels[i])
            else:
                ax1.loglog(r, ngr, color=cols[i], linewidth=1.7)
            if (i == 0):
                ax2.loglog(r, mlr/Mpluto, color=cols[i], \
                           ls=lins[j], linewidth=1.7, label=denslabels[j])
            else:
                ax2.loglog(r, mlr/Mpluto, color=cols[i], \
                           ls=lins[j], linewidth=1.7)

    # Finish the figure
    ax1.xaxis.grid(True)
    ax1.set_xlim([0.5,5.0e2])
    ax1.set_ylim([1,1.0e8])
    ax1.set_ylabel('$N_>(R)$', fontsize='x-large')
    ax1.set_title('Centaur Radius and Mass Distribution', \
                   fontsize='x-large')
    ax1.legend(loc='best', fontsize='medium', title='Largest Centaur')
    ax2.xaxis.grid(True)
    ax2.set_ylim([1.0e-4,1.0])
    ax2.set_xlabel('Radius [km]', fontsize='x-large')
    ax2.set_ylabel('$M_<(R)$ [$M_{Pluto}$]', fontsize='x-large')
    ax2.legend(loc='best', fontsize='medium', title='Density')
    f.show()

    # Another 2 figures, just to check proportion in medium size
    # bin and large bin:
    plt.figure()
    testrad = np.where(r >= 10.)
    ml10km = mlr[(testrad[(0)][0])]
    print 'm < 10 km', ml10km
    mlr = mlr - ml10km
    plt.loglog(r, mlr/Mpluto, 'k-')
    plt.xlim([10.,30.])
    # ylim([0.1, 0.13])
    plt.xlabel('Radius')
    plt.ylabel('$M_<(R)$')
    plt.show()

    plt.figure()
    testrad = np.where(r >= 30.)
    ml30km = mlr[(testrad[(0)][0])]
    print 'm < 30 km', ml30km
    mlr = mlr - ml30km
    plt.loglog(r, mlr/Mpluto, 'k-')
    plt.xlim([30.,300.])
    # ylim([0.1, 0.13])
    plt.xlabel('Radius')
    plt.ylabel('$M_<(R)$')
    plt.show()
