import matplotlib.pyplot as plt
from quicktools.bin_spec import bin_spec
import numpy as np

def plot_SpecAll(wave, flux, err, continuum=None,
                 targetname='TargetAll', filedir='.'):

    '''
    Plot the overall spectrum and continuum.
    '''

def plot_OneLine(data, redchi2=-1, figname='test.pdf',
                 targetname=None, linename=None, nbin=3, 
                 velwidth=1000., pltrange=[-200, 200]):

    '''
    Plot and (evaluate) the continuum fits.
    '''

    flux, sig, conti, nflux, nsig, vel = data[1], data[2], data[3], data[4], data[5], data[6]
    fig = plt.figure(figsize=(6, 6))
    ax1 = fig.add_axes([0.12, 0.48, 0.82, 0.41])   # original spec, +/- velwidth
    ax2 = fig.add_axes([0.12, 0.27, 0.82, 0.2])   # normalized, +/- velwidth
    ax3 = fig.add_axes([0.12, 0.08, 0.82, 0.15])

    fig.text(0.94, 0.95, '%s'%(targetname), fontsize=12, horizontalalignment='right')
    fig.text(0.94, 0.905, '%s'%(linename), fontsize=12, horizontalalignment='right', color='r', fontweight='bold')

    # ax1, original spectra, with continuum on top
    # flux array
    if nbin > 1:
        x1, y1 = bin_spec(vel, flux, nbin)
        xsig1, ysig1 = bin_spec(vel, sig, nbin)
    else:
        x1, y1 = vel, flux
        xsig1, ysig1 = vel, sig

    xx1, yy1 = np.repeat(x1, 2)[1:], np.repeat(y1, 2)[:-1]
    xxsig1, yysig1 = np.repeat(xsig1, 2)[1:], np.repeat(ysig1, 2)[:-1]
    # plot
    ax1.plot(xx1, yy1, color='k', lw=0.7)
    ax1.plot(xxsig1, yysig1, color='b', lw=0.7)
    ax1.plot(vel, conti, color='r', lw=1.2)

    ax1.set_xlim([-velwidth, velwidth])
    ax1.minorticks_on()
    ax1.set_xticklabels([])

    #### ax2, normalized spec
    # flux array
    if nbin > 1:
        x2, y2 = bin_spec(vel, nflux, nbin)
        xsig2, ysig2 = bin_spec(vel, nsig, nbin)
    else:
        x2, y2 = vel, nflux
        xsig2, ysig2 = vel, nsig

    xx2, yy2 = np.repeat(x2, 2)[1:], np.repeat(y2, 2)[:-1]
    xxsig2, yysig2 = np.repeat(xsig2, 2)[1:], np.repeat(ysig2, 2)[:-1]

    ax2.hlines(1.0, -velwidth, velwidth, linestyle=':')
    ax2.plot(xx2, yy2, color='k', lw=0.7)
    ax2.plot(xxsig2, yysig2, color='b', lw=0.7)

    ax2.set_xlim([-velwidth, velwidth])
    ax2.set_ylim(0., 1.8)
    ax2.vlines(0, 0, 1.8, linestyle='--')
    ax2.set_xticks(np.mgrid[-velwidth:velwidth+1:500])
    ax2.set_yticks(np.mgrid[0:2.:0.5])
    ax2.minorticks_on()

    ##### ax3, zoom in normalized spec
    ax3.hlines(1.0, pltrange[0], pltrange[1], linestyle=':')
    ax3.plot(xx2, yy2, color='k')
    ax3.plot(xxsig2, yysig2, color='b')
    ax3.set_xlim(pltrange)
    ax3.set_ylim(0, 1.8)
    ax3.vlines(0, 0, 1.8, linestyle='--')
    minr = (pltrange[0]//100)*100
    maxr = (pltrange[1]//100+1)*100
    ax3.set_xticks(np.mgrid[minr:maxr:200][1:])
    ax3.set_yticks(np.mgrid[0:2.:0.5])
    ax3.minorticks_on()
    ax3.set_xlabel('Velocity (km/s)')
    ax3.set_ylabel('Norm. Flux')
    fig.savefig(figname)
    plt.close()
    return figname

