import os
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def build_axes(line_number, pltrange=[-400, 400]):
    '''
    Setup the figure axes for the stack plotting. If line_number>10, use the bigger canvas.
    '''

    # the maximum is 40 sub figures. 
    if line_number > 11: 
        axpos40 = np.asarray([[0.08, 0.8], [0.08, 0.7], [0.08, 0.6], [0.08, 0.5],
                              [0.08, 0.4], [0.08, 0.3], [0.08, 0.2], [0.08, 0.1],
                              [0.26, 0.8], [0.26, 0.7], [0.26, 0.6], [0.26, 0.5],
                              [0.26, 0.4], [0.26, 0.3], [0.26, 0.2], [0.26, 0.1],
                              [0.44, 0.8], [0.44, 0.7], [0.44, 0.6], [0.44, 0.5],
                              [0.44, 0.4], [0.44, 0.3], [0.44, 0.2], [0.44, 0.1],
                              [0.62, 0.8], [0.62, 0.7], [0.62, 0.6], [0.62, 0.5],
                              [0.62, 0.4], [0.62, 0.3], [0.62, 0.2], [0.62, 0.1],
                              [0.80, 0.8], [0.80, 0.7], [0.80, 0.6], [0.80, 0.5],
                              [0.80, 0.4], [0.80, 0.3], [0.80, 0.2], [0.80, 0.1]])
        do_xlabel=[0, 0, 0, 0, 0, 0, 0, 1,
                   0, 0, 0, 0, 0, 0, 0, 1,
                   0, 0, 0, 0, 0, 0, 0, 1,
                   0, 0, 0, 0, 0, 0, 0, 1,
                   0, 0, 0, 0, 0, 0, 0, 1]
        do_ylabel=[1, 1, 1, 1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0]
    
        if line_number > axpos40.size:
            logger.info("Too many lines. Only plot the first 40.")
            pltax = axpos40
            do_xlabel = do_xlabel
            do_ylabel = do_ylabel
        else:
            pltax = axpos40[0:line_number]
            do_xlabel = do_xlabel[0:line_number]
            do_xlabel[-1] = 1
            do_ylabel = do_ylabel[0:line_number]
        
        fig = plt.figure(figsize=(10, 8))
        axwd, axht = 0.17, 0.09
    else:   # for only a small set of lines
        axpos11 = np.asarray([[0.2, 0.780], [0.2, 0.705], [0.2, 0.630], [0.2, 0.555], [0.2, 0.480], [0.2, 0.405],
                              [0.2, 0.330], [0.2, 0.255], [0.2, 0.180], [0.2, 0.105], [0.2, 0.030]])
        do_xlabel = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        do_ylabel = [1]*11
        fig = plt.figure(figsize=(2.5, 9))
        axwd, axht = 0.72, 0.07
        pltax = axpos11[0:line_number]
        do_xlabel = do_xlabel[0:line_number]
        do_xlabel[-1] = 1
        do_ylabel = do_ylabel[0:line_number]

    axes = []
    for i in range(len(pltax)):
        iaxpos, dx, dy = pltax[i], do_xlabel[i], do_ylabel[i]
        iax = fig.add_axes([iaxpos[0], iaxpos[1], axwd, axht])
        minr = (pltrange[0]//100)*100
        maxr = (pltrange[1]//100+1)*100
        iax.set_xticks(np.mgrid[minr:maxr:200][1:])
        iax.minorticks_on()
        if dx == 0: iax.set_xticklabels([])
        if dy == 0: iax.set_yticklabels([])
        if i>0: iax.set_yticks([0, 0.5, 1.0, 1.5])
        iax.minorticks_on()
        axes.append(iax)

    return axes, fig

def stack_spec(target_info, filedir, lines='All', pltrange=[-400, 400], vline=0., nbin=1, savedir='./'):
    '''
    Find all the available sliced lines in filedir, and stack them together. 
    '''

    linefiles = os.listdir(filedir)
    newfiles = []
    for iline in linefiles:
        line_ele = iline.split('_')
        if 'HI21cm' not in line_ele:
            newfiles.append(iline)
    newfiles.append('%s_HI21cm_HI4PI_Beam16arcmin.fits.gz'%(target_info['NAME']))
    linefiles = newfiles
    linefiles.sort()

    ## only stack a few lines, but always have HI
    if lines != 'All':
        if 'HI21cm' not in lines: lines.append('HI21cm')
        customized_files = []
        for iline in lines:
            for ifile in linefiles:
                if iline in ifile: 
                    customized_files.append(ifile)

        linefiles = customized_files
        
    if 'HI21cm' in '-'.join(linefiles): axnumber = len(linefiles)
    else: axnumber = len(linefiles)+1

    axes, fig = build_axes(axnumber)
    vmin, vmax = pltrange[0], pltrange[1]
    # plot the HI 21cm line if exist
    if 'HI21cm' in '-'.join(linefiles):
        for ifile in linefiles: 
            if 'HI21cm' in ifile: 
                break
        hitb = fits.open('%s/%s'%(filedir, ifile))
        hivel, hispec = hitb[1].data.field('VLSR'), hitb[1].data.field('FLUX')
        ind = np.where(np.all([hivel>=vmin, hivel<=vmax], axis=0) == True)
        iax = axes[0]
        xx, yy = np.repeat(hivel[ind], 2)[1:], np.repeat(hispec[ind], 2)[:-1]
        iax.plot(xx, yy, color='k', lw=0.8)
        tmax, tmin = np.nanmax(hispec[ind]), np.nanmin(hispec[ind])
        iax.hlines(0, vmin, vmax, linestyle=':')
        iax.vlines(vline, tmin, tmax*1.4, linestyle='--')
        iax.set_xlim(vmin, vmax)
        iax.set_ylim(tmin, tmax*1.4)
        iax.text(vmin+0.05*np.fabs(vmax-vmin), tmax, 'HI-21cm', color='r')
        fig.text(0.05, 0.9, '(%s)'%(ifile), fontsize=8)
        hitb.close()
    else:
        iax = axes[0]
        iax.set_xlim(vmin, vmax)
        iax.set_ylim(-1, 1)

    # now plot other ion lines
    nion = 0
    for ifile in linefiles:
        if ifile[0] == '.': continue
        if 'HI21cm' in ifile: continue
        nion = nion+1
        iontb = fits.open('%s/%s'%(filedir, ifile))
        from quicktools.vhelio2vlsr import vhelio2vlsr_Westmeier
        vhel2vlsr = vhelio2vlsr_Westmeier(0, target_info['l'], target_info['b'], doradec=False)
        ivel = iontb[1].data.field('VELOCITY')+vhel2vlsr
        iflux = iontb[1].data.field('NORMFLUX')
        isig = iontb[1].data.field('NORMERR')
        iax = axes[nion]
        ind = np.where(np.all([ivel>=vmin, ivel<=vmax], axis=0) == True)[0]
        if len(ivel[ind])<=3:
            iax.set_xlim(vmin, vmax)
            iax.set_ylim(0., 1.8)
            iax.text(vmin+0.05*np.fabs(vmax-vmin), 1.4, ifile.split('_')[-1][:-4], color='r')
        else:
            if nbin>1:
                x, y = bin_spec(ivel[ind], iflux[ind], nbin)
                z = bin_spec(ivel[ind], isig[ind], nbin)[1]
            else:
                x, y, z = ivel[ind], iflux[ind], isig[ind]

            xx, yy = np.repeat(x, 2)[1:], np.repeat(y, 2)[:-1]
            zz = np.repeat(z, 2)[:-1]
            iax.plot(xx, yy, color='k', lw=0.8)
            iax.plot(xx, zz, color='b', lw=0.8)
            iax.hlines(1., vmin, vmax, linestyle=':')
            iax.vlines(vline, 0., 1.8, linestyle='--')
            iax.set_xlim(vmin, vmax)
            iax.set_ylim(0, 1.8)
            from yzSpec.read_linelibrary import read_linelibrary
            thisline = read_linelibrary(lines=ifile.split('_')[-1][:-4], doprint=False)
            iax.text(vmin+0.05*np.fabs(vmax-vmin), 1.3, '%s        f%.4f'%(iontb[0].header['LINE'],
                                                                           iontb[0].header['FVAL']),
                                                                           color='r', fontsize=10)
            iontb.close()

    fig.text(0.05, 0.96, '%s'%(target_info['NAME']), fontsize=12, horizontalalignment='left')
    fig.text(0.05, 0.94, '(RA%.2f, DEC%.2f, gl%.2f, gb%.2f)'%(target_info['RA'], target_info['DEC'],
                                                              target_info['l'], target_info['b']),
                                                              fontsize=8)
    fig.text(0.05, 0.92, '(SN%d, z%.3f)'%(target_info['S/N'], target_info['z']), fontsize=8)
    # figname = '%s/%s_stackspec_bin%d.pdf'%(savedir, target_info['NAME'], nbin)
    figname = '%s/%s_stackspec.pdf'%(savedir, target_info['NAME'])  # for the final product
    fig.savefig(figname)
    plt.close()
    return 

def plot_OneLine(target_info, linefile, pltrange=[-400, 400], nbin=1, filedir='.', velwidth=1000):
    '''
    Plot and (evaluate) the continuum fits.
    '''

    import astropy.io.fits as fits
    spechdu = fits.open(linefile)
    
    line = spechdu[0].header['LINE']
    fval = spechdu[0].header['FVAL']

    flux = spechdu[1].data.field('FLUX')
    sig = spechdu[1].data.field('ERROR')
    conti = spechdu[1].data.field('CONTINUUM')
    nflux = spechdu[1].data.field('NORMFLUX')
    nsig = spechdu[1].data.field('NORMERR')
    vel = spechdu[1].data.field('VELOCITY')
    spechdu.close()

    from quicktools.vhelio2vlsr import vhelio2vlsr_Westmeier
    vhel2vlsr = vhelio2vlsr_Westmeier(0, target_info['l'], target_info['b'], doradec=False)
    vel = vel+vhel2vlsr

    fig = plt.figure(figsize=(6, 6))
    ax1 = fig.add_axes([0.12, 0.48, 0.82, 0.41])   # original spec
    ax2 = fig.add_axes([0.12, 0.27, 0.82, 0.2])   # normalized
    ax3 = fig.add_axes([0.12, 0.08, 0.82, 0.15])

    fig.text(0.94, 0.95, '%s'%(target_info['NAME']), fontsize=12, horizontalalignment='right')
    fig.text(0.94, 0.905, '%s  fval=%.4f'%(line, fval), fontsize=12, 
             horizontalalignment='right', color='r', fontweight='bold')

    # ax1, original spectra, with continuum on top
    ax1.plot(vel, flux, color='k', lw=0.7)
    ax1.plot(vel, sig, color='b', lw=0.7)
    ax1.plot(vel, conti, color='r', lw=1.2)

    ax1.set_xlim([-velwidth, velwidth])
    ax1.minorticks_on()
    ax1.set_xticklabels([])
    ax1.set_ylabel('Flux (erg/s/cm2/A)')

    ax2.hlines(1.0, -velwidth, velwidth, linestyle=':')
    ax2.plot(vel, nflux, color='k', lw=0.7)
    ax2.plot(vel, nsig, color='b', lw=0.7)
    ax2.set_xlim([-velwidth, velwidth])
    ax2.set_ylim(0., 1.8)
    ax2.vlines(0, 0, 1.8, linestyle='--')
    ax2.set_xticks(np.mgrid[-velwidth:velwidth+1:500])
    ax2.set_yticks(np.mgrid[0:2.:0.5])
    ax2.set_ylabel('Norm. Flux')
    ax2.minorticks_on()

    ##### ax3, zoom in normalized spec
    ax3.hlines(1.0, pltrange[0], pltrange[1], linestyle=':')
    ax3.plot(vel, nflux, color='k', lw=0.7)
    ax3.plot(vel, nsig, color='b', lw=0.7)
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

    figname = '%s/%s_%s.pdf'%(filedir, target_info['NAME'], line.replace(' ', ''))
    fig.savefig(figname)
    plt.close()
