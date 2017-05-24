import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from quicktools.bin_spec import bin_spec
from quicktools.vhelio2vlsr import vhelio2vlsr_Westmeier
from yzSpec.read_linelibrary import read_linelibrary


def build_axes(line_number, pltrange=[-400, 400]):
    # the maximum is 40 sub figures. 
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
        print("Too many lines. Only plot the first 20.")
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
    
def stack_spec(filedir='.', targetname=None, lines='All', pltrange=[-400, 400], hivline=0., 
               redshift=0.0, nbin=3, infoHI=None, coords=[0, 0, 0, 0, 0], figname='none'):
     
    linedir = filedir + '/lines'
    figdir = filedir + '/figs'
    if os.path.isdir(linedir) is False: 
        print(linedir + ' does not exist, please check')
        return np.nan
    if os.path.isdir(figdir) is False: os.makedirs(figdir)

    linefiles = os.listdir(linedir)
    linefiles.sort()

    if 'HI21cm' in '-'.join(linefiles): axnumber = len(linefiles)
    else: axnumber = len(linefiles)+1

    axes, fig = build_axes(axnumber)
    vmin, vmax = pltrange[0], pltrange[1]
    # plot the HI 21cm line if exist
    if 'HI21cm' in '-'.join(linefiles): 
        hitb = Table.read('%s/%s_HI21cm_%s.dat'%(linedir, targetname, infoHI[0]), format='ascii')
        hivel, hispec = hitb.columns[0], hitb.columns[1]
        ind = np.where(np.all([hivel>=vmin, hivel<=vmax], axis=0) == True)
        iax = axes[0]
        xx, yy = np.repeat(hivel[ind], 2)[1:], np.repeat(hispec[ind], 2)[:-1]
        iax.plot(xx, yy, color='k', lw=0.8)
        tmax, tmin = np.nanmax(hispec[ind]), np.nanmin(hispec[ind])
        iax.hlines(0, vmin, vmax, linestyle=':')
        iax.vlines(hivline, tmin, tmax*1.4, linestyle='--')
        iax.set_xlim(vmin, vmax)
        iax.set_ylim(tmin, tmax*1.4)
        iax.text(vmin+0.05*np.fabs(vmax-vmin), tmax, 'HI-21cm', color='r')
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
        iontb = Table.read(linedir+'/'+ifile, format='ascii')
        vhel2vlsr = vhelio2vlsr_Westmeier(0, coords[2], coords[3], doradec=False)
        ivel, iflux, isig = iontb['VEL_z']+vhel2vlsr, iontb['NORMFLUX'], iontb['NORM_SIG']
        iax = axes[nion]
        ind = np.where(np.all([ivel>=vmin, ivel<=vmax], axis=0) == True)[0]
        if len(ivel[ind])<=3: 
            iax.set_xlim(vmin, vmax)
            iax.set_ylim(0., 1.8)
            iax.text(vmin+0.05*np.fabs(vmax-vmin), 1.4, ifile.split('_')[-1][:-4], color='r')
            continue
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
            iax.vlines(hivline, 0., 1.8, linestyle='--')
            iax.set_xlim(vmin, vmax)
            iax.set_ylim(0, 1.8)
            thisline = read_linelibrary(lines=ifile.split('_')[-1][:-4], doprint=False)
            iax.text(vmin+0.05*np.fabs(vmax-vmin), 1.3, '%s        f%.4f'%(ifile.split('_')[-1][:-4], 
                                                                           thisline[-1][0]), 
                                                                           color='r', fontsize=10)

    fig.text(0.08, 0.91, '%s'%(targetname), fontsize=12, horizontalalignment='left')
    fig.text(0.96, 0.94, '(RA%.2f, DEC%.2f, gl%.2f, gb%.2f, SN%d)'%(coords[0], coords[1], 
                                                                   coords[2], coords[3], 
                                                                   coords[4]), fontsize=9, 
                                                                   horizontalalignment='right')
    if infoHI != None:
        fig.text(0.96, 0.91, '(%s, %s, RA%.2f, DEC%.2f)'%(infoHI[0], infoHI[1], infoHI[2], infoHI[3]), 
                 fontsize=9, horizontalalignment='right')
    if figname == 'none': fig.savefig('%s/%s_stackspec_bin%d.pdf'%(filedir, targetname, nbin))
    else: fig.savefig(figname)
    plt.close()
    return 0  

def stack_spec_customized(filedir='.', targetname=None, lines='All', pltrange=[-400, 400], hivline=0.,
                          redshift=0.0, nbin=3, infoHI=None, coords=[0, 0, 0, 0, 0], figname='none'):

    linedir = filedir + '/lines'
    figdir = filedir + '/figs'
    if os.path.isdir(linedir) is False:
        print(linedir + ' does not exist, please check')
        return np.nan
    if os.path.isdir(figdir) is False: os.makedirs(figdir)

    linefiles = os.listdir(linedir)
    linefiles.sort()

    if 'HI21cm' in '-'.join(linefiles): axnumber = len(linefiles)
    else: axnumber = len(linefiles)+1

    axes, fig = build_axes(axnumber)
    vmin, vmax = pltrange[0], pltrange[1]
    # plot the HI 21cm line if exist
    if 'HI21cm' in '-'.join(linefiles):
        hitb = Table.read('%s/%s_HI21cm_%s.dat'%(linedir, targetname, infoHI[0]), format='ascii')
        hivel, hispec = hitb.columns[0], hitb.columns[1]
        ind = np.where(np.all([hivel>=vmin, hivel<=vmax], axis=0) == True)
        iax = axes[0]
        xx, yy = np.repeat(hivel[ind], 2)[1:], np.repeat(hispec[ind], 2)[:-1]
        iax.plot(xx, yy, color='k', lw=0.8)
        tmax, tmin = np.nanmax(hispec[ind]), np.nanmin(hispec[ind])
        iax.hlines(0, vmin, vmax, linestyle=':')
        iax.vlines(hivline, tmin, tmax*1.4, linestyle='--')
        iax.set_xlim(vmin, vmax)
        iax.set_ylim(tmin, tmax*1.4)
        iax.text(vmin+0.05*np.fabs(vmax-vmin), tmax, 'HI-21cm', color='r')
    else:
        iax = axes[0]
        iax.set_xlim(vmin, vmax)
        iax.set_ylim(-1, 1)

    # now plot other ion lines
    nion = 0
    for ifile in linefiles:
        if 'HI21cm' in ifile: continue
        nion = nion+1
        iontb = Table.read(linedir+'/'+ifile, format='ascii')
        vhel2vlsr = vhelio2vlsr_Westmeier(0, coords[2], coords[3], doradec=False)
        ivel, iflux, isig = iontb['VEL_z']+vhel2vlsr, iontb['NORMFLUX'], iontb['NORM_SIG']
        iax = axes[nion]
        ind = np.where(np.all([ivel>=vmin, ivel<=vmax], axis=0) == True)[0]
        if len(ivel[ind])<=3:
            iax.set_xlim(vmin, vmax)
            iax.set_ylim(0., 1.8)
            iax.text(vmin+0.05*np.fabs(vmax-vmin), 1.4, ifile.split('_')[-1][:-4], color='r')
            continue
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
            iax.vlines(hivline, 0., 1.8, linestyle='--')
            iax.set_xlim(vmin, vmax)
            iax.set_ylim(0, 1.8)
            thisline = read_linelibrary(lines=ifile.split('_')[-1][:-4], doprint=False)
            iax.text(vmin+0.05*np.fabs(vmax-vmin), 1.3, '%s        f%.4f'%(ifile.split('_')[-1][:-4],
                                                                           thisline[-1][0]),
                                                                           color='r', fontsize=10)

    fig.text(0.08, 0.91, '%s'%(targetname), fontsize=12, horizontalalignment='left')
    fig.text(0.96, 0.94, '(RA%.2f, DEC%.2f, gl%.2f, gb%.2f, SN%d)'%(coords[0], coords[1],
                                                                   coords[2], coords[3],
                                                                   coords[4]), fontsize=9,
                                                                   horizontalalignment='right')
    if infoHI != None:
        fig.text(0.96, 0.91, '(%s, %s, RA%.2f, DEC%.2f)'%(infoHI[0], infoHI[1], infoHI[2], infoHI[3]),
                 fontsize=9, horizontalalignment='right')
    if figname == 'none': fig.savefig('%s/%s_stackspec_bin%d.pdf'%(filedir, targetname, nbin))
    else: fig.savefig(figname)
    plt.close()
    return 0

