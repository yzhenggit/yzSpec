import numpy as np
import h5py as hp
import re, os
from astropy.table import Table
import astropy.io.fits as fits
from astropy.constants import c as speed_of_light_ms # speed of light, in m/s
import matplotlib.pyplot as plt
from yzSpec.plot_spec import plot_OneLine

def stat_OneLine(dataset):
    '''
    Evaluate line continuuum fitting with reduced Chi2 value.
    '''
    
    #wave[ind], flux[ind], err[ind], continuum[ind], normflux[ind], normerr[ind], velocity[ind]

def save_SpecAll(wave, flux, sig, continuum=None, save_ascii=False, save_fits=True,
                 targetname='TargetAll', filedir=''):
    
    '''
    Save the Continuum fits into ascii files.
    '''

    if continuum != None:
        # to avoid 0 in the denominator
        indnan = np.where(continuum == 0.)
        continuum[indnan] = 1.

        # normlized flux
        normflux = flux / continuum
        normflux[indnan] = np.nan

        # normlized error
        normsig = sig / continuum
        normsig[indnan] = np.nan

    else:
        continuum = np.empty(wave.size)+np.nan
        normflux = np.empty(wave.size)+np.nan
        normsig = np.empty(wave.size)+np.nan
       
    # save the data
    if save_ascii == True:
        dataset = np.c_[wave, flux, sig, continuum, normflux, normsig]
        filename = filedir+'/'+targetname+'.dat'
        np.savetxt(filename, dataset, fmt='%15.4f  %15.3e  %15.3e  %15.3e  %15.4f  %15.4f',
                   header="WAVE           FLUX             SIG            CONTINUUM        NORMFLUX         NORMSIG")
    if save_fits == True:
        col_wave = fits.Column(name='WAVE', format='D', array=wave)
        col_flux = fits.Column(name='FLUX', format='D', array=flux)
        col_sig = fits.Column(name='SIG', format='D', array=sig)
        col_cont = fits.Column(name='CONTINUUM', format='D', array=continuum)
        col_normflux = fits.Column(name='NORMFLUX', format='D', array=normflux)
        col_normsig = fits.Column(name='NORMSIG', format='D', array=normsig)
        
        cols = fits.ColDefs([col_wave, col_flux, col_sig, col_cont, col_normflux, col_normsig])
        tbhdu = fits.BinTableHDU.from_columns(cols)

        import datetime
        todate = str(datetime.datetime.now())
        head0 = hdulist[0].header
        head0['HISTORY'] = 'YZ added continuum with linetools. %s'%(todate)
        prihdu = fits.PrimaryHDU(header=head0)
        thdulist = fits.HDUList([prihdu, tbhdu])
        filename = filedir+'/'+targetname+'.fits.gz'
        thdulist.writeto(filename)

    return filename

def slice_OneLine(line, rest_lambda, wave, flux, err, continuum=None,
                  redshift=0.0, targetname=None, filedir='.', 
                  velwidth=1000., pltrange=[-200, 200]):

    # bulid/check the folder that contain lines and figures
    linedir = filedir + '/lines'
    figdir = filedir + '/figs'

    if targetname == None: 
        filename = '%s/%s.dat'%(linedir, line.replace(' ', ''))
        figname = '%s/%s.pdf'%(figdir, line.replace(' ', ''))
    else: 
        filename = '%s/%s_%s.dat'%(linedir, targetname, line.replace(' ', ''))
        figname = '%s/%s_%s.pdf'%(figdir, targetname, line.replace(' ', ''))

    # find the corresponding pixel of the line in the spectra
    obs_lambda = rest_lambda * (1 + redshift)
    obs_ind = np.argmin(np.fabs(wave - obs_lambda))

    # check if enough pixels for IVC
    lll = obs_lambda - 100/(speed_of_light_ms.value/1e3)*obs_lambda
    rrr = obs_lambda + 100/(speed_of_light_ms.value/1e3)*obs_lambda
    indlr = np.where(np.all([wave>=lll, wave<=rrr], axis=0)==True)[0]
    if indlr.size <=3.: # no data, or very few data, then don't plot
        print("Not enough pixels to save/plot. Skip this line: ", line)
        return False

    if continuum != None:
        # to avoid 0 in the denominator
        indnan = np.where(continuum == 0.)
        continuum[indnan] = 1.
        
        # normlized flux and error
        normflux = flux / continuum
        normflux[indnan] = np.nan
        normerr = err / continuum
        normerr[indnan] = np.nan
        
    else: 
        continuum = np.empty(wave.size)+np.nan
        normflux = np.empty(wave.size)+np.nan
        normerr = np.empty(wave.size)+np.nan

    # slice about +/- 1000 km/s from the line center, obs_lambda
    lf_lambda = obs_lambda - velwidth/(speed_of_light_ms.value/1e3)*obs_lambda
    rt_lambda = obs_lambda + velwidth/(speed_of_light_ms.value/1e3)*obs_lambda
    ind = np.where(np.all([wave>=lf_lambda, wave<=rt_lambda], axis=0)==True)[0]

    # velocity vector 
    velocity = (wave - obs_lambda)/obs_lambda*(speed_of_light_ms.value/1e3)      

    # save the data
    dataset = np.c_[wave[ind], flux[ind], err[ind], continuum[ind], 
                    normflux[ind], normerr[ind], velocity[ind]]    
    np.savetxt(filename, dataset, fmt='%15.4f  %15.3e  %15.3e  %15.3e  %15.4f  %15.4f  %15.4f', 
               header="WAVE           FLUX             SIGMA            CONTINUUM        NORMFLUX         NORM_SIG         VEL_z")
    
    redchi2 = stat_OneLine(dataset.transpose())
    plot_OneLine(dataset.transpose(), redchi2=redchi2, figname=figname, 
                 targetname=targetname, linename=line, velwidth=velwidth, 
                 pltrange=pltrange)

    return figname

def slice_SpecLine(wave, flux, err, continuum=None, lines='All', coords=[0, 0, 0, 0, 0], 
                   redshift=0.0, targetname=None, filedir='.', doHI=False, 
                   pltrange=[-200, 200], velwidth=1000):

    from yzSpec.read_linelibrary import read_linelibrary
    # read the line library, the accepted line name, wavelength, (and oscillator strength)
    line_library = read_linelibrary(lines=lines)
    liblines, line_lambda = line_library[0], line_library[1]  

    # Do spec slicing for all the lines.
    totfile = []
    for il in range(len(liblines)):
        rfile = slice_OneLine(liblines[il], line_lambda[il], wave, flux, err,
                              continuum=continuum, redshift=redshift, 
                              targetname=targetname, filedir=filedir, 
                              velwidth=velwidth, pltrange=pltrange)
        totfile.append(rfile)

    from yzSpec.combine_figpdfs import combine_figpdfs
    outname = targetname+'_stackcont.pdf'
    combine_figpdfs(filedir+'/figs', outputdir=filedir, outpdfname=outname)
    return outname
