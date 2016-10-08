import numpy as np
import h5py as hp
import re
from astropy.constants import c as speed_of_light_ms # speed of light, in m/s

def save_SpecAll(wave, flux, err, continuum=None, filename='specall.h5'):
    ispecname = filename+'.h5'
    ispec = hp.File(filename, 'w')
    ispec.create_dataset('WAVE', dataset=wave)
    ispec.create_dataset('FLUX', dataset=flux)
    ispec.create_dataset('SIGMA', dataset=err)
    ispec.create_dataset('CONTINUUM', dataset=continuum)

    # add normalized flux file 
    if continuum != None:
        # to avoid 0 in the denominator
        ind0 = np.where(continuum == 0.)[0]
        continuum[ind0] = 1.
        normflux = flux / continuum
        normflux[ind0] = np.nan
        ispec.create_dataset('NORMFLUX', data=normflux)
    else: 
        ispec.create_dataset('NORMFLUX', data=None)

    ispec.close()
    return ispecname

def slice_OneLine(line, rest_lambda, wave, flux, err, continuum=None, 
                  redshift=0.0, targetname=None):
    
    if targetname == None: filename = '%s.h5'%(line.replace(' ', ''))
    else: filename = '%s_%s.h5'%(targetname, line.replace(' ', ''))

    # find the corresponding pixel of the line in the spectra
    obs_lambda = rest_lambda * (1 + redshift)
    obs_ind = np.argmin(np.fabs(wave - obs_lambda))    

    # slice about +/- 1000 km/s from the line center, obs_lambda
    velstep = 1000.
    lf_lambda = obs_lambda - velstep/(speed_of_light_ms/1e3)
    rt_lambda = obs_lambda + velstep/(speed_of_light_ms/1e3)

    ind = np.where(np.all([wave>=lf_lambda, wave<=rt_lambda], axis=0)==True)[0]
    
    ispec = hp.File(filename, 'w')
    ispec.create_dataset('WAVE', dataset=wave[ind])
    ispec.create_dataset('FLUX', dataset=flux[ind])
    ispec.create_dataset('SIGMA', dataset=err[ind])

    # add normalized flux file 
    if continuum != None:
        ispec.create_dataset('CONTINUUM', dataset=continuum[ind])
        # to avoid 0 in the denominator
        ind0 = np.where(continuum == 0.)[0]
        continuum[ind0] = 1.
        normflux = flux / continuum
        normflux[ind0] = np.nan
        ispec.create_dataset('NORMFLUX', data=normflux)

    else:
        ispec.create_dataset('CONTINUUM', data=None)
        ispec.create_dataset('NORMFLUX', data=None)

    ispec.close()
    return filename

def slice_SpecLine(wave, flux, err, continuum=None, lines='All', redshift=0.0, 
                   targetname=None):

    # read in the line library 
    all_lin = fits.getdata('files/xidl_all_lin.fits')

    # Lines that covered in G130 M
    if lines == 'All': 
        liblines = ['FeII 1142', 'FeII 1143', 'FeII 1144', 'PII 1152', 'OI 1302',  
                    'SII 1250', 'SII 1253', 'SII 1259', 'SiIV 1393', 'SiIV 1402', 
                    'CII 1334', 'SiII 1190', 'SiII 1193', 'SiII 1260', 'SiIII 1206'] 
    else: 
        # Do spec slicing for some particular lines
        # lines could be just one: lines = ['FeII 1142']
        #              or several: lines = ['SiII 1250', 'SiII 1253', 'SiII 1259']
        
        # reformat the lines, in case it is not in forms, like, 'FeII 1142', but 'Fe II 1142'
        liblines
        for il in lines:
            eles = re.split('(\d+)', il.replace(' ', '')) 
            newline = '%s %s'%(eles[0], eles[1])
            liblines.append(newline)

    # read line wavelength and oscillator strength
    line_lambda = []
    line_fval = []
    for i in range(len(liblines)):
        for j in range(len(all_lin)):
            if liblines[i] == all_lin[j][0]: 
                print all_lin[j]
                line_lambda.append(all_lin[j][1])
                line_fval.append(all_lin[j][2])
    line_fval[2] = 0.083 # FeII 1144. Check emails exchange with X.   

    # Do spec slicing for all the lines.
    totfile = []
    for il in range(len(liblines)): 
        rfile = slice_OneLine(liblines[il], line_lambda[il], wave, flux, err, 
                              continuum=continuum, redshift=redshift, targetname=targetname)        
        totfile.append(rfile)
   
    return totfile
