import datetime
import numpy as np
import astropy.io.fits as fits
from astropy.constants import c as speed_of_light_ms # speed of light, in m/s
from yzSpec.read_linelibrary import read_linelibrary

def save_spec(lt_spec, target_info, has_continuum=False, line='none', filedir='', do_redshift=False, velwidth=1000):
    '''
    Save the spectra into fits. Could save the whole spectra, or a segment near a specific line. 
   
    lt_spec: linetools.spectra.xspectrum1d.XSpectrum1D object
    target_info: a dict object that contains ('NAME', 'z', 'RA', 'DEC', 'l', 'b', 'S/N', 'DATAFILE')
    has_continuum: if the XSpectrum1D.fit_continuum is performed, this should be True.
    line: default to 'none' if want to save the whole lt_spec; 
          otherwise set to, e.g., line='SII1250', to save a segment of the spec within 
          +/-1000 km/s of the line centroid. In this case, the fits header primary hdu
          would have the line info, i.e., line name, rest wavelength, and f value. 
    '''

    wave = lt_spec.wavelength.value
    flux = lt_spec.flux.value
    sig = lt_spec.sig.value

    if has_continuum == True:
        # to avoid 0 in the denominator
 
        continuum = lt_spec.co.value
        indnan = np.where(continuum == 0.)
        continuum[indnan] = 1.

        # normlized flux
        normflux = flux / continuum
        normflux[indnan] = np.nan

        # normlized error
        normsig = sig / continuum
        normsig[indnan] = np.nan  # still keep the nan since we want the velocity array to be continous

    else:
        continuum = np.zeros(wave.size)+np.nan
        normflux = np.zeros(wave.size)+np.nan
        normsig = np.zeros(wave.size)+np.nan

    if line != 'none': 
        # only save a specific line
        # read the line library, the accepted line name, wavelength, (and oscillator strength)
        line_info = read_linelibrary(lines=line, doprint=False)
        rest_lambda = line_info[1][0]
        if do_redshift == True:
            obs_lambda = rest_lambda*(1+target_info['z']) 
        else: 
            obs_lambda = rest_lambda*(1+0.0) # MW line
  
        # velocity vector 
        velocity = (wave-obs_lambda)/obs_lambda*(speed_of_light_ms.value/1e3)
    
        # slice +/- 1000 km/s from the line center, obs_lambda
        lf_lambda = obs_lambda-velwidth/(speed_of_light_ms.value/1e3)*obs_lambda
        rt_lambda = obs_lambda+velwidth/(speed_of_light_ms.value/1e3)*obs_lambda
        slice_ind = np.where(np.all([wave>=lf_lambda, wave<=rt_lambda], axis=0)==True)[0]
    else: 
        slice_ind = np.arange(wave.size) 

    # check if the data exist. 
    if wave[slice_ind].size <= 3: 
        return 'none'  
    else:  
        # save the data; first generate the header
        todate = str(datetime.datetime.now())
        hdulist = fits.open(target_info['DATAFILE'])
        head0 = hdulist[0].header
        head0['HISTORY'] = 'YZ added continuum using linetools. %s'%(todate)
   
        col_wave = fits.Column(name='WAVE', format='D', array=wave[slice_ind])
        col_flux = fits.Column(name='FLUX', format='D', array=flux[slice_ind])
        col_sig = fits.Column(name='ERROR', format='D', array=sig[slice_ind])
        col_cont = fits.Column(name='CONTINUUM', format='D', array=continuum[slice_ind])
        col_normflux = fits.Column(name='NORMFLUX', format='D', array=normflux[slice_ind])
        col_normsig = fits.Column(name='NORMERR', format='D', array=normsig[slice_ind])
        col_arrs = [col_wave, col_flux, col_sig, col_cont, col_normflux, col_normsig]
       
        # consider the case if just one line is needed 
        if line != 'none': 
            head0['LINE'] = line_info[0][0]
            head0['LAMBDA'] = line_info[1][0]
            head0['FVAL'] = line_info[2][0]
            head0['LINEREF'] = line_info[3]
            filename = '%s/%s_%s.fits.gz'%(filedir, target_info['NAME'], line.replace(' ', ''))
    
            col_velo = fits.Column(name='VELOCITY', format='D', array=velocity[slice_ind])
            col_arrs.append(col_velo)
        else: 
            filename = '%s/%s.fits.gz'%(filedir, target_info['NAME'])
    
        cols = fits.ColDefs(col_arrs)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        prihdu = fits.PrimaryHDU(header=head0)
        thdulist = fits.HDUList([prihdu, tbhdu])
        thdulist.writeto(filename, clobber=True)
        
        return filename
