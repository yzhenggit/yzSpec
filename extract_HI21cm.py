import os
import numpy as np
from astropy.table import Table
import astropy.io.fits as fits
from yzGALFAHI.get_cubeinfo import get_cubeinfo
from astropy.coordinates import SkyCoord

import warnings
warnings.filterwarnings('ignore')

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def find_the_cube(ra, dec, observation):
    '''
    Use the input (ra, dec) to decide which cube the data locates.

    observation: most of the time, it is HI4PI. 
    '''

    if observation == 'HI4PI': 
        datadir = '/Volumes/YongData2TB/'+observation
    elif observation == 'GALFA_HI': 
        datadir = '/Volumes/YongData2TB/GALFAHI_DR2/DR2W_RC5/DR2W_RC5/Wide'

    clt = Table.read('%s/%s_RADEC.dat'%(datadir, observation), format='ascii')
    cramin, cramax = clt['min_ra'], clt['max_ra']
    cdcmin, cdcmax = clt['min_dec'], clt['max_dec']
    indra = np.all([ra>cramin, ra<cramax], axis=0)
    inddc = np.all([dec>cdcmin, dec<cdcmax], axis=0)
    indall = np.where(np.all([indra, inddc], axis=0) == True)

    cubename = clt['cubename'][indall]
    if len(cubename)==0:
        logger.info("No corresponding HI cube in %s"%(observation))
        return '', ''
    else:
        cubename = cubename[0]
        if observation == 'EBHIS': cubedir = datadir+'/'+cubename+'.fit'
        else: cubedir = datadir+'/'+cubename+'.fits'

        return cubename, cubedir


def extract_HI21cm(target_info, filedir='.', observation='HI4PI', beam=1.):
    '''
    To obtain the corresponding HI data for the QSO sightlines. Can be used 
    to obtain from HI4PI (EBHIS+GASS) cubes. HI4PI has res of 10.8 arcmin, 
    each pixel has 3.25 arcmin. 

    beam: to decide within what diameter (in deg) the HI spec is averaged. 
    '''

    target = target_info['NAME']
    beam_radius = beam/2.

    if observation == 'LAB':
        labfile = '/Users/Yong/Dropbox/databucket/LAB/labh_glue.fits'
        labdata = fits.getdata(labfile)
        labhdr = fits.getheader(labfile)
        
        gl, gb, cvel = get_cubeinfo(labhdr)
        tar_coord = SkyCoord(l=target_info['l'], b=target_info['b'], unit='deg', frame='galactic')
        cube_coord = SkyCoord(l=gl, b=gb, unit='deg', frame='galactic')
        dist = tar_coord.separation(cube_coord)
        
        dist_copy = dist.value.copy()
        within_beam_2d = dist_copy<=beam_radius
        within_beam_3d = np.asarray([within_beam_2d]*cvel.size)
                 
        labdata_copy = labdata.copy()
        labdata_copy[np.logical_not(within_beam_3d)] = np.nan

        ispec = np.nanmean(np.nanmean(labdata_copy, axis=2), axis=1)

        # save the spectrum 
        prihdr = fits.Header()
        prihdr['OBS'] = observation
        prihdr.comments['OBS'] = 'See %s publication.'%(observation)
        prihdr['CREATOR'] = "YZ"
        prihdr['COMMENT'] = "HI 21cm spectrum averaged within beam size of %.2f deg. "%(beam)
        import datetime as dt
        prihdr['DATE'] = str(dt.datetime.now())
        prihdu = fits.PrimaryHDU(header=prihdr)

        ## table 
        col1 = fits.Column(name='VLSR', format='D', array=cvel)
        col2 = fits.Column(name='FLUX', format='D', array=ispec)
        cols = fits.ColDefs([col1, col2])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        thdulist = fits.HDUList([prihdu, tbhdu])

        if os.path.isdir(filedir) is False: os.makedirs(filedir)
        if beam >= 1.:
            hifile = '%s/%s_HI21cm_%s_Beam%ddeg.fits.gz'%(filedir, target, observation, beam)
        else:
            hifile = '%s/%s_HI21cm_%s_Beam%darcmin.fits.gz'%(filedir, target, observation, beam*60)

        thdulist.writeto(hifile, clobber=True)
          
    elif observation in ['HI4PI', 'GALFA_HI']:
        tar_coord = SkyCoord(ra=target_info['RA'], dec=target_info['DEC'], unit='deg')

        ### use the input ra/dec to decide which cube to explore. 
        if observation == 'HI4PI': 
            datadir = '/Volumes/YongData2TB/'+observation
        else: 
            datadir = '/Volumes/YongData2TB/GALFAHI_DR2/DR2W_RC5/DR2W_RC5/Wide'

        clt = Table.read('%s/%s_RADEC.dat'%(datadir, observation), format='ascii')

        cubefiles = []
        for ic in range(len(clt)):
            cubefile = datadir+'/'+clt['cubename'][ic]+'.fits'
            cubehdr = fits.getheader(cubefile)
            cra, cdec, cvel = get_cubeinfo(cubehdr)
            cube_coord = SkyCoord(ra=cra, dec=cdec, unit='deg')
            dist_coord = tar_coord.separation(cube_coord)
        
            dist = dist_coord.value
            within_beam_2d = dist<=beam_radius
            if dist[within_beam_2d].size > 0:
                cubefiles.append(cubefile)
 
        specs = []
        for cubefile in cubefiles:
            cubehdr = fits.getheader(cubefile)
            cra, cdec, cvel = get_cubeinfo(cubehdr)
            cube_coord = SkyCoord(ra=cra, dec=cdec, unit='deg')
            dist_coord = tar_coord.separation(cube_coord)
        
            dist = dist_coord.value
            within_beam_2d = dist<=beam_radius
            within_beam_3d = np.asarray([within_beam_2d]*cvel.size)
        
            cubedata = fits.getdata(cubefile)
            cubedata[np.logical_not(within_beam_3d)] = np.nan
            
            ispec = np.nanmean(np.nanmean(cubedata, axis=2), axis=1)
            specs.append(ispec)
            
        ispec = np.mean(np.asarray(specs), axis=0)

        # save the spectrum 
        prihdr = fits.Header()
        prihdr['OBS'] = observation
        prihdr.comments['OBS'] = 'See %s publication.'%(observation)
        prihdr['CREATOR'] = "YZ"
        prihdr['COMMENT'] = "HI 21cm spectrum averaged within beam size of %.2f deg. "%(beam)
        import datetime as dt
        prihdr['DATE'] = str(dt.datetime.now())
        prihdu = fits.PrimaryHDU(header=prihdr)
        
        ## table 
        col1 = fits.Column(name='VLSR', format='D', array=cvel)
        col2 = fits.Column(name='FLUX', format='D', array=ispec)
        cols = fits.ColDefs([col1, col2])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        thdulist = fits.HDUList([prihdu, tbhdu])
         
        if os.path.isdir(filedir) is False: os.makedirs(filedir)
       
        if beam >= 1.: 
            hifile = '%s/%s_HI21cm_%s_Beam%ddeg.fits.gz'%(filedir, target, observation, beam)
        else:
            hifile = '%s/%s_HI21cm_%s_Beam%darcmin.fits.gz'%(filedir, target, observation, beam*60)
        thdulist.writeto(hifile, clobber=True)

    else: 
        logger.info('%s are not in [LAB, HI4PI]'%(observation)) 
        hifile = ''
    return hifile
