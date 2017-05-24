import sys, os
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from yzGALFAHI.get_cubeinfo import get_cubeinfo
homedir = os.path.expanduser('~')

def find_the_cube(ra, dec, observation):
    '''
    Use the input (ra, dec) to decide which cube the data locates.
    '''
    
    if homedir == '/Users/Yong': datadir = '/Volumes/YongData2TB/'+observation
    else: datadir = '/home/yzheng/Work/'+observation

    clt = Table.read('%s/%s_RADEC.dat'%(datadir, observation), format='ascii') 
    cramin, cramax = clt['min_ra'], clt['max_ra']
    cdcmin, cdcmax = clt['min_dec'], clt['max_dec']
    indra = np.all([ra>cramin, ra<cramax], axis=0)
    inddc = np.all([dec>cdcmin, dec<cdcmax], axis=0)
    indall = np.where(np.all([indra, inddc], axis=0) == True)
    
    cubename = clt['cubename'][indall]
    if len(cubename)==0: 
        print("No corresponding HI cube in %s"%(observation))
        return '', ''
    else:
        cubename = cubename[0]
        if observation == 'EBHIS': cubedir = datadir+'/'+cubename+'.fit'
        else: cubedir = datadir+'/'+cubename+'.fits'

        return cubename, cubedir

def locate_category(targetname):
    ## find the target in either QSOALS, GALAXy, or STAR_EARLY
    listdir = homedir+'/Dropbox/HSLA_Feb16/targetlists'
    qsolist = Table.read(listdir+'/QSOALS_sample.txt', format='ascii')
    if targetname in qsolist['ID']: return 'QSOALS'

    gallist = Table.read(listdir+'/GALAXY_sample.txt', format='ascii')
    if targetname in gallist['ID']: return 'GALAXY'

    strlist = Table.read(listdir+'/STAR_EARLY_sample.txt', format='ascii')
    if targetname in strlist['ID']: return 'STAR_EARLY'

def find_targetinfo(targetname):
    category = locate_category(targetname)
    objs = Table.read('%s/Dropbox/HSLA_Feb16/targetlists/%s_sample.txt'%(homedir, category), format='ascii')
    ind = np.where(targetname == objs['ID'])[0]

    if len(ind) == 0: qra, qdc = np.nan, np.nan
    else: qra, qdc = objs['RA'][ind][0], objs['DEC'][ind][0]
    return qra, qdc

def extract_HI(ra=-999, dec=-999, targetname=None, filedir='.', observation='HI4PI'):
    '''
    To obtain the corresponding HI data for the QSO sightlines. 
    
    Can be used to obtain from the EBHIS cube. EBHIS has res of 10.8 arcmin, 
    each pixel has 3.25 arcmin. 
    
    If qsoname is set, use the name to find the ra, dec
    Otherwise use the provided ra, dec. 
    '''

    # find the ra, dec for the qso; othewise use the provided ra, dec
    if ra == -999 or dec == -999: ra, dec = find_targetinfo(targetname)
    
    ### use the input ra/dec to decide which cube to explore. 
    cubename, cubefile = find_the_cube(ra, dec, observation)
    if cubename is np.nan: return []
    else:
        cubehdr = fits.open(cubefile)[0]
        cra, cdec, cvel = get_cubeinfo(cubehdr.header)
        cra, cdec = cra[100, :], cdec[:, 100]
        cubedata = cubehdr.data

        minra  = np.nanmin(np.fabs(ra-cra))
        mindc = np.nanmin(np.fabs(dec-cdec))
        
        indra = np.where(np.fabs(ra-cra)==minra)[0][0]
        inddc = np.where(np.fabs(dec-cdec)==mindc)[0][0]
        
        ispec = cubedata[:, inddc, indra] 
        infoHI = (observation, cubename, cra[indra], cdec[inddc])
        
        # save the spectrum 
        savedir = filedir + '/lines'
        if os.path.isdir(savedir) is False: os.makedirs(savedir)
        if targetname == None: 
            filename = '%s/%s_RA%.2f_DEC%.2f.dat'%(savedir, infoHI[0], infoHI[2], infoHI[3])
        else:
            filename = '%s/%s_HI21cm_%s.dat'%(savedir, targetname, observation)
        dataset = np.c_[cvel, ispec]
        np.savetxt(filename, dataset, fmt='%15.4f  %15.4f', header="VEL         FLUX")  

        # add additional information to the end
        f = open(filename, 'a')
        f.write('### %s cube: %s\n'%(observation, cubename))
        f.write('### RA :%.2f degree    pixel index :%d\n'%(cra[indra], indra))
        f.write('### DEC:%.2f degree    pixel index:%d\n'%(cdec[inddc], inddc))
        f.close()

        return infoHI
