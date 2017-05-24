import os, sys
import numpy as np
from astropy.table import Table
from shutil import copyfile
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('TkAgg')  # backend for graphical interfaces
import linetools.spectra.xspectrum1d as lsx
from yzSpec.slice_spec import slice_SpecLine, save_SpecAll
from yzSpec.extract_HI import extract_HI
from yzSpec.stack_spec import stack_spec, stack_spec_customized
homedir = os.path.expanduser('~')

def question_yes_no(question):
    while True:
        sys.stdout.write(question + '[y/n]\n')
        answer = input()
        if answer == 'y':
            return True
        elif answer == 'n':
            return False
        else:
            sys.stdout.write("Please respond with y or n \n")

def locate_category(targetname):
    '''
    With the input of targetname (provided it is in HSLA), this func
    helps to find the category for the target. Locally, for now, 
    it will search through QSOALS, GALAXY, and STAR_EARLY.
    Return: category name (if founded)
    '''
    qsolist = homedir+'/Dropbox/HSLA_Feb16/targetlists/QSOALS_sample.txt'
    qsotb = Table.read(qsolist, format='ascii')
    if targetname in qsotb['ID']: return 'QSOALS'

    gallist = homedir+'/Dropbox/HSLA_Feb16/targetlists/GALAXY_sample.txt'
    galtb = Table.read(gallist, format='ascii')
    if targetname in galtb['ID']: return 'GALAXY'

    strlist = homedir+'/Dropbox/HSLA_Feb16/targetlists/STAR_EARLY_sample.txt'
    strtb = Table.read(strlist, format='ascii')
    if targetname in strtb['ID']: return 'STAR_EARLY'

def find_targetinfo(targetname, category):
    '''
    Provided the target name, this func would find which category (QSOALS, 
    GALAXY, or STAR_EARLY) it belongs to, and the gratings/coordinates info
    for this target. 
    '''
    objlist = '%s/Dropbox/HSLA_Feb16/targetlists/%s_sample.txt'%(homedir, category)
    objs = Table.read(objlist, format='ascii')
    ind = np.where(targetname == objs['ID'])[0]
    if len(ind) == 0: 
        return 0, 'No such target on the list', [np.nan]
    else:
        qra, qdc = objs['RA'][ind][0], objs['DEC'][ind][0]
        ql, qb = objs['l'][ind][0], objs['b'][ind][0]
        qsn = objs['S/N'][ind][0]
        g130, g160 = objs['G130M'][ind][0], objs['G160M'][ind][0]
        
        objcoords = [qra, qdc, ql, qb, qsn] 
        if g130+g160 == 2: yes_do, gtag = 1, 'FUVM'
        elif g130+g160 == 1:
            if g130 == 1: gtag = 'G130M'
            else: gtag = 'G160M'
            yes_do = 1
        else: yes_do, gtag = 0, 'No G130M or G160M'

        return yes_do, gtag, objcoords

def make_customized_stackspec(targetname, lines='All', pltrange=[-400, 400]):
    '''
    This function is used to make some customized stacked spectra, provided 
    that the line spectra have been processed already.
    The corresponding function stack_spec_customized can be changed accordingly 
    to best suit the purpose of the plot. (So that there is no need to change
    the stack_spec funtion which is mainly used for MWDH.)
    '''
    category = locate_category(targetname)
    datadir = homedir+'/Dropbox/HSLA_Feb16/datapile/'+category
    savedir = homedir+'/Dropbox/HSLA_Feb16/QSOSpec_YZ/'+category

    ## check if QSO exist, and if observed with G130M and/or G160M
    yes_do, gtag, objcoords = find_targetinfo(targetname, category)
    if yes_do == 0: return gtag

    ######## If everthing if fine, let's start! #########
    filedir = savedir+'/'+targetname
    hifile = '%s/lines/%s_HI21cm_HI4PI.dat'%(filedir, targetname)
    if os.path.isfile(hifile) == True:
        print("HI 21cm line has been extracted from HI4PI. ")
        fHI = Table.read('tables/QSOs_HIextracted.dat', format='ascii')
        ind = np.where(targetname == fHI['ID'])[0][0]
        infoHI = [fHI['obsers.'][ind], fHI['cube'][ind],
                  fHI['ra'][ind], fHI['dec'][ind], np.nan, np.nan]

    ## stacking the spectra
    print("Stacking Spectra ... ")
    for nbin in [1, 3]:
        stack_spec_customized(filedir=filedir, targetname=targetname,
                              pltrange=pltrange, lines=lines,
                              coords=objcoords, infoHI=infoHI,
                              nbin=nbin, observation=observation)
    return targetname

def read_knots(knotname):
    '''
    This is to read the **_knots.jsn file for the continuum fitting func.
    '''
    knotfile = open(knotname).readline().split('],')
    knotlist = []
    for iknot in knotfile: 
        a = float(iknot.split(',')[0].replace('[', ''))
        b = float(iknot.split(',')[1].replace(']', ''))
        knotlist.append([a, b])    
    return knotlist

def find_HIspec(targetname, observation='HI4PI'):
    '''
    If the target does not have HI spectrum, will extract it from 
    HI4PI or other designated observations. If has, then the sightline
    coordinate info should be recorded in QSOs_HIextracted.dat; will 
    read the sightline HI info from that file.  
    '''
    category = locate_category(targetname)
    datadir = homedir+'/Dropbox/HSLA_Feb16/datapile/'+category
    savedir = homedir+'/Dropbox/HSLA_Feb16/QSOSpec_YZ/'+category

    ## check if QSO exist, and if observed with G130M and/or G160M
    objcoords = find_targetinfo(targetname, category)[2]
    
    ## save the whole continuum 
    filedir = savedir+'/'+targetname
    if os.path.isdir(filedir) == False:
        os.makedirs(filedir)
        os.makedirs(filedir+'/lines')
        os.makedirs(filedir+'/figs')

    ## check if HI data exist; if not, extract it from HI4PI
    hilist = homedir+'/Dropbox/HSLA_Feb16/code/tables/QSOs_HIextracted.dat'
    hifile = '%s/lines/%s_HI21cm_%s.dat'%(filedir, targetname, observation)
    if os.path.isfile(hifile) == False:
        print('No HI data found. Now extracting from %s ...'%(observation))
        infoHI = extract_HI(targetname=targetname, filedir=filedir, observation=observation)
        qextractHI = open(hilist, 'a')
        qextractHI.write('%25s    %5s   %5s  %7s  %.2f\n'%(targetname,
                                                           infoHI[0], infoHI[1],
                                                           infoHI[2], infoHI[3]))
        qextractHI.close()
    else:
        print("HI 21cm line has been extracted. ")
        fHI = Table.read(hilist, format='ascii')
        ind = np.where(targetname == fHI['ID'])[0][0]
        infoHI = [fHI['obsers.'][ind], fHI['cube'][ind], fHI['ra'][ind], fHI['dec'][ind]]
    return infoHI    
 
def continuumfitting(targetname, lines='All', ask_if_final=True):
    '''
    This is a complicated func... It would check whether the target
    needs to be continuum normalized, then check if HI substraion
    (from HI4PI) is needed, and then save all the data, and make
    the stacked spectra plot for this target. 
    '''
    category = locate_category(targetname)
    datadir = homedir+'/Dropbox/HSLA_Feb16/datapile/'+category
    savedir = homedir+'/Dropbox/HSLA_Feb16/QSOSpec_YZ/'+category

    ## check if QSO exist, and if observed with G130M and/or G160M
    yes_do, gtag, objcoords = find_targetinfo(targetname, category)
    if yes_do == 0: return gtag

    gts = Table.read('%s/%s/all_exposures.txt'%(datadir, targetname), format='ascii')
    if ('G130M' not in gts['Grating']) and ('G160M' not in gts['Grating']):
        return 'No G130M or G160M observed'

    ######## If everthing is fine, let's start! #########
    cosspecname = '%s/%s/%s_coadd_%s_final_all.fits.gz'%(datadir, targetname, targetname, gtag)

    ## ask if use existing knots 
    ## (if the QSO has been processed before, the knots should be saved in knots/) 
    knotname = '%s/Dropbox/HSLA_Feb16/code/knots/%s_coadd_%s_knots.jsn'%(homedir, targetname, gtag)
    if os.path.isfile(knotname) is True:
        useknot = question_yes_no("Yong: found %s. Use it? "%(knotname))
        if useknot == True: knotlist = read_knots(knotname)
        else: knotlist = None
    else: knotlist = None  # the default

    ## Fit continuum
    spec = lsx.XSpectrum1D.from_file(cosspecname)
    try:
        spec.fit_continuum(knots=knotlist)
    except RuntimeError:
        print('RuntimeError: Problem generating continuum spline knots. Use random knots.')
        tempknots = '%s/Dropbox/HSLA_Feb16/code/knots/tempknots.jsn'%(homedir)
        knotlist = read_knots(tempknots)
        spec.fit_continuum(knots=knotlist)

    ## save the fitting knots
    copyfile('_knots.jsn', knotname)
    os.remove('_knots.jsn')

    ## save the whole continuum 
    filedir = savedir+'/'+targetname
    if os.path.isdir(filedir) == False:
        os.makedirs(filedir)
        os.makedirs(filedir+'/lines')
        os.makedirs(filedir+'/figs')
    print("Saving the total spectra in %s ... "%(filedir))
    save_SpecAll(spec.wavelength.value, spec.flux.value,
                 spec.sig.value, continuum=spec.co.value,
                 targetname=targetname, filedir=filedir, save_fits=True)

    ## slice the spec for each lines
    slice_SpecLine(spec.wavelength.value, spec.flux.value,
                   spec.sig.value, continuum=spec.co.value,
                   targetname=targetname, lines=lines,
                   coords=objcoords, velwidth=2000, 
                   filedir=filedir, pltrange=[-400, 400])

    ## check if HI data exist; if not, extract it from HI4PI
    observation = 'HI4PI'
    infoHI = find_HIspec(targetname, observation=observation) 

    ## stacking the spectra
    print("Stacking Spectra ... ")
    for nbin in [1, 3]:
        stack_spec(filedir=filedir, targetname=targetname, pltrange=[-400, 400],  
                   lines=lines, coords=objcoords, infoHI=infoHI, nbin=nbin)

    ## record if this one is good for analyses, i.e., has been double checked 
    if ask_if_final == True:
        if_final_yz = question_yes_no("Yong: write %s to QSOs_good2go.txt?"%(targetname))
        if if_final_yz == True:
            good2go = open(homedir+'/Dropbox/HSLA_Feb16/code/tables/QSOs_good2go.txt', 'a')
            good2go.write('%25s \n'%(targetname))
            good2go.close()
    return targetname 
