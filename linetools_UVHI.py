import numpy as np
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('TkAgg')  # backend for graphical interfaces
from yzSpec.yes_or_no import yes_or_no
from yzSpec.read_knots import read_knots

import os
homedir = os.path.expanduser('~')

import logging 
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def linetools_UVHI(target, lines='All', ask_if_final=True):
    '''
    Call the linetools continuum fitting function to fit the spectra 
    for the input target. It will also find the corresponding HI 21cm spectra.
   
    ask_if_final: this is a yes or no question; if yes or True, then 
                  the target will be recorded in good2go. 
    '''
   
    from yzSpec.find_targetinfo import find_targetinfo
    category, target_info = find_targetinfo(target)
    if category == 'none':
        logger.info(target+' is not in QSOALS, GALAXY, STAR_EARLY.')
        return 'none'
    elif len(target_info) == 0: 
        logger.info(target+' does not have either G130M or G160M observation.')
        return 'none'
    else: 
        # OK, now fit continuum 
        ## (if the QSO has been processed before, the knots should be saved in knots/) 
        knotdir = '%s/Dropbox/HSLA_Feb16/code/knots'%(homedir)
        knotname = '%s/%s_coadd_%s_knots.jsn'%(knotdir, target_info['NAME'], target_info['Grating'])
        knotlist = read_knots(knotname)
 
        print(target_info['DATAFILE']) 
        logger.info('Start fitting %s ...'%(target)) 
        import linetools.spectra.xspectrum1d as lsx 
        lt_spec = lsx.XSpectrum1D.from_file(target_info['DATAFILE'])
        try:
            lt_spec.fit_continuum(knots=knotlist)
        except RuntimeError:
            logger.info('RuntimeError: Problem generating continuum spline knots. Use random knots.')
            tempknots = '%s/Dropbox/HSLA_Feb16/code/knots/tempknots.jsn'%(homedir)
            knotlist = read_knots(tempknots)
            lt_spec.fit_continuum(knots=knotlist)
   
        from shutil import copyfile 
        copyfile('_knots.jsn', knotname)
        os.remove('_knots.jsn')
    
        ## save the whole continuum 
        savedir = homedir+'/Dropbox/HSLA_Feb16/QSOSpec_YZ/'+category
        filedir = savedir+'/'+target_info['NAME']
        if os.path.isdir(filedir) == False:
            os.makedirs(filedir)
            os.makedirs(filedir+'/lines')
            os.makedirs(filedir+'/figs')

        ## add on 02/20/2018, to add HLSP-formatted product 
        hlsp_savedir = homedir+'/Dropbox/HSLA_Feb16/QSOSpec_YZ/hlsp_qsoals'
        hlsp_filedir = hlsp_savedir+'/'+target_info['NAME'].lower()
        if os.path.isdir(hlsp_filedir) == False:
            os.makedirs(hlsp_filedir)
            os.makedirs(hlsp_filedir+'/linedata_uv')
            os.makedirs(hlsp_filedir+'/linedata_21cm')
            os.makedirs(hlsp_filedir+'/figs')
  
        logger.info('Saving, ploting, stacking spec for %s ...'%(target)) 
        import time
        t1 = time.time()
        from yzSpec.save_spec import save_spec
        save_spec(lt_spec, target_info, has_continuum=True, filedir=filedir)
        save_spec(lt_spec, target_info, has_continuum=True, filedir=hlsp_filedir) # add on 02/20/2018 to add HLSP
        t2 = time.time()
        logger.info('Saving: %.1f seconds.'%(t2-t1))

        ## slice the spec, and save the data
        from yzSpec.plot_spec import plot_OneLine
        from yzSpec.read_linelibrary import import_lines
        if lines == 'All': lines = import_lines()

        for iline in lines:
            savename = save_spec(lt_spec, target_info, has_continuum=True, line=iline, filedir=filedir+'/lines')
            if savename != 'none':
                plot_OneLine(target_info, savename, filedir=filedir+'/figs')
        t3 = time.time()
        logger.info('Slicing & plotting: %.1f seconds.'%(t3-t2))
 
        ## check if HI data exist; if not, extract it from HI4PI
        #from yzSpec.extract_HI21cm import extract_HI21cm
        #hifile = extract_HI21cm(target_info, observation='HI4PI', filedir=filedir+'/lines') 
        #hifile = '/Users/Yong/Dropbox/HSLA_Feb16/HI4PI_HI21cm/four_arcmin/%s_HI21cm_HI4PI.fits.gz'%(target_info['NAME'])
        #tothisplace = '/Users/Yong/Dropbox/HSLA_Feb16/QSOSpec_YZ/QSOALS/%s/lines/%s_HI21cm_HI4PI.fits.gz'%(target_info['NAME'], target_info['NAME'])
        #copyfile(hifile, tothisplace)
        t4 = time.time()
        #logger.info('Finding HI 21cm: %.1f seconds.'%(t4-t3))

        ## stacking the spectra
        from yzSpec.plot_spec import stack_spec
        stack_spec(target_info, filedir+'/lines', savedir=filedir)
        t5 = time.time()
        logger.info('Stacking: %.1f seconds.'%(t5-t4))
       
        ## record if this one is good for analyses, i.e., has been double checked 
        if ask_if_final == True:
            if_final_yz = yes_or_no("Yong: write %s to QSOs_good2go.txt?"%(target_info['NAME']))
            if if_final_yz == True:
                good2go = open(homedir+'/Dropbox/HSLA_Feb16/code/tables/QSOs_good2go.txt', 'a')
                good2go.write('%25s \n'%(target_info['NAME']))
                good2go.close()
        return target_info['NAME']
