import astropy.io.fits as fits
import sys, re
import numpy as np

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def import_lines():
    '''
    Load the interested lines.
    '''

    defaultlines = ['PII 1152', 'PII 1301',  
                    'FeII 1142', 'FeII 1143', 'FeII 1144', 'FeII 1608',
                    'CII 1334', 'CII* 1335', 'CIV 1548', 'CIV 1550', 
                    'SiII 1304', 'SiII 1190', 'SiII 1193', 'SiII 1260', 'SiII 1526', 
                    'SiII* 1533', 'SiII* 1264', 'SiIII 1206', 'SiIV 1393', 'SiIV 1402', 
                    'SII 1250', 'SII 1253', 'SII 1259', 'OI 1302', 'NV 1238', 'NV 1242']
                    #'SI 1296', 'SI 1316', 'SI 1425', 'SI 1295',
                    #'HI 1215', 'PII 1532','SiI 1595', 'SiI 1589',
		    #'SiII* 1309', 'SiII* 1197'

    return defaultlines

def read_linelibrary(lines='All', doprint=True):
    '''
    For the input lines, find the wavelength and fval. 
    '''

    # find the path, read in the line library 
    import sys
    for ipath in sys.path:
        if 'GitRepo' in ipath:
            linepath = ipath
            break
    all_lin = fits.getdata(linepath+'/yzSpec/files/xidl_all_lin.fits')
    # Lines that covered in G130 M, should really expand this later on. 
    if lines == 'All':
        dolines = import_lines()
    else:
        # Do spec slicing for some particular lines
        # lines could be just one: lines = ['FeII 1142']
        #              or several: lines = ['SiII 1250', 'SiII 1253', 'SiII 1259']

        # reformat the lines, in case it is not in forms, like, 'FeII 1142', but 'Fe II 1142'
        dolines = []
        defaultlines = import_lines()
        if type(lines) is str: lines = [lines]  # make sure the for loop will work
        for il in range(len(lines)):
            eles = re.split('(\d+)', lines[il].replace(' ', ''))
            newline = '%s %s'%(eles[0], eles[1])
            if newline in defaultlines: dolines.append(newline)
            else: logger.info(lines[il]+' not in the right format. skip.')

    # read line wavelength and oscillator strength
    line_lambda = []
    line_fval = []
    liblines = []
    for i in range(len(dolines)):
        for j in range(len(all_lin)):
            if dolines[i] == all_lin[j][0]:
                line_lambda.append(all_lin[j][1])
                liblines.append(dolines[i])

                # # FeII 1144. Check emails exchange with X.  
                if liblines[-1] == 'FeII 1144': line_fval.append(0.083)
                else: line_fval.append(all_lin[j][2])
                
                break

    if doprint == True: logger.info("Found these in our library: ",liblines)
    line_ref = 'Morton (2003)'
    return liblines, line_lambda, line_fval, line_ref

