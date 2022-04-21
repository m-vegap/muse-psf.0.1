# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 11:26:08 2022

@author: macarena
"""
import sys
import os
import subprocess
from astropy.table import Table, hstack
from astropy.io import fits, ascii 
from astropy.stats import sigma_clipped_stats
import numpy as np
from photutils import DAOStarFinder
import json
import paths as path

#########################################
############ Initial Message   ##########
#########################################

class Initial_Message():
    
    def __init__(self):
        self.cubedir = sys.argv[1]
        self.nfiles = sys.argv[2]
        print("")
    def show_message(self):
        print('\n\t \t --------------------------------------------\n')
        print('\t \t               MUSE PSF Fitting              ') 
        print('\t \t                                     ')
        print('\t \t                Author: M.Vega')
        print('\n\t \t --------------------------------------------\n')
        print ("")
        print ('Datacubes directory: '+os.path.realpath(self.cubedir))
        print ('Number of files to process: '+self.nfiles)
        print ("")
        
class psf_fitter():
    
    def __init__(self):
        pass
        
    def scan(self, directory, check_mode = True, save_file = True):
        
        cubedir = sys.argv[1]
        
        fitsfiles = [f for f in os.listdir(cubedir) if f.endswith('.fits')] 
        if len(fitsfiles) != 0:
            pass
        else:
            while len(fitsfiles) == 0:
                print('No fits files found.')
                yes = ['y', 'Y', 'yes', 'Yes', 'YES']
                no = ['n', 'N', 'no', 'No', 'NO']
                print('Do you want to scan another directory?: [y/n]')
                answer = str(input())
                continue_yes = answer in yes
                if continue_yes:
                    print('Choose your data directory:')
                    directory = str(input())
                    fitsfiles = [f for f in os.listdir(directory) if f.endswith('.fits')] 
                elif answer in no:
                    print('Process stopped.')
                    exit()
            
            
        # set up the housekeeping table
        iband = np.zeros_like(fitsfiles, dtype=int)
        sources = np.zeros_like(fitsfiles, dtype=int)
        prep = np.zeros_like(fitsfiles, dtype=int)
        ppmuse = np.zeros_like(fitsfiles, dtype=int)
        nsources = np.zeros_like(fitsfiles, dtype=int)
        
        tab = Table()
        tab['file'] = fitsfiles
        tab['iband'] = iband
        tab['sources'] = sources
        tab['prep'] = prep
        tab['ppmuse'] = ppmuse
        tab['nsources'] = nsources
        
    
        self.monitoring_file = tab


        if check_mode == True:
            
            modes = []
            
            for file in self.monitoring_file['file']:
                with fits.open(os.path.join(cubedir,file)) as hdul:
            
                ### agregar un bloque try
                
                    CTYPE1 = hdul[1].header['CTYPE1'] 
                    CTYPE2 = hdul[1].header['CTYPE2'] 
                    CONFIG = hdul[0].header['HIERARCH ESO INS AO FOCU1 CONFIG']
                    
                    if (type(CTYPE1) != str or type(CTYPE2) != str or type(CONFIG) != str):  ###check rápido de los header keys
                        
                        print('Header keywords must be strings')
                    else:
                    
                    
                        '''mode 0 : WFM Science Data
                           mode 1 : NFM Standar Stars
                           mode 2 : NFM Science Data'''
                        
                        
                        if (CTYPE1 == 'RA---TAN' and CTYPE2 == 'DEC--TAN' and CONFIG == 'WFM'):
                            mode = 'WFM-WCS'
                            
                        elif (CTYPE1 == 'PIXEL' and CTYPE2 == 'PIXEL' and CONFIG == 'NFM'):
                            mode = 'NFM-PIXEL'
                        
                        elif (CTYPE1 == 'RA---TAN' and CTYPE2 == 'DEC--TAN' and CONFIG == 'NFM'):
                            mode = 'NFM-WCS'
                        
                        modes.append(mode)
                        
            self.monitoring_file['modes'] = modes       
            
        if save_file == True: 
            
            tab.write('monitoring_file.txt', format='ascii.tab', overwrite=True)
            print('Monitoring file created.')
        
            
            
            
    def mk_filter_images(self, band = path.passband):
        
        cubedir = sys.argv[1]
        nfiles = sys.argv[2]
        
        if(os.path.exists(path.filterdir)==False):
            print ('Creating the filter images directory')
            subprocess.call("mkdir "+ path.filterdir, shell=True)
        inputlist = self.monitoring_file
        print('Making filter images...')
        ncount = 0
        for i in range(len(inputlist)):
            if ncount == nfiles:
                print('No more files to process')
                break
            else:
                file = inputlist['file'][i]
                if not inputlist['iband'][i]:
                    fname, fformat = os.path.splitext(file)
                    nname = fname + '_' + band + fformat
                    print(nname)
                    args_esorex = {'image': os.path.join(os.path.realpath(cubedir), file),
                                   'passband': band,
                                   'filterlist': path.filterlist}
                    args_mv = {'input': os.path.join(os.path.realpath(cubedir), nname),
                               'output': os.path.join(path.filterdir, nname)}
                    os.system('muse_cube_filter -f {passband} {image} {filterlist}'.format(**args_esorex))
                    os.system('mv {input} {output}'.format(**args_mv))
                    # write to the housekeeping table
                    inputlist['iband'][i] = 1
                    ncount += 1
                    print('Progress: {0}/{1}'.format(ncount, nfiles))
        # save the updated housekeeping table
        inputlist.write('monitoring_file.txt', format='ascii.tab', overwrite=True)
        print('Filter images created.')
        print('Updated monitoring file.')
                
                
                
                
    def match(self, x1, y1, x2, y2, tol=1.5):
        """
        Coordinate matching of (x1, y1) and (x2, y2) with tolerance tol.
        """
        cat1 = []
        cat2 = []
        for i in range(len(x1)):
            for j in range(len(x2)):
                dx_tmp = x1[i] - x2[j]
                dy_tmp = y1[i] - y2[j]
                if dx_tmp**2 + dy_tmp**2 <= tol**2:
                    cat1.append(i)
                    cat2.append(j)
        return np.array(cat1), np.array(cat2)
                
                
    def SExtraction(self, band = path.passband):
        
        nfiles = sys.argv[2]
        
        if(os.path.exists(path.catdir)==False):
            print ('Creating the catalogue directory')
            subprocess.call("mkdir "+ path.catdir, shell=True)
        # load housekeeping table
        inputlist = self.monitoring_file
        # zeropoint list
        zpt_list = ascii.read(path.zero_point_list)
        # get the correct zeropoint
        ind = zpt_list['Filter'] == band.upper()
        ZPT = zpt_list['ZPT'][ind]
        
        
        
        # %% loop over everything
        ncount = 0
        
        # %% loop
        for i in range(len(inputlist)):
            if ncount == nfiles:
                print('No more files to process')
                break
            else:
                if not inputlist['sources'][i]:
                    mode = inputlist['modes'][i]
                    basename = inputlist['file'][i]
                    fname, ext = os.path.splitext(basename)
                    file = fname + '_' + band + ext
                    try:
                        HDUL = fits.open(os.path.join(path.filterdir, file))
                        data = HDUL['DATA'].data
                        header_amb = HDUL['PRIMARY'].header # ambiental data stored here
                        exptime = header_amb['EXPTIME']
                        
                        
                        # first guess on FWHM from seeing and airmass
                        avg_seeing = np.mean([header_amb['HIERARCH ESO TEL AMBI FWHM END'],
                                              header_amb['HIERARCH ESO TEL AMBI FWHM START']])
                        airm_start = header_amb['ESO TEL AIRM START']
                        airm_end = header_amb['ESO TEL AIRM END']
                        z_start = np.arccos(1./airm_start)
                        z_end = np.arccos(1./airm_end)
                        avg_z = np.mean([z_start, z_end])
                        avg_airm = 1./np.cos(avg_z)
                        if mode == 'WFM-WCS':
                            FWHM_guess = avg_seeing * avg_airm**(3./5.)  # in arcsec
                            
                        elif (mode == 'NFM-PIXEL' or mode == 'NFM-WCS'):
                            FWHM_guess = 0.1  # in arcsec
                            
        
                        # close the file
                        HDUL.close()
        
                        # how do we want to call the output?
                        catname = os.path.join(path.catdir,
                                               fname + '_' + band + '.se.cat')
        
                        # arguments to call SExtractor
                        args_se = {'IMAGE': os.path.join(path.filterdir, file),
                                   'NAME': catname,
                                   'ZPT': ZPT[0] + 2.5*np.log10(exptime),
                                   'SEEING_FWHM': FWHM_guess,
                                   'GAIN': path.muse_gain}
        
                        # call SExtractor
                        print('Running SExtractor on File:', file)
                        print('Parameters')
                        print('ZPT: ', ZPT[0] + 2.5*np.log10(exptime))
                        print('EXPTIME: ', exptime)
                        print('FWHM: ', FWHM_guess)
                        os.system('sex {IMAGE} -c config/se.config -CATALOG_NAME {NAME} -MAG_ZEROPOINT {ZPT} -SEEING_FWHM {SEEING_FWHM} -GAIN {GAIN}'.format(**args_se))
        
                            
                        if (mode == 'WFM-WCS' or mode == 'NFM-WCS'):
                            # read the input catalogue
                            cat_se = ascii.read(catname)
            
                            # calculate image statistics and run daofind
                            mean, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
                            daofind = DAOStarFinder(fwhm=FWHM_guess, threshold=5*std)
                            cat_dao = daofind(data - median)
                            # adjust zpt
                            cat_dao['mag'] += ZPT[0]
                        
        
                            # merge the two tables
                            id_se, id_dao = match(x1=cat_se['X_IMAGE_DBL'],
                                                  y1=cat_se['Y_IMAGE_DBL'],
                                                  x2=cat_dao['xcentroid'],
                                                  y2=cat_dao['ycentroid'])
            
                            # output table
                            tab = hstack([cat_dao[id_dao], cat_se[id_se]])
                            tab.write(os.path.join(path.catdir,
                                                   file.replace('.fits', '.detections.cat')),
                                      format='ascii.ecsv',
                                      overwrite=True)
            
                            # delete sextractor table
                            os.system('rm {se_cat}'.format(se_cat=catname))
                        
                            
        
                        # write to the housekeeping table
                        inputlist['sources'][i] = 1
        
                    except FileNotFoundError:
                        inputlist['iband'][i] = -1
        
                    ncount += 1
                    print('Progress: {0}/{1}'.format(ncount, nfiles))
        
        # save the updated monitoring table
        inputlist.write('monitoring_file.txt', format='ascii.tab', overwrite=True)
        print('Updated monitoring file.')            
               
    def incat_ppmuse(self, band = path.passband):
        
        nfiles = sys.argv[2]
        
        # %% read from housekeeping table
        inputlist = self.monitoring_file
        datadir = path.catdir
        
        # %% cut levels
        #mag_crowd = 25
        mag_lim = 24
        stellarity_lim = 0.5
        flux_lim = 1000
        #passband = "Cousins_I"
        
        # %% loop over everything
        ncount = 0
        
        # %% loop over all files
        for i in range(len(inputlist)):
            mode = inputlist['modes'][i]
            if (mode == 'WFM-WCS' or mode == 'NFM-WCS'):
                if ncount == nfiles:
                    print('No more files to process')
                    break
                else:
                    if not inputlist['prep'][i] and inputlist['iband'][i] > 0:
                        basename = inputlist['file'][i]
                        fname, ext = os.path.splitext(basename)
                        name = fname + '_' + band + '.detections.cat'
                        print('Reading file: ', name)
                        cat = ascii.read(os.path.join(datadir, name))
                        cut_mag = cat['mag'] < mag_lim
                        cut_stellarity = cat['CLASS_STAR'] > stellarity_lim
                        cut = np.logical_and(cut_mag, cut_stellarity)
                        print('Selecting {n} objects for analysis with PampelMuse'.format(n=sum(cut)))
                        tab = Table()
                        tab['id'] = cat['id'][cut]
                        tab['ra'] = cat['ALPHA_J2000'][cut]
                        tab['dec'] = cat['DELTA_J2000'][cut]
                        tab[band.upper()] = cat['mag'][cut]
                        tab.write(os.path.join(datadir, name.replace('.detections.cat', '.ppmuse_in.dat')),
                                  format='ascii.csv', overwrite=True)
                        print('Wrote to: ', name.replace('.detections.cat', '.ppmuse_in.dat'))
                        inputlist['nsources'][i] = len(tab)
            
                        # write to the housekeeping table
                        inputlist['prep'][i] = 1
                        ncount += 1
                        print('Progress: {0}/{1}'.format(ncount, nfiles))
            
            elif mode == 'NFM-PIXEL':
                if ncount == nfiles:
                    print('No more files to process')
                    break
                else:
                    if not inputlist['prep'][i] and inputlist['iband'][i] > 0:
                        
                        basename = inputlist['file'][i]
                        fname, ext = os.path.splitext(basename)
                        name = fname + '_' + band + '.se.cat'
                        print('Reading file: ', name)
                        cat = ascii.read(os.path.join(datadir, name))
                        cut_mag = cat['FLUX_AUTO'] > flux_lim
                        cut_stellarity = cat['CLASS_STAR'] > stellarity_lim
                        cut = np.logical_and(cut_mag, cut_stellarity)
                        print('Selecting {n} objects for analysis with PampelMuse'.format(n=sum(cut)))
                        tab = Table()
                        tab['id'] = cat['NUMBER'][cut]
                        tab['id'] = tab['id'].astype(int)
                        tab['x'] = cat['X_IMAGE_DBL'][cut]
                        tab['y'] = cat['Y_IMAGE_DBL'][cut]
                        tab[band.upper()] = cat['MAG_AUTO'][cut]
                        tab['status'] = '3'
                        tab.write(os.path.join(datadir, name.replace('.se.cat', '.ppmuse_in.dat')),
                                  format='ascii.csv', overwrite=True)
            
                        print('Wrote to: ', name.replace('.se.cat', '.ppmuse_in.dat'))
            
                        # write to the monitoring table
                        inputlist['nsources'][i] = len(tab)
                        inputlist['prep'][i] = 1
            
                        ncount += 1
                        print('Progress: {0}/{1}'.format(ncount, nfiles))
                        
                
                
    def run_ppmuse(self):
        
        cubedir = sys.argv[1]
        nfiles = sys.argv[2]
        
        if(os.path.exists(path.analysisdir)==False):
            print ('Creating the analysis directory')
            subprocess.call("mkdir "+ path.analysisdir, shell=True)
        
        # open sample ppmuse config file
        config = json.load(open(path.pampelmuse_file, 'r'))
        # housekeeping table and directories
        inputlist = self.monitoring_file
        

        # %% loop over everything
        ncount = 0
        for i in range(len(inputlist)):
            mode = inputlist['modes'][i]
            if ncount == nfiles:
                print('No more files to process')
                break
            else:
                if not inputlist['ppmuse'][i]:
                    if inputlist['nsources'][i] > 0 and inputlist['iband'][i] > 0:
                        basename = inputlist['file'][i]
                        fname, ext = os.path.splitext(basename)
                        catname = fname + '_' + path.passband + '.ppmuse_in.dat'
                        out_config = open('tmp.pampelmuse.json', 'w')
                        config['global']['prefix'] = os.path.join(cubedir, fname)
                        config['catalog']['name'] = os.path.join(path.catdir, catname)
                        s = json.dumps(config)
                        out_config.write(s)
                        out_config.close()
        
                        try:
                            
                            if (mode == 'WFM-WCS' or mode == 'NFM-WCS'):
                                subprocess.run(['PampelMuse', 'tmp.pampelmuse.json',
                                                'INITFIT'],
                                               check=True)
                            elif mode == 'NFM-PIXEL':
                                subprocess.run(['PampelMuse', 'tmp.pampelmuse.json',
                                                'SINGLESRC'],
                                               check=True)
                                
                            subprocess.run(['PampelMuse', 'tmp.pampelmuse.json',
                                            'CUBEFIT'],
                                           check=True)
                            subprocess.run(['PampelMuse', 'tmp.pampelmuse.json',
                                            'POLYFIT'],
                                           check=True)
                            subprocess.run('mv {0}.prm.fits {1}'.format(os.path.join(cubedir, fname), path.analysisdir),
                                        shell=True, check=True)
                            subprocess.run('mv {0}.psf.fits {1}'.format(os.path.join(cubedir, fname), path.analysisdir),
                                        shell=True, check=True)
                            # write to the housekeeping table
                            inputlist['ppmuse'][i] = 1
        
                        except subprocess.CalledProcessError:
                            inputlist['ppmuse'][i] = -1
        
                    ncount += 1
        
        # save the updated monitoring table
        inputlist.write('monitoring_file.txt', format='ascii.tab', overwrite=True)
        print('Updated monitoring file.')


###########################################################################################
############################### MAIN FUNCTION #############################################
###########################################################################################
def main():

    ####################### 0.- Initial Configuration ###############################

    #### A.- Check the parameters ####
    if (len(sys.argv)<3):
        print ('Usage: '+sys.argv[0]+' <Datacubes Directory> <n° of files>')
        sys.exit(0)

    #### B.- Initializating ####
    process = psf_fitter()
    #Current Directory
    cubedir = sys.argv[1]
    #Target Light Curve
    nfiles = sys.argv[2]


    #### C.- Display initial message ####
    init=Initial_Message()
    init.show_message()
    

    
    
    ###################### 1.- Scan the input directory ###############################
    print('Initializing...')
    print('Preparing scan...')

    process.scan(cubedir)
        
    print('Scan finished without errors.')
    
    
    
    
    ####################### 2.- Make the filter images ################################
    print("Do You Want To Continue? [y/n]")
    answer = str(input())
    yes = ['y', 'Y', 'yes', 'Yes', 'YES']
    no = ['n', 'N', 'no', 'No', 'NO']
    continue_yes = answer in yes
    if  continue_yes:
          #nfiles = int(input('How many files to process: '))
          process.mk_filter_images()
          print('Filter images finished without errors.')
    elif answer in no:
        continue_yes = False
        print('Process stopped.')
        exit()
    
    
    
    
    ####################### 3.- Source Extraction ################################
    print("Do You Want To Continue? [y/n]")
    answer = str(input())
    yes = ['y', 'Y', 'yes', 'Yes', 'YES']
    no = ['n', 'N', 'no', 'No', 'NO']
    continue_yes = answer in yes
    if  continue_yes:
          process.SExtraction()
          print('Source extraction finished without errors.')
    elif answer in no:
        continue_yes = False
        print('Process stopped.')
        exit()
        
    
    
    
    ####################### 4.- Input catalogue for PampelMuse ################################
    print("Do You Want To Continue? [y/n]")
    answer = str(input())
    yes = ['y', 'Y', 'yes', 'Yes', 'YES']
    no = ['n', 'N', 'no', 'No', 'NO']
    continue_yes = answer in yes
    if  continue_yes:
          process.incat_ppmuse()
          print('PampelMuse input catalogue finished without errors.')
    elif answer in no:
        continue_yes = False
        print('Process stopped.')
        exit()
    
    
    
    
    ####################### 5.- Run PampelMuse ################################
    print("Do You Want To Continue? [y/n]")
    answer = str(input())
    yes = ['y', 'Y', 'yes', 'Yes', 'YES']
    no = ['n', 'N', 'no', 'No', 'NO']
    continue_yes = answer in yes
    if  continue_yes:
          process.run_ppmuse()
          print('PampelMuse routines finished without errors.')
    elif answer in no:
        continue_yes = False
        print('Process stopped.')
        exit()
    

#####################################
###### CALL THE MAIN FUNCTION #######
if __name__ == "__main__":
    main()
#####################################
