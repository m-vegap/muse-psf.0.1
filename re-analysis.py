# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 15:40:28 2022

@author: macarena
"""
import sys
import os
import subprocess
import json
import paths as path
from astropy.io import ascii

###########################################################################################
############################### MAIN FUNCTION #############################################
###########################################################################################
def main():

    ####################### 0.- Initial Configuration ###############################

    #### A.- Check the parameters ####
    if (len(sys.argv)<4):
        print ('Usage: '+sys.argv[0]+' <Datacubes Directory> <Output Directory> <n° of files>')
        sys.exit(0)

    #### B.- Housekeeping table and directories ####

    cubedir = sys.argv[1]
    outdir = sys.argv[2]
    nfiles = sys.argv[3]
    inputlist = ascii.read('monitoring_file.txt')
            
    
    if(os.path.exists(outdir)==False):
        print ('Creating the analysis directory')
        subprocess.call("mkdir "+ outdir, shell=True)
    
    # open sample ppmuse config file
    config = json.load(open(path.pampelmuse_file, 'r'))
    
    
    
    # %% loop over everything
    ncount = 0
    for i in range(len(inputlist)):
        mode = inputlist['modes'][i]
        if ncount == nfiles:
            print('No more files to process')
            break
        else:
            if inputlist['ppmuse'][i] > 0:
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
                        subprocess.run('mv {0}.prm.fits {1}'.format(os.path.join(cubedir, fname), outdir),
                                    shell=True, check=True)
                        subprocess.run('mv {0}.psf.fits {1}'.format(os.path.join(cubedir, fname), outdir),
                                    shell=True, check=True)
                        # write to the housekeeping table
                        inputlist['ppmuse'][i] = 1
    
                    except subprocess.CalledProcessError:
                        inputlist['ppmuse'][i] = -1
    
                ncount += 1
    
    # save the updated monitoring table
    inputlist.write('monitoring_file.txt', format='ascii.tab', overwrite=True)
    print('Updated monitoring file.')
    
#####################################
###### CALL THE MAIN FUNCTION #######
if __name__ == "__main__":
    main()
#####################################
