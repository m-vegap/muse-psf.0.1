# muse-psf.0.1

This repository stores scripts for PSF fitting from IFS data using PampelMuse. (https://gitlab.gwdg.de/skamann/pampelmuse).

The main script is utils.py, which contains all the necessary steps to prepare the input catalogue necessary to run PampelMuse. This steps are described below: 

1. Scan directory: Scan the directory that contains your single datacubes and create a monitoring file to keep the track of the progress.
  
2. Make filter images from single datacubes: Calls esorex to create I-band images. Other filters are available in the "filters" folder.

3. Source extraction from filter images: Creates input catalogues for PampelMuse from the filter images.

4. Input catalogue for PampelMuse: Clip catalogue extracted in the previous step. 

5. Run PampelMuse: Calls the PampelMuse routines INITFIT/SINGLESRC, CUBEFIT, and POLYFIT.


DEPENDENCIES
------------

This script makes use of the following important libraries and softwares:

	+ Numpy.
	+ Astropy.
	+ Photutils.
	+ Esorex.
	+ SExtractor.
	+ Pampelmuse.

USAGE
-----

Before running the code, you must change the hardcoded.py file, which contains the paths of the folders where you want to store the filter images, catalogues and psf files. If this folders is not in the library folder, the script creates them in the correspondent process. Also contains the instrument gain and the passband that you want to use to create the filter images.

After this is done, you must change the configuration files stored in the "config" folder. The default files and se.params are not necessary to change unless you want to use other configuration and obtain more or less parameters in the first catalogue.  
The files pampelmuse.json and se.config contain the information for the PSF fit, so you must change it with the parameters of the instrument and initial parameters of the PSF. 

After all this is done, you simply run:

		python psf_fit.py

And the code will do all the steps sequentially one by one. When the script finish one of the process, it will ask you if you want to continue with the next one. If the answer is "no", then the process is stopped. 

OUTPUTS
-------

The code will generate the following outputs: 

1. Monitoring file: monitoring_file.txt stored in the main folder. 

2. Filter images: By default, the output images are stored in the "single_iband" folder.

3. SExtractor catalogue: Table in ASCII.ecsv format with the name “filename_filter.se.cat”.
The output table contains a lot of parameters like flux, magnitude, x and y positions, Alpha_J2000 (RA), Delta_J2000 (DEC), ellipticity, class_star (stellarity), etc. You can add or remove parameters changing the config/se.params file.

4. Input catalog for PampelMuse: These catalogues have the suffix ‘ppmuse_in.dat'.

5. psf.fits and prm.fits files: The parameter file (prm.fits) contains the parameters of the PSF and the psf file (psf.fits) contains the image of the PSF and the residuals.
