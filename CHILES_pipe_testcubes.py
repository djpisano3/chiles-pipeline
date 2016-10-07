# CHILES_pipe_testcubes.py
# This module is run after the target module in the CHILES pipeline and makes some
# small test cubes to verify that the pipeline calibration and flagging is of
# sufficient quality.  
# 3/3/16 DJP
# 6/28/16, DJP:  Reducing number of iterations to 100 from 1000.  
#                Testing splitting data first. 
# 9/21/16 DJP:  Reduced number of iterations to 1.  As with original code, but running on each spw separately.
 

logprint ("Starting CHILES_pipe_testcubes.py", logfileout='logs/testcubes.log')
time_list=runtiming('testcubes', 'start')

# If running this module right after starting casa, must run CHILES_pipe_restore first.
# Step 1, load needed functions.
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

# Test if imaging is quicker after  averaging by 2x in time & 4x in frequency
outputms=ms_active[:-3]+'_calibrated_deepfield'


seq=range(0,15)
for ii in seq:

    logprint ('Split calibrated deepfield uv data, Spw='+str(ii), logfileout='logs/testcubes.log')

    targetfile=outputms+'_spw'+str(ii)+'.ms'
    
    default('split')
    vis=ms_active
    datacolumn='corrected'
    outputvis=targetfile
    field='deepfield'
    spw =str(ii)
    width=4
    timebin='16s'
    split()

# Make test cube of 10 mJy source, all channels
    logprint ('Make cube of 10 mJy source, all channels, Spw='+str(ii), logfileout='logs/testcubes.log')

    default('clean')

    image_name='targetcube_10mJy_spw'+str(ii)
    fieldid='deepfield'
    phasecenter='J2000 10h01m31.4 +02d26m40'
    grid_mode=''
    number_w=1
    image_size=[64,64]
    iteration=1  # Changed to 1.  
    mask_name=['']

    vis=targetfile
    imagename=image_name
    selectdata=False
    field=fieldid
    spw=''
    mode='frequency'
    nterms=1
    niter=iteration
    gain=0.1
    gridmode=grid_mode
    wprojplanes=number_w
    threshold='0.0mJy'
    psfmode='clark'
    imagermode='csclean'
    cyclefactor=1.5
    cyclespeedup=-1
    multiscale=[]
    interactive=False
    mask=mask_name
    imsize=image_size
    cell=['1.5arcsec','1.5arcsec']
    stokes='I'
    weighting='briggs'
    robust=0.8
    uvtaper=False
    modelimage=''
    restoringbeam=['']
    pbcor=False
    usescratch=False
    allowchunk=False
    async=False
    clean()

# Make small cube of lowest z detection
    logprint ('Make cube of strongest pilot detection, Spw='+str(ii), logfileout='logs/testcubes.log')

    default('clean')

    image_name='targetcube_HIdet_spw'+str(ii)
    fieldid='deepfield'
    phasecenter='J2000 10h01m15.2 +02d18m24'
    grid_mode=''
    number_w=1
    image_size=[64,64]
    iteration=1  # Changed to 1 to save time
    mask_name=['']

    vis=targetfile
    imagename=image_name
    selectdata=False
    field=fieldid
    spw=''
    mode='frequency'
    nterms=1
    niter=iteration
    gain=0.1
    gridmode=grid_mode
    wprojplanes=number_w
    threshold='0.0mJy'
    psfmode='clark'
    imagermode='csclean'
    cyclefactor=1.5
    cyclespeedup=-1
    multiscale=[]
    interactive=False
    mask=mask_name
    imsize=image_size
    cell=['1.5arcsec','1.5arcsec']
    stokes='I'
    weighting='briggs'
    robust=0.8
    uvtaper=False
    modelimage=''
    restoringbeam=['']
    pbcor=False
    usescratch=False
    allowchunk=False
    async=False
    clean()
    
#Extract spectrum for both cubes
    logprint ('Extract Spectra, Spw='+str(ii), logfileout='logs/testcubes.log')

    default('imval')
    stokes='I'
    imagename='targetcube_10mJy_spw'+str(ii)+'.image'
    results=imval()

    freq1=results['coords'][:,3]
    freq1/=1.e6
    chan1=range(len(freq1))
    flux1=results['data']

    fig=pylab.figure()
    pylab.plot(freq1,flux1,'k-')
    pylab.xlabel('Frequency [MHz]')
    pylab.ylabel('Flux [Jy/bm]')
    pylab.ylim(-0.002,0.025)
    pylab.title('10 mJy source, Spw '+str(ii))
    pylab.savefig('target_10mJy_spectrum_spw'+str(ii)+'.png')
    pylab.close(fig)

    default('imval')
    stokes='I'
    box='16,16,48,48'
    imagename='targetcube_HIdet_spw'+str(ii)+'.image'
    results=imval()

    freq2=results['coords'][3,3,:,3]
    freq2/=1.e6

    chan2=range(len(freq2))
    flux2=np.mean(np.mean(results['data'],axis=0),axis=0)
    
    fig=pylab.figure()
    pylab.plot(freq2,flux2,'k-')
    pylab.xlabel('Frequency [MHz]')
    pylab.ylabel('Flux [Jy/bm]')
    pylab.ylim(-0.002,0.002)
    pylab.title('HI detection, Spw '+str(ii))
    pylab.savefig('target_HI_spectrum_spw'+str(ii)+'.png')
    pylab.close(fig)

#Move plots, images to sub-directory
os.system("mv *.png plots")
os.system("mv targetcube*.* images")

#Make output webpage
if os.path.exists("testcubes.html"):
    os.system("rm testcubes.html")
wlog = open("testcubes.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_testcubes results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/testcubes.log">Test Cubes Log</a></li>\n')
wlog.write('<hr>\n')
wlog.write('<li> Spectrum towards 10mJy source & HI detection at z=0.02 source: \n')
wlog.write('<table>\n')
for ii in seq:
    wlog.write('<tr>\n')
    wlog.write('<td><img src="plots/target_10mJy_spectrum_spw'+str(ii)+'.png"></td>')
    wlog.write('<td><img src="plots/target_HI_spectrum_spw'+str(ii)+'.png"></td>\n')
    wlog.write('</tr>\n')
    wlog.write('<hr>\n')
wlog.write('</table>\n')
wlog.write('</li>')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()


logprint ("Finished CHILES_pipe_testcubes.py", logfileout='logs/testcubes.log')
time_list=runtiming('testcubes', 'end')

pipeline_save()
