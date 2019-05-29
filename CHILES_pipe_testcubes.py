# CHILES_pipe_testcubes.py
# This module is run after the target module in the CHILES pipeline and makes some
# small test cubes to verify that the pipeline calibration and flagging is of
# sufficient quality.  
# 3/3/16 DJP
# 6/28/16, DJP:  Reducing number of iterations to 100 from 1000.  
#                Testing splitting data first. 
# 9/21/16 DJP:  Reduced number of iterations to 1.  As with original code, but running on each spw separately.
# 12/8/16 DJP: Added try, except to spw loop in case one is entirely flagged.
# 12/8/16 DJP: Calculated, print noise for each spectrum.  Changed output layout.
# 3/5/18 DJP: Changed clean to tclean, imaging in topocentric frequency to help with RFI identification.
# 4/22/18 DJP: Changing split to oldsplit
# 5/15/19 DJP: Removing flux and phase calibrator cubes, limiting low-z source to spw='13'

 

logprint ("Starting CHILES_pipe_testcubes.py", logfileout='logs/testcubes.log')
time_list=runtiming('testcubes', 'start')

# If running this module right after starting casa, must run CHILES_pipe_restore first.
# Step 1, load needed functions.
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

# Split data with no averaging (will slow down imaging, but avoids other issues
# fluxms=ms_active[:-3]+'_calibrated_fluxcal'
# phasems=ms_active[:-3]+'_calibrated_phasecal'
targetms=ms_active[:-3]+'_calibrated_deepfield'

seq=range(0,15)

for ii in seq:

    logprint ('Split calibrated uv data, Spw='+str(ii), logfileout='logs/testcubes.log')

    # fluxfile=fluxms+'_spw'+str(ii)+'.ms'
    # phasefile=phasems+'_spw'+str(ii)+'.ms'
    targetfile=targetms+'_spw'+str(ii)+'.ms'
    # if os.path.exists(fluxfile):
    #     os.system("rm -rf "+fluxfile)
    #     os.system("rm -rf fluxcube_spw"+str(ii)+".*")
    #     os.system("rm -rf images/fluxcube_spw"+str(ii)+".*")
    # if os.path.exists(phasefile):
    #     os.system("rm -rf "+phasefile)
    #     os.system("rm -rf phasecube_spw"+str(ii)+".*")
    #     os.system("rm -rf images/phasecube_spw"+str(ii)+".*")
    if os.path.exists(targetfile):
        os.system("rm -rf "+targetfile)
        os.system("rm -rf targetcube_10mJy_spw"+str(ii)+".*")
        os.system("rm -rf targetcube_HIdet_spw"+str(ii)+".*")
        os.system("rm -rf images/targetcube_10mJy_spw"+str(ii)+".*")
        os.system("rm -rf images/targetcube_HIdet_spw"+str(ii)+".*")

    try:

# # Split off flux calibrator
#         default('oldsplit')
#         vis=ms_active
#         datacolumn='corrected'
#         outputvis=fluxfile
#         field='1331+305=3C286'
#         spw =str(ii)
#         width=1
#         timebin=''
#         oldsplit()
 
# # Split off phase calibrator
#         default('oldsplit')
#         vis=ms_active
#         datacolumn='corrected'
#         outputvis=phasefile
#         field='J0943-0819'
#         spw =str(ii)
#         width=1
#         timebin=''
#         oldsplit()
               
# Split off deepfield                
        default('oldsplit')
        vis=ms_active
        datacolumn='corrected'
        outputvis=targetfile
        field='deepfield'
        spw =str(ii)
        width=1
        timebin=''
        oldsplit()
    except:
        logprint ('Unable to split uv data, Spw='+str(ii), logfileout='logs/testcubes.log')


# Now make cubes
for ii in seq:
    logprint ('Image calibrated uv data, Spw='+str(ii), logfileout='logs/testcubes.log')

    try:
# # Make test cube of flux calibrator, all channels
#         logprint ('Make cube of 3C286, all channels, Spw='+str(ii), logfileout='logs/testcubes.log')
# 
#         default('tclean')
#     
#         fluxfile=fluxms+'_spw'+str(ii)+'.ms'
#         image_name='fluxcube_spw'+str(ii)
#         fieldid='1331+305=3C286'
#         grid_mode=''
#         number_w=1
#         image_size=[64,64]
#         iteration=10000   # Total number of iterations for the entire cube (per spw), ~4 per channel 
#         mask_name=['']    
# 
#         vis=fluxfile
#         imagename=image_name
#         selectdata=False
#         field=fieldid
#         spw=''
#         specmode='cubedata'
#         nterms=1
#         niter=iteration
#         gain=0.1
#         gridmode=grid_mode
#         wprojplanes=number_w
#         threshold='10.0mJy'  # Adding threshold
#         deconvolver='clark'
#         imagermode='csclean'
#         cyclefactor=1.5
#         cyclespeedup=-1
#         multiscale=[]
#         interactive=False
#         mask=mask_name
#         imsize=image_size
#         cell=['1.5arcsec','1.5arcsec']
#         stokes='I'
#         weighting='briggs'
#         robust=0.8
#         uvtaper=[]
#         modelimage=''
#         restoringbeam=['']
#         pblimit=-0.2
#         pbcor=False
#         usescratch=False
#         allowchunk=False
#         async=False
#         tclean()
# 
# # Make test cube of phase calibrator, all channels
#         logprint ('Make cube of J0943-0819, all channels, Spw='+str(ii), logfileout='logs/testcubes.log')
# 
#         default('tclean')
#     
#         phasefile=phasems+'_spw'+str(ii)+'.ms'
#         image_name='phasecube_spw'+str(ii)
#         fieldid='J0943-0819'
#         iteration=10000   # Total number of iterations for the entire cube (per spw), ~4 per channel 
# 
#         vis=phasefile
#         imagename=image_name
#         selectdata=False
#         field=fieldid
#         spw=''
#         specmode='cubedata'
#         nterms=1
#         niter=iteration
#         gain=0.1
#         gridmode=grid_mode
#         wprojplanes=number_w
#         threshold='10.0mJy'  # Adding threshold
#         deconvolver='clark'
#         imagermode='csclean'
#         cyclefactor=1.5
#         cyclespeedup=-1
#         multiscale=[]
#         interactive=False
#         mask=mask_name
#         imsize=image_size
#         cell=['1.5arcsec','1.5arcsec']
#         stokes='I'
#         weighting='briggs'
#         robust=0.8
#         uvtaper=[]
#         modelimage=''
#         restoringbeam=['']
#         pblimit=-0.2
#         pbcor=False
#         usescratch=False
#         allowchunk=False
#         async=False
#         tclean()


# Make test cube of 10 mJy source, all channels
        logprint ('Make cube of 10 mJy source, all channels, Spw='+str(ii), logfileout='logs/testcubes.log')

        default('tclean')
    
        targetfile=targetms+'_spw'+str(ii)+'.ms'
        image_name='targetcube_10mJy_spw'+str(ii)
        fieldid='deepfield'
        phasecenter='J2000 10h01m31.4 +02d26m40'
        iteration=2048   # Total number of iterations for the entire cube (per spw), ~1 per channel 

        vis=targetfile
        imagename=image_name
        selectdata=False
        field=fieldid
        spw=''
        specmode='cubedata'
        nterms=1
        niter=iteration
        gain=0.1
        gridmode=grid_mode
        wprojplanes=number_w
        threshold='10.0mJy'  # Adding threshold
        deconvolver='clark'
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
        uvtaper=[]
        modelimage=''
        restoringbeam=['']
        pblimit=-0.2
        pbcor=False
        usescratch=False
        allowchunk=False
        async=False
        tclean()

# Make small cube of lowest z detection
        if ii==13:
            logprint ('Make cube of strongest pilot detection, Spw='+str(ii), logfileout='logs/testcubes.log')

            default('tclean')
        
            image_name='targetcube_HIdet_spw'+str(ii)
            fieldid='deepfield'
            phasecenter='J2000 10h01m15.2 +02d18m24'
            iteration=10000   # Total number of iterations for the entire spw cube, ~4 per channel, only used for HI source spw 
        
            vis=targetfile
            imagename=image_name
            selectdata=False
            field=fieldid
            spw=''
            specmode='cubedata'
            nterms=1
            niter=iteration
            gain=0.1
            gridmode=grid_mode
            wprojplanes=number_w
            threshold='10.0mJy'
            deconvolver='clark'
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
            uvtaper=[]
            modelimage=''
            restoringbeam=['']
            pblimit=-0.2
            pbcor=False
            usescratch=False
            allowchunk=False
            async=False
            tclean()
                        
    except:
        logprint("Unable to image "+str(field)+", spw "+str(ii), logfileout='logs/testcubes.log')

logprint("Finished making all cubes", logfileout='logs/testcubes.log')


# Prepare to extract information about cubes.
# Initialize variables

sigma=np.zeros(15,float)

# bmaj
# bmaj_f=[]
# bmaj_p=[]
bmaj_t10=[]
bmaj_tHI=[]

# bmin
# bmin_f=[]
# bmin_p=[]
bmin_t10=[]
bmin_tHI=[]

# rms
# fsigma=[]
# psigma=[]
t10sigma=[]
tHIsigma=[]

# flux
# bp_flux=[]
# p_flux=[]
t10_flux=[]
tHI_flux=[]

# Frequencies
# freq_f=[]
# freq_p=[]
freq_t10=[]
freq_tHI=[]

for ii in seq:       
    logprint ('Extract Data from Cubes, Spw='+str(ii), logfileout='logs/testcubes.log')
    
# #Extract Flux calibrator data
#     logprint ('Extract Flux Calibrator data, Spw='+str(ii), logfileout='logs/testcubes.log')
# 
#     try:
#         image_name='fluxcube_spw'+str(ii)+'.image'
#         header=imhead(image_name)
#         default('imval')
#         stokes='I'
#         imagename=image_name
#         results=imval()
#         for i in range(2048):
#             bmaj_f.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['major']['value'])
#             bmin_f.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['minor']['value'])
#             freq_f.append(results['coords'][i,3]/1e6)
#             bp_flux.append(results['data'][i])
#                     
#         default('imval')
#         stokes='I'
#         box='0,0,16,16'
#         imagename=image_name
#         results=imval()
#         for i in range(2048):
#             fsigma.append(np.std(results['data'][:,:,i]))
#         
#     except:
#         logprint("Unable to extract data from flux calibrator cube, spw "+str(ii), logfileout='logs/testcubes.log')
#         
#  
# #Extract Phase calibrator data
#     logprint ('Extract Phase Calibrator data, Spw='+str(ii), logfileout='logs/testcubes.log')
#     try:
#         image_name='phasecube_spw'+str(ii)+'.image'
#         header=imhead(image_name)
#         default('imval')
#         stokes='I'
#         imagename=image_name
#         results=imval()
#         for i in range(2048):
#             bmaj_p.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['major']['value'])
#             bmin_p.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['minor']['value'])
#             freq_p.append(results['coords'][i,3]/1e6)
#             p_flux.append(results['data'][i])
#                     
#         default('imval')
#         stokes='I'
#         box='0,0,16,16'
#         imagename=image_name
#         results=imval()
#         for i in range(2048):
#             psigma.append(np.std(results['data'][:,:,i]))
#     except:
#         logprint("Unable to extract data from phase calibrator cube, spw "+str(ii), logfileout='logs/testcubes.log')

#Extract spectra for 10 mJy cube
    try:
        logprint ('Extract Spectra from 10 mJy cube, Spw='+str(ii), logfileout='logs/testcubes.log')
        image_name='targetcube_10mJy_spw'+str(ii)+'.image'
        header=imhead(image_name)
        default('imval')
        stokes='I'
        imagename=image_name
        results=imval()
        for i in range(2048):
            bmaj_t10.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['major']['value'])
            bmin_t10.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['minor']['value'])
            freq_t10.append(results['coords'][i,3]/1e6)
            t10_flux.append(results['data'][i])
            
        default('imval')
        stokes='I'
        box='0,0,16,16'
        imagename=image_name
        results=imval()
        for i in range(2048):
            t10sigma.append(np.std(results['data'][:,:,i]))

# Measure noise averaged over entire spw
        sigma[ii]=np.std(results['data'][:,:,:])
        
    except:
        logprint("Unable to extract data from 10 mJy cube, spw "+str(ii), logfileout='logs/testcubes.log')

#Extract spectrum for deepfield cubes
    if ii==13:
        try:    
            logprint ('Extract Spectra from HI source cube, Spw='+str(ii), logfileout='logs/testcubes.log')    
            image_name='targetcube_HIdet_spw'+str(ii)+'.image'
            header=imhead(image_name)
            default('imval')
            stokes='I'
            box='16,16,48,48'
            imagename=image_name
            results=imval()
            for i in range(2048):
                bmaj_tHI.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['major']['value'])
                bmin_tHI.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['minor']['value'])
                freq_tHI.append(results['coords'][0,0,i,3]/1e6)
                tHI_flux.append(np.mean(results['data'][:,:,i]))
                
            default('imval')
            stokes='I'
            box='0,0,16,16'
            imagename=image_name
            results=imval()
            for i in range(2048):
                tHIsigma.append(np.std(results['data'][:,:,i]))
                    
        except:
            logprint("Unable to extract data from HI source cube, spw "+str(ii), logfileout='logs/testcubes.log')

# Convert lists to arrays and replace zeros with nan (for plotting)
# bmaj_f=np.array(bmaj_f,float)
# bmaj_p=np.array(bmaj_p,float)
bmaj_t10=np.array(bmaj_t10,float)
bmaj_tHI=np.array(bmaj_tHI,float)
# bmin_f=np.array(bmin_f,float)
# bmin_p=np.array(bmin_p,float)
bmin_t10=np.array(bmin_t10,float)
bmin_tHI=np.array(bmin_tHI,float)
# fsigma=np.array(fsigma,float)
# psigma=np.array(psigma,float)
t10sigma=np.array(t10sigma,float)
tHIsigma=np.array(tHIsigma,float)
# bp_flux=np.array(bp_flux,float)
# p_flux=np.array(p_flux,float)
t10_flux=np.array(t10_flux,float)
tHI_flux=np.array(tHI_flux,float)

# bmaj_f[bmaj_f==0.0]=np.nan
# bmaj_p[bmaj_p==0.0]=np.nan
bmaj_t10[bmaj_t10==0.0]=np.nan
bmaj_tHI[bmaj_tHI==0.0]=np.nan
# bmin_f[bmin_f==0.0]=np.nan
# bmin_p[bmin_p==0.0]=np.nan
bmin_t10[bmin_t10==0.0]=np.nan
bmin_tHI[bmin_t10==0.0]=np.nan
# fsigma[fsigma==0.0]=np.nan
# psigma[psigma==0.0]=np.nan
t10sigma[t10sigma==0.0]=np.nan
tHIsigma[tHIsigma==0.0]=np.nan
# bp_flux[bp_flux==0.0]=np.nan
# p_flux[p_flux==0.0]=np.nan
t10_flux[t10_flux==0.0]=np.nan
tHI_flux[tHI_flux==0.0]=np.nan

# Plot spectra for flux cal
logprint("Plotting Results", logfileout='logs/testcubes.log')
# fig=pylab.figure()
# pylab.plot(freq_f,bmaj_f,'b-')
# pylab.plot(freq_f,bmin_f,'r-')
# pylab.xlim(min(freq_f),max(freq_f))
# pylab.ylim(4.,25.)
# pylab.xlabel('Frequency [MHz]')
# pylab.ylabel('Beam size ["]')
# pylab.title('Beamsize for 3C286')
# pylab.savefig('image_fluxcal_beamsize.png')
# pylab.close(fig)
# 
# fig=pylab.figure()
# pylab.plot(freq_f,bp_flux,'k-')
# pylab.xlim(min(freq_f),max(freq_f))
# pylab.xlabel('Frequency [MHz]')
# pylab.ylabel('Flux [Jy/bm]')
# pylab.title('Flux of 3C286')
# pylab.savefig('image_fluxcal_spectrum.png')
# pylab.close(fig)
# 
# fig=pylab.figure()
# pylab.plot(freq_f,fsigma,'k-')
# pylab.xlim(min(freq_f),max(freq_f))
# pylab.xlabel('Frequency [MHz]')
# pylab.ylabel('RMS [Jy/bm]')
# pylab.title('Noise of 3C286 Image')
# pylab.savefig('image_fluxcal_rms.png')
# pylab.close(fig)
# 
# # Plot spectra for phase cal
# fig=pylab.figure()
# pylab.plot(freq_p,bmaj_p,'b-')
# pylab.plot(freq_p,bmin_p,'r-')
# pylab.xlim(min(freq_p),max(freq_p))
# pylab.ylim(4.,12.)
# pylab.xlabel('Frequency [MHz]')
# pylab.ylabel('Beam size ["]')
# pylab.title('Beamsize for J0943-0819')
# pylab.savefig('image_phasecal_beamsize.png')
# pylab.close(fig)
# 
# fig=pylab.figure()
# pylab.plot(freq_p,p_flux,'k-')
# pylab.xlim(min(freq_p),max(freq_p))
# pylab.xlabel('Frequency [MHz]')
# pylab.ylabel('Flux [Jy/bm]')
# pylab.title('Flux of J0943-0819')
# pylab.savefig('image_phasecal_spectrum.png')
# pylab.close(fig)
# 
# fig=pylab.figure()
# pylab.plot(freq_p,psigma,'k-')
# pylab.xlim(min(freq_p),max(freq_p))
# pylab.xlabel('Frequency [MHz]')
# pylab.ylabel('RMS [Jy/bm]')
# pylab.title('Noise of J0943-0819 Image')
# pylab.savefig('image_phasecal_rms.png')
# pylab.close(fig)

# Plot spectra for 10 mJy source
fig=pylab.figure()
pylab.plot(freq_t10,bmaj_t10,'b-')
pylab.plot(freq_t10,bmin_t10,'r-')
pylab.xlim(min(freq_t10),max(freq_t10))
pylab.ylim(4.,12.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Beam size ["]')
pylab.title('Beamsize for 10 mJy Source')
pylab.savefig('image_10mJy_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(freq_t10,t10_flux,'k-')
pylab.xlim(min(freq_t10),max(freq_t10))
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Flux [Jy/bm]')
pylab.title('Flux of 10 mJy source')
pylab.savefig('image_10mJy_spectrum.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(freq_t10,t10sigma,'k-')
pylab.xlim(min(freq_t10),max(freq_t10))
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('RMS [Jy/bm]')
pylab.title('Noise of 10 mJy source')
pylab.savefig('image_10mJy_rms.png')
pylab.close(fig)

# Plot spectra for HI field
fig=pylab.figure()
pylab.plot(freq_tHI,bmaj_tHI,'b-')
pylab.plot(freq_tHI,bmin_tHI,'r-')
pylab.xlim(min(freq_tHI),max(freq_tHI))
pylab.ylim(4.,12.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Beam size ["]')
pylab.title('Beamsize for HI field')
pylab.savefig('image_HIdet_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(freq_tHI,tHI_flux,'k-')
pylab.xlim(min(freq_tHI),max(freq_tHI))
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Flux [Jy/bm]')
pylab.title('Flux of HI field')
pylab.savefig('image_HIdet_spectrum.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(freq_tHI,tHIsigma,'k-')
pylab.xlim(min(freq_tHI),max(freq_tHI))
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('RMS [Jy/bm]')
pylab.title('Noise of HI field')
pylab.savefig('image_HIdet_rms.png')
pylab.close(fig)

#Move plots, images to sub-directory
os.system("mv *.png plots")
os.system("mv *cube*.* images")

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
# wlog.write('<li> Beam Size vs. Frequency for 3C286</li>\n')
# wlog.write('<li><img src="plots/image_fluxcal_beamsize.png"></li><br>\n')
# wlog.write('<li> Flux vs. Frequency for 3C286</li>\n')
# wlog.write('<li><img src="plots/image_fluxcal_spectrum.png"></li><br>\n')
# wlog.write('<li> RMS Noise vs. Frequency for 3C286</li>\n')
# wlog.write('<li><img src="plots/image_fluxcal_rms.png"></li><br>\n')
# wlog.write('<li> Beam Size vs. Frequency for J0943-0819</li>\n')
# wlog.write('<li><img src="plots/image_phasecal_beamsize.png"></li><br>\n')
# wlog.write('<li> Flux vs. Frequency for J0943-0819</li>\n')
# wlog.write('<li><img src="plots/image_phasecal_spectrum.png"></li><br>\n')
# wlog.write('<li> RMS Noise vs. Frequency for J0943-0819</li>\n')
# wlog.write('<li><img src="plots/image_phasecal_rms.png"></li><br>\n')
wlog.write('<li> Beam Size vs. Frequency for 10 mJy source</li>\n')
wlog.write('<li><img src="plots/image_10mJy_beamsize.png"></li><br>\n')
wlog.write('<li> Flux vs. Frequency for 10 mJy source</li>\n')
wlog.write('<li><img src="plots/image_10mJy_spectrum.png"></li><br>\n')
wlog.write('<li> RMS Noise vs. Frequency for 10 mJy source</li>\n')
wlog.write('<li><img src="plots/image_10mJy_rms.png"></li><br>\n')
wlog.write('<li> Beam Size vs. Frequency for HI source</li>\n')
wlog.write('<li><img src="plots/image_HIdet_beamsize.png"></li><br>\n')
wlog.write('<li> Flux vs. Frequency for HI source</li>\n')
wlog.write('<li><img src="plots/image_HIdet_spectrum.png"></li><br>\n')
wlog.write('<li> RMS Noise vs. Frequency for HI source</li>\n')
wlog.write('<li><img src="plots/image_HIdet_rms.png"></li><br>\n')
wlog.write('<br>\n')
wlog.write('<li> RMS Noise per Spw for HI source</li>\n')
wlog.write('<table>\n')
for ii in seq:
    wlog.write('<tr>\n')
    wlog.write('<td> Noise in spw '+str(ii)+'= '+str(sigma[ii]*1e3)+' mJy</td>\n')
    wlog.write('</tr>\n')
wlog.write('</table>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()

logprint ("Finished CHILES_pipe_testcubes.py", logfileout='logs/testcubes.log')
time_list=runtiming('testcubes', 'end')

pipeline_save()
