# CHILES_pipe_phasecal_plots.py
# This module of the CHILES pipeline re-applies phase calibration and makes
# the diagnostic plots again.    
# 9/23/16 DJP


logprint ("Starting CHILES_pipe_phasecal_plots.py", logfileout='logs/phasecal.log')
time_list=runtiming('phase', 'start')

# Step 1:  Set default parameters (to match previous modules), and clean up old files:
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

# Remove old images of phase calibrator
os.system("rm -rf images/phasecalibrator_spw*.*")


# Step 10: Makes diagnostic plots for assessment
# Look at calibration tables for phase cal (amp, phase vs. time, frequency)
# Make images of phase cal and look at flux,beam vs. spw

logprint ("Making diagnostic plots", logfileout='logs/phasecal.log')

#Plot UV spectrum (averaged over all baselines & time) of phase calibrator
default('plotms')
vis=ms_active   
field='0'        # Only for J0943
xaxis='freq'
yaxis='amp'
xdatacolumn='corrected'
ydatacolumn='corrected'
averagedata=True
avgtime='1e5'
avgscan=True
avgbaseline=True
scalar=False
spw=''
avgspw=False
showlegend=False
coloraxis='spw'
plotrange=[0.95,1.43,0,0]
showgui=False
clearplots=True
plotfile='phasecal_spectrum_full.png'
plotms()

plotrange=[0.95,1.43,2.7,3.5]
plotfile='phasecal_spectrum_zoom.png'
plotms()


ms_name=ms_active[:-3]
output_ms=ms_name+'_phasecal_flux_averaged.ms'

# Remove old averaged MS if it exists.
if os.path.exists(output_ms):
    os.system("rm -rf "+output_ms)
    
#Split the phase cal.
default('split')
vis=ms_active
outputvis=output_ms
datacolumn='corrected'
field='0'
spw='0~14:128~1920'  # Average over all channels, except the very edges
width='1792'
antenna=''
timebin=''
timerange=''
scan=''
intent=''
array=''
uvrange=''
correlation=''
observation=''
keepflags=False
keepmms=False
split()

# Use plotms to make plot of amplitude vs. phase, divided by spw
seq=range(0,15)
for ii in seq:
    default('plotms')
    vis=output_ms  # File only contains J0943, phase calibrator data
    field=''       
    xaxis='phase'
    yaxis='amp'
    xdatacolumn='corrected'
    ydatacolumn='corrected'
    averagedata=False  # Data already averaged
    spw=str(ii)
    gridrows=1
    showlegend=False
    iteraxis='spw'
    coloraxis='corr'
    showgui=False
    clearplots=True
    plotfile='phasecal_ampphase.png'
    plotms()

# Use plotms to make plot of amplitude and phase vs. time, for each spw.  
for ii in seq:
    default('plotms')
    vis=output_ms  # File only contains J0943, phase calibrator data
    field=''       
    xaxis='time'
    yaxis='amp'
    xdatacolumn='corrected'
    ydatacolumn='corrected'
    averagedata=False  # Data already averaged
    spw=str(ii)
    gridrows=1
    showlegend=False
    iteraxis='spw'
    coloraxis='corr'
    showgui=False
    clearplots=True
    plotfile='phasecal_amptime.png'
    plotms()
    yaxis='phase'
    plotfile='phasecal_phasetime.png'
    plotms()

seq=range(0,15)
#Image phase cal: 
for ii in seq:
    print 'STARTS IMAGING PHASE CALIBRATOR OF SPW='+str(ii)
    default('clean')
    image_name='phasecalibrator_spw'+str(ii)
    fieldid='J0943*'
    grid_mode=''
    number_w=1
    image_size=[512,512]
    iteration=1000
    mask_name=['']
    vis=output_ms
    imagename=image_name
    selectdata=False
    field=fieldid
    spw=str(ii)
    mode='mfs'
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

#     Calculate bmaj, bmin, peak, and rms for images of flux calibrator in each spw.
#     Output results as text file and plots of these values.
box_phase='300,30,460,200'

bmaj_phase=[]
bmin_phase=[]
max_phase=[]
rms_phase=[]

for ii in seq:
    
    image_phasecal='phasecalibrator_spw'+str(ii)+'.image'
    bmaj1=imhead(image_phasecal,mode='get',hdkey='beammajor')
    if bmaj1==None:
        bmaj_phase.append(0.0)
        bmin_phase.append(0.0)
        max_phase.append(0.0)
        rms_phase.append(0.0)
    else:
        bmaj11=bmaj1['value']
        bmaj_phase.append(bmaj11)
        bmin1=imhead(image_phasecal,mode='get',hdkey='beamminor')
        bmin11=bmin1['value']
        bmin_phase.append(bmin11)
        imst1=imstat(image_phasecal)
        max1=imst1['max']
        max_phase.append(max1[0])
        imst1=imstat(image_phasecal,box=box_phase)
        rms1=imst1['rms']
        rms_phase.append(rms1[0]*1e3)

# Output Image statistics in text file
f=open('statisticsPhase.txt','w')

print >> f, "Phase calibrator"
print >> f, "major axis [\"] \t", "".join(["%12.3f \t" %x for x in bmaj_phase])
print >> f, "minor axis [\"] \t", "".join(["%12.3f \t" %x for x in bmin_phase])
print >> f, "peak value [Jy] \t", "".join(["%12.3f \t" %x for x in max_phase])
print >> f, "RMS noise [mJy] \t", "".join(["%12.3f \t" %x for x in rms_phase])
f.close()

#Make plots of bmaj, bmin, peak, and rms
fig=pylab.figure()
pylab.plot(seq,bmaj_phase,'r--x',label='Bmaj')
pylab.plot(seq,bmin_phase,'b--x',label='Bmin')
pylab.xlabel('Spectral Window')
pylab.ylabel('Beam Size ["]')
pylab.legend()
pylab.savefig('phasecal_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,max_phase,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('Peak Flux [Jy]')
#pylab.yscale('log')
pylab.savefig('phasecal_peak.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,rms_phase,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('RMS [mJy]')
pylab.yscale('log')
pylab.savefig('phasecal_rms.png')
pylab.close(fig)

#Want to plot image of phase calibrator in each spw.  Use "imview"
# Dropping figures of phase calibrator, leaving in diagnostic plots.
#for ii in seq:
#    image_phasecal='phasecalibrator_spw'+str(ii)+'.image'
#    kntr_levels=[-2*rms_phase[ii]/1000.,2*rms_phase[ii]/1000.,0.1*max_phase[ii],0.3*max_phase[ii],0.5*max_phase[ii],0.7*max_phase[ii],0.9*max_phase[ii]]
#    imview(raster={'file':image_phasecal,'scaling':-3, 'colorwedge':True},contour={'file':image_phasecal,'levels':kntr_levels},out='phasecal_spw'+str(ii)+'.png')

#Plot calibration tables
default('plotcal')
caltable='finalamp.gcal'
xaxis='time'
yaxis='amp'
showgui=False
figfile='caltable_finalamp_amp.png'
plotcal()

yaxis='phase'
figfile='caltable_finalamp_phase.png'
plotcal()

# Removed these plots since I don't find them useful.  DJP 2/19/16
#default('plotcal')
#caltable='finalflux.gcal'
#xaxis='time'
#yaxis='amp'
#showgui=False
#figfile='caltable_finalflux_amp.png'
#plotcal()
#
#yaxis='phase'
#figfile='caltable_finalflux_phase.png'
#plotcal()

#Move plots, images to sub-directory

os.system("mv *.png plots")
os.system("mv phasecalibrator_spw*.* images")

#Create webpage with results

if os.path.exists('phasecal.html'):
    os.system("rm phasecal.html")
wlog = open("phasecal.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_phasecal results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/phasecal.log">Phasecal Log</a></li>\n')
wlog.write('<li>finalamp.gcal Amp vs. Time: \n')
wlog.write('<br><img src="plots/caltable_finalamp_amp.png"></li>\n')
wlog.write('<li>finalamp.gcal Phase vs. Time: \n')
wlog.write('<br><img src="plots/caltable_finalamp_phase.png"></li>\n')
#wlog.write('<li>finalflux.gcal Amp vs. Time: \n')
#wlog.write('<br><img src="plots/caltable_finalflux_amp.png"></li>\n')
#wlog.write('<li>finalflux.gcal Phase vs. Time: \n')
#wlog.write('<br><img src="plots/caltable_finalflux_phase.png"></li>\n')
wlog.write('<li> Amp vs. Phase: \n')
for ii in seq:
    wlog.write('<br><img src="plots/phasecal_ampphase_Spw'+str(ii)+'.png">\n')
wlog.write('<li> Spectrum of Phase calibrator (both LL & RR, averaged over all time & baselines): \n')
wlog.write('<br><img src="plots/phasecal_spectrum_full.png">\n')
wlog.write('<br><img src="plots/phasecal_spectrum_zoom.png">\n')
wlog.write('<li> Amp. & Phase vs. time for Phase Calibrator (averaged over frequency): \n')
for ii in seq:
    wlog.write('<br><img src="plots/phasecal_amptime_Spw'+str(ii)+'.png">\n')
    wlog.write('<br><img src="plots/phasecal_phasetime_Spw'+str(ii)+'.png">\n')
wlog.write('</li>')
wlog.write('<li> Measured properties of phase calibrator: \n')
wlog.write('<br><img src="plots/phasecal_beamsize.png">\n')
wlog.write('<br><img src="plots/phasecal_peak.png">\n')
wlog.write('<br><img src="plots/phasecal_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Percentage of J0943 data flagged: '+str(phase_flag*100)+'\n')
wlog.write('<br>')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()


logprint ("Finished CHILES_pipe_phasecal_plots.py", logfileout='logs/phasecal.log')
time_list=runtiming('phase', 'end')

pipeline_save()
