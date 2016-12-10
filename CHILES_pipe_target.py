# -*- coding: utf-8 -*-
# CHILES_pipe_target.py
# This task runs after the "phasecal" task so all calibrations have been applied
# to the target field.  This task runs RFLAG and extend on the target field and
# makes some diagnostic plots.  
# 9/9/15 DJP
# Removed clipping due to overflagging 12/7/15 DJP
# Set extendpols=False in RFLAG 1/27/16 DJP
# Included flag summary for target, and applycal for target moved here, 2/9/16 DJP
# Removed averaging from plots, 2/11/16 DJP
# 2/15/16 DJP: Changed ff value to a float and re-implemented "extend".
# 2/18/16 DJP: Plot UV spectrum of field (averaging over time, baseline)
# 2/18/16 DJP: Make diagnostic plots using averaged data.
# 2/19/16 DJP: Make 2 UVSPEC plots (one with full range, one with zoom).  Changed averaging.  
# 6/28/16 DJP: Including time/frequency averaging on the target after initial rflag and extend.
# 7/1/16 DJP:  Including time-averaged RFLAG and using parameters determined by XF.  Setting rflag scale=3, growtime=60.  
# 9/21/16 DJP:  Fixed variable names for time-averaged RFLAG. 
# 12/7/16 DJP:  Changed plots to frequency on x-axis. Including amp v. time plots. Saving flagging statistics to file.
# 12/9/16 DJP:  Added additional flagging to completely flag any channel that is already more than 90% flagged (code from XF).

logprint ("Starting CHILES_pipe_target.py", logfileout='logs/target.log')
time_list=runtiming('target', 'start')

# If running this module right after starting casa, must run CHILES_pipe_restore first.
# Step 1, load needed functions.
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

# Clean up data from past runs of module:
os.system("rm -rf images/target_spw*.*")
os.system("rm -rf plots/target_*.png")
os.system("rm -rf *_target_flux_averaged.ms")

logprint ("Apply calibration to target", logfileout='logs/target.log')

TargetTables=copy.copy(priorcals)
TargetTables.append('finalphase_scan.gcal')
TargetTables.append('finalamp.gcal')
TargetTables.append('finalflux.gcal')
TargetSpwMapValues=copy.copy(priorspwmap)
TargetSpwMapValues.append([])
TargetSpwMapValues.append([])
TargetSpwMapValues.append([])

if os.path.exists('antposcal.p')==True:
  TargetFields=['','','2','2','0','0','0']

if os.path.exists('antposcal.p')==False:
  TargetFields=['','2','2','0','0','0']

default('applycal')
vis=ms_active
field='1'            # Apply calibration to target
spw=''
intent=''
selectdata=True
gaintable=TargetTables    # HG : Change this from AllCalTables to TargetTables
gainfield=TargetFields
interp=['']
spwmap=TargetSpwMapValues # Was [] in previous version, now corresponds to TargetTables
gaincurve=False
opacity=[]
parang=False
calwt=False
flagbackup=False
async=False
applycal()


# Step 2: RFLAG with extend on deepfield
logprint("Flag deepfield with RFLAG",logfileout='logs/target.log')
 
# Using parameters determined by XF.  Set sigma clip level at 3.0x the noise level

ff=float(1)
freqnoavg=0.4
timenoavg=0.5
scaling1=[2.3,1.8,1.5,1.3,1.3,1.2,1.2,1.2,1.1,1.0,1.0,1.0,1.0,1.0,1.0]
scaling=np.asarray(scaling1)
sigmacut=4.0
freqd1noavg=scaling*freqnoavg*sigmacut
timed1noavg=scaling*timenoavg*sigmacut

default('flagdata')
vis=ms_active
mode='rflag'
field='deepfield'
spw='0~14'
correlation=''
ntime='scan'
combinescans=False
datacolumn='corrected'
extendflags=False       
extendpols=False      
winsize=3
freqdev=[[ff,0.0,freqd1noavg[0]],[ff,1.0,freqd1noavg[1]],[ff,2.0,freqd1noavg[2]],[ff,3.0,freqd1noavg[3]],[ff,4.0,freqd1noavg[4]],[ff,5.0,freqd1noavg[5]],[ff,6.0,freqd1noavg[6]],[ff,7.0,freqd1noavg[7]],[ff,8.0,freqd1noavg[8]],[ff,9.0,freqd1noavg[9]],[ff,10.0,freqd1noavg[10]],[ff,11.0,freqd1noavg[11]],[ff,12.0,freqd1noavg[12]],[ff,13.0,freqd1noavg[13]],[ff,14.0,freqd1noavg[14]]]
timedev=[[ff,0.0,timed1noavg[0]],[ff,1.0,timed1noavg[1]],[ff,2.0,timed1noavg[2]],[ff,3.0,timed1noavg[3]],[ff,4.0,timed1noavg[4]],[ff,5.0,timed1noavg[5]],[ff,6.0,timed1noavg[6]],[ff,7.0,timed1noavg[7]],[ff,8.0,timed1noavg[8]],[ff,9.0,timed1noavg[9]],[ff,10.0,timed1noavg[10]],[ff,11.0,timed1noavg[11]],[ff,12.0,timed1noavg[12]],[ff,13.0,timed1noavg[13]],[ff,14.0,timed1noavg[14]]]
timedevscale=1.0
freqdevscale=1.0
action='apply'
display=''
flagbackup=False
savepars=True
flagdata()

clearstat()

logprint ("Extending flags on target", logfileout='logs/target.log')


default('flagdata')
vis=ms_active
mode='extend'
field='deepfield'
datacolumn='corrected'
correlation=''
combinescans=False
extendpols=False      
growtime=60       
growfreq=90       
growaround=True       
flagneartime=True      
flagnearfreq=True       
action='apply'
display=''
flagbackup=False    
savepars=False
flagdata()

clearstat()

# time-averaged rflag for target

logprint ("Time-averaged flagging on target", logfileout='logs/target.log')

# Using parameters determined by XF.  Set sigma clip level at 3.0x the noise level

freqavgval=0.08
timeavgval=9.2e-6


scaling1=[2.3,1.8,1.5,1.3,1.3,1.2,1.2,1.2,1.1,1.0,1.0,1.0,1.0,1.0,1.0]
scaling=np.asarray(scaling1)
sigmacut=3.0
freqd1avg=scaling*freqavgval*sigmacut
timed1avg=scaling*timeavgval*sigmacut

default('flagdata')
vis=ms_active
mode='rflag'
field='deepfield'
correlation=''
ntime='scan'
combinescans=False
datacolumn='corrected'
extendpols=False      
extendflags=False        
freqdev=[[ff,0.0,freqd1avg[0]],[ff,1.0,freqd1avg[1]],[ff,2.0,freqd1avg[2]],[ff,3.0,freqd1avg[3]],[ff,4.0,freqd1avg[4]],[ff,5.0,freqd1avg[5]],[ff,6.0,freqd1avg[6]],[ff,7.0,freqd1avg[7]],[ff,8.0,freqd1avg[8]],[ff,9.0,freqd1avg[9]],[ff,10.0,freqd1avg[10]],[ff,11.0,freqd1avg[11]],[ff,12.0,freqd1avg[12]],[ff,13.0,freqd1avg[13]],[ff,14.0,freqd1avg[14]]]
timedev=[[ff,0.0,timed1avg[0]],[ff,1.0,timed1avg[1]],[ff,2.0,timed1avg[2]],[ff,3.0,timed1avg[3]],[ff,4.0,timed1avg[4]],[ff,5.0,timed1avg[5]],[ff,6.0,timed1avg[6]],[ff,7.0,timed1avg[7]],[ff,8.0,timed1avg[8]],[ff,9.0,timed1avg[9]],[ff,10.0,timed1avg[10]],[ff,11.0,timed1avg[11]],[ff,12.0,timed1avg[12]],[ff,13.0,timed1avg[13]],[ff,14.0,timed1avg[14]]]
timedevscale=1.0
freqdevscale=1.0
channelavg=False
chanbin=1
timeavg=True
timebin='1000s'
action='apply'
flagbackup=False
savepars=True
flagdata()

clearstat()

#EM: back to normal logger output

#casalog.filter("INFO")

# Summary of flagging, after final flagging (for testing purposes only)
logprint ("Summary of flags after flagging target", logfileout='logs/target.log')
default('flagdata')
vis=ms_active
mode='summary'
spw='0~14'
field='1'           # Only flagged deepfield, no need to check others here.
correlation='RR,LL'
spwchan=True
spwcorr=True
action='calculate'
s_t=flagdata()


target_flag=s_t['flagged']/s_t['total']
logprint ("Percentage of all data flagged after RFLAG+Extend+time-avg RFLAG: "+str(target_flag*100)+"%", logfileout='logs/target.log')

# Flagging any channel that is already more than 90% flagged.  
limFlag=0.9
flagChannels=[]
for a in s_t['spw:channel']:
    fp=s_t['spw:channel'][a]['flagged']/s_t['spw:channel'][a]['total']
    if  fp>limFlag:
        flagChannels.append(a)
strChan = ','.join(flagChannels)

flagdata(vis=ms_active,field='1',mode="manual",spw=strChan, autocorr=False)


# Save final version of flags
logprint("Saving flags",logfileout='logs/target.log')

default('flagmanager')
vis=ms_active
mode='save'
versionname='targetflags'
comment='Final flags saved after calibrations and target rflag+extend'
merge='replace'
flagmanager()
logprint ("Flag column saved to "+versionname, logfileout='logs/target.log')


# Step 3: Diagnostic Plots
logprint("Making diagnostic plots",logfileout='logs/target.log')

# Make plot of percentage of data flagged as a function of channel (for both correlations combined):
flag_frac=[]
ct=-1
chan=[]
freq=[]
flagged=[]
totals=[]
# Extract frequency of first channel of spw=0 from listobs output
nu0=reference_frequencies[0]/1.e6 #get reference frequency in MHz
dnu=0.015625 # channel width in MHz
freq=[]

for s in range(15):
    for c in range(2048):
        ct+=1
        chan.append(ct)
        freq.append(nu0+dnu*ct)
        flagged.append(s_t['spw:channel'][str(s)+':'+str(c)]['flagged'])
        totals.append(s_t['spw:channel'][str(s)+':'+str(c)]['total'])
        flag_frac.append(flagged[ct]/totals[ct])

fig=pylab.figure()
pylab.plot(freq,flag_frac,'k-')
pylab.xlim(940.,1445.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Fraction of data flagged')
pylab.savefig("target_flag.png")
pylab.close(fig)

#Plot UV spectrum (averaged over all baselines & time) of target
default('plotms')
vis=ms_active   
field='1'        # Only for deepfield
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
plotfile='target_spectrum_full.png'
plotms()

plotrange=[0.95,1.43,-0.0001,0.015]
plotfile='target_spectrum_zoom.png'
plotms()

# Use plotms to make plot of amplitude vs. frequency, divided by spw
seq=range(0,15)
for ii in seq:
    default('plotms')
    vis=ms_active  # Use standard MS
    field='1'       # Only for deepfield
    xaxis='frequency'
    yaxis='amp'
    xdatacolumn='corrected'
    ydatacolumn='corrected'
    averagedata=True
    avgtime='1e5'
    avgscan=True   
    spw=str(ii)
    gridrows=1
    showlegend=False
    iteraxis='spw'
    coloraxis='corr'
    showgui=False
    clearplots=True
    plotfile='target_ampchannel.png'
    plotms()
    
# Make amplitude vs. time plot
    xaxis='time'
    avgtime=''
    avgscan=False
    avgchannel='2048'
    plotfile='target_amptime.png'
    plotms()

#Split the target.
ms_name=ms_active[:-3]
output_ms=ms_name+'_target_flux_averaged.ms'

default('split')
vis=ms_active
outputvis=output_ms
datacolumn='corrected'
field='deepfield'
spw='0~14:128~1920'  # Extend averaging to include all but 128 edge channels
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

seq=range(0,15)
#Image target: 
for ii in seq:
    print 'STARTS IMAGING Deepfield OF SPW='+str(ii)
    default('clean')
    image_name='target_spw'+str(ii)
    fieldid='deepfield'
    grid_mode=''
    number_w=1
    image_size=[2048,2048]
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

# Measure statistics of deepfield image:
box_target='1300,1100,1900,1600'
bmaj_target=[]
bmin_target=[]
max_target=[]
rms_target=[]

seq=range(0,15)
logprint("START READING STATISTICS FROM IMHEAD() AND IMSTAT()",logfileout='logs/target.log')
for ii in seq:
	image_target='target_spw'+str(ii)+'.image'
	bmaj1=imhead(image_target,mode='get',hdkey='beammajor')
	if bmaj1==None :
		bmaj_target.append(0.0)
		bmin_target.append(0.0)
		max_target.append(0.0)
		rms_target.append(0.0)
	else:
		bmaj11=bmaj1['value']
		bmaj_target.append(bmaj11)
		bmin1=imhead(image_target,mode='get',hdkey='beamminor')
		bmin11=bmin1['value']
		bmin_target.append(bmin11)
		imst1=imstat(image_target)
		max1=imst1['max']
		max_target.append(max1[0]*1e3)
		imst1=imstat(image_target,box=box_target)
		rms1=imst1['rms']
		rms_target.append(rms1[0]*1e6)

logprint("FINISH READING STATISTICS FROM IMHEAD() AND IMSTAT()",logfileout='logs/target.log')

f=open('statisticsTarget.txt','w')
print >> f, "Deep field"
print >> f, "major axis [\"] \t", "".join(["%12.3f \t" %x for x in bmaj_target])
print >> f, "minor axis [\"] \t", "".join(["%12.3f \t" %x for x in bmin_target])
print >> f, "peak value [mJy] \t", "".join(["%12.3f \t" %x for x in max_target])
print >> f, "RMS noise [uJy] \t", "".join(["%12.3f \t" %x for x in rms_target])

f.close()

#Make plots of bmaj, bmin, peak, and rms
fig=pylab.figure()
pylab.plot(seq,bmaj_target,'r--x',label='Bmaj')
pylab.plot(seq,bmin_target,'b--x',label='Bmin')
pylab.xlabel('Spectral Window')
pylab.ylabel('Beam Size ["]')
pylab.legend()
pylab.savefig('target_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,max_target,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('Peak Flux [mJy]')
#pylab.yscale('log')
pylab.savefig('target_peak.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,rms_target,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('RMS [uJy]')
pylab.yscale('log')
pylab.savefig('target_rms.png')
pylab.close(fig)

#Want to plot image of flux calibrator in each spw.  Use "imview"

for ii in seq:
    image_target='target_spw'+str(ii)+'.image'
    kntr_levels=[-2*rms_target[ii]/1000.,2*rms_target[ii]/1000.,0.1*max_target[ii],0.3*max_target[ii],0.5*max_target[ii],0.7*max_target[ii],0.9*max_target[ii]]
    imview(raster={'file':image_target, 'colorwedge':True},contour={'file':image_target,'levels':kntr_levels},out='target_spw'+str(ii)+'.png')


#Move plots, images to sub-directory

os.system("mv *.png plots")
os.system("mv target_spw*.* images")

#Create webpage with results

if os.path.exists("target.html"):
    os.system("rm target.html")
wlog = open("target.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_target results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/target.log">Target Log</a></li>\n')
wlog.write('<li> Amp vs. Frequency (time-averaged, all channels shown) and Amp vs. Time (all channels averaged): \n')
wlog.write('<table> \n')
for ii in seq:
    wlog.write('<tr>\n')
    wlog.write('<td><img src="plots/target_ampchannel_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('<td><img src="plots/target_amptime_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('</tr>\n')
wlog.write('</table> \n')
wlog.write('<li> Spectrum of Deepfield (both LL & RR, averaged over all time & baselines): \n')
wlog.write('<br><img src="plots/target_spectrum_full.png">\n')
wlog.write('<br><img src="plots/target_spectrum_zoom.png">\n')
wlog.write('<li> Images of Deepfield: \n')
for ii in seq:
    wlog.write('<br><img src="plots/target_spw'+str(ii)+'.png">\n')
wlog.write('</li>')
wlog.write('<li> Measured properties of deepfield: \n')
wlog.write('<br><img src="plots/target_beamsize.png">\n')
wlog.write('<br><img src="plots/target_peak.png">\n')
wlog.write('<br><img src="plots/target_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Percentage of deepfield data flagged: '+str(target_flag*100)+'\n')
wlog.write('<br> Flagging percentage vs. frequency (before removing channels that are more than 90% flagged)\n')
wlog.write('<br><img src="plots/target_flag.png">\n')
wlog.write('<br>')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()


logprint ("Finished CHILES_pipe_target.py", logfileout='logs/target.log')
time_list=runtiming('target', 'end')

pipeline_save()
