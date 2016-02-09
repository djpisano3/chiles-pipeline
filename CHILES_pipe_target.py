# CHILES_pipe_target.py
# This task runs after the "phasecal" task so all calibrations have been applied
# to the target field.  This task runs RFLAG and extend on the target field and
# makes some diagnostic plots.  
# 9/9/15 DJP
# Removed clipping due to overflagging 12/7/15 DJP
# Set extendpols=False in RFLAG 1/27/16 DJP
# Included flag summary for target, and applycal for target moved here, 2/9/16 DJP

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
field='1'            # Apply same calibration to target
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
 

#EM: clip the target data by pre-set values determined by XF
#DJP, v0.7: re-instituted clipping since it does not appear to result in overflagging.  

#default(flagcmd)
#vis      = ms_active
#inpmode  ='list'
#inpfile  = pipepath+'flagcmd_clip_target.txt'
#action   ='apply'
#flagbackup= True
#savepars=True
#flagcmd()
#
#clearstat()

#EM: rflag on the target source using pre set values
ff='1'

default('flagdata')
vis=ms_active
mode='rflag'
field='deepfield'
spw='0~14'
correlation=''
ntime='scan'
combinescans=False
datacolumn='corrected'
extendflags=False        # Added (KH)
extendpols=False     # Default is True.  May allow some weak RFI through, but try it.   
winsize=3
freqdev=[[ff,0.0,6.1],[ff,1.0,4.7],[ff,2.0,3.9],[ff,3.0,3.5],[ff,4.0,3.4],[ff,5.0,3.2],[ff,6.0,3.1],[ff,7.0,3.2],[ff,8.0,2.9],[ff,9.0,2.6],[ff,10.0,2.6],[ff,11.0,2.6],[ff,12.0,2.5],[ff,13.0,2.6],[ff,14.0,2.6]]
timedev=[[ff,0.0,8.0],[ff,1.0,6.2],[ff,2.0,5.1],[ff,3.0,4.6],[ff,4.0,4.4],[ff,5.0,4.2],[ff,6.0,4.1],[ff,7.0,4.2],[ff,8.0,3.7],[ff,9.0,3.4],[ff,10.0,3.3],[ff,11.0,3.4],[ff,12.0,3.3],[ff,13.0,3.4],[ff,14.0,3.4]]
timedevscale=1.0
freqdevscale=1.0
action='apply'
display=''
flagbackup=False
savepars=True
async=False
flagdata()

clearstat()

#EM: extend the flags on the target only
#DJP: commented out in v0.7

#default('flagdata')
#vis=ms_active
#mode='extend'
#field='deepfield'
#correlation=''
#combinescans=False
#datacolumn='corrected'
#extendpols=False      
#growtime=80       
#growfreq=90       
#growaround=True      
#flagneartime=True       
#flagnearfreq=True       
#action='apply'
#display=''
#flagbackup=True    # Saving flags here will allow us to easily fix overflagging by EXTEND, if needed.  
#savepars=True
#async=False
#flagdata()
#
#clearstat()

logprint ("Not extending flags on target", logfileout='logs/target.log')

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

logprint ("Percentage of all data flagged after target module: "+s_i['flagged']/s_i['total']*100+'%', logfileout='logs/target.log')


# Save final version of flags
logprint("Saving flags",logfileout='logs/target.log')

default('flagmanager')
vis=ms_active
mode='save'
versionname='finalflags'
comment='Final flags saved after calibrations and target rflag+extend'
merge='replace'
async=False
flagmanager()
logprint ("Flag column saved to "+versionname, logfileout='logs/target.log')


# Step 3: Diagnostic Plots
logprint("Making diagnostic plots",logfileout='logs/target.log')

ms_name=ms_active[:-3]
output_ms=ms_name+'_target_flux_averaged.ms'

#Split the target.
default('split')
vis=ms_active
outputvis=output_ms
datacolumn='corrected'
field='deepfield'
spw='0~14:200~1850'
width='1651'
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
async=False
split()

# Use plotms to make plot of amplitude vs. uvdist, divided by spw
seq=range(0,15)
for ii in seq:
    default('plotms')
    vis=ms_active
    field='1'       # Hard coded to deepfield
    xaxis='uvdist'
    yaxis='amp'
    xdatacolumn='corrected'
    ydatacolumn='corrected'
    averagedata=True
    avgtime='1e5'
    avgscan=True
    #avgchannel='1e6'
    avgspw=False
    spw=str(ii)
    gridrows=1
    #gridcolumns=
    showlegend=False
    iteraxis='spw'
    coloraxis='corr'
    showgui=False
    clearplots=True
    plotfile='target_ampuvdist.png'
    plotms()



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
pylab.yscale('log')
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
wlog.write('<li> Amp vs. UVdist (time-averaged): \n')
for ii in seq:
    wlog.write('<br><img src="plots/target_ampuvdist_Spw'+str(ii)+'.png">\n')
wlog.write('<li> Images of Deepfield: \n')
for ii in seq:
    wlog.write('<br><img src="plots/target_spw'+str(ii)+'.png">\n')
wlog.write('</li>')
wlog.write('<li> Measured properties of target field: \n')
wlog.write('<br><img src="plots/target_beamsize.png">\n')
wlog.write('<br><img src="plots/target_peak.png">\n')
wlog.write('<br><img src="plots/target_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()


logprint ("Finished CHILES_pipe_target.py", logfileout='logs/target.log')
time_list=runtiming('target', 'end')

pipeline_save()
