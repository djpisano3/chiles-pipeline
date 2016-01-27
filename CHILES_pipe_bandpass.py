# CHILES_pipe_bandpass.py
# This module of the CHILES pipeline does the delay and bandpass correction
# on the data as well as the flux calibration.  It is heavily based on the code 
# written by Ximena Fernandez, and not on the NRAO continuum pipeline. 
# 8/27/15 DJP

#Need to run CHILES_pipe_restore.py first, if not immediately following CHILES_pipeline_initial.py
#8/31/15 DJP

# HG made a number of changes to make sure that the caltable lists and spw selection
# were done correctly.  I explicitly removed all prior caltables.  12/7/15 DJP

# 1/14/16, DJP: Set fillgaps=10 with interp='linear,spline' in bandpass (probably no effect)
#               Set interp=',spline' for BPcal tables for frequency interpolation.

# 1/27/16 DJP: Set extendpols=False in RFLAG step.

#Part I: define some variables that will be used later
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

logprint ("Starting CHILES_pipe_bandpass.py", logfileout='logs/bandpass.log')
time_list=runtiming('bandpass', 'start')

#Clean up calibration tables from past runs:
if os.path.exists('initialBPcal.b'):
    rmtables(tablenames='initialBPcal.b')
if os.path.exists('finalBPcal.b'):
    rmtables(tablenames='finalBPcal.b')
if os.path.exists('initialdelay.k'):
    rmtables(tablenames='initialdelay.k')
if os.path.exists('finaldelay.k'):
    rmtables(tablenames='finaldelay.k')
if os.path.exists('finalBPinitialgain.g'):
    rmtables(tablenames='finalBPinitialgain.g')
if os.path.exists('finaldelayinitialgain.g'):
    rmtables(tablenames='finaldelayinitialgain.g')
if os.path.exists('gain_curves.g'):
    rmtables(tablenames='gain_curves.g')
if os.path.exists('initialBPinitialgain.g'):
    rmtables(tablenames='initialBPinitialgain.g')
if os.path.exists('testdelayintiialgain.g'):
    rmtables(tablenames='testdelayintiialgain.g')
# Remove old images of flux calibrator
os.system("rm -rf images/fluxcalibrator_spw*.*")


#This creates an array that defines a channel range for each spectral window
tb.open(ms_active+'/SPECTRAL_WINDOW')
channels = tb.getcol('NUM_CHAN')
numSpws = len(channels)
tst_delay_spw=''
all_spw=''

for ispw in range(numSpws):
    endch1=int(channels[ispw]/3.0)
    endch2=int(2.0*channels[ispw]/3.0)+1
    if (ispw<max(range(numSpws))):
        tst_delay_spw=tst_delay_spw+str(ispw)+':'+str(endch1)+'~'+str(endch2)+','
        all_spw=all_spw+str(ispw)+','
    else:
        tst_delay_spw=tst_delay_spw+str(ispw)+':'+str(endch1)+'~'+str(endch2)
        all_spw=all_spw+str(ispw)

tst_bpass_spw=tst_delay_spw

#Set minimum # of baselines need for a solution to 8, based on experience
minBL_for_cal=8

#Set uvrange to apply in order to optimally exclude RFI without flagging:
uvr_cal='>1500m'

#Reference Antenna should have already been selected:
while bool(refAnt)==False:
    refAnt=raw_input("Enter name of one or more reference antenna(s) (e.g. ea01, ea02, etc.;this is needed to proceed): ")


#--------------------------------------------------------------------
# Part II:  Setup initial calibration tables.
logprint ("Setup initial calibration tables", logfileout='logs/bandpass.log')

#Clear previous models of flux calibrator
default('delmod')
vis=ms_active
otf=True
field='2'
scr=False
delmod()

#1 setjy:
default('setjy')
vis=ms_active
field='2'
spw=''
selectdata=False
model='3C286_L.im'
listmodimages=False
scalebychan=True
fluxdensity=-1
standard='Perley-Butler 2013'
usescratch=False         # DJP: Okay to be False in version 4.5.  
async=False
setjy()


#2: Initial gain curve, TABLE: gain_curves.g

if os.path.exists('gain_curves.g')==False:
    default(gencal)
    vis=ms_active
    caltable='gain_curves.g'
    caltype='gc'
    spw=''
    antenna=''
    pol=''
    parameter=[]
    async=False
    gencal()

priorcals=['gain_curves.g']
priorspwmap=[[]]

#3: Antenna position calibration:
if os.path.exists('antposcal.p')==False:
    default(gencal)
    vis=ms_active
    caltable='antposcal.p'
    caltype='antpos'
    spw=''
    antenna=''
    pol=''
    parameter=[]
    async=False
    gencal()

if os.path.exists('antposcal.p')==True:  # HG (12/04/2015) : include antposcal.p in the priorcals if it exists
  priorcals.append('antposcal.p')
  priorspwmap.append([])


#--------------------------------------------------------------------
#Part III: initial bandpass calibration + rflag
logprint ("Initial delay, gain, BP calibration", logfileout='logs/bandpass.log')

#4: Gain calibration on delay calibrator, TABLE: testdelayinitialgain.g
default('gaincal')
vis=ms_active
caltable='testdelayinitialgain.g'
field='2'            # Hard-coded to field 2 as this is always 3C286
spw=tst_delay_spw
intent=''
selectdata=False
solint='int'
combine='scan'
preavg=-1.0
refant=refAnt
uvrange=uvr_cal      # Set uvrange to exclude worst of RFI
minblperant= minBL_for_cal
minsnr=3.0
solnorm=False
gaintype='G'
smodel=[]
calmode='p'
append=False
docallib=False
gaintable=priorcals
gainfield=['']
interp=['']
spwmap=[]
parang=False
async=False
gaincal()

GainTables=copy.copy(priorcals)
GainTables.append('testdelayinitialgain.g')

#Delay calibration to all spws. TABLE: initialdelay.k
default('gaincal')
vis=ms_active
caltable='initialdelay.k'
field='2'            # Hard-coded to field 2 as this is always 3C286
spw=''
intent=''
selectdata=False
uvrange=''           # No uvrange for delay calibration, since it exclude antennas.
solint='inf'
combine='scan'
preavg=-1.0
refant=refAnt
minblperant=minBL_for_cal
minsnr=3.0
solnorm=False
gaintype='K'
smodel=[]
calmode='p'
append=False
gaintable=GainTables
gainfield=['']
interp=['']
spwmap=[]
opacity=[]
parang=False
async=False
gaincal()

#specify spw14 to apply delays from there

spwac=14
spwmapdelayarray=pl.zeros(15,int)
spwmapdelayarray[0:15]=spwac
spwmapdelay=list(spwmapdelayarray)

GainTables=copy.copy(priorcals)
GainTables.append('initialdelay.k')
SpwMapValues=copy.copy(priorspwmap)
SpwMapValues.append(spwmapdelay)


#5: This second gaincal applies the delays of spw 14? TABLE:initialBPinitialgain.g
default('gaincal')
vis=ms_active
caltable='initialBPinitialgain.g'
field='2'            # Hard-coded to field 2 as this is always 3C286
spw=tst_bpass_spw
selectdata=False
solint='int'
combine='scan'
preavg=-1.0
refant=refAnt
uvrange=uvr_cal      # Set uvrange to exclude worst of RFI
minblperant=minBL_for_cal
minsnr=3.0
solnorm=False
gaintype='G'
smodel=[]
calmode='p'
append=False
gaintable=GainTables
gainfield=['']
interp=['']
spwmap=SpwMapValues
opacity=[]
parang=False
async=False
gaincal()

BPGainTables=copy.copy(GainTables)            # HG (12/04/2015) : simplifying steps for gaintables and spwmaps
BPGainTables.append('initialBPinitialgain.g')
BPSpwMapValues=copy.copy(SpwMapValues)
BPSpwMapValues.append([])

#6: Bandpass solution: TABLE: initialBPcal.b
default('bandpass')
vis=ms_active
caltable='initialBPcal.b'
field='2'            # Hard-coded to field 2 as this is always 3C286
spw='0~14:64~1982'
selectdata=True
uvrange=uvr_cal      # Set uvrange to exclude worst of RFI
solint='inf'
combine='scan'
refant=refAnt
minblperant=minBL_for_cal
minsnr=5.0
solnorm=False       # Must be set to False otherwise solution doesn't apply BP to flagged channels
bandtype='B'
fillgaps=10          # Interpolation appears to be occurring regardless so set this to 10 channels
interp=['linear,spline']          # Use linear interpolation in time, spline in frequency.
append=False
gaintable=BPGainTables
gainfield=['']
spwmap=BPSpwMapValues
gaincurve=False
opacity=[]
parang=False
async=False
bandpass()


AllCalTables=copy.copy(GainTables)           # HG (12/04/2015) : simplifying steps for gaintables and spwmaps
AllCalTables.append('initialBPcal.b')
AllCalTables.append('initialBPinitialgain.g')
AllSpwMapValues=copy.copy(SpwMapValues)
AllSpwMapValues.append([])
AllSpwMapValues.append([])

#7: #Apply calibrations to the flux calibrator:
default('applycal')
vis=ms_active
field='2'            # Hard-coded to field 2 as this is always 3C286
spw=''
intent=''
selectdata=True
gaintable=AllCalTables
interp=['','','',',spline','']         # Default is to interpolate linearly in time & frequency; interpolating frequency with spline for BP.
spwmap=AllSpwMapValues
gaincurve=False
opacity=[]
parang=False
calwt=False
flagbackup=False
async=False
applycal()

# Summary of flagging, after initial calib (for testing purposes only)
logprint ("Summary of flags after initial BP cal", logfileout='logs/bandpass.log')
default('flagdata')
vis=ms_active
mode='summary'
spw='0~14'
correlation='RR,LL'
spwchan=True
spwcorr=True
action='calculate'
flagdata()


#8: First rflag run 
logprint ("RFLAG+Extend on flux calibrator", logfileout='logs/bandpass.log')

f=ms.msseltoindex(vis=ms_active,field='1331*')['field']
ff= float(f[0])


default('flagdata')
vis=ms_active
mode='rflag'
field='2'            # Hard-coded to field 2 as this is always 3C286
spw='0~14'
correlation=''
ntime='scan'
combinescans=False
datacolumn='corrected'
extendflags=False    # Explicitly set to False since default is True.  Extending on next step. 
extendpols=False     # Default is True.  May allow some weak RFI through, but try it.   
winsize=3
freqdev=[[ff,0.0,7.1],[ff,1.0,5.5],[ff,2.0,4.5],[ff,3.0,4.1],[ff,4.0,3.9],[ff,5.0,3.7],[ff,6.0,3.6],[ff,7.0,3.7],[ff,8.0,3.3],[ff,9.0,3.0],[ff,10.0,3.0],[ff,11.0,3.0],[ff,12.0,2.9],[ff,13.0,3.0],[ff,14.0,3.0]]
timedev=[[ff,0.0,9.4],[ff,1.0,7.3],[ff,2.0,5.9],[ff,3.0,5.4],[ff,4.0,5.2],[ff,5.0,4.9],[ff,6.0,4.8],[ff,7.0,4.9],[ff,8.0,4.4],[ff,9.0,4.0],[ff,10.0,3.9],[ff,11.0,3.9],[ff,12.0,3.9],[ff,13.0,4.0],[ff,14.0,4.0]]
timedevscale=1.0
freqdevscale=1.0
action='apply'
display=''
flagbackup=False
savepars=True
async=False
flagdata()

# Summary of flagging, after initial RFLAG (for testing purposes only)
logprint ("Summary of flags after initial RFLAG", logfileout='logs/bandpass.log')
default('flagdata')
vis=ms_active
mode='summary'
spw='0~14'
correlation='RR,LL'
spwchan=True
spwcorr=True
action='calculate'
flagdata()

#9:
default('flagdata')
vis=ms_active
mode='extend'
field='2'            # Hard-coded to field 2 as this is always 3C286
correlation=''
ntime='scan'
combinescans=False
datacolumn='corrected'
extendpols=False      
growtime=70       
growfreq=80       
growaround=True      
flagneartime=True       
flagnearfreq=True       
action='apply'
display=''
flagbackup=False   # Will backup final flagging manually.
savepars=True
async=False
flagdata()

# Summary of flagging, after RFLAG+extend (for testing purposes only)
logprint ("Summary of flags after RFLAG+extend", logfileout='logs/bandpass.log')
default('flagdata')
vis=ms_active
mode='summary'
spw='0~14'
correlation='RR,LL'
spwchan=True
spwcorr=True
action='calculate'
flagdata()

# Save flags
logprint ("Saving flags", logfileout='logs/bandpass.log')


default('flagmanager')
vis=ms_active
mode='save'
versionname='fluxcal_flags'
comment='Flux calibrator flags saved after application'
merge='replace'
async=False
flagmanager()
logprint ("Flag column saved to "+versionname, logfileout='logs/bandpass.log')

#--------------------------------------------------------------------
#Part IV: second round of bandpass calibration
logprint ("Final delay, gain, BP calibration", logfileout='logs/bandpass.log')

#10:Gain calibration on delay calibrator, TABLE: finaldelayinitialgain.g
default('gaincal')
vis=ms_active
caltable='finaldelayinitialgain.g'
field='2'            # Hard-coded to field 2 as this is always 3C286
spw=tst_delay_spw
intent=''
selectdata=False
solint='int'
combine='scan'
preavg=-1.0
refant=refAnt
uvrange=uvr_cal      # Set uvrange to exclude worst of RFI
minblperant= minBL_for_cal
minsnr=3.0
solnorm=False
gaintype='G'
smodel=[]
calmode='p'
append=False
docallib=False
gaintable=priorcals
gainfield=['']
interp=['']
spwmap=[]
parang=False
async=False
gaincal()

DelayTables=copy.copy(priorcals)
DelayTables.append('finaldelayinitialgain.g')

#11: Delay calibration to all spws. TABLE: finaldelay.k
default('gaincal')
vis=ms_active
caltable='finaldelay.k'
field='2'            # Hard-coded to field 2 as this is always 3C286
spw=''
intent=''
selectdata=False
uvrange=''           # No uvrange for delay calibration, since it exclude antennas.
solint='inf'
combine='scan'
preavg=-1.0
refant=refAnt
minblperant=minBL_for_cal
minsnr=3.0
solnorm=False
gaintype='K'
smodel=[]
calmode='p'
append=False
gaintable=DelayTables
gainfield=['']
interp=['']
spwmap=[]
opacity=[]
parang=False
async=False
gaincal()

#specify spw14 to apply delays from there

spwac=14
spwmapdelayarray=pl.zeros(15,int)
spwmapdelayarray[0:15]=spwac
spwmapdelay=list(spwmapdelayarray)

GainTables=copy.copy(priorcals)
GainTables.append('finaldelay.k')
SpwMapValues=copy.copy(priorspwmap)
SpwMapValues.append(spwmapdelay)


#12:This second gaincal applies the delays of spw 14. TABLE:finalBPinitialgain.g
default('gaincal')
vis=ms_active
caltable='finalBPinitialgain.g'
field='2'            # Hard-coded to field 2 as this is always 3C286
spw=tst_bpass_spw
selectdata=False
solint='int'
combine='scan'
preavg=-1.0
refant=refAnt
uvrange=uvr_cal      # Set uvrange to exclude worst of RFI
minblperant=minBL_for_cal
minsnr=3.0
solnorm=False
gaintype='G'
smodel=[]
calmode='p'
append=False
gaintable=GainTables
gainfield=['']
interp=['']
spwmap=SpwMapValues
opacity=[]
parang=False
async=False
gaincal()

BPGainTables=copy.copy(GainTables)          # HG (12/04/2015) : simplifying steps for caltables and spwmap
BPGainTables.append('finalBPinitialgain.g')
BPSpwMapValues=copy.copy(SpwMapValues)
BPSpwMapValues.append([])

#13: Bandpass solution: TABLE: finalBPcal.b
default('bandpass')
vis=ms_active
caltable='finalBPcal.b'
field='2'            # Hard-coded to field 2 as this is always 3C286
spw='0~14:64~1982'
selectdata=True
uvrange=uvr_cal      # Set uvrange to exclude worst of RFI
solint='inf'
combine='scan'
refant=refAnt
minblperant=minBL_for_cal
minsnr=5.0
solnorm=False       # Must be set to False otherwise solution doesn't apply BP to flagged channels
bandtype='B'
fillgaps=10          # Interpolation appears to be occurring regardless so set this to 10 channels
interp=['linear,spline']          # Use linear interpolation in time, spline in frequency.
append=False
gaintable=BPGainTables
gainfield=['']
spwmap=BPSpwMapValues
gaincurve=False
opacity=[]
parang=False
async=False
bandpass()


AllCalTables=copy.copy(GainTables)
AllCalTables.append('finalBPcal.b')
AllCalTables.append('finalBPinitialgain.g')
AllSpwMapValues=copy.copy(SpwMapValues)
AllSpwMapValues.append([])
AllSpwMapValues.append([])

#15:
#Apply calibrations to the flux calibrator:
default('applycal')
vis=ms_active
field='2'            # Hard-coded to field 2 as this is always 3C286
spw=''
intent=''
selectdata=True
gaintable=AllCalTables
interp=['','','',',spline','']         # Default is to interpolate linearly in time & frequency; interpolating frequency with spline for BP.
spwmap=AllSpwMapValues
gaincurve=False
opacity=[]
parang=False
calwt=False
flagbackup=False
async=False
applycal()

# Summary of flagging, after final BP calib (for testing purposes only)
logprint ("Summary of flags after final BP cal", logfileout='logs/bandpass.log')
default('flagdata')
vis=ms_active
mode='summary'
spw='0~14'
correlation='RR,LL'
spwchan=True
spwcorr=True
action='calculate'
flagdata()

#--------------------------------------------------------------------
#Part V: Data inspection
logprint ("Create data inspection plots", logfileout='logs/bandpass.log')


#16: plot delays
nplots=int(numAntenna/3)

if ((numAntenna%3)>0):
    nplots = nplots + 1

for ii in range(nplots):
    filename='finaldelay'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)
    
    antPlot=str(ii*3)+'~'+str(ii*3+2)
    
    default('plotcal')
    caltable='finaldelay.k'
    xaxis='freq'
    yaxis='delay'
    poln=''
    field=''
    antenna=antPlot
    spw=''
    timerange=''
    subplot=311
    overplot=False
    clearpanel='Auto'
    iteration='antenna'
    plotrange=[]
    showflags=False
    plotsymbol='o'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    async=False
    plotcal()

#17: Bandpass solutions
default('plotbandpass')
caltable='finalBPcal.b'
yaxis='both'
xaxis='freq'
figfile='bandpass'
buildpdf=True
convert='/usr/bin/convert'
interactive=False
subplot='42'
plotbandpass()



ms_name=ms_active[:-3]
output_ms=ms_name+'_flux_averaged.ms'

# Remove old averaged MS if it exists.
if os.path.exists(output_ms):
    os.system("rm -rf "+output_ms)
    
#18: average the data
default('split')
vis=ms_active
outputvis=output_ms
datacolumn='corrected'
field='2'
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

# Use plotms to make plot of amplitude vs. phase, divided by spw
seq=range(0,15)
for ii in seq:
    default('plotms')
    vis=ms_active
    field='2'       # Hard coded to 3C286, flux/bandpass calibrator
    xaxis='phase'
    yaxis='amp'
    xdatacolumn='corrected'
    ydatacolumn='corrected'
    averagedata=True
    #avgtime='1e5'
    #avgscan=True
    avgchannel='2048'    # Number of channels in all spws for CHILES
    spw=str(ii)
    avgspw=False
    #gridcolumns=
    showlegend=False
    iteraxis='spw'
    coloraxis='corr'
    showgui=False
    clearplots=True
    plotfile='fluxcal_ampphase.png'
    plotms()

seq=range(0,15)
#19: 
for ii in seq:
    print 'STARTS IMAGING FLUX CALIBRATOR OF SPW='+str(ii)
    default('clean')
    image_name='fluxcalibrator_spw'+str(ii)
    fieldid='1331*'
    grid_mode=''
    number_w=1
    image_size=[512,512]
    iteration=1000
#    mask_name=['fluxcal_mask.txt']
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

   
#20:  Calculate bmaj, bmin, peak, and rms for images of flux calibrator in each spw.
#     Output results as text file and plots of these values.
box_flux='300,30,460,200'

bmaj_flux=[]
bmin_flux=[]
max_flux=[]
rms_flux=[]

for ii in seq:
    
    image_fluxcal='fluxcalibrator_spw'+str(ii)+'.image'
    bmaj1=imhead(image_fluxcal,mode='get',hdkey='beammajor')
    if bmaj1==None:
        bmaj_flux.append(0.0)
        bmin_flux.append(0.0)
        max_flux.append(0.0)
        rms_flux.append(0.0)
    else:
        bmaj11=bmaj1['value']
        bmaj_flux.append(bmaj11)
        bmin1=imhead(image_fluxcal,mode='get',hdkey='beamminor')
        bmin11=bmin1['value']
        bmin_flux.append(bmin11)
        imst1=imstat(image_fluxcal)
        max1=imst1['max']
        max_flux.append(max1[0])
        imst1=imstat(image_fluxcal,box=box_flux)
        rms1=imst1['rms']
        rms_flux.append(rms1[0]*1e3)

# Output Image statistics in text file
f=open('statisticsFlux.txt','w')

print >> f, "Flux calibrator"
print >> f, "major axis [\"] \t", "".join(["%12.3f \t" %x for x in bmaj_flux])
print >> f, "minor axis [\"] \t", "".join(["%12.3f \t" %x for x in bmin_flux])
print >> f, "peak value [Jy] \t", "".join(["%12.3f \t" %x for x in max_flux])
print >> f, "RMS noise [mJy] \t", "".join(["%12.3f \t" %x for x in rms_flux])
f.close()

#Make plots of bmaj, bmin, peak, and rms
fig=pylab.figure()
pylab.plot(seq,bmaj_flux,'r--x',label='Bmaj')
pylab.plot(seq,bmin_flux,'b--x',label='Bmin')
pylab.xlabel('Spectral Window')
pylab.ylabel('Beam Size ["]')
pylab.legend()
pylab.savefig('fluxcal_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,max_flux,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('Peak Flux [Jy]')
pylab.yscale('log')
pylab.savefig('fluxcal_peak.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,rms_flux,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('RMS [mJy]')
pylab.yscale('log')
pylab.savefig('fluxcal_rms.png')
pylab.close(fig)

#Want to plot image of flux calibrator in each spw.  Use "imview"

for ii in seq:
    image_fluxcal='fluxcalibrator_spw'+str(ii)+'.image'
    kntr_levels=[-2*rms_flux[ii]/1000.,2*rms_flux[ii]/1000.,0.1*max_flux[ii],0.3*max_flux[ii],0.5*max_flux[ii],0.7*max_flux[ii],0.9*max_flux[ii]]
    imview(raster={'file':image_fluxcal,'scaling':-4, 'colorwedge':True},contour={'file':image_fluxcal,'levels':kntr_levels},out='fluxcal_spw'+str(ii)+'.png')


#Move plots, images to sub-directory

os.system("mv *.png plots")
os.system("mv bandpass.pdf plots")
if os.path.exists('images')==False:
    os.system("mkdir images")
os.system("mv fluxcalibrator_spw*.* images")

#Create webpage with results

if os.path.exists('bandpass.html'):
    os.system("rm bandpass.html")
wlog = open("bandpass.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_bandpass results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/bandpass.log">Bandpass Log</a></li>\n')
wlog.write('<li>Delay Solutions: \n')
wlog.write('<br><img src="plots/finaldelay0.png">\n')
wlog.write('<br><img src="plots/finaldelay1.png">\n')
wlog.write('<br><img src="plots/finaldelay2.png">\n')
wlog.write('<br><img src="plots/finaldelay3.png">\n')
wlog.write('<br><img src="plots/finaldelay4.png">\n')
wlog.write('<br><img src="plots/finaldelay5.png">\n')
wlog.write('<br><img src="plots/finaldelay6.png">\n')
wlog.write('<br><img src="plots/finaldelay7.png">\n')
wlog.write('<br><img src="plots/finaldelay8.png"></li>\n')
wlog.write('<li> Bandpass solutions: \n')
#wlog.write('<br><embed src="plots/bandpass.pdf"></li>\n')
wlog.write('<br><object data="plots/bandpass.pdf" type="application/pdf" width="100%" height="100%"><p>Download <a href="plots/bandpass.pdf">bandpass.pdf</a></p></object></li>\n')
wlog.write('<li> Amp vs. Phase (averaged over all channels): \n')
for ii in seq:
    wlog.write('<br><img src="plots/fluxcal_ampphase_Spw'+str(ii)+'.png">\n')
#    wlog.write('<br><img src="plots/fluxcal_phaseuvdist_Spw'+str(ii)+'.png">\n')
wlog.write('<li> Images of Flux Calibrator: \n')
for ii in seq:
    wlog.write('<br><img src="plots/fluxcal_spw'+str(ii)+'.png">\n')
wlog.write('</li>')
wlog.write('<li> Measured properties of flux calibrator: \n')
wlog.write('<br><img src="plots/fluxcal_beamsize.png">\n')
wlog.write('<br><img src="plots/fluxcal_peak.png">\n')
wlog.write('<br><img src="plots/fluxcal_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()



logprint ("Finished CHILES_pipe_bandpass.py", logfileout='logs/bandpass.log')
time_list=runtiming('bandpass', 'end')

pipeline_save()
