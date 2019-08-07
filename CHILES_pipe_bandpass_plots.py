# CHILES_pipe_bandpass_plots.py
# This module of the CHILES pipeline remakes the diagnostic plots
# 9/23/16 DJP

#Need to run CHILES_pipe_restore.py first
#8/31/15 DJP

# 9/6/17 DJP, backs up old plots/html file before making new ones.
# 8/29/18 DJP: Changed field='2' to field='1331+305=3C286'.  

#Part I: define some variables that will be used later
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

logprint ("Starting CHILES_pipe_bandpass_plots.py", logfileout='logs/bandpass.log')
time_list=runtiming('bandpass', 'start')

#Put old plots & html file in Backup directory
if os.path.exists('BACKUP')!=True:
    os.system('mkdir BACKUP')
os.system('mv bandpass.html BACKUP/.')
os.system('mv plots/fluxcal*.png BACKUP/.')
os.system('mv plots/bpamp*.png BACKUP/.')
os.system('mv plots/bpphase*.png BACKUP/.')
os.system('mv plots/finaldelay*.png BACKUP/.')
# Remove old images of flux calibrator
os.system("rm -rf images/fluxcalibrator_spw*.*")

priorcals=['gain_curves.g']
priorspwmap=[[]]

if os.path.exists('antposcal.p')==True:  # HG (12/04/2015) : include antposcal.p in the priorcals if it exists
  priorcals.append('antposcal.p')
  priorspwmap.append([])

GainTables=copy.copy(priorcals)
GainTables.append('finaldelay.k')
SpwMapValues=copy.copy(priorspwmap)
SpwMapValues.append(spwmapdelay)

AllCalTables=copy.copy(GainTables)
AllCalTables.append('finalBPcal.b')
AllCalTables.append('finalBPinitialgain.g')
AllSpwMapValues=copy.copy(SpwMapValues)
AllSpwMapValues.append([])
AllSpwMapValues.append([])


#--------------------------------------------------------------------
#Part V: Data inspection
logprint ("Create data inspection plots", logfileout='logs/bandpass.log')

# Make plot of flagging statistics
# s_b is the output of flagdata run (above)
# Get information for flagging percentage vs. uvdistance
#gantdata = get_antenna_data(ms_active)
#create adictionary with flagging info
#base_dict = create_baseline_dict(ms_active, gantdata)
#gantdata and base_dict are already in the initial module so no need to retrieve that information again.
#match flagging data to dictionary entry
datamatch = flag_match_baseline(s_b['baseline'], base_dict)
#bin the statistics
binned_stats = bin_statistics(datamatch, 'B', 25)  # 25 is the number of uvdist bins such that there is minimal error in uvdist.

#Plot flagging % vs. uvdist
### Plot the Data
barwidth = binned_stats[0][1]
totflagged = 'BP Flagging: '+ str(flux_flag*100) + '% Data Flagged'
pylab.close()
pylab.bar(binned_stats[0],binned_stats[1], width=barwidth, color='grey', align='edge')
pylab.title(totflagged)
pylab.grid()
pylab.ylabel('flagged data [%]')
pylab.xlabel('average UV distance [m]')
pylab.savefig('bp_flag_uvdist.png')
os.system("mv bp_flag_uvdist.png plots/.") 

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
for s in range(15):
    for c in range(2048):
        ct+=1
        chan.append(ct)
        freq.append(nu0+dnu*ct)
        flagged.append(s_b['spw:channel'][str(s)+':'+str(c)]['flagged'])
        totals.append(s_b['spw:channel'][str(s)+':'+str(c)]['total'])
        flag_frac.append(flagged[ct]/totals[ct])

fig=pylab.figure()
pylab.plot(freq,flag_frac,'k-')
pylab.xlim(940.,1445.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Fraction of data flagged')
pylab.savefig("bp_flag.png")
pylab.close(fig)

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
    plotcal()

#17: Bandpass solutions
seq=range(0,15)
for ii in seq:
    default('plotms')
    vis='finalBPcal.b'
    field='1331+305=3C286'        # Only for 3C286
    xaxis='freq'
    yaxis='amp'
    xdatacolumn=''
    ydatacolumn=''
    averagedata=False
    spw=str(ii)
    iteraxis='spw'
    showlegend=False
    coloraxis='Antenna1'
    title='Spw'+str(ii)+': Amp v. Freq'
    showgui=False
    clearplots=True
    plotfile='bpamp.png'
    plotms()
    
    yaxis='phase'
    title='Spw'+str(ii)+': Phase v. Freq'
    plotfile='bpphase.png'
    plotms()

#Plot UV spectrum (averaged over all baselines & time) of flux calibrator
default('plotms')
vis=ms_active   
field='1331+305=3C286'        # Only for 3C286
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
plotfile='fluxcal_spectrum_full.png'
plotms()

plotrange=[0.95,1.43,14.5,18.5]
plotfile='fluxcal_spectrum_zoom.png'
plotms()

# Plot calibrated phase vs. frequency
xaxis='freq'
yaxis='phase'
plotrange=[0.95,1.43,0,0]
clearplots=True
plotfile='fluxcal_phasespectrum.png'
plotms()


ms_name=ms_active[:-3]
output_ms=ms_name+'_flux_averaged.ms'

# Remove old averaged MS if it exists.
if os.path.exists(output_ms):
    os.system("rm -rf "+output_ms)
    
#18: average the data
default('oldsplit')
vis=ms_active
outputvis=output_ms
datacolumn='corrected'
field='1331+305=3C286'
spw='0~14:128~1920'    # Extend the region used for imaging/plots, only excluding the edges.
width='1793'
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
oldsplit()

# Use plotms to make plot of amplitude vs. phase, divided by spw
seq=range(0,15)
for ii in seq:
    default('plotms')
    vis=output_ms   # Take averaged 3C286 data for diagnostic plots
    field=''        # Only 3C286 in this MS.
    xaxis='phase'
    yaxis='amp'
    xdatacolumn='corrected'
    ydatacolumn='corrected'
    averagedata=False     
    spw=str(ii)
    avgspw=False
    showlegend=False
    iteraxis='spw'
    coloraxis='corr'
    showgui=False
    clearplots=True
    plotfile='fluxcal_ampphase.png'
    plotms()

# Use plotms to make plot of amplitude and phase vs. time, for each spw.  
for ii in seq:
    default('plotms')
    vis=output_ms  # File only contains 3C286, flux calibrator data
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
    plotfile='fluxcal_amptime.png'
    plotms()
    yaxis='phase'
    plotfile='fluxcal_phasetime.png'
    plotms()

seq=range(0,15)
#19: 
for ii in seq:
    #19: Imaging flux calibrator
    print 'STARTS IMAGING FLUX CALIBRATOR OF SPW='+str(ii)
    default('tclean')
    image_name='fluxcalibrator_spw'+str(ii)
    fieldid='1331*'
    grid_mode=''
    number_w=1
    image_size=[512,512]
    iteration=1000
    mask_name=['']
    
    vis=output_ms
    imagename=image_name
    selectdata=True
    datacolumn='data'
    field=fieldid
    spw=str(ii)
    specmode='mfs'
    nterms=1
    niter=iteration
    gain=0.1
    gridmode=grid_mode
    wprojplanes=number_w
    threshold='0.0mJy'
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
#pylab.yscale('log')
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
#Dropping this from diagnostic plots

#for ii in seq:
#    image_fluxcal='fluxcalibrator_spw'+str(ii)+'.image'
#    kntr_levels=[-2*rms_flux[ii]/1000.,2*rms_flux[ii]/1000.,0.1*max_flux[ii],0.3*max_flux[ii],0.5*max_flux[ii],0.7*max_flux[ii],0.9*max_flux[ii]]
#    imview(raster={'file':image_fluxcal,'scaling':-4, 'colorwedge':True},contour={'file':image_fluxcal,'levels':kntr_levels},out='fluxcal_spw'+str(ii)+'.png')


#Move plots, images to sub-directory

os.system("mv *.png plots")
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
wlog.write('<style>table tr {page-break-inside: avoid}</style>\n')
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
wlog.write('<li> Bandpass solutions (amplitude and phase) for reference antenna: \n')
wlog.write('<li> Color coded by antenna, both polarizations shown \n')
wlog.write('<table> \n')
for ii in seq:
    wlog.write('<tr><td><img src="plots/bpamp_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('<td><img src="plots/bpphase_Spw'+str(ii)+'.png"></td></tr>\n')
wlog.write('</table> \n')
wlog.write('<br>')
wlog.write('<li> Amp vs. Phase (averaged over all channels in a spw): \n')
for ii in seq:
    wlog.write('<br><img src="plots/fluxcal_ampphase_Spw'+str(ii)+'.png">\n')
wlog.write('<li> Spectrum of Flux calibrator (both LL & RR, averaged over all time & baselines): \n')
wlog.write('<br><img src="plots/fluxcal_spectrum_full.png">\n')
wlog.write('<br><img src="plots/fluxcal_spectrum_zoom.png">\n')
wlog.write('<br><img src="plots/fluxcal_phasespectrum.png">\n')
wlog.write('<li> Amp. & Phase vs. time for Flux Calibrator (averaged over frequency): \n')
wlog.write('<table> \n')
for ii in seq:
    wlog.write('<tr><td><img src="plots/fluxcal_amptime_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('<td><img src="plots/fluxcal_phasetime_Spw'+str(ii)+'.png"></td></tr>\n')
wlog.write('</table> \n')
wlog.write('</li>')
wlog.write('<li> Measured properties of flux calibrator: \n')
wlog.write('<br><img src="plots/fluxcal_beamsize.png">\n')
wlog.write('<br><img src="plots/fluxcal_peak.png">\n')
wlog.write('<br><img src="plots/fluxcal_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Flagging percentage vs. frequency:\n')
wlog.write('<li><br><img src="./plots/bp_flag.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Flagging percentage vs. uvdist\n')
wlog.write('<li><br><img src="./plots/bp_flag_uvdist.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Percentage of 3C286 data flagged: '+str(flux_flag*100)+'\n')
wlog.write('<br>')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()

logprint ("Finished CHILES_pipe_bandpass.py", logfileout='logs/bandpass.log')
time_list=runtiming('bandpass', 'end')

pipeline_save()
