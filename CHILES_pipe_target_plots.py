# CHILES_pipe_target_plots.py
# This task runs after the "phasecal" task so all calibrations have been applied
# to the target field.  This task runs RFLAG and extend on the target field and
# makes some diagnostic plots.  
# 9/9/15 DJP
# 4/22/18 DJP: Changing flagging and split to oldsplit
# 8/29/18 DJP:  Updated references from field='1' to field='deepfield'
# 05/15/19 DJP:  Removed images from plots.



logprint ("Starting CHILES_pipe_target_plots.py", logfileout='logs/target.log')
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

# If target module crashed after splits:
if os.path.exists('target_flux_averaged.ms'):
    os.system("rm -rf target_flux_averaged.ms'))

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
field='deepfield'        # Only for deepfield
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
    field='deepfield'       # Only for deepfield
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

default('oldsplit')
vis=ms_active
outputvis=output_ms
datacolumn='corrected'
field='deepfield'
spw='0~14:128~1920'  # Extend averaging to include all but 128 edge channels
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

seq=range(0,15)
#Image target: 
for ii in seq:
    print 'STARTS IMAGING Deepfield OF SPW='+str(ii)
    default('tclean')
    image_name='target_spw'+str(ii)
    fieldid='deepfield'
    grid_mode=''
    number_w=1
    image_size=[2048,2048]
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

# for ii in seq:
#     image_target='target_spw'+str(ii)+'.image'
#     kntr_levels=[-2*rms_target[ii]/1000.,2*rms_target[ii]/1000.,0.1*max_target[ii],0.3*max_target[ii],0.5*max_target[ii],0.7*max_target[ii],0.9*max_target[ii]]
#     imview(raster={'file':image_target, 'colorwedge':True},contour={'file':image_target,'levels':kntr_levels},out='target_spw'+str(ii)+'.png')


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
wlog.write('<style>table tr {page-break-inside: avoid}</style>\n')
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
# wlog.write('<li> Images of Deepfield: \n')
# for ii in seq:
#     wlog.write('<br><img src="plots/target_spw'+str(ii)+'.png">\n')
# wlog.write('</li>')
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


logprint ("Finished CHILES_pipe_target_plots.py", logfileout='logs/target.log')
time_list=runtiming('target', 'end')

pipeline_save()
