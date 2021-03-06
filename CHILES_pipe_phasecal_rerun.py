# CHILES_pipe_phasecal_rerun.py
# This module of the CHILES pipeline reruns the phasecal module
# 9/21/16 DJP
# 8/29/18 DJP: Changing field='0' to field='J0943-0819' (and similarly for other fields)
# 10/10/18 DJP: Updated plots
# 05/15/19 DJP:  No longer doing phase calibration on flux calibrator, make UVSPEC of phase v. frequency

logprint ("Starting CHILES_pipe_phasecal_rerun.py", logfileout='logs/phasecal.log')
time_list=runtiming('phase', 'start')

# Step 1:  Set default parameters (to match previous modules), and clean up old files:
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

#Put old plots & html file in Backup directory
if os.path.exists('BACKUP')!=True:
    os.system('mkdir BACKUP')
os.system('mv plots/phase*.png BACKUP/.')
os.system('mv phasecal.html BACKUP/.')

#Remove old files from previous attempts using rmtables
#Clear model of phase calibrator if module previously run.
if os.path.exists('initialphase.gcal'):
    rmtables(tablenames='initialphase.gcal')
if os.path.exists('initialamp.gcal'):
    rmtables(tablenames='initialamp.gcal')
if os.path.exists('initialflux.gcal'):
    rmtables(tablenames='initialflux.gcal')
if os.path.exists('finalphase_int.gcal'):
    rmtables(tablenames='finalphase_int.gcal')
if os.path.exists('finalphase_scan.gcal'):
    rmtables(tablenames='finalphase_scan.gcal')
if os.path.exists('finalamp.gcal'):
    rmtables(tablenames='finalamp.gcal')
if os.path.exists('finalflux.gcal'):
    rmtables(tablenames='finalflux.gcal')
    # If caltables exist, then delmod must be run.
    default('delmod')
    vis=ms_active
    otf=True
    field='J0943-0819'  # Hard-coded for phase calibrator
    scr=False
    delmod()
os.system("rm -rf images/phasecalibrator_spw*.*")


# The following parameters are set in the initial pipeline script by default or user input.
##Set minimum # of baselines need for a solution to 8, based on experience
minBL_for_cal=8

##Set uvrange to apply in order to optimally exclude RFI without flagging:
uvr_cal='>1500m'

#Set prior cals including gain, delay, and BP calibration
# HG : adding the spwmap for prior cals. & adding the conditional statements for antposcal.p

spwac=14
spwmapdelayarray=pl.zeros(15,int)
spwmapdelayarray[0:15]=spwac
spwmapdelay=list(spwmapdelayarray)

# HG: if no antenna position changes, then table not generated.
if os.path.exists('antposcal.p')==True:
  priorcals=['antposcal.p','gain_curves.g','finaldelay.k','finalBPcal.b']
  priorspwmap=[[],[],spwmapdelay,[]]

if os.path.exists('antposcal.p')==False:
  priorcals=['gain_curves.g','finaldelay.k','finalBPcal.b']
  priorspwmap=[[],spwmapdelay,[]]


#Set channel/spw range to use for calibrations (based on BP code)
tst_gain_spw=tst_bpass_spw

# This module is hard-coded for CHILES, so:
# field='J0943-0819' is J0943-0819, phase calibrator
# Field='2' is 1331+305 (3C286), flux & bandpass calibrator
# Field='1' is "deepfield", the target


  
# Step 2:  gaincal on both calibrators
# Need short solint = "int" (for calibrators)
# Long solint = "inf" (to apply to target) not needed until final iteration.

# Step 7: repeat gaincal 
# Short solint = "int" (length of individual calibrator scan)

logprint ("Running final gaincal", logfileout='logs/phasecal.log')

default('gaincal')
vis=ms_active
caltable='finalphase_int.gcal'
field='J0943-0819,1331+305=3C286'
spw=tst_gain_spw
intent=''
selectdata=False
solint='int'
combine=''
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
spwmap=priorspwmap # In previous version this was [].  These spw's correspond to priorcals.
parang=False
async=False
gaincal()

default('gaincal')
vis=ms_active
caltable='finalphase_scan.gcal'
field='J0943-0819,1331+305=3C286'
spw=tst_gain_spw
intent=''
selectdata=False
solint='inf'
combine=''
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
spwmap=priorspwmap  # In previous version this was [].  These spw's correspond to priorcals.
parang=False
async=False
gaincal()

logprint ("Deriving final amp & phase solutions", logfileout='logs/phasecal.log')

GainTables=copy.copy(priorcals)
GainTables.append('finalphase_int.gcal')
SpwMapValues=copy.copy(priorspwmap)
SpwMapValues.append([])

default('gaincal')
vis=ms_active
caltable='finalamp.gcal'
field='J0943-0819,1331+305=3C286'
spw=tst_gain_spw
intent=''
selectdata=False
solint='inf'
combine=''
preavg=-1.0
refant=refAnt
uvrange=uvr_cal      # Set uvrange to exclude worst of RFI
minblperant= minBL_for_cal
minsnr=3.0
solnorm=False
gaintype='G'
smodel=[]
calmode='ap'         # Solve for both amp & phase
append=False
docallib=False
gaintable=GainTables
gainfield=['']
interp=['']
spwmap=SpwMapValues  # In previous version this was [].  These spw's correspond to GainTables.
parang=False
async=False
gaincal()

# Step 8: repeat bootstrapping
os.system('rm -rf '+fluxscale_output)
logprint ("Final flux bootstrapping", logfileout='logs/phasecal.log')
logprint ("Flux densities will be written to "+fluxscale_output, logfileout='logs/phasecal.log')

#Clear previous models of phase calibrator
default('delmod')
vis=ms_active
otf=True
field='J0943-0819'  # Hard-coded for phase calibrator
scr=False
delmod()

casalog.setlogfile(fluxscale_output)

default('fluxscale')
vis=ms_active
caltable='finalamp.gcal'
fluxtable='finalflux.gcal'
reference='1331+305=3C286'
transfer=['J0943-0819']
listfile=''
append=False
refspwmap=[-1]        
incremental=True
fitorder=1

try:
    fluxscale_result=fluxscale()

    casalog.setlogfile(maincasalog)


    logprint ("Fitting data with power law", logfileout='logs/phasecal.log')


#
# the variable center_frequencies should already have been filled out
# with the reference frequencies of the spectral window table
#

    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    try:
        ff = open(fluxscale_output, 'r')
    except IOError as err:
        logprint (fluxscale_output+" doesn't exist, error: "+err.filename, logfileout='logs/phasecal.log')

# looking for lines like:
#2012-03-09 21:30:23     INFO    fluxscale::::    Flux density for J1717-3342 in SpW=3 is: 1.94158 +/- 0.0123058 (SNR = 157.777, N= 34)
# sometimes they look like:
#2012-03-09 21:30:23     INFO    fluxscale::::    Flux density for J1717-3342 in SpW=0 is:  INSUFFICIENT DATA 
# so watch for that.

    sources = []
    flux_densities = []
    spws = []

#Find the field_ids in the dictionary returned from the CASA task fluxscale
    dictkeys = fluxscale_result.keys()
    keys_to_remove = ['freq', 'spwName', 'spwID']
    dictkeys = [field_id for field_id in dictkeys if field_id not in keys_to_remove]

    logprint("Finding field_ids returned from fluxscale.", logfileout='logs/phasecal.log')

    for field_id in dictkeys:
        sourcename = fluxscale_result[field_id]['fieldName']
        secondary_keys = fluxscale_result[field_id].keys()
        secondary_keys_to_remove=['fitRefFreq', 'spidxerr', 'spidx', 'fitFluxd', 'fieldName', 'fitFluxdErr']
        spwkeys = [spw_id for spw_id in secondary_keys if spw_id not in secondary_keys_to_remove]

        for spw_id in spwkeys:
            flux_d = list(fluxscale_result[field_id][spw_id]['fluxd'])
            flux_d_err = list(fluxscale_result[field_id][spw_id]['fluxdErr'])
            #spwslist  = list(int(spw_id))

            #flux_d = list(fluxscale_result[field_id]['fluxd'])
            #flux_d_err = list(fluxscale_result[field_id]['fluxdErr'])
            #spwslist  = list(fluxscale_result['spwID'])

            for i in range(0,len(flux_d)):
                if (flux_d[i] != -1.0 and flux_d[i] != 0.0):
                    sources.append(sourcename)
                    flux_densities.append([float(flux_d[i]), float(flux_d_err[i])])
                    spws.append(int(spw_id))

    logprint("Finding sources returned from fluxscale.", logfileout='logs/phasecal.log')

    ii = 0
    unique_sources = list(np.unique(sources))
    results = []
    for source in unique_sources:
        indices = []
        for ii in range(len(sources)):
            if (sources[ii] == source):
                indices.append(ii)
        unique_bands = ['L']
        band=unique_bands[0]  # Needs to be defined for a logprint statement.  Meaningless otherwise
        lfreqs = []
        lfds = []
        lerrs = []
        uspws = []
        for ii in range(len(indices)):
                lfreqs.append(log10(center_frequencies[spws[indices[ii]]]))
                lfds.append(log10(flux_densities[indices[ii]][0]))
                lerrs.append(log10(e) * flux_densities[indices[ii]][1]/flux_densities[indices[ii]][0])
                uspws.append(spws[indices[ii]])
# if we didn't care about the errors on the data or the fit coefficients, just:
#       coefficients = np.polyfit(lfreqs, lfds, 1)
# or, if we ever get to numpy 1.7.x, for weighted fit, and returning
# covariance matrix, do:
#       ...
#       weights = []
#       weight_sum = 0.0
#       for ii in range(len(lfreqs)):
#           weights.append(1.0 / (lerrs[ii]*lerrs[ii]))
#           weight_sum += weights[ii]
#       for ii in range(len(weights)):
#           weights[ii] /= weight_sum
#       coefficients = np.polyfit(lfreqs, lfds, 1, w=weights, cov=True)
# but, for now, use the full scipy.optimize.leastsq route...
#
# actually, after a lot of testing, np.polyfit does not return a global
# minimum solution.  sticking with leastsq (modified as below to get the
# proper errors), or once we get a modern enough version of scipy, moving
# to curve_fit, is better.
#

        if len(lfds) < 2:
            aa = lfds[0]
            bb = 0.0
            SNR = 0.0
        else:
            alfds = scp.array(lfds)
            alerrs = scp.array(lerrs)
            alfreqs = scp.array(lfreqs)
            pinit = [0.0, 0.0]
            fit_out = scpo.leastsq(errfunc, pinit, args=(alfreqs, alfds, alerrs), full_output=1)
            pfinal = fit_out[0]
            covar = fit_out[1]
            aa = pfinal[0]
            bb = pfinal[1]
#
# the fit is of the form:
#     log(S) = a + b * log(f)
# with a = pfinal[0] and b = pfinal[1].  the errors on the coefficients are
# sqrt(covar[i][i]*residual_variance) with the residual covariance calculated
# as below (it's like the reduced chi squared without dividing out the errors).
# see the scipy.optimize.leastsq documentation and 
# http://stackoverflow.com/questions/14854339/in-scipy-how-and-why-does-curve-fit-calculate-the-covariance-of-the-parameter-es
#
            summed_error = 0.0
            for ii in range(len(alfds)):
                model = aa + bb*alfreqs[ii]
                residual = (model - alfds[ii]) * (model - alfds[ii])
                summed_error += residual
            residual_variance = summed_error / (len(alfds) - 2)
            SNR = fabs(bb) / sqrt(covar[1][1] * residual_variance)

#
# take as the reference frequency the lowest one.  (this shouldn't matter,
# in principle).
#

        logprint("Setting flux density, ref. freq, and spix", logfileout='logs/phasecal.log')
        
        reffreq = 10.0**lfreqs[0]/1.0e9
        fluxdensity = 10.0**(aa + bb*lfreqs[0])
        spix = bb
        results.append([ source, uspws, fluxdensity, spix, SNR, reffreq ])
        logprint(source + ' ' + ', '+ band + ' fitted spectral index & SNR = ' + str(spix) + ' ' + str(SNR), logfileout='logs/phasecal.log')
        logprint("Frequency, data, error, and fitted data:", logfileout='logs/phasecal.log')
        for ii in range(len(lfreqs)):
            SS = fluxdensity * (10.0**lfreqs[ii]/reffreq/1.0e9)**spix
            fderr = lerrs[ii]*(10**lfds[ii])/log10(e)
            logprint('    '+str(10.0**lfreqs[ii]/1.0e9)+'  '+ str(10.0**lfds[ii])+'  '+str(fderr)+'  '+str(SS), logfileout='logs/phasecal.log')
    
    
    logprint ("Setting power-law fit in the model column", logfileout='logs/phasecal.log')

    for result in results:
        for spw_i in result[1]:
#
# here, check on SNR, but don't do this yet, until we know what typical SNRs are
#
#           if result[4] > SNRlimit:
            logprint('Running setjy on spw '+str(spw_i), logfileout='logs/phasecal.log')
            default('setjy')
            vis=ms_active
            field = str(result[0])
            #spw = ','.join(["%s" % ii for ii in result[1]])
            spw = str(spw_i)
            selectdata=False
            scalebychan=True
            standard='manual'
            fluxdensity = [ result[2], 0, 0, 0 ]
            spix = result[3]
            reffreq = str(result[5])+'GHz'
            usescratch=False
            try:
                setjy()
                if (abs(spix) > 5.0):
                    QA2_fluxboot='Fail'
            except:
                logprint("Unable to complete flux scaling operation for field "+str(field)+", spw "+str(spw), logfileout='logs/phasecal.log')
except:
    logprint("A problem was detected while running fluxscale.  Please review the CASA log.", logfileout='logs/phasecal.log')


# Step 9: do final applycal to all fields

logprint ("Running final applycal", logfileout='logs/phasecal.log')

AllCalTables=copy.copy(GainTables)      # using both amp.gcal & flux.gcal tables
AllCalTables.append('finalamp.gcal')
AllCalTables.append('finalflux.gcal')
AllSpwMapValues=copy.copy(SpwMapValues)
AllSpwMapValues.append([])
AllSpwMapValues.append([])

# Backup finalamp.gcal table
os.system('cp finalamp.gcal finalamp_backup.gcal')

# HG: set gainfield appropriately depending if antposcal.p exists or not.  
if os.path.exists('antposcal.p')==True:  
  PhaseFields=['','','1331+305=3C286','1331+305=3C286','J0943-0819','J0943-0819','J0943-0819']
  FluxFields =['','','1331+305=3C286','1331+305=3C286','1331+305=3C286','1331+305=3C286','1331+305=3C286']

if os.path.exists('antposcal.p')==False:
  PhaseFields=['','1331+305=3C286','1331+305=3C286','J0943-0819','J0943-0819','J0943-0819']
  FluxFields =['','1331+305=3C286','1331+305=3C286','1331+305=3C286','1331+305=3C286','1331+305=3C286']


default('applycal')
vis=ms_active
field='J0943-0819'            # Apply final calibration to phase cal
spw=''
intent=''
selectdata=True
gaintable=AllCalTables
gainfield=PhaseFields
interp=['']
spwmap=AllSpwMapValues   # Was [] in previous version, now corresponds to AllCalTables
gaincurve=False
opacity=[]
parang=False
calwt=False
flagbackup=False
async=False
applycal()

default('applycal')
vis=ms_active
field='1331+305=3C286'            # Apply same calibration to phase cal
spw=''
intent=''
selectdata=True
gaintable=AllCalTables
gainfield=FluxFields
interp=['']
spwmap=AllSpwMapValues  # Was [] in previous version, now corresponds to AllCalTables
gaincurve=False
opacity=[]
parang=False
calwt=False
flagbackup=False
async=False
applycal()

# Save flags
logprint ("Saving flags", logfileout='logs/phasecal.log')


default('flagmanager')
vis=ms_active
mode='save'
versionname='phasecal_flags'
comment='Phase Calibrator flags saved after application'
merge='replace'
flagmanager()
logprint ("Flag column saved to "+versionname, logfileout='logs/phasecal.log')

# Step 10: Makes diagnostic plots for assessment
# Look at calibration tables for phase cal (amp, phase vs. time, frequency)
# Make images of phase cal and look at flux,beam vs. spw

logprint ("Making diagnostic plots", logfileout='logs/phasecal.log')

# Make plot of flagging statistics
# s_p is the output of flagdata run (above)
# Get information for flagging percentage vs. uvdistance
#gantdata = get_antenna_data(ms_active)
#create adictionary with flagging info
#base_dict = create_baseline_dict(ms_active, gantdata)
#gantdata and base_dict are already in the initial module so no need to retrieve that information again.
#match flagging data to dictionary entry
datamatch = flag_match_baseline(s_p['baseline'], base_dict)
#bin the statistics
binned_stats = bin_statistics(datamatch, 'B', 25)  # 25 is the number of uvdist bins such that there is minimal error in uvdist.

#Plot flagging % vs. uvdist
### Plot the Data
barwidth = binned_stats[0][1]
totflagged = 'Phase Cal Flagging: '+ str(phase_flag*100) + '% Data Flagged'
pylab.close()
pylab.bar(binned_stats[0],binned_stats[1], width=barwidth, color='grey', align='edge')
pylab.title(totflagged)
pylab.grid()
pylab.ylabel('flagged data [%]')
pylab.xlabel('average UV distance [m]')
pylab.savefig('phase_flag_uvdist.png')
pylab.close()

os.system("mv phase_flag_uvdist.png plots/.") 

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
        flagged.append(s_p['spw:channel'][str(s)+':'+str(c)]['flagged'])
        totals.append(s_p['spw:channel'][str(s)+':'+str(c)]['total'])
        flag_frac.append(flagged[ct]/totals[ct])

fig=pylab.figure()
pylab.plot(freq,flag_frac,'k-')
pylab.xlim(940.,1445.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Fraction of data flagged')
pylab.savefig("phase_flag.png")
pylab.close(fig)

#Plot UV spectrum (averaged over all baselines & time) of phase calibrator
default('plotms')
vis=ms_active   
field='J0943-0819'        # Only for J0943
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

# Plot calibrated phase vs. frequency
xaxis='freq'
yaxis='phase'
plotrange=[0.95,1.43,0,0]
clearplots=True
plotfile='phasecal_phasespectrum.png'
plotms()

ms_name=ms_active[:-3]
output_ms=ms_name+'_phasecal_flux_averaged.ms'

# Remove old averaged MS if it exists.
if os.path.exists(output_ms):
    os.system("rm -rf "+output_ms)
    
#Split the phase cal.
default('oldsplit')
vis=ms_active
outputvis=output_ms
datacolumn='corrected'
field='J0943-0819'
spw='0~14:128~1920'  # Average over all channels, except the very edges
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
    default('tclean')
    image_name='phasecalibrator_spw'+str(ii)
    fieldid='J0943*'
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
if os.path.exists('images')==False:
    os.system("mkdir images")
os.system("mv phasecalibrator_spw*.* images")

#Create webpage with results

if os.path.exists('phasecal.html'):
    os.system("rm phasecal.html")
wlog = open("phasecal.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('<style>table tr {page-break-inside: avoid}</style>\n')
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
wlog.write('<li> Amp vs. Phase: \n')
for ii in seq:
    wlog.write('<br><img src="plots/phasecal_ampphase_Spw'+str(ii)+'.png">\n')
wlog.write('<li> Spectrum of Phase calibrator (both LL & RR, averaged over all time & baselines): \n')
wlog.write('<br><img src="plots/phasecal_spectrum_full.png">\n')
wlog.write('<br><img src="plots/phasecal_spectrum_zoom.png">\n')
wlog.write('<br><img src="plots/phasecal_phasespectrum.png">\n')
wlog.write('<li> Amp. & Phase vs. time for Phase Calibrator (averaged over frequency): \n')
wlog.write('<table> \n')
for ii in seq:
    wlog.write('<tr>\n')
    wlog.write('<td><img src="plots/phasecal_amptime_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('<td><img src="plots/phasecal_phasetime_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('</tr>\n')
wlog.write('</table> \n')
wlog.write('</li>')
wlog.write('<li> Measured properties of phase calibrator: \n')
wlog.write('<br><img src="plots/phasecal_beamsize.png">\n')
wlog.write('<br><img src="plots/phasecal_peak.png">\n')
wlog.write('<br><img src="plots/phasecal_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Flagging percentage vs. frequency:\n')
wlog.write('<li><br><img src="./plots/phase_flag.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Flagging percentage vs. uvdist\n')
wlog.write('<li><br><img src="./plots/phase_flag_uvdist.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Percentage of J0943 data flagged: '+str(phase_flag*100)+'\n')
wlog.write('<br>')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()


logprint ("Finished CHILES_pipe_phasecal.py", logfileout='logs/phasecal.log')
time_list=runtiming('phase', 'end')

pipeline_save()
