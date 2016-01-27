# CHILES_pipe_phasecal.py
# This module of the CHILES pipeline follows the "bandpass" module using the BP 
# and flux calibration from that task to solve for the complex gain (amp & phase)
# as a function of time using the phase calibrator.  This will involve doing an
# initial solution, running an RFLAG with extend, then getting a final solution
# ready to apply to the observations of the deepfield itself.
# 9/7/15 DJP
# HG fixed a number of bugs and explicitly sets caltable and spw lists correctly.
# 12/7/15 DJP
# 1/27/16 DJP: Set extendpols=False in RFLAG step.


logprint ("Starting CHILES_pipe_phasecal.py", logfileout='logs/phasecal.log')
time_list=runtiming('phase', 'start')

# Step 1:  Set default parameters (to match previous modules), and clean up old files:
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

#Remove old files from previous attempts using rmtables
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
os.system("rm -rf images/phasecalibrator_spw*.*")


#Set minimum # of baselines need for a solution to 8, based on experience
minBL_for_cal=8

#Set uvrange to apply in order to optimally exclude RFI without flagging:
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
# Field='0' is J0943-0819, phase calibrator
# Field='2' is 1331+305 (3C286), flux & bandpass calibrator
# Field='1' is "deepfield", the target


  
# Step 2:  gaincal on both calibrators
# Need short solint = "int" (for calibrators)
# Long solint = "inf" (to apply to target) not needed until final iteration.

logprint ("Running initial gaincal", logfileout='logs/phasecal.log')

default('gaincal')
vis=ms_active
caltable='initialphase.gcal'
field='0,2'
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
spwmap=priorspwmap     # In previous version this was [].  These spw's correspond to priorcals.
parang=False
async=False
gaincal()
  
# Step 3: Now use to derive amp & phase solutions
# This step is used to derive solutions with solint='inf' (averaging over scans only)

logprint ("Deriving initial amp & phase solutions", logfileout='logs/phasecal.log')

GainTables=copy.copy(priorcals)
GainTables.append('initialphase.gcal')
SpwMapValues=copy.copy(priorspwmap)
SpwMapValues.append([])


default('gaincal')
vis=ms_active
caltable='initialamp.gcal'
field='0,2'
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

# Step 4:  Bootstrap fluxes using code from EVLA_pipe_fluxboot.
fluxscale_output=msname.rstrip('ms')+'fluxdensities'
fluxcalfields=flux_field_select_string
if os.path.exists(fluxscale_output):
    os.system('rm -rf '+fluxscale_output)

logprint ("Initial flux bootstrapping", logfileout='logs/phasecal.log')
logprint ("Flux densities will be written to "+fluxscale_output, logfileout='logs/phasecal.log')

casalog.setlogfile(fluxscale_output)

default('fluxscale')
vis=ms_active
caltable='initialamp.gcal'
fluxtable='initialflux.gcal'
reference='2'      # Reference field is always 3C286
transfer=['0']
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

    ii = 0
    unique_sources = list(np.unique(sources))
    results = []
    for source in unique_sources:
        indices = []
        for ii in range(len(sources)):
            if (sources[ii] == source):
                indices.append(ii)
        unique_bands = ['L']
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
        reffreq = 10.0**lfreqs[0]/1.0e9
        fluxdensity = 10.0**(aa + bb*lfreqs[0])
        spix = bb
        results.append([ source, uspws, fluxdensity, spix, SNR, reffreq ])
        logprint(source + ' ' + band + ' fitted spectral index & SNR = ' + str(spix) + ' ' + str(SNR), logfileout='logs/phasecal.log')
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


# Step 5: applycal to flux & phase cal

logprint ("Initial applycal", logfileout='logs/phasecal.log')

AllCalTables=copy.copy(GainTables)        # Using both amp.gcal and flux.gcal tables.
AllCalTables.append('initialamp.gcal')
AllCalTables.append('initialflux.gcal')
AllSpwMapValues=copy.copy(SpwMapValues)
AllSpwMapValues.append([])
AllSpwMapValues.append([])

# HG: set gainfield appropriately depending if antposcal.p exists or not.  
if os.path.exists('antposcal.p')==True:  
  PhaseFields=['','','2','2','0','0','0']
  FluxFields =['','','2','2','2','2','2']

if os.path.exists('antposcal.p')==False:
  PhaseFields=['','2','2','0','0','0']
  FluxFields =['','2','2','2','2','2']


default('applycal')
vis=ms_active
field='0'            # Apply one calibration to phase cal
spw=''
intent=''
selectdata=True
gaintable=AllCalTables
gainfield=PhaseFields
interp=['']             # Default interpolation is linear in time & frequency
spwmap=AllSpwMapValues  # In previous version this was [].  These spw's correspond to AllCalTables.
gaincurve=False
opacity=[]
parang=False
calwt=False
flagbackup=False
async=False
applycal()

default('applycal')
vis=ms_active
field='2'            # Apply same calibration to flux cal
spw=''
intent=''
selectdata=True
gaintable=AllCalTables
gainfield=FluxFields
interp=['']             # Default interpolation is linear in time & frequency
spwmap=AllSpwMapValues  # In previous version this was [].  These spw's correspond to AllCalTables.
gaincurve=False
opacity=[]
parang=False
calwt=False
flagbackup=False
async=False
applycal()

# Summary of flagging, after initial phase calib (for testing purposes only)
logprint ("Summary of flags after initial phase cal", logfileout='logs/phasecal.log')
default('flagdata')
vis=ms_active
mode='summary'
spw='0~14'
correlation='RR,LL'
spwchan=True
spwcorr=True
action='calculate'
flagdata()

# Step 6: run RFLAG and extend

logprint ("Initial RFLAG", logfileout='logs/phasecal.log')

ff= 0   # This sets the field being flagged; here it is the phase cal, field 0.

default('flagdata')
vis=ms_active
mode='rflag'
field='0'            # Hard-coded to field 0 as this is always phase cal
spw='0~14'
correlation=''
ntime='scan'
combinescans=False
datacolumn='corrected'
extendflags=False    # Explicitly set to False as default is True.  Extend is explicitly done later.  
extendpols=False     # Default is True.  May allow some weak RFI through, but try it.   
winsize=3
# The following noise levels are taken from tests by XF.
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

# Summary of flagging, after RFLAG (for testing purposes only)
logprint ("Summary of flags after RFLAG", logfileout='logs/phasecal.log')
default('flagdata')
vis=ms_active
mode='summary'
spw='0~14'
correlation='RR,LL'
spwchan=True
spwcorr=True
action='calculate'
flagdata()

# Extend flags
default('flagdata')
vis=ms_active
mode='extend'
field='0'            # Hard-coded to field 0 as this is always phase cal
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
flagbackup=True
savepars=True
async=False
flagdata()

# Summary of flagging, after RFLAG+extend (for testing purposes only)
logprint ("Summary of flags after RFLAG+extend", logfileout='logs/phasecal.log')
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
logprint ("Saving flags", logfileout='logs/phasecal.log')


default('flagmanager')
vis=ms_active
mode='save'
versionname='phasecal_flags'
comment='Phase Calibrator flags saved after application'
merge='replace'
async=False
flagmanager()
logprint ("Flag column saved to "+versionname, logfileout='logs/phasecal.log')

# Step 7: repeat gaincal 
# Short solint = "int" (length of individual calibrator scan)

logprint ("Running final gaincal", logfileout='logs/phasecal.log')

default('gaincal')
vis=ms_active
caltable='finalphase_int.gcal'
field='0,2'
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
field='0,2'
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
field='0,2'
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
field='0'  # Hard-coded for phase calibrator
scr=False
delmod()

casalog.setlogfile(fluxscale_output)

default('fluxscale')
vis=ms_active
caltable='finalamp.gcal'
fluxtable='finalflux.gcal'
reference='2'
transfer=['0']
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

    ii = 0
    unique_sources = list(np.unique(sources))
    results = []
    for source in unique_sources:
        indices = []
        for ii in range(len(sources)):
            if (sources[ii] == source):
                indices.append(ii)
        unique_bands = ['L']
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
        reffreq = 10.0**lfreqs[0]/1.0e9
        fluxdensity = 10.0**(aa + bb*lfreqs[0])
        spix = bb
        results.append([ source, uspws, fluxdensity, spix, SNR, reffreq ])
        logprint(source + ' ' + band + ' fitted spectral index & SNR = ' + str(spix) + ' ' + str(SNR), logfileout='logs/phasecal.log')
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

# HG: set gainfield appropriately depending if antposcal.p exists or not.  
if os.path.exists('antposcal.p')==True:  
  PhaseFields=['','','2','2','0','0','0']
  FluxFields =['','','2','2','2','2','2']

if os.path.exists('antposcal.p')==False:
  PhaseFields=['','2','2','0','0','0']
  FluxFields =['','2','2','2','2','2']


default('applycal')
vis=ms_active
field='0'            # Apply final calibration to phase cal
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
field='2'            # Apply same calibration to phase cal
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

# Summary of flagging, after final phase calib (for testing purposes only)
logprint ("Summary of flags after final phase cal", logfileout='logs/phasecal.log')
default('flagdata')
vis=ms_active
mode='summary'
spw='0~14'
correlation='RR,LL'
spwchan=True
spwcorr=True
action='calculate'
flagdata()

#Clear model of phase calibrator for future runs of phasecal; has no effect on corrected column if all is okay already.
default('delmod')
vis=ms_active
otf=True
field='0'  # Hard-coded for phase calibrator
scr=False
delmod()

# Step 10: Makes diagnostic plots for assessment
# Look at calibration tables for phase cal (amp, phase vs. time, frequency)
# Make images of phase cal and look at flux,beam vs. spw

logprint ("Making diagnostic plots", logfileout='logs/phasecal.log')

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
    field='0'       # Hard coded to J0943, phase calibrator
    xaxis='phase'
    yaxis='amp'
    xdatacolumn='corrected'
    ydatacolumn='corrected'
    averagedata=True
    #avgtime='1e5'
    #avgscan=True
    avgchannel='2048'  # Hard-coded for channels in each spw for CHILES 
    avgspw=False
    spw=str(ii)
    gridrows=1
    #gridcolumns=
    showlegend=False
    iteraxis='spw'
    coloraxis='corr'
    showgui=False
    clearplots=True
    plotfile='phasecal_ampphase.png'
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
pylab.yscale('log')
pylab.savefig('phasecal_peak.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,rms_phase,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('RMS [mJy]')
pylab.yscale('log')
pylab.savefig('phasecal_rms.png')
pylab.close(fig)

#Want to plot image of flux calibrator in each spw.  Use "imview"

for ii in seq:
    image_phasecal='phasecalibrator_spw'+str(ii)+'.image'
    kntr_levels=[-2*rms_phase[ii]/1000.,2*rms_phase[ii]/1000.,0.1*max_phase[ii],0.3*max_phase[ii],0.5*max_phase[ii],0.7*max_phase[ii],0.9*max_phase[ii]]
    imview(raster={'file':image_phasecal,'scaling':-3, 'colorwedge':True},contour={'file':image_phasecal,'levels':kntr_levels},out='phasecal_spw'+str(ii)+'.png')

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

default('plotcal')
caltable='finalflux.gcal'
xaxis='time'
yaxis='amp'
showgui=False
figfile='caltable_finalflux_amp.png'
plotcal()

yaxis='phase'
figfile='caltable_finalflux_phase.png'
plotcal()

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
wlog.write('<li>finalflux.gcal Amp vs. Time: \n')
wlog.write('<br><img src="plots/caltable_finalflux_amp.png"></li>\n')
wlog.write('<li>finalflux.gcal Phase vs. Time: \n')
wlog.write('<br><img src="plots/caltable_finalflux_phase.png"></li>\n')
wlog.write('<li> Amp vs. Phase: \n')
for ii in seq:
    wlog.write('<br><img src="plots/phasecal_ampphase_Spw'+str(ii)+'.png">\n')
#    wlog.write('<br><img src="plots/phasecal_phaseuvdist_Spw'+str(ii)+'.png">\n')
wlog.write('<li> Images of Phase Calibrator: \n')
for ii in seq:
    wlog.write('<br><img src="plots/phasecal_spw'+str(ii)+'.png">\n')
wlog.write('</li>')
wlog.write('<li> Measured properties of phase calibrator: \n')
wlog.write('<br><img src="plots/phasecal_beamsize.png">\n')
wlog.write('<br><img src="plots/phasecal_peak.png">\n')
wlog.write('<br><img src="plots/phasecal_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()


logprint ("Finished CHILES_pipe_phasecal.py", logfileout='logs/phasecal.log')
time_list=runtiming('phase', 'end')

pipeline_save()
