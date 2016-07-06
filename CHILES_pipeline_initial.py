# CHILES_pipeline_initial.py
# CHILES pipeline, first module
# based on the EVLA pipeline, the past CHILES pipeline, the CHILES AIPS 
# pipeline, and work done by Ximena, this pipeline is designed to fully 
# calibrate CHILES data taken in 2015.  
# This first module imports data into a MS, gets summary information, identifies
# possible reference antennas, applies default flags, hanning smooths the data,
# and removes the continuum spectral windows.  This module requires minimal 
# human interaction, but prepares the user for running subsequent modules.
# v0.1: 08/25/15 DJP
# v0.2: changed to CASA 4.5 11/13/15 DJP
# v0.3: Using rmtables to remove old calibration tables, updated plots. 11/22/15 DJP
# v0.5: The entire pipeline now appears to work and produce reasonable results.
#       Target flagging may still need to be tweaked.  11/23/15 DJP
# v0.6: Removed clipping from target, 
#       removed Shadow flagging since unnecessary in B array, and fixed small
#       bugs in calibration routines (mostly done by HG). 12/7/15 DJP
# v0.7: Fixed more bugs, removed extend from target (but included clipping), 
#       and changed incremental back to True in fluxscale
#       12/10/15 DJP
# v0.8: Set extendflags=False for all rflag runs.  Removed clipping and extend from target module.
#       12/21/15 DJP
# v0.9: Fixed typos in code, run delmod before setjy (to clear past models), 
#       added flagdata, mode='summary' runs after all applycals and flagdatas 
#       in order to ID excessive flagging (for testing purposes only).  Testing 
#       bandpass with fillgaps=10, interp='linear','spline' since BP seems to be interpolating anyways.  
#       1/13/16 DJP
# v0.10:  Fixed code to apply online flags & zero flags when importing SDM with importevla
#         Also extendpols=False set for all RFLAG runs.
#         1/27/16 DJP
# v0.11:  Added check to apply zero flags (and remind user about online flags) when working with MS.
#         removed multiple flagdata, mode='summary' runs and limit to specific sources.
#         remove explicit interpolation in bandpass, applycal runs (unclear if it is useful to vary)
# v0.12:  Updated to run under CASA 4.5.1
# v1.0:  Fixed ff parameter for rflag runs, setting to float instead of string or integer; reimplement extend on target.
#        Extend BP solution to all channels, but only plot average values in diagnostics.
#        Make plot of UV spectrum for fluxcal, phasecal, and target.
#        Make plot of fraction of flagged data on deepfield as a function of channel.
# v1.1:  Fixed some small bugs and tweaked plots.
# v1.2:  Updated code to work on CASA 4.5.2, validated results to be the same.
# v1.3:  Code includes minimum SNR setting for calibration solutions, includes the "testcubes" module, and is updated to run on CASA 4.6
# v1.3.1:  Fixed code to account for new mode of hanningsmooth using fix from Emmanuel Momjian.
# v1.4:    Update path to work on Epoch 3 data.  
# v1.4.1:  Fix error with snrval, and allow uvmin, minBL_for_cal, and snrval to be set by user. Updated diagnostic plots.
# v1.5:  Revert minsnr to original pipeline settings.  Remove user input of these values.  Including amp/phase vs. time plots for phase/flux calib.
# v1.6:  Including new extend values for target plus time averaged rflag.  testcubes module does a split with time/frequency averaging before imaging.

version = "1.6"
svnrevision = '11nnn'
date = "2016Jul06"

print "Pipeline version "+version+" for use with CASA 4.6"
import sys
import pylab as pylab
# include additional packages for hanningsmooth
import shutil
import glob
import os


# Check that we are using the correct version of CASA
[major,minor,revision] = casadef.casa_version.split('.')
casa_version = 100*int(major)+10*int(minor)+int(revision)
if casa_version != 460:
    sys.exit("Your CASA version is "+casadef.casa_version+", please re-start using CASA 4.6")

# Define location of pipeline
#pipepath='/lustre/aoc/cluster/pipeline/script/prod/'
pipepath='/data/dpisano/CHILES/chiles_pipeline/'
#pipepath='/lustre/aoc/projects/chiles/chiles_pipeline/'

# Define location of output files from NRAO continuum pipeline
nrao_weblog_path='/data/dpisano/CHILES/weblogs/'
#nrao_weblog_path='/lustre/aoc/projects/chiles/weblogs/'

# To find the path to the weblog index.html file for the NRAO pipeline run:
def find(name,path):
    for root,dirs,files in os.walk(path):
        if name in files:
            return os.path.join(root,name)

#This is the default time-stamped casa log file, in case we
#    need to return to it at any point in the script
log_dir='logs'
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

maincasalog = casalogger.func_globals['thelogfile']

def logprint(msg, logfileout=maincasalog):
    print (msg)
    casalog.setlogfile(logfileout)
    casalog.post(msg)
    casalog.setlogfile(maincasalog)
    casalog.post(msg)
    return

#Create timing profile list and file if they don't already exist
if 'time_list' not in globals():
    time_list = []

timing_file='logs/timing.log'

if not os.path.exists(timing_file):
    timelog=open(timing_file,'w')
else:
    timelog=open(timing_file,'a')
    
def runtiming(pipestate, status):
    '''Determine profile for a given state/stage of the pipeline
    '''
    time_list.append({'pipestate':pipestate, 'time':time.time(), 'status':status})
#    
    if (status == "end"):
        timelog=open(timing_file,'a')
        timelog.write(pipestate+': '+str(time_list[-1]['time'] - time_list[-2]['time'])+' sec \n')
        timelog.flush()
        timelog.close()
        #with open(maincasalog, 'a') as casalogfile:
        #    tempfile = open('logs/'+pipestate+'.log','r')
        #    casalogfile.write(tempfile.read())
        #    tempfile.close()
        #casalogfile.close()
#        
    return time_list

# Now to get things setup for the pipeline

# The following script includes all the definitions and functions and
# prior inputs needed by a run of the pipeline.

time_list=runtiming('startup', 'start')
execfile(pipepath+'CHILES_pipe_startup.py')
time_list=runtiming('startup', 'end')
pipeline_save()

weblog_file=find('index.html',nrao_weblog_path+SDM_name)

try:

######################################################################

    time_list=runtiming('initial', 'start')

# Query user to set snrval, minBL_for_cal, and uvmin (all for calibration solutions)
    
#Set minimum # of baselines need for a solution to 8, based on experience
#    minBL_for_cal=raw_input("What is the minimum number of baselines required for a calibration solution (default=8; hit return for default): ")
#    if minBL_for_cal=='':
#        minBL_for_cal=int(8)
#    else:
#        minBL_for_cal=int(minBL_for_cal)
#
##Set uvrange to apply in order to optimally exclude RFI without flagging:
#    uvr_cal=raw_input("What should the minimum uvrange be for a calibration solution (default='>1500m'; hit return for default): ")
#    if uvr_cal=='':
#        uvr_cal='>1500m'
#
##Set minsnr value to use in calibration tasks:
#    snrval=raw_input("What is the minimum SNR needed to solve for calibration (default=8; hit return for default): ")
#    if snrval=='':
#        snrval=float(8.)
#    else:
#        snrval=float(snrval)

# IMPORT THE DATA TO CASA

    if (os.path.exists(msname) == False):
        logprint ("Creating measurement set", logfileout='logs/initial.log')

        scanlist=raw_input("If desired, enter list of subset of scans (using CASA nomenclature; leave blank for all scans): ")

        default('importevla')
        asdm=SDM_name
        vis=msname
        ocorr_mode='co'
        compression=False
        asis=''
        scans=scanlist
        verbose=True
        overwrite=False
        online=True         #KMH, calculate online flags when importing. 
        flagzero=True       # Calculate zero flags when importing.
        flagpol=False
        shadow=False
        tolerance=0.0
        addantenna=''
        applyflags=True     #Apply online, zero flags on import.
        savecmds=False
        flagbackup=False
        importevla()

        logprint ("Measurement set "+msname+" created", logfileout='logs/initial.log')
        ms_create_flag=True
    else:
        logprint ("Measurement set already exists, will use "+msname, logfileout='logs/initial.log')
        ms_create_flag=False
# Online flagging should be done when importing ASDM or creating MS.  Remind user of this.
        logprint ("If MS created outside this module, be sure online flags were applied then.", logfileout='logs/initial.log')
        
######################################################################

# HANNING SMOOTH (OPTIONAL, MAY BE IMPORTANT IF THERE IS NARROWBAND RFI)

    if myHanning.lower() == "y":
        logprint ("Hanning smoothing the data", logfileout='logs/initial.log')

        default('hanningsmooth')
        vis=msname
        datacolumn='data'
        outputvis='temphanning.ms'
        hanningsmooth()
        myHanning="n"

        logprint ("Copying xml files to the output ms")
        for file in glob.glob(msname+'/*.xml'):
                shutil.copy2(file , 'temphanning.ms/')
        logprint ("Removing original VIS '+msname, logfileout='logs/initial.log")
        shutil.rmtree(msname)
        logprint('Renaming temphanning.ms to '+msname, logfileout='logs/initial.log')
        os.rename('temphanning.ms', msname)
        logprint ("Hanning smoothing finished, myHanning parameter reset to 'n' to avoid further smoothing on restarts", logfileout='logs/initial.log')
    else:
        logprint ("NOT Hanning smoothing the data, keeping original ms unchanged", logfileout='logs/initial.log')

######################################################################

# GET SOME INFORMATION FROM THE MS THAT WILL BE NEEDED LATER, LIST
# THE DATA, AND MAKE SOME PLOTS

# Run listobs

    logprint ("Listing ms contents", logfileout='logs/initial.log')

    listname=msname.rstrip('ms') + 'listobs'
    syscommand='rm -rf '+listname
    os.system(syscommand)

    default('listobs')
    vis=ms_active
    selectdata=False
    verbose=True
    listfile=listname
    listobs()

# Identify spw information

    tb.open(ms_active+'/SPECTRAL_WINDOW')
    channels = tb.getcol('NUM_CHAN')
    originalBBClist = tb.getcol('BBC_NO')
    spw_bandwidths=tb.getcol('TOTAL_BANDWIDTH')
    reference_frequencies = tb.getcol('REF_FREQUENCY')
    center_frequencies = []
    for ii in range(len(reference_frequencies)):
        center_frequencies.append(reference_frequencies[ii]+spw_bandwidths[ii]/2)
    tb.close()
    unique_bands = 'L'				# KMH
    unique_bands_string = ','.join(["%s" % ii for ii in unique_bands])
    logprint("unique band string = " + unique_bands_string, logfileout='logs/initial.log')

    numSpws = len(channels)

# Set up spw selection for initial gain solutions

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

# Identify number of fields, positions, and source IDs

    tb.open(ms_active+'/FIELD')
    numFields = tb.nrows()
    field_positions = tb.getcol('PHASE_DIR')
    field_ids=range(numFields)
    field_names=tb.getcol('NAME')
    tb.close()

# Map field IDs to spws

    field_spws = []
    for ii in range(numFields):
        field_spws.append(spwsforfield(ms_active,ii))
    
# Identify scan numbers, map scans to field ID, and run scan summary
# (needed for figuring out integration time later)

    tb.open(ms_active)
    scanNums = sorted(np.unique(tb.getcol('SCAN_NUMBER')))
    field_scans = []
    for ii in range(0,numFields):
        subtable = tb.query('FIELD_ID==%s'%ii)
        field_scans.append(list(np.unique(subtable.getcol('SCAN_NUMBER'))))
        subtable.close()
    tb.close()

# Identify intents

    tb.open(ms_active+'/STATE')
    intents=tb.getcol('OBS_MODE')
    tb.close()
    
# Figure out integration time used

    ms.open(ms_active)
    scan_summary = ms.getscansummary()
    ms_summary = ms.summary()
    ms.close()
    startdate=float(ms_summary['BeginTime'])
#
# scan list
#
    integ_scan_list = []
    for scan in scan_summary:
        integ_scan_list.append(int(scan))
    sorted_scan_list = sorted(integ_scan_list)
    
# Set the integration time:
    int_time=8.0                      #KMH
    logprint ("Maximum integration time is "+str(int_time)+"s", logfileout='logs/initial.log')
    
# Find scans for quacking

    scan_list = [1]
    old_scan = scan_summary[str(sorted_scan_list[0])]['0']
    old_field = old_scan['FieldId']
    old_spws = old_scan['SpwIds']
    for ii in range(1,len(sorted_scan_list)):
        new_scan = scan_summary[str(sorted_scan_list[ii])]['0']
        new_field = new_scan['FieldId']
        new_spws = new_scan['SpwIds']
        if ((new_field != old_field) or (set(new_spws) != set(old_spws))):
            scan_list.append(sorted_scan_list[ii])
            old_field = new_field
            old_spws = new_spws
    quack_scan_string = ','.join(["%s" % ii for ii in scan_list])
    
# Identify scans and fields associated with different calibrator intents

    bandpass_state_IDs = []
    delay_state_IDs = []
    flux_state_IDs = []
    polarization_state_IDs = []
    phase_state_IDs = []
    amp_state_IDs = []
    calibrator_state_IDs = []
    pointing_state_IDs = []
    for state_ID in range(0,len(intents)):
        state_intents = intents[state_ID].rsplit(',')
        for intent in range(0,len(state_intents)):
            scan_intent = state_intents[intent].rsplit('#')[0]
            subscan_intent = state_intents[intent].rsplit('#')[1]
            if (scan_intent == 'CALIBRATE_BANDPASS'):
                bandpass_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_DELAY'):
                delay_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_FLUX'):
                flux_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_POLARIZATION'):
                polarization_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_AMPLI'):
                amp_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_PHASE'):
                phase_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_POINTING'):
                pointing_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
    
    tb.open(ms_active)

    if (len(flux_state_IDs) == 0):
        logprint("ERROR: No flux density calibration scans found", logfileout='logs/initial.log')
        raise Exception("No flux density calibration scans found")
    else:
        flux_state_select_string = ('STATE_ID in [%s'%flux_state_IDs[0])
        for state_ID in range(1,len(flux_state_IDs)):
            flux_state_select_string += (',%s')%flux_state_IDs[state_ID]
        flux_state_select_string += ']'
        subtable = tb.query(flux_state_select_string)
        flux_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        flux_scan_select_string = ','.join(["%s" % ii for ii in flux_scan_list])
        logprint ("Flux density calibrator(s) scans are "+flux_scan_select_string, logfileout='logs/initial.log')
        flux_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        flux_field_select_string = ','.join(["%s" % ii for ii in flux_field_list])
        logprint ("Flux density calibrator(s) are fields "+flux_field_select_string, logfileout='logs/initial.log')
    
    if (len(bandpass_state_IDs) == 0):
        logprint ("No bandpass calibration scans defined, using flux density calibrator", logfileout='logs/initial.log')
        bandpass_scan_select_string=flux_scan_select_string
        logprint ("Bandpass calibrator(s) scans are "+bandpass_scan_select_string, logfileout='logs/initial.log')
        bandpass_field_select_string=flux_field_select_string
        logprint ("Bandpass calibrator(s) are fields "+bandpass_field_select_string, logfileout='logs/initial.log')
    else:
        bandpass_state_select_string = ('STATE_ID in [%s'%bandpass_state_IDs[0])
        for state_ID in range(1,len(bandpass_state_IDs)):
            bandpass_state_select_string += (',%s')%bandpass_state_IDs[state_ID]
        bandpass_state_select_string += ']'
        subtable = tb.query(bandpass_state_select_string)
        bandpass_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        bandpass_scan_select_string = ','.join(["%s" % ii for ii in bandpass_scan_list])
        logprint ("Bandpass calibrator(s) scans are "+bandpass_scan_select_string, logfileout='logs/initial.log')
        bandpass_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        bandpass_field_select_string = ','.join(["%s" % ii for ii in bandpass_field_list])
        logprint ("Bandpass calibrator(s) are fields "+bandpass_field_select_string, logfileout='logs/initial.log')
        if (len(bandpass_field_list) > 1):
            logprint ("WARNING: More than one field is defined as the bandpass calibrator.", logfileout='logs/initial.log')
            logprint ("WARNING: Models are required for all BP calibrators if multiple fields", logfileout='logs/initial.log')
            logprint ("WARNING: are to be used, not yet implemented; the pipeline will use", logfileout='logs/initial.log')
            logprint ("WARNING: only the first field.", logfileout='logs/initial.log')
            bandpass_field_select_string = str(bandpass_field_list[0])
    
    if (len(delay_state_IDs) == 0):
        logprint ("No delay calibration scans defined, using bandpass calibrator", logfileout='logs/initial.log')
        delay_scan_select_string=bandpass_scan_select_string
        logprint ("Delay calibrator(s) scans are "+delay_scan_select_string, logfileout='logs/initial.log')
        delay_field_select_string=bandpass_field_select_string
        logprint ("Delay calibrator(s) are fields "+delay_field_select_string, logfileout='logs/initial.log')
    else:
        delay_state_select_string = ('STATE_ID in [%s'%delay_state_IDs[0])
        for state_ID in range(1,len(delay_state_IDs)):
            delay_state_select_string += (',%s')%delay_state_IDs[state_ID]
        delay_state_select_string += ']'
        subtable = tb.query(delay_state_select_string)
        delay_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        delay_scan_select_string = ','.join(["%s" % ii for ii in delay_scan_list])
        logprint ("Delay calibrator(s) scans are "+delay_scan_select_string, logfileout='logs/initial.log')
        delay_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        delay_field_select_string = ','.join(["%s" % ii for ii in delay_field_list])
        logprint ("Delay calibrator(s) are fields "+delay_field_select_string, logfileout='logs/initial.log')
    
    if (len(polarization_state_IDs) == 0):
        logprint ("No polarization calibration scans defined, no polarization calibration possible", logfileout='logs/initial.log')
        polarization_scan_select_string=''
        polarization_field_select_string=''
    else:
        logprint ("Warning: polarization calibration scans found, but polarization calibration not yet implemented", logfileout='logs/initial.log')
        polarization_state_select_string = ('STATE_ID in [%s'%polarization_state_IDs[0])
        for state_ID in range(1,len(polarization_state_IDs)):
            polarization_state_select_string += (',%s')%polarization_state_IDs[state_ID]
        polarization_state_select_string += ']'
        subtable = tb.query(polarization_state_select_string)
        polarization_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        polarization_scan_select_string = ','.join(["%s" % ii for ii in polarization_scan_list])
        logprint ("Polarization calibrator(s) scans are "+polarization_scan_select_string, logfileout='logs/initial.log')
        polarization_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        polarization_field_select_string = ','.join(["%s" % ii for ii in polarization_field_list])
        logprint ("Polarization calibrator(s) are fields "+polarization_field_select_string, logfileout='logs/initial.log')
    
    if (len(phase_state_IDs) == 0):
        QA2_msinfo='Fail'
        logprint("ERROR: No gain calibration scans found", logfileout='logs/initial.log')
        raise Exception("No gain calibration scans found")
    else:
        phase_state_select_string = ('STATE_ID in [%s'%phase_state_IDs[0])
        for state_ID in range(1,len(phase_state_IDs)):
            phase_state_select_string += (',%s')%phase_state_IDs[state_ID]
        phase_state_select_string += ']'
        subtable = tb.query(phase_state_select_string)
        phase_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        phase_scan_select_string = ','.join(["%s" % ii for ii in phase_scan_list])
        logprint ("Phase calibrator(s) scans are "+phase_scan_select_string, logfileout='logs/initial.log')
        phase_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        phase_field_select_string = ','.join(["%s" % ii for ii in phase_field_list])
        logprint ("Phase calibrator(s) are fields "+phase_field_select_string, logfileout='logs/initial.log')
    
    if (len(amp_state_IDs) == 0):
        logprint ("No amplitude calibration scans defined, will use phase calibrator", logfileout='logs/initial.log')
        amp_scan_select_string=phase_scan_select_string
        logprint ("Amplitude calibrator(s) scans are "+amp_scan_select_string, logfileout='logs/initial.log')
        amp_field_select_string=phase_scan_select_string
        logprint ("Amplitude calibrator(s) are fields "+amp_field_select_string, logfileout='logs/initial.log')
    else:
        amp_state_select_string = ('STATE_ID in [%s'%amp_state_IDs[0])
        for state_ID in range(1,len(amp_state_IDs)):
            amp_state_select_string += (',%s')%amp_state_IDs[state_ID]
        amp_state_select_string += ']'
        subtable = tb.query(amp_state_select_string)
        amp_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        amp_scan_select_string = ','.join(["%s" % ii for ii in amp_scan_list])
        logprint ("Amplitude calibrator(s) scans are "+amp_scan_select_string, logfileout='logs/initial.log')
        amp_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        amp_field_select_string = ','.join(["%s" % ii for ii in amp_field_list])
        logprint ("Amplitude calibrator(s) are fields "+amp_field_select_string, logfileout='logs/initial.log')
    
# Find all calibrator scans and fields

    calibrator_state_select_string = ('STATE_ID in [%s'%calibrator_state_IDs[0])
    for state_ID in range(1,len(calibrator_state_IDs)):
        calibrator_state_select_string += (',%s')%calibrator_state_IDs[state_ID]
        
    calibrator_state_select_string += ']' 
    subtable = tb.query(calibrator_state_select_string)
    calibrator_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
    calibrator_scan_select_string = ','.join(["%s" % ii for ii in calibrator_scan_list])
    calibrator_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
    subtable.close()
    calibrator_field_select_string = ','.join(["%s" % ii for ii in calibrator_field_list])
    
    tb.close()
    
#Prep string listing of correlations from dictionary created by method buildscans
#For now, only use the parallel hands.  Cross hands will be implemented later.
    scandict = buildscans(ms_active)
    corrstring_list = scandict['DataDescription'][0]['corrdesc']
    removal_list = ['RL', 'LR', 'XY', 'YX']
    corrstring_list = list(set(corrstring_list).difference(set(removal_list)))
    corrstring = string.join(corrstring_list,',')
    logprint ("Correlations shown in plotms will be "+corrstring, logfileout='logs/initial.log')
    
#Get number of antennas, store in numAntenna
    tbLoc = casac.table()
    tbLoc.open( '%s/ANTENNA' % ms_active)
    nameAntenna = tbLoc.getcol( 'NAME' )
    numAntenna = len(nameAntenna)
    tbLoc.close()

#    minBL_for_cal=max(3,int(numAntenna/2.0))

# Set 3C84 variables so the pipeline doesn't complain later:
    cal3C84_d = False
    cal3C84_bp = False
    cal3C84 = False
    #uvrange3C84 = '0~1800klambda'
    uvrange3C84 = ''

# Identify bands/basebands/spws

    tb.open(ms_active+'/SPECTRAL_WINDOW')
    spw_names = tb.getcol('NAME')
    tb.close()
    
# If the dataset is too old to have the bandname in it, assume that
# either there are 8 spws per baseband (and allow for one or two for
# pointing), or that this is a dataset with one spw per baseband

    if (len(spw_names)>=8):
        critfrac=0.9/int(len(spw_names)/8.0)
    else:
        critfrac=0.9/float(len(spw_names))
    
    if '#' in spw_names[0]:
#
# i assume that if any of the spw_names have '#', they all do...
#
        bands_basebands_subbands = []
        for spw_name in spw_names:
            receiver_name, baseband, subband = spw_name.split('#')
            receiver_band = (receiver_name.split('_'))[1]
            bands_basebands_subbands.append([receiver_band, baseband, int(subband)])
        spws_info = [[bands_basebands_subbands[0][0], bands_basebands_subbands[0][1], [], []]]
        bands = [bands_basebands_subbands[0][0]]
        for ii in range(len(bands_basebands_subbands)):
            band,baseband,subband = bands_basebands_subbands[ii]
            found = -1
            for jj in range(len(spws_info)):
                oband,obaseband,osubband,ospw_list = spws_info[jj]
                if band==oband and baseband==obaseband:
                    osubband.append(subband)
                    ospw_list.append(ii)
                    found = jj
                    break
            if found >= 0:
                spws_info[found] = [oband,obaseband,osubband,ospw_list]
            else:
                spws_info.append([band,baseband,[subband],[ii]])
                bands.append(band)
        logprint("Bands/basebands/spws are:", logfileout='logs/initial.log')
        for spw_info in spws_info:
            spw_info_string = spw_info[0] + '   ' + spw_info[1] + '   [' + ','.join(["%d" % ii for ii in spw_info[2]]) + ']   [' + ','.join(["%d" % ii for ii in spw_info[3]]) + ']'
            logprint(spw_info_string, logfileout='logs/initial.log')
# Critical fraction of flagged solutions in delay cal to avoid an
# entire baseband being flagged on all antennas
        critfrac=0.9/float(len(spws_info))
    elif ':' in spw_names[0]:
        logprint("old spw names with :", logfileout='logs/initial.log')
    else:
        logprint("unknown spw names", logfileout='logs/initial.log')
    
# Check for missing scans

    missingScans = 0
    missingScanStr = ''
    
    for i in range(max(scanNums)):
        if scanNums.count(i+1) == 1: pass
        else:
            logprint ("WARNING: Scan "+str(i+1)+" is not present", logfileout='logs/initial.log')
            missingScans += 1
            missingScanStr = missingScanStr+str(i+1)+', '
    
    if (missingScans > 0):
        logprint ("WARNING: There were "+str(missingScans)+" missing scans in this MS", logfileout='logs/initial.log')
    else:
        logprint ("No missing scans found.", logfileout='logs/initial.log')
    
#PRTAN
#This program mimics the behavior of PRTAN in AIPS by plotting the locations
#of VLA antennas based on their station location, not their physical position,
#so that they are reasonably spaced on a plot.
#3/24/15 DJP

    tb.open(ms_active+"/ANTENNA")
    antname=tb.getcol("NAME")
    station=tb.getcol("STATION")
    x=np.zeros(len(antname))
    y=np.zeros(len(antname))
    
#For N antennas plot at x=0, y=Station #/4
#For W antennas plot at x=Station #/4*-cos(-45), y=Station #/4*sin(-45)
#For E antennas plot at x=Station #/4*cos(45), y=Station #/4*sin(-45)

    for i in range(len(antname)):
        if station[i].find("N") >= 0:
            y[i]=float(re.findall('\d+',station[i])[0])/4.
        if station[i].find("E") >= 0:
            x[i]=float(re.findall('\d+',station[i])[0])/4.*(np.cos(-np.pi/4.))
            y[i]=float(re.findall('\d+',station[i])[0])/4.*(np.sin(-np.pi/4.))
        if station[i].find("W") >= 0:
            x[i]=float(re.findall('\d+',station[i])[0])/4.*(-np.cos(-np.pi/4.))
            y[i]=float(re.findall('\d+',station[i])[0])/4.*(np.sin(-np.pi/4.))
    pylab.figure()
    pylab.ioff()
    pylab.plot(x,y,'w.')
    for i in range(len(antname)):
        pylab.text(x[i],y[i],antname[i])
    
#Compare list of antennas in array to those not in use, print those not in use at bottom of plot
    antlist=['ea01','ea02','ea03','ea04','ea05','ea06','ea07','ea08','ea09','ea10','ea11','ea12','ea13','ea14','ea15','ea16','ea17','ea18','ea19','ea20','ea21','ea22','ea23','ea24','ea25','ea26','ea27','ea28']

    badant=list(set(antlist)-set(antname))
    badname=''
    for i in range(len(badant)):
        badname+=badant[i]
        badname+=' '
    xb=-0.05*len(badname)
    yb=min(y)
    pylab.text(-xb,yb,badname)
    pylab.tick_params(reset=True,axis='both',which='both',top='off',labeltop='off',bottom='off',labelbottom='off',left='off',labelleft='off',right='off',labelright='off')
    pylab.savefig('antpos.png')
#END PRTAN
    if os.path.exists('plots')==False:
        os.system("mkdir plots")
    os.system("mv antpos.png plots/.")

######################################################################

# Flagging bad antennas

    logprint("Flagging user-selected bad antennas", logfileout='logs/initial.log')
    if badants=='':
        logprint ("No antenna flagging performed", logfileout='logs/initial.log')
    else:
        default('flagdata')
        vis=ms_active
        mode='manual'
        antenna=badants
        action='apply'
        flagbackup=False
        savepars=False
        async=False
        flagdata()
        clearstat()
        logprint ("Bad antenna Flagging completed", logfileout='logs/initial.log')

# Flagging Continuum spws

    logprint ("Flagging the continuum spws in the spectral line data set", logfileout='logs/initial.log')

    contspw='15~18'

    default ('flagdata')
    vis=msname
    mode ='manual'
    spw =contspw
    action ='apply'
    async=False
    flagbackup=False
    flagdata()
    
    logprint ("The continuum spws have been flagged in the spectral line data set", logfileout='logs/initial.log')



# DETERMINISTIC FLAGGING:
# TIME-BASED: online flags, shadowed data, zeroes, pointing scans, quacking
# CHANNEL-BASED: end 5% of channels of each spw, 10 end channels at
# edges of basebands
# Online flags already applied, DJP 8/26/15

    logprint ("Deterministic flagging", logfileout='logs/initial.log')

    outputflagfile = 'flagging_commands1.txt'
    syscommand='rm -rf '+outputflagfile
    os.system(syscommand)
        
# Zero flagging done on import in version 0.10 and following, but in case it wasn't.
# First do zero flagging (reason='CLIP_ZERO_ALL')
    if ms_create_flag==False:
        default('flagdata')
        vis=ms_active
        mode='clip'
        clipzeros=True
        correlation='ABS_ALL'
        action='apply'
        flagbackup=False
        savepars=False
        async=False
        outfile=outputflagfile
        myzeroflags = flagdata()
        clearstat()
        logprint ("Zero flags carried out", logfileout='logs/initial.log')
    else:
        logprint ("Zero flagging already applied on import", logfileout='logs/initial.log')
         
# Now shadow flagging
# Not needed for B configuration observations
    #default('flagdata')
    #vis=ms_active
    #mode='shadow'
    #tolerance=0.0
    #action='apply'
    #flagbackup=False
    #savepars=False
    #async=False
    #flagdata()
    #clearstat()
    #logprint ("Shadow flags carried out", logfileout='logs/initial.log')
    
#Define list of flagdata parameters to use in 'list' mode
    flagdata_list=[]
    cmdreason_list=[]


# Quack the data
    logprint ("Quack the data", logfileout='logs/initial.log')
    flagdata_list.append("mode='quack' scan=" + quack_scan_string +
        " quackinterval=" + str(1.5*int_time) + " quackmode='beg' " +
        "quackincrement=False")
    
#Write out list for use in flagdata mode 'list'
    f = open(outputflagfile, 'a')
    for line in flagdata_list:
        f.write(line+"\n")
    f.close()
    
# Apply all flags
    logprint ("Applying all flags to data", logfileout='logs/initial.log')

    default('flagdata')
    vis=ms_active
    mode='list'
    inpfile=outputflagfile
    correlation=corrstring
    action='apply'
    flagbackup=False
    savepars=True
    cmdreason=string.join(cmdreason_list, ',')
    async=False
    flagdata()
    clearstat()
    
    logprint ("Flagging completed ", logfileout='logs/initial.log')
    logprint ("Flag commands saved in file "+outputflagfile, logfileout='logs/initial.log')
    
# Placeholder comment for pre-flagging phase calibrator based on template of worst RFI frequencies
# waiting for JvG to provide template.

# Save flags
    logprint ("Saving flags", logfileout='logs/initial.log')


    default('flagmanager')
    vis=ms_active
    mode='save'
    versionname='initialflags'
    comment='Deterministic flags saved after application'
    merge='replace'
    async=False
    flagmanager()
    logprint ("Flag column saved to "+versionname, logfileout='logs/initial.log')
    
# Summary of flagging, after initial flags (for testing purposes only)
    logprint ("Summary of flags at end of initial flagging", logfileout='logs/initial.log')
    default('flagdata')
    vis=ms_active
    mode='summary'
    spw='0~14'
    correlation='RR,LL'
    spwchan=True
    spwcorr=True
    action='calculate'
    s_i=flagdata() # Save results to dictionary
    initial_flag=s_i['flagged']/s_i['total']
    logprint ("Percentage of all data flagged after initial module: "+str(initial_flag*100)+'%', logfileout='logs/initial.log')
    

######################################################################

# Assemble graphical output for this stage.
# Need to include:  links to NRAO continuum pipeline, prtan output, logs
# Anything else?

    syscommand="rm -rf *.html"
    os.system(syscommand)

    wlog = open("initial.html","w")
    wlog.write('<html>\n')
    wlog.write('<head>\n')
    wlog.write('<title>CHILES Pipeline Web Log</title>\n')
    wlog.write('</head>\n')
    wlog.write('<body>\n')
    wlog.write('<br>\n')
    wlog.write('<hr>\n')
    wlog.write('<li> Session: '+SDM_name+'\n')
    wlog.write('<li><a href="'+weblog_file+'">NRAO continuum pipeline weblog</a></li>\n')
    wlog.write('<li>Antenna positions: \n')
    wlog.write('<br><img src="./plots/antpos.png"></li>\n')
    wlog.write('<li><a href="logs/initial.log">Initial Module Log</a></li>\n')
    wlog.write('<br>\n')
    wlog.write('<br> Percentage of data flagged: '+str(initial_flag*100)+'\n')
    wlog.write('<br>')
    wlog.write('<hr>\n')
    wlog.write('</body>\n')
    wlog.write('</html>\n')
    wlog.close()


# Quit if there have been any exceptions caught:

except KeyboardInterrupt, keyboardException:
    logprint ("Keyboard Interrupt: " + str(keyboardException))
except Exception, generalException:
    logprint ("Exiting script: " + str(generalException))

time_list=runtiming('initial', 'end')
pipeline_save()