# CHILES_pipe_split.py
# This module is run after the user is satisfied with the results of the CHILES 
# pipeline.
# 9/21/16 DJP
# 1/10/17 DJP: Include html, logs, and plots in FINAL directory

logprint ("Starting CHILES_pipe_split.py", logfileout='logs/split.log')
time_list=runtiming('split', 'start')

# If running this module right after starting casa, must run CHILES_pipe_restore first.
# Step 1, load needed functions.
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

# Set up directory to save results
if os.path.exists("FINAL")==False:
    os.system("mkdir FINAL")

logprint ('Split calibrated deepfield uv data', logfileout='logs/split.log')

# Do final flag summary and compile statistics.
# Summary of flagging, after final flagging (for testing purposes only)
logprint ("Summary of flags after flagging target", logfileout='logs/split.log')
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

logprint ("Final percentage of all data flagged: "+str(target_flag*100)+'%', logfileout='logs/split.log')

# Save final version of flags
logprint("Saving flags before split",logfileout='logs/target.log')

default('flagmanager')
vis=ms_active
mode='save'
versionname='finalflags'
comment='Final flags saved before split'
merge='replace'
flagmanager()
logprint ("Flag column saved to "+versionname, logfileout='logs/split.log')


# Step 3: Diagnostic Plots
logprint("Making diagnostic plots",logfileout='logs/split.log')

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

# Make updated plot of fraction of data flagged.
fig=pylab.figure()
pylab.plot(freq,flag_frac,'k-')
pylab.xlim(940.,1445.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Fraction of data flagged')
pylab.savefig("finaltarget_flag.png")
pylab.close(fig)


# Write flagging statistics to file
flagstatname=msname.rstrip('ms') + 'flag_stats.txt'

np.savetxt(flagstatname,np.column_stack((freq,flagged,totals,flag_frac)),fmt='%s',header='freq [MHz] #flags #total %flagged')


# Split with averaging by 2x in time & 4x in frequency
outputms=ms_active[:-3]+'_calibrated_deepfield.ms'
targetfile=outputms

# Delete final MS if already present
if os.path.exists(targetfile):
    os.system("rm -rf "+targetfile)
    os.system("rm -rf FINAL/"+targetfile)
        
default('split')
vis=ms_active
datacolumn='corrected'
outputvis=targetfile
field='deepfield'
spw ='0~14'
width=4
timebin='16s'
split()

os.system("mv "+targetfile+" FINAL/")

# Save calibration tables
if os.path.exists('antposcal.p')==True:
    os.system("cp -r antposcal.p FINAL/.")

os.system("cp -r gain_curves.g FINAL/.")
os.system("cp -r finaldelay.k FINAL/.")
os.system("cp -r finalBPcal.b FINAL/.")
os.system("cp -r finalphase_scan.gcal FINAL/.")
os.system("cp -r finalamp.gcal FINAL/.")
os.system("cp -r finalflux.gcal FINAL/.")

# Save flagging commands & statistics

os.system("cp manual_flags_*.txt FINAL/.")
os.system("cp *flag_stats.txt FINAL/.")

# Save Flag tables

os.system("cp -r "+ms_active+".flagversions/flags.finalflags FINAL/.")
os.system("cp -r "+ms_active+".flagversions/flags.finalfinalflags FINAL/.")

#Move plots, images to sub-directory

os.system("mv *.png plots")

#Create webpage with results

if os.path.exists("split.html"):
    os.system("rm split.html")
wlog = open("split.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_split results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/split.log">Split Log</a></li>\n')
wlog.write('<li> Percentage of data flagged as a function of frequency: \n')
wlog.write('<br><img src="plots/finaltarget_flag.png">\n')
wlog.write('<br>')
wlog.write('<li><a href="logs/timing.log">Timing Log</a></li>\n')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()

# Copy html files and plots to FINAL directory
os.system("cp -r *html FINAL/.")
os.system("cp -r plots FINAL/.")
# Save variable values
os.system("cp -r CHILES_pipe_restore.list FINAL/.")

logprint ("Finished CHILES_pipe_split.py", logfileout='logs/split.log')
time_list=runtiming('split', 'end')

pipeline_save()

# Copy logs to FINAL directory (needs to be after final "save" to preserve all information
os.system("cp -r logs FINAL/.")
