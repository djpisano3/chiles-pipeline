# CHILES_pipe_split.py
# This module is run after the user is satisfied with the results of the CHILES 
# pipeline.
# 9/21/16 DJP

logprint ("Starting CHILES_pipe_split.py", logfileout='logs/split.log')
time_list=runtiming('split', 'start')

# If running this module right after starting casa, must run CHILES_pipe_restore first.
# Step 1, load needed functions.
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

# Set up to save results
if os.path.exists("FINAL")==False:
    os.system("mkdir FINAL")

logprint ('Split calibrated deepfield uv data', logfileout='logs/split.log')

# Split with averaging by 2x in time & 4x in frequency
outputms=ms_active[:-3]+'_calibrated_deepfield.ms'
targetfile=outputms
    
default('split')
vis=ms_active
datacolumn='corrected'
outputvis=targetfile
field='deepfield'
spw ='0~14'
width=4
timebin='16s'
split()

os.system("mv "+targetfile+" FINAL/.")

# Save calibration tables
if os.path.exists('antposcal.p')==True:
    os.system("cp antposcal.p FINAL/.")

os.system("cp gain_curves.g FINAL/.")
os.system("cp finaldelay.k FINAL/.")
os.system("cp finalBPcal.b FINAL/.")
os.system("cp finalphase_scan.gcal FINAL/.")
os.system("cp finalamp.gcal FINAL/.")
os.system("cp finalflux.gcal FINAL/.")

# Save flagging commands

os.system("cp manual_flags_*.txt FINAL/.")

# Save Flag tables

os.system("cp "+ms_active+".flagversions/flags.finalflags FINAL/.")
os.system("cp "+ms_active+".flagversions/flags.finalfinalflags FINAL/.")



logprint ("Finished CHILES_pipe_split.py", logfileout='logs/split.log')
time_list=runtiming('split', 'end')

pipeline_save()
