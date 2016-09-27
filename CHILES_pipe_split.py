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

# Test if imaging is quicker after  averaging by 2x in time & 4x in frequency
outputms=ms_active[:-3]+'_calibrated_deepfield.ms'
targetfile=outputms

logprint ('Split calibrated deepfield uv data', logfileout='logs/split.log')

    
default('split')
vis=ms_active
datacolumn='corrected'
outputvis=targetfile
field='deepfield'
spw ='0~14'
width=4
timebin='16s'
split()

# Save calibration tables

# Save Flag tables

# Save flagging commands


logprint ("Finished CHILES_pipe_split.py", logfileout='logs/split.log')
time_list=runtiming('split', 'end')

pipeline_save()
