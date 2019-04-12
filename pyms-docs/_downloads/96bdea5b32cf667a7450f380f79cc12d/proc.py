"""proc.py

This example uses the matplotlib library for plotting in Python

To use the Display class as shown in this script, the matplotlib
library must be installed. To check for this type "import matplotlib" 
in an interactive python session. 

The authors suggest the use of the TKAgg backend, which can be set in matplotlib
by changing the backend parameter in your matplotlibrc file.

The location of your matplotlibrc file can be found by typing
>>>import matplotlib
>>>matplotlib.matplotlib_fname()

in an interactive session
"""


from pyms.GCMS.IO.ANDI.Function import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix


from pyms.Display.Function import plot_ic

# read the raw data as a GCMS_data object
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)


# IntensityMatrix
# must build intensity matrix before accessing any intensity matrix methods.

# default, float masses with interval (bin interval) of one from min mass
print("default intensity matrix, bin interval = 1, boundary +/- 0.5")
im = build_intensity_matrix(data)


#
# IonChromatogram
#

# TIC from raw data
tic = data.get_tic()
# save TIC to a file

#Call function to store a plot of the TIC
plot_ic(tic, line_label = 'TIC', plot_title = 'TIC for gc01_0812_066')


