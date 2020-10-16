Example: Comparing two GC-MS data sets
======================================

Occasionally it is useful to compare two data sets. For example, one may
want to check the consistency between the data set exported in netCDF
format from the manufacturer’s software, and the JCAMP format exported
from a third party software.

First, setup the paths to the datafiles and the output directory, then
import JCAMP_reader and AMDI_reader.

.. nbinput:: ipython3
    :execution-count: 1

    import pathlib
    data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
    # Change this if the data files are stored in a different location

    output_directory = pathlib.Path(".").resolve() / "output"

    from pyms.GCMS.IO.JCAMP import JCAMP_reader
    from pyms.GCMS.IO.ANDI import ANDI_reader

Then the raw data is read as before.

.. nbinput:: ipython3
    :execution-count: 2

    andi_file = data_directory / "gc01_0812_066.cdf"
    data1 = ANDI_reader(andi_file)
    data1


.. parsed-literal::

     -> Reading netCDF file '/home/vagrant/PyMassSpec/pyms-data/gc01_0812_066.cdf'




.. parsed-literal::

    GCMS_data(rt range 305.582 - 4007.721, time_step 0.37531822789943226, length 9865)



.. nbinput:: ipython3
    :execution-count: 3

    jcamp_file = data_directory / "gc01_0812_066.jdx"
    data2 = JCAMP_reader(jcamp_file)
    data2


.. parsed-literal::

     -> Reading JCAMP file '/home/vagrant/PyMassSpec/pyms-data/gc01_0812_066.jdx'




.. parsed-literal::

    GCMS_data(rt range 305.582 - 4007.722, time_step 0.3753183292781833, length 9865)



To compare the two data sets, use the function :py:meth:`diff() <pyms.GCMS.Function.diff()>`

.. nbinput:: ipython3
    :execution-count: 4

    from pyms.GCMS.Function import diff

    diff(data1, data2)



.. parsed-literal::

     Data sets have the same number of time points.
       Time RMSD: 3.54e-04
     Checking for consistency in scan lengths ...OK
     Calculating maximum RMSD for m/z values and intensities ...
       Max m/z RMSD: 1.03e-05
       Max intensity RMSD: 0.00e+00


If the data cannot be compared, for example because of different number
of scans, or inconsistent number of m/z values in between two scans,
:py:meth:`diff() <pyms.GCMS.Function.diff()>` will report the difference. For example:

.. nbinput:: ipython3
    :execution-count: 5

    data2.trim(begin=1000, end=2000)



.. parsed-literal::

    Trimming data to between 1000 and 2001 scans


.. nbinput:: ipython3
    :execution-count: 6

    diff(data1, data2)




.. parsed-literal::

     The number of retention time points differ.
    	First data set: 9865 time points
    	Second data set: 1002 time points
     Data sets are different.
