Example: Reading JCAMP GC-MS data
=================================

The PyMS package :mod:`pyms.GCMS.IO.JCAMP` provides capabilities to read
the raw GC-MS data stored in the JCAMP-DX format.

First, setup the paths to the datafile and the output directory, then
import JCAMP_reader.

.. nbinput:: ipython3
    :execution-count: 1

    import pathlib
    data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
    # Change this if the data files are stored in a different location

    output_directory = pathlib.Path(".").resolve() / "output"

    from pyms.GCMS.IO.JCAMP import JCAMP_reader

Read the raw JCAMP-dx data.

.. nbinput:: ipython3
    :execution-count: 2

    jcamp_file = data_directory / "gc01_0812_066.jdx"
    data = JCAMP_reader(jcamp_file)
    data


.. parsed-literal::

     -> Reading JCAMP file '/home/vagrant/PyMassSpec/pyms-data/gc01_0812_066.jdx'




.. parsed-literal::

    <GCMS_data(305.582 - 4007.722 seconds, time step 0.3753183292781833, 9865 scans)>



A GCMS_data Object
~~~~~~~~~~~~~~~~~~

The object ``data`` (from the two previous examples) stores the raw data
as a :class:`pyms.GCMS.Class.GCMS_data` object. Within the :class:`~pyms.GCMS.Class.GCMS_data`
object, raw data are stored as a list of :class:`pyms.Spectrum.Scan` objects
and a list of retention times. There are several methods available to
access data and attributes of the :class:`~pyms.GCMS.Class.GCMS_data` and :class:`~pyms.Spectrum.Scan` objects.

The :class:`~pyms.GCMS.Class.GCMS_data` object’s methods relate to the raw data. The main
properties relate to the masses, retention times and scans. For example,
the minimum and maximum mass from all of the raw data can be returned by
the following:

.. nbinput:: ipython3
    :execution-count: 3

    data.min_mass





.. parsed-literal::

    50.0



.. nbinput:: ipython3
    :execution-count: 4

    data.max_mass




.. parsed-literal::

    599.9



A list of the first 10 retention times can be returned with:

.. nbinput:: ipython3
    :execution-count: 5

    data.time_list[:10]




.. parsed-literal::

    [305.582,
     305.958,
     306.333,
     306.708,
     307.084,
     307.459,
     307.834,
     308.21,
     308.585,
     308.96]



The index of a specific retention time (in seconds) can be returned
with:

.. nbinput:: ipython3
    :execution-count: 6

    data.get_index_at_time(400.0)




.. parsed-literal::

    252



Note that this returns the index of the retention time in the data
closest to the given retention time of 400.0 seconds.

The :attr:`GCMS_data.tic <pyms.GCMS.Class.GCMS_data.tic>` attribute returns a total ion chromatogram (TIC)
of the data as an :class:`~pyms.IonChromatogram.IonChromatogram` object:

.. nbinput:: ipython3
    :execution-count: 7

    data.tic




.. parsed-literal::

    <pyms.IonChromatogram.IonChromatogram at 0x7f6b22ff9d68>



The :class:`~pyms.IonChromatogram.IonChromatogram` object is explained in a later example.

A Scan Object
~~~~~~~~~~~~~

A :class:`pyms.Spectrum.Scan` object contains a list of masses and a
corresponding list of intensity values from a single mass-spectrum scan
in the raw data. Typically only non-zero (or non-threshold) intensities
and corresponding masses are stored in the raw data.

A list of the first 10 :class:`pyms.Spectrum.Scan` objects can be returned
with:

.. nbinput:: ipython3
    :execution-count: 8

    scans = data.scan_list
    scans[:10]




.. parsed-literal::

    [<pyms.Spectrum.Scan at 0x7f6b4117a518>,
     <pyms.Spectrum.Scan at 0x7f6b22ff9400>,
     <pyms.Spectrum.Scan at 0x7f6b22ff9dd8>,
     <pyms.Spectrum.Scan at 0x7f6b22ff9e80>,
     <pyms.Spectrum.Scan at 0x7f6b22ff9f28>,
     <pyms.Spectrum.Scan at 0x7f6b22ff9fd0>,
     <pyms.Spectrum.Scan at 0x7f6b22ff9e48>,
     <pyms.Spectrum.Scan at 0x7f6b22ff9668>,
     <pyms.Spectrum.Scan at 0x7f6b22ff9d30>,
     <pyms.Spectrum.Scan at 0x7f6b22ff9cf8>]



A list of the first 10 masses in a scan (e.g. the 1st scan) is returned
with:

.. nbinput:: ipython3
    :execution-count: 9

    scans[0].mass_list[:10]




.. parsed-literal::

    [50.1, 51.1, 53.1, 54.2, 55.1, 56.2, 57.2, 58.2, 59.1, 60.1]



A list of the first 10 corresponding intensities in a scan is returned
with:

.. nbinput:: ipython3
    :execution-count: 10

    scans[0].intensity_list[:10]




.. parsed-literal::

    [22128.0,
     10221.0,
     31400.0,
     27352.0,
     65688.0,
     55416.0,
     75192.0,
     112688.0,
     152256.0,
     21896.0]



The minimum and maximum mass in an individual scan (e.g. the 1st scan)
are returned with:

.. nbinput:: ipython3
    :execution-count: 11

    scans[0].min_mass




.. parsed-literal::

    50.1



.. nbinput:: ipython3
    :execution-count: 12

    scans[0].max_mass




.. parsed-literal::

    599.4



Exporting data and obtaining information about a data set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often it is of interest to find out some basic information about the
data set, e.g. the number of scans, the retention time range, and m/z
range and so on. The :class:`~pyms.GCMS.Class.GCMS_data` class provides a method :py:meth:`info() <pyms.GCMS.Class.GCMS_data.info()>`
that can be used for this purpose.

.. nbinput:: ipython3
    :execution-count: 13

    data.info()


.. parsed-literal::

     Data retention time range: 5.093 min -- 66.795 min
     Time step: 0.375 s (std=0.000 s)
     Number of scans: 9865
     Minimum m/z measured: 50.000
     Maximum m/z measured: 599.900
     Mean number of m/z values per scan: 56
     Median number of m/z values per scan: 40


The entire raw data of a :class:`~pyms.GCMS.Class.GCMS_data` object can be exported to a file
with the method :py:meth:`write() <pyms.GCMS.Class.GCMS_data.write()>`:

.. nbinput:: ipython3
    :execution-count: 14

    data.write(output_directory / "data")


.. parsed-literal::

     -> Writing intensities to '/home/vagrant/PyMassSpec/pyms-demo/jupyter/output/data.I.csv'
     -> Writing m/z values to '/home/vagrant/PyMassSpec/pyms-demo/jupyter/output/data.mz.csv'


This method takes the filename (“output/data”, in this example) and
writes two CSV files. One has extension “.I.csv” and contains the
intensities (“output/data.I.csv” in this example), and the other has the
extension “.mz” and contains the corresponding table of m/z value
(“output/data.mz.csv” in this example). In general, these are not
two-dimensional matrices, because different scans may have different
number of m/z values recorded.
