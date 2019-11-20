*********************
GC-MS Raw Data Model
*********************

.. contents:: :local: Table of Contents

Introduction
=============

PyMassSpec can read gas chromatography-mass spectrometry (GC-MS) data stored in
Analytical Data Interchange for Mass Spectrometry (ANDI-MS), [#ANDI-MS]_
and Joint Committee on Atomic and Molecular Physical Data (JCAMP-DX) [#JCAMP-DX]_
formats. The information contained in the data files can vary significantly
depending on the instrument, vendor's software, or conversion utility.
PyMassSpec makes the following assumptions about the information contained in the data file:

* The data contain the m/z and intensity value pairs across a scan.
* Each scan has a retention time.

Internally, PyMassSpec stores the raw data from ANDI files or JCAMP files as a
:class:`~pyms.GCMS.Class.GCMS_data` object.

Creating a "GCMS_data" object
================================

Reading JCAMP GC-MS data
----------------------------

.. note::This example is in :ref:`pyms-demo/20a <demo-20a>`

The PyMS package :mod:`pyms.GCMS.IO.JCAMP` provides capabilities to read the raw
GC-MS data stored in the JCAMP-DX format.

The file ``gc01_0812_066.jdx`` (located in "data") is a GC-MS experiment
converted from Agilent ChemStation format to JCAMP format using File
Translator Pro [#ChemSW]_. This file can be loaded in Python
as follows:

.. code-block:: python

    >>> from pyms.GCMS.IO.JCAMP import JCAMP_reader
    >>> jcamp_file = "data/gc01_0812_066.jdx"
    >>> data = JCAMP_reader(jcamp_file)
     -> Reading JCAMP file 'data/gc01_0812_066.jdx'
    >>> data
    <pyms.GCMS.Class.GCMS_data at 0x7f3ec77da0b8>

The above command creates the object ``data`` which is an instance
of the class :class:`~pyms.GCMS.Class.GCMS_data`.

Reading ANDI GC-MS data
------------------------

.. note:: This example is in :ref:`pyms-demo/20b <demo-20b>`

The PyMS package :mod:`pyms.GCMS.IO.ANDI` provides capabilities to read the raw
GC-MS data stored in the ANDI-MS format.

The file ``gc01_0812_066.cdf`` is a GC-MS experiment converted to ANDI-MS
format from Agilent ChemStation (the same data as in the example 20a above).
This file can be loaded as follows:

.. code-block:: python

    >>> from pyms.GCMS.IO.ANDI import ANDI_reader
    >>> ANDI_file = "data/gc01_0812_066.cdf"
    >>> data = ANDI_reader(ANDI_file)
     -> Reading netCDF file 'data/gc01_0812_066.cdf'
    >>> data
    <pyms.GCMS.Class.GCMS_data at 0x7f3ec77da0b8>

The above command creates the object ``data`` which is an instance
of the class :class:`~pyms.GCMS.Class.GCMS_data`.

A GCMS_data object
======================

Methods
---------

.. note:: The examples below are in :ref:`pyms-demo/20a <demo-20a>` and :ref:`pyms-demo/20b <demo-20b>`

The object ``data`` (from the two previous examples) stores the raw data as a
:class:`~pyms.GCMS.Class.GCMS_data` object. Within the
:class:`~pyms.GCMS.Class.GCMS_data` object, raw data are stored as a list
of :class:`~pyms.Spectrum.Scan` objects and a list of retention times.
There are several methods available to access data and attributes of the
:class:`~pyms.GCMS.Class.GCMS_data` and
:class:`~pyms.Spectrum.Scan` objects.

The :class:`~pyms.GCMS.Class.GCMS_data` object's methods relate to the raw data. The main properties
relate to the masses, retention times and scans. For example, the
minimum and maximum mass from all of the raw data can be returned by the
following:

    >>> data.min_mass
    50.5
    >>> data.max_mass
    599.9000244140625

A list of all retention times can be returned by:

    >>> data.time_list
    [305.582, 305.958, 306.333, 306.70799999999997, 307.084, ...]

The index of a specific retention time (in seconds) can be returned by:

    >>> data.get_index_at_time(400.0)
    252

Note that this returns the index of the retention time in the
data closest to the given retention time of 400.0 seconds.

The :attr:`GCMS_data.tic <pyms.GCMS.Class.GCMS_data.tic>` attribute
returns a total ion chromatogram (TIC) of the data
as an :class:`~pyms.IonChromatogram.IonChromatogram` object:


    >>> data.tic
    <pyms.IonChromatogram.IonChromatogram at 0x7f99a27bd320>

The :class:`~pyms.IonChromatogram.IonChromatogram`
object is explained in a later chapter.

A Scan data object
----------------------

A :class:`~pyms.Spectrum.Scan`object contains a list of masses and a corresponding list of intensity values from a single mass-spectrum scan in the raw data. Typically only non-zero (or non-threshold) intensities and corresponding masses are stored in the raw data.

.. note:: The following examples are the same in :ref:`pyms-demo/20a <demo-20a>` and :ref:`pyms-demo/20b <demo-20b>`

A list of all the raw :class:`~pyms.Spectrum.Scan` objects can be returned with:

    >>> scans = data.scan_list
    >>> scans
    [<pyms.Spectrum.Scan at 0x7f999d2209b0>, <pyms.Spectrum.Scan at 0x7f999d220828>, <pyms.Spectrum.Scan at 0x7f999d220390>, ...]

A list of all masses in a scan (e.g. the 1st scan) is returned with:

    >>> scans[0].mass_list
    [50.099998474121094, 51.099998474121094, 53.099998474121094, ...]

A list of all corresponding intensities in a scan is returned with:

    >>> scans[0].intensity_list
    [22128.0, 10221.0, 31400.0, 27352.0, 65688.0, ...]

The minimum and maximum mass in an individual scan (e.g. the 1st scan) are
returned with:

    >>> scans[0].min_mass
    50.099998474121094
    >>> scans[0].max_mass
    599.4000244140625

Exporting data and obtaining information about a data set
----------------------------------------------------------

.. note:: This example is in :ref:`pyms-demo/20c <demo-20c>`

Often it is of interest to find out some basic information about the
data set, e.g. the number of scans, the retention time range, and
m/z range and so on. The :class:`~pyms.GCMS.Class.GCMS_data`
class provides a method :py:meth:`info() <pyms.GCMS.Class.GCMS_data.info()>`
that can be used for this purpose.

.. code-block:: python

    >>> from pyms.GCMS.IO.ANDI import ANDI_reader
    >>> andi_file = "data/gc01_0812_066.cdf"
    >>> data = ANDI_reader(andi_file)
     -> Reading netCDF file 'data/gc01_0812_066.cdf'
    >>> data.info()
     Data retention time range: 5.093 min -- 66.795 min
     Time step: 0.375 s (std=0.000 s)
     Number of scans: 9865
     Minimum m/z measured: 50.000
     Maximum m/z measured: 599.900
     Mean number of m/z values per scan: 56
     Median number of m/z values per scan: 40
    >>>

The entire raw data can be exported to a file with the method
:py:meth:`write() <pyms.GCMS.Class.GCMS_data.write()>` :

.. code-block:: python

    >>> data.write("output/data")
     -> Writing intensities to 'output/data.I.csv'
     -> Writing m/z values to 'output/data.mz.csv'

This method takes the string ("output/data", in this example)
and writes two CSV files. One has extension ".I.csv" and
contains the intensities ("output/data.I.csv" in this example),
and the other has the extension ".mz" and contains the
corresponding table of m/z value ("output/data.mz.csv" in
this example). In general, these are not two-dimensional matrices,
because different scans may have different number of m/z
values recorded.

Comparing two GC-MS data sets
----------------------------------

.. note:: This example is in :ref:`pyms-demo/20d <demo-20d>`

Occasionally it is useful to compare two data sets. For example,
one may want to check the consistency between the data set
exported in netCDF format from the manufacturer's software, and
the JCAMP format exported from a third party software.

For example:

.. code-block:: python

    >>> from pyms.GCMS.IO.JCAMP import JCAMP_reader
    >>> from pyms.GCMS.IO.ANDI import ANDI_reader
    >>> andi_file = "data/gc01_0812_066.cdf"
    >>> jcamp_file = "data/gc01_0812_066.jdx"
    >>> data1 = ANDI_reader(andi_file)
     -> Reading netCDF file 'data/gc01_0812_066.cdf'
    >>> data2 = JCAMP_reader(jcamp_file)
     -> Reading JCAMP file 'data/gc01_0812_066.jdx'

To compare the two data sets:

.. code-block:: python

    >>> from pyms.GCMS.Function import diff
    >>> diff(data1, data2)
     Data sets have the same number of time points.
       Time RMSD: 1.80e-13
     Checking for consistency in scan lengths ... OK
     Calculating maximum RMSD for m/z values and intensities ...
       Max m/z RMSD: 1.03e-05
       Max intensity RMSD: 0.00e+00

If the data cannot be compared, for example because of
different number of scans, or inconsistent number of m/z values
in between two scans, :py:meth:`diff() <pyms.GCMS.Function.diff>`
will report the difference. For example:

.. code-block:: python

    >>> data2.trim(begin=1000, end=2000)
    Trimming data to between 1000 and 2000 scans
    >>> diff(data1, data2)
     -> The number of retention time points different.
     First data set: 9865 time points
     Second data set: 1001 time points
     Data sets are different.

.. rubric:: Footnotes

.. [#ANDI-MS] ANDI-MS was developed by the Analytical Instrument Association
.. [#JCAMP-DX] JCAMP-DX is maintained by the International Union of Pure and Applied Chemistry
.. [#ChemSW] ChemSW, Inc.