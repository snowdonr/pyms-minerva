*********************
GC-MS Raw Data Model
*********************

.. contents:: Table of Contents

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
:py:class:`GCMS_data <pyms.GCMS.Class.GCMS_data>` object.

Creating a "GCMS_data" object
================================

Reading JCAMP GC-MS data
----------------------------

.. note::This example is in `pyms-test/20a <../pyms-test/20a/20a.html>`__

The PyMS package :mod:`pyms.GCMS.IO.JCAMP` provides capabilities to read the raw
GC-MS data stored in the JCAMP-DX format.

The file ``gc01_0812_066.jdx`` (located in "data") is a GC-MS experiment
converted from Agilent ChemStation format to JCAMP format using File
Translator Pro [#ChemSW]_. This file can be loaded in Python
as follows:

.. code-block:: python

    >>> from pyms.GCMS.IO.JCAMP.Function import JCAMP_reader
    >>> jcamp_file = "data/gc01_0812_066.jdx"
    >>> data = JCAMP_reader(jcamp_file)
     -> Reading JCAMP file 'data/gc01_0812_066.jdx'
    >>>

The above command creates the object ``data`` which is an instance
of the class :py:class:`GCMS_data <pyms.GCMS.Class.GCMS_data>`.

Reading ANDI GC-MS data
------------------------

.. note:: This example is in `pyms-test/20b <../pyms-test/20b/20b.html>`__

The PyMS package :mod:`pyms.GCMS.IO.ANDI` provides capabilities to read the raw
GC-MS data stored in the ANDI-MS format.

The file ``gc01_0812_066.cdf`` is a GC-MS experiment converted to ANDI-MS
format from Agilent ChemStation (the same data as in the example 20a above).
This file can be loaded as follows:

.. code-block:: python

    >>> from pyms.GCMS.IO.ANDI.Function import ANDI_reader
    >>> ANDI_file = "data/gc01_0812_066.cdf"
    >>> data = ANDI_reader(ANDI_file)
     -> Reading netCDF file 'data/gc01_0812_066.cdf'
    >>>

The above command creates the object ``data`` which is an instance
of the class :py:class:`GCMS_data <pyms.GCMS.Class.GCMS_data>`.

A GCMS_data object
======================

Methods
---------

.. note:: The examples below are in `pyms-test/20a <../pyms-test/20a/20a.html>`__ and `pyms-test/20b <../pyms-test/20b/20b.html>`__

The object ``data`` (from the two previous examples) stores the raw data as a
:py:class:`GCMS_data <pyms.GCMS.Class.GCMS_data>` object. Within the
:py:class:`GCMS_data <pyms.GCMS.Class.GCMS_data>` object, raw data are stored as a list
of :py:class:`Scan <pyms.Scan.Scan>` objects and a list of retention times.
There are several methods available to access data and attributes of the
:py:class:`GCMS_data <pyms.GCMS.Class.GCMS_data>` and
:py:class:`Scan <pyms.Scan.Scan>` objects.

The :py:class:`GCMS_data <pyms.GCMS.Class.GCMS_data>` object's methods relate to the raw data. The main properties
relate to the masses, retention times and scans. For example, the
minimum and maximum mass from all of the raw data can be returned by the
following:

    >>> data.min_mass
    >>> data.max_mass


A list of all retention times can be returned by:

    >>> time = data.time_list

The index of a specific retention time (in seconds) can be returned by:

    >>> data.get_index_at_time(400.0)

Note that this returns the index of the retention time in the
data closest to the given retention time of 400.0 seconds.

The :py:attr:`GCMS_data.tic <pyms.GCMS.Class.GCMS_data.tic>` attribute
returns a total ion chromatogram (TIC) of the data
as an :py:class:`IonChromatogram <pyms.IonChromatogram.IonChromatogram>` object:


    >>> tic = data.tic

The :py:class:`IonChromatogram <pyms.IonChromatogram.IonChromatogram>`
object is explained in a later chapter.

A Scan data object
----------------------

A Scan object contains a list of masses and a corresponding list of intensity
values from a single mass-spectrum scan in the raw data. Typically only
non-zero (or non-threshold) intensities and corresponding masses are stored in
the raw data.

.. note:: The following examples are the same in `pyms-test/20a <../pyms-test/20a/20a.html>`__ and `pyms-test/20b <../pyms-test/20b/20b.html>`__

A list of all the raw Scan objects can be returned by:

    >>> scans = data.scan_list

A list of all masses in a scan (e.g. the 1st scan) is returned by:

    >>> scans[0].mass_list

A list of all corresponding intensities in a scan is returned by:

    >>> scans[0].intensity_list

The minimum and maximum mass in an individual scan (e.g. the 1st scan) are
returned by:

    >>> scans[0].min_mass
    >>> scans[0].max_mass

Exporting data and obtaining information about a data set
----------------------------------------------------------

.. note:: This example is in `pyms-test/20c <../pyms-test/20c/20c.html>`__

Often it is of interest to find out some basic information about the
data set, e.g. the number of scans, the retention time range, and
m/z range and so on. The :py:class:`GCMS_data <pyms.GCMS.Class.GCMS_data>`
class provides a method :py:meth:`info() <pyms.GCMS.Class.GCMS_data.info()>`
that can be used for this purpose.

.. code-block:: python

    >>> from pyms.GCMS.IO.ANDI.Function import ANDI_reader
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

.. note:: This example is in `pyms-test/20d <../pyms-test/20d/20d.html>`__

Occasionally it is useful to compare two data sets. For example,
one may want to check the consistency between the data set
exported in netCDF format from the manufacturer's software, and
the JCAMP format exported from a third party software.

For example:

.. code-block:: python

    >>> from pyms.GCMS.IO.JCAMP.Function import JCAMP_reader
    >>> from pyms.GCMS.IO.ANDI.Function import ANDI_reader
    >>> andi_file = "data/gc01_0812_066.cdf"
    >>> jcamp_file = "data/gc01_0812_066.jdx"
    >>> data1 = ANDI_reader(andi_file)
     -> Reading netCDF file 'data/gc01_0812_066.cdf'
    >>> data2 = JCAMP_reader(jcamp_file)
     -> Reading JCAMP file 'data/gc01_0812_066.jdx'

To compare the two data sets:

.. code-block:: python

    >>> from pyms.GCMS.Function import diff
    >>> diff(data1,data2)
     Data sets have the same number of time points.
       Time RMSD: 1.80e-13
     Checking for consistency in scan lengths ... OK
     Calculating maximum RMSD for m/z values and intensities ...
       Max m/z RMSD: 1.03e-05
       Max intensity RMSD: 0.00e+00

If the data is not possible to compare, for example because of
different number of scans, or inconsistent number of m/z values
in between two scans, :py:meth:`diff() <pyms.GCMS.Function.diff>`
will report the difference. For example:

.. code-block:: python

    >>> data2.trim(begin=1000,end=2000)
    Trimming data to between 1000 and 2000 scans
    >>> diff(data1,data2)
     -> The number of retention time points different.
     First data set: 9865 time points
     Second data set: 1001 time points
     Data sets are different.

.. rubric:: Footnotes

.. [#ANDI-MS] ANDI-MS was developed by the Analytical Instrument Association
.. [#JCAMP-DX] JCAMP-DX is maintained by the International Union of Pure and Applied Chemistry
.. [#ChemSW] ChemSW, Inc.