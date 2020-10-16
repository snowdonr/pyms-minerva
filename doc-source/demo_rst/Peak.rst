Example: Peak Objects
=====================

Fundamental to GC-MS analysis is the identification of individual
components of the sample mix. The basic component unit is represented as
a signal peak. In PyMassSpec a signal peak is represented as :class:`~pyms.Peak.Class.Peak`
object. PyMassSpec provides functions to detect peaks and create peaks
(discussed at the end of the chapter).

A peak object stores a minimal set of information about a signal peak,
namely, the retention time at which the peak apex occurs and the mass
spectra at the apex. Additional information, such as, peak width, TIC
and individual ion areas can be filtered from the GC-MS data and added
to the Peak object information.

Creating a Peak Object
----------------------

A peak object can be created for a scan at a given retention time by
providing the retention time (in minutes or seconds) and the
:class:`~pyms.Spectrum.MassSpectrum` object of the scan.

First, setup the paths to the datafiles and the output directory, then
import JCAMP_reader.

.. nbinput:: ipython3
    :execution-count: 1

    import pathlib
    data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
    # Change this if the data files are stored in a different location

    output_directory = pathlib.Path(".").resolve() / "output"

    from pyms.GCMS.IO.JCAMP import JCAMP_reader

Read the raw data files.

.. nbinput:: ipython3
    :execution-count: 2

    jcamp_file = data_directory / "gc01_0812_066.jdx"
    data = JCAMP_reader(jcamp_file)


.. parsed-literal::

     -> Reading JCAMP file '/home/vagrant/PyMassSpec/pyms-data/gc01_0812_066.jdx'


Build the :class:`~pyms.IntensityMatrix.IntensityMatrix`.

.. nbinput:: ipython3
    :execution-count: 3

    from pyms.IntensityMatrix import build_intensity_matrix_i

    im = build_intensity_matrix_i(data)

Extract the :class:`~pyms.Spectrum.MassSpectrum` at 31.17 minutes in this example.

.. nbinput:: ipython3
    :execution-count: 4

    index = im.get_index_at_time(31.17*60.0)
    ms = im.get_ms_at_index(index)

Create a :class:`~pyms.Peak.Class.Peak` object for the given retention time.

.. nbinput:: ipython3
    :execution-count: 5

    from pyms.Peak.Class import Peak
    peak = Peak(31.17, ms, minutes=True)


By default the retention time is assumed to be in seconds. The parameter
``minutes`` can be set to ``True`` if the retention time is given in
minutes. Internally, PyMassSpec stores retention times in seconds, so
the ``minutes`` parameter ensures the input and output of the retention
time are in the same units.

Peak Object properties
----------------------

The retention time of the peak, in seconds, can be returned with
:attr:`pyms.Peak.Class.Peak.rt`. The mass spectrum can be returned with
:attr:`pyms.Peak.Class.Peak.mass_spectrum`.

The :class:`~pyms.Peak.Class.Peak` object constructs a unique identification (UID) based on
the spectrum and retention time. This helps in managing lists of peaks
(covered in the next chapter). The UID can be returned with
:attr:`pyms.Peak.Class.Peak.UID`. The format of the UID is the masses of the
two most abundant ions in the spectrum, the ratio of the abundances of
the two ions, and the retention time (in the same units as given when
the Peak object was created). The format is:

Mass1-Mass2-Ratio-RT

For example:

.. nbinput:: ipython3
    :execution-count: 6

    peak.rt




.. parsed-literal::

    1870.2



.. nbinput:: ipython3
    :execution-count: 7

    peak.UID




.. parsed-literal::

    '319-73-74-1870.20'



.. nbinput:: ipython3
    :execution-count: 8

    index = im.get_index_of_mass(73.3)

    index




.. parsed-literal::

    23



Modifying a Peak Object
-----------------------

The Peak object has methods for modifying the mass spectrum. The mass
range can be cropped to a smaller range with :meth:`crop_mass() <pyms.Peak.Class.Peak.crop_mass>`, and the
intensity values for a single ion can be set to zero with
:meth:`null_mass() <pyms.Peak.Class.Peak.null_mass>`. For example, the mass range can be set from 60 to 450
:math:`m/z`, and the ions related to sample preparation can be ignored by
setting their intensities to zero as follows:

.. nbinput:: ipython3
    :execution-count: 9

    peak.crop_mass(60, 450)
    peak.null_mass(73)
    peak.null_mass(147)

The UID is automatically updated to reflect the changes:

.. nbinput:: ipython3
    :execution-count: 10

    peak.UID




.. parsed-literal::

    '319-205-54-1870.20'



It is also possible to change the peak mass spectrum by setting the
attribute :attr:`pyms.Peak.Class.Peak.mass_spectrum`.
