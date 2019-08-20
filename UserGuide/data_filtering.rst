****************
Data filtering
****************

.. contents:: Table of Contents

Introduction
=============

In this chapter filtering techniques that allow pre-processing of GC-MS data for analysis and comparison to other pre-processed GC-MS data are covered.

Time strings
==============

Before considering the filtering techniques, the mechanism for representing retention times is outlined here.

A time string is the specification of a time interval, that takes the format ``NUMBERs`` or ``NUMBERm`` for time interval in seconds or minutes. For example, these are valid time strings: ``10s`` (10 seconds) and ``0.2m`` (0.2 minutes).

Intensity Matrix resizing
==============================

Once an IntensityMatrix has been constructed from the raw GC-MS data, the entries of the matrix can be modified. These modifications can operate on the entire matrix, or individual masses or scans.

Retention time range
-----------------------

.. note:: This example is in :ref:`demo-40a`

A basic operation on the GC-MS data is to select a specific time range for processing. In |pkgname|, any data outside the chosen time range is discarded. The :meth:`trim() <pyms.GCMS.Class.GCMS_data.trim>` method operates on the raw data, so any subsequent processing only refers to the trimmed data.

Given a previously loaded raw GC-MS data file, ``data``, the data can be
trimmed to specific scans;

    >>> data.trim(1000, 2000)
    >>> data.info()

or specific retention times (in ``seconds`` or ``minutes``);

    >>> data.trim("6.5m", "21m")
    >>> data.info()

Mass spectrum range and entries
---------------------------------

.. note:: This example is in `pyms-demo/40b <../pyms-demo/40b/40b.html>`__

An :class:`~pyms.IntensityMatrix.IntensityMatrix` object has a set mass range and interval that is derived
from the data at the time of building the intensity matrix. The range of mass
values can be cropped. This is done, primarily, to ensure that the range of
masses used are consistent when comparing samples.

Given a previously loaded raw GC-MS data file that has been converted into an
IntensityMatrix, `im`, the mass range can be "cropped" to a new (smaller)
range;

    >>> im.crop_mass(60, 400)
    >>> print(im.min_mass, im.max_mass)

It is also possible to set all intensities for a given mass to zero. This is
useful for ignoring masses associated with sample preparation. The mass can be
"nulled" via;

    >>> data.null_mass(73)
    >>> print(sum(im.get_ic_at_mass(73).get_intensity_array()))


Noise smoothing
=================

The purpose of noise smoothing is to remove high-frequency noise from
data, and thereby increase the contribution of the signal relative to
the contribution of the noise.

Window averaging
-----------------

.. note:: This example is in `pyms-demo/41a <../pyms-demo/41a/41a.html>`__

A simple approach to noise smoothing is moving average window smoothing.
In this approach the window of a fixed size (:math:`2N+1` points) is moved
across the ion chromatogram, and the intensity value at each point is
replaced with the mean intensity calculated over the window size.
The example below illustrates smoothing of TIC by window averaging.

Load the data and get the TIC:

    >>> andi_file = "data/gc01_0812_066.cdf"
    >>> data = ANDI_reader(andi_file)
     -> Reading netCDF file 'data/gc01_0812_066.cdf'
    >>> tic = data.get_tic()

Apply the mean window smoothing with the 5-point window:

.. code-block:: python

    from pyms.Noise.Window import window_smooth
    tic1 = window_smooth(tic, window=5)
     -> Window smoothing (mean): the wing is 2 point(s)

Apply the median window smoothing with the 5-point window:

    >>> tic2 = window_smooth(tic, window=5, median=True)
     -> Window smoothing (median): the wing is 2 point(s)

Apply the mean windows smoothing, but specify the window as
a time string (in this example, 7 seconds):

    >>> tic3 = window_smooth(tic, window='7s')
    -> Window smoothing (mean): the wing is 9 point(s)

Time strings are explained in the Section `Time Strings`_.

Window Averaging on Intensity Matrix
------------------------------------
.. note:: This example is in `pyms-demo/41b <../pyms-demo/41b/41b.html>`__

In the previous section, window averaging was applied to an
Ion Chromatogram object (in that case a TIC). Where filtering
is to be performed on all Ion Chromatograms, the
:py:meth:`window_smooth_im() <pyms.Noise.Window.window_smooth_im>`
function may be used instead.

The use of this function is identical to the Ion Chromatogram
:py:meth:`window_smooth() <pyms.Noise.Window.window_smooth>`
function, except that an Intensity Matrix
is passed to it.

For example, to perform window smoothing on an
:py:meth:`IntensityMatrix <pyms.GCMS.Class.IntensityMatrix>`
object with a 5 point window and mean window smoothing:

    >>> from pyms.Noise.Window import window_smooth_im()
    ... im is a PyMS IntensityMatrix object
    >>> im_smooth = window_smooth_im(im, window = 5, median = False)

Savitzky--Golay noise filter
------------------------------

.. note:: This example is in `pyms-demo/41c <../pyms-demo/41c/41c.html>`__

A more sophisticated noise filter is the Savitzky-Golay filter.
Given the data loaded as above, this filter can be applied as
follows:

    >>> from pyms.Noise.SavitzkyGolay import savitzky_golay
    >>> tic1 = savitzky_golay(tic)
     -> Applying Savitzky-Golay filter
          Window width (points): 7
          Polynomial degree: 2

In this example the default parameters were used.

Savitzky-Golay Noise filtering of Intensity Matrix Object
----------------------------------------------------------

.. note:: This example is in `pyms-demo/41d <../pyms-demo/41d/41d.html>`__

The :py:meth:`savitzky_golay() <pyms.Noise.SavitzkyGolay.savitzky_golay>`
function described in the previous section acts on a single
Ion Chromatogram. Where it is desired to perform Savitzky Golay
filtering on the whole Intensity matrix the function
:py:meth:`savitzky_golay_im() <pyms.Noise.SavitzkyGolay.savitzky_golay_im>`
may be used as follows:

    >>> from pyms.Noise.SavitzkyGolay import savitzky_golay_im
    ... im is a PyMS IntensityMatrix object
    >>> im_smooth = savitzky_golay(im)


Baseline correction
====================
.. note:: This example is in `pyms-demo/62a <../pyms-demo/62a/62a.html>`__

Baseline distortion originating from instrument imperfections and
experimental setup is often observed in mass spectrometry data,
and off-line baseline correction is often an important step in
data pre-processing. There are many approaches for baseline
correction. One advanced approach is based top-hat transform
developed in mathematical morphology [1]_, and used
extensively in digital image processing for tasks such as image
enhancement. Top-hat baseline correction was previously applied
in proteomics based mass spectrometry [2]_.

PyMS currently implements only top-hat baseline corrector, using
the SciPy package ``ndimage``.

Application of the top-hat baseline corrector requires the size
of the structural element to be specified. The structural element
needs to be larger than the features one wants to retain in the
spectrum after the top-hat transform. In the example below, the
top-hat baseline corrector is applied to the TIC of the data set
``gc01_0812_066.cdf``, with the structural element of 1.5 minutes:


    >>> from pyms.GCMS.IO.ANDI.Function import ANDI_reader
    >>> andi_file = "data/gc01_0812_066.cdf"
    >>> data = ANDI_reader(andi_file)
     -> Reading netCDF file 'data/gc01_0812_066.cdf'
    >>> tic = data.get_tic()
    >>> from pyms.Noise.SavitzkyGolay import savitzky_golay
    >>> tic1 = savitzky_golay(tic)
     -> Applying Savitzky-Golay filter
          Window width (points): 7
          Polynomial degree: 2
    >>> from pyms.Baseline.TopHat import tophat
    >>> tic2 = tophat(tic1, struct="1.5m")
     -> Top-hat: structural element is 239 point(s)
    >>> tic.write("output/tic.dat",minutes=True)
    >>> tic1.write("output/tic_smooth.dat",minutes=True)
    >>> tic2.write("output/tic_smooth_bc.dat",minutes=True)

In the interactive session shown above, the data set if first loaded,
Savitzky-Golay smoothing was applied, followed by baseline correction.
Finally the original, smoothed, and smoothed and baseline corrected
TIC were saved in the directory ``output/``.

Tophat Baseline correction on an Intensity Matrix object
-----------------------------------------------------------

.. note:: This example is in `pyms-demo/42b <../pyms-demo/42b/42b.html>`__

The :py:meth:`tophat() <pyms.Baseline.TopHat.tophat>` function
outlined in the instructions above, acts on a single
:py:meth:`IonChromatogram <pyms.GCMS.Class.IonChromatogram>`.
To perform baseline correction on an
:py:meth:`IntensityMatrix <pyms.GCMS.Class.IntensityMatrix>`
object (i.e. on all `Ion Chromatograms`) the
:py:meth:`tophat_im() <pyms.Baseline.TopHat.tophat_im>`
function may be used.

Using the same definition for "`struct`" as above, use of the
:py:meth:`tophat_im() <pyms.Baseline.TopHat.tophat_im>`
function is as follows:

    >>> from pyms.Baseline.TopHat import tophat_im()
    ... im is an Intensity Matrix object
    >>> im_base_corr = tophat(im, struct="1.5m")


Pre-processing the IntensityMatrix
====================================

.. note:: This example is in `pyms-demo/43 <../pyms-demo/43/43.html>`__

The entire noise smoothing and baseline correction can be applied to each ion
chromatogram in the intensity matrix;

    >>> jcamp_file = "data/gc01_0812_066.jdx"
    >>> data = JCAMP_reader(jcamp_file)
    >>> im = build_intensity_matrix(data)
    >>> n_scan, n_mz = im.get_size()
    >>> for ii in range(n_mz):
    ...     print "Working on IC#", ii+1
    ...     ic = im.get_ic_at_index(ii)
    ...     ic_smooth = savitzky_golay(ic)
    ...     ic_bc = tophat(ic_smooth, struct="1.5m")
    ...     im.set_ic_at_index(ii, ic_bc)
    ...

Alternatively, the filtering may be performed on the Intensity Matrix without
using a ``for`` loop, as outlined in the sections above. However filtering by
Ion Chromatogram in a ``for`` loop as described here is much faster.

The resulting IntensityMatrix object can be ``dumped'' to a file for later
retrieval. There are general perpose object file handling methods in
:py:meth:`pyms.Utils.IO <pyms.Utils.IO>`. For example;

    >>> from pyms.Utils.IO import dump_object
    >>> dump_object(im, "output/im-proc.dump")


References
============

.. [1] Serra J. `Image Analysis and Mathematical Morphology`. Academic Press, Inc, Orlando, 1983. ISBN 0126372403

.. [2] Sauve AC and Speed TP. Normalization, baseline correction and alignment of high-throughput mass spectrometry data. `Procedings Gensips`, 2004