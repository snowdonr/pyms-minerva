***********************************
Peak detection and representation
***********************************

.. contents:: Table of Contents

Peak Object
=============

Fundamental to GC-MS analysis is the identification of individual components of the sample mix.
The basic component unit is represented as a signal peak.
In |pkgname| a signal peak is represented as :meth:`~pyms.Peak.Class.Peak` object (the class defined in :mod:`pyms.Peak`.
|pkgname| provides functions to detect peaks and create peaks (discussed at the end of the chapter).

A peak object stores a minimal set of information about a signal peak, namely, the retention time at which the peak apex occurs and the mass spectra at the apex.
Additional information, such as, peak width, TIC and individual ion areas can be filtered from the GC-MS data and added to the Peak object information.

Creating a Peak Object
-------------------------

.. note:: This example is in :ref:`pyms-demo/50 <demo-50>`

A peak object can be created for a scan at a given retention time by providing the retention time (in minutes or seconds) and the :class:`~pyms.MassSpectrum.MassSpectrum` object of the
scan. In the example below, first a file is loaded and an :class:`~pyms.IntensityMatrix.IntensityMatrix`, ``im``, built, then a :class:`~pyms.MassSpectrum.MassSpectrum`, ``ms``, can be selected at a given time (31.17 minutes in this example).

    >>> from pyms.IntensityMatrix import build_intensity_matrix_i
    >>> from pyms.GCMS.IO.ANDI import ANDI_reader
    >>> andi_file = "data/gc01_0812_066.cdf"
    >>> data = ANDI_reader(andi_file)
    >>> im = build_intensity_matrix_i(data)
    >>> index = im.get_index_at_time(31.17*60.0)
    >>> ms = im.get_ms_at_index(index)

Now a :meth:`~pyms.Peak.Class.Peak` object can be created for the given retention time and :class:`~pyms.MassSpectrum.MassSpectrum`.

    >>> from pyms.Peak.Class import Peak
    >>> peak = Peak(31.17, ms, minutes=True)


By default the retention time is assumed to be in seconds. The parameter ``minutes`` can be set to ``True`` if the retention time is given in minutes.
As a matter of convention, |pkgname| internally stores retention times in seconds, so the ``minutes`` parameter ensures the input and output of the retention time are in the same units.

Peak Object properties
------------------------

.. note:: This example is in :ref:`pyms-demo/50 <demo-50>`

The retention time of the peak can be returned with :attr:`~pyms.Peak.Class.Peak.rt`.
The retention time is returned in seconds with this method.
The mass spectrum can be returned with :attr:`~pyms.Peak.Class.Peak.mass_spectrum`.

The :class:`~pyms.Peak.Class.Peak` object constructs a unique identification (UID) based on the spectrum
and retention time. This helps in managing lists of peaks (covered in the next chapter).
The UID can be returned with :arre:`pyms.Peak.Class.Peak.UID`.
The format of the UID is the masses of the two most abundant ions in the spectrum, the ratio of the abundances of the two ions, and the retention time (in the same units as given when the Peak object was created).
The format is:

.. code-block:: text

    Mass1-Mass2-Ratio-RT

For example,

    >>> peak.rt
    1870.2
    >>> peak.UID
    319-73-74-31.17


Modifying a Peak Object
-------------------------

.. note:: This example is in :ref:`pyms-demo/51 <demo-51>`

The Peak object has methods for modifying the mass spectrum. The mass range can be cropped to a smaller range with :meth:`crop_mass() <pyms.Peak.Class.Peak.crop_mass>`, and the intensity values for a single ion can be set to zero with :meth:`null_mass() <pyms.Peak.Class.Peak.null_mass>`.
For example, the mass range can be set from 60 to 450 :math:`m/z`, and the ions related to
sample preparation can be ignored by setting their intensities to zero as follows:

    >>> peak.crop_mass(60, 450)
    >>> peak.null_mass(73)
    >>> peak.null_mass(147)

The UID is automatically updated to reflect the changes;

    >>> peak.UID
    319-205-54-31.17

It is also possible to change the peak mass spectrum by setting the attribute :attr:`~pyms.Peak.Class.Peak.mass_spectrum`.

Peak Detection
================

The general use of a :class:`~pyms.Peak.Class.Peak` object is to extract them from the GC-MS data and build a list of peaks. In |pkgname|, the function for peak detection is based on the method of Biller and Biemann (1974) [1]_.
The basic process is to find all maximising ions in a pre-set window of scans, for a given scan.
The ions that maximise at a given scan are taken to belong to the same peak.

The function is :py:meth:`BillerBiemann() <pyms.BillerBiemann.BillerBiemann>`. in :mod:`pyms.BillerBiemann`.
The function has parameters for the window width for detecting the local maxima (``points``), and the number of ``scans`` across which neighbouring, apexing, ions are combined and considered as belonging to the same peak.
The number of neighbouring scans to combine is related to the likelihood of detecting a peak apex at a single scan or several neighbouring scans.
This is more likely when there are many scans across the peak.
It is also possible, however, when there are very few scans across the peak.
The scans are combined by taking all apexing ions to have occurred at the scan that had to greatest TIC prior to combining scans.

Sample processing and Peak detection
-------------------------------------

.. note:: This example is in :ref:`pyms-demo/52 <demo-52>`

The process for detecting peaks is to pre-process the data by performing noise smoothing and baseline correction on each ion (as in :ref:`pyms-demo/51 <demo-51>`).
The first steps then are:

    >>> from pyms.GCMS.IO.ANDI import ANDI_reader
    >>> from pyms.IntensityMatrix import build_intensity_matrix
    >>> from pyms.Noise.SavitzkyGolay import savitzky_golay
    >>> from pyms.TopHat import tophat
    >>>
    >>> andi_file = "/x/PyMS/data/gc01_0812_066.cdf"
    >>> data = ANDI_reader(andi_file)
    >>>
    >>> im = build_intensity_matrix(data)
    >>> n_scan, n_mz = im.size
    >>>
    >>> for ii in range(n_mz):
    ...     ic = im.get_ic_at_index(ii)
    ...     ic_smooth = savitzky_golay(ic)
    ...     ic_bc = tophat(ic_smooth, struct="1.5m")
    ...     im.set_ic_at_index(ii, ic_bc)
    ...

Now the Biller and Biemann based technique can be applied to detect peaks.

    >>> from pyms.BillerBiemann import BillerBiemann
    >>> peak_list = BillerBiemann(im)
    >>> len(peak_list)
    9845

Note that this is nearly as many peaks as there are scans in the data (9865 scans).
This is due to noise and the simplicity of the technique.

The number of detected peaks can be constrained by the selection of better parameters.
Parameters can be determined by counting the number of points across a peak, and examining where peaks are found.
For example, the peak list can be found with the parameters of a window of 9 points and by combining 2 neighbouring scans if they apex next to each other:

    >>> peak_list = BillerBiemann(im, points=9, scans=2)
    >>> len(peak_list)
    3698

The number of detected peaks has been reduced, but there are still many more than would be expected from the sample. Functions to filter the peak list are covered in the next section.

Filtering Peak Lists
====================

.. note:: This example is in :ref:`pyms-demo/53 <demo-53>`

There are two functions to filter the list of Peak objects.
The first, :meth:`rel_threshold() <pyms.BillerBiemann.rel_threshold>`, modifies the mass spectrum stored in each peak so any intensity that is less than a given percentage of the maximum intensity for the peak is removed.
The second, :meth:`num_ions_threshold() <pyms.BillerBiemann.num_ions_threshold>`, removes any peak that has less than a given number of ions above a given threshold.
Once the peak list has been constructed, the filters can be applied by:

    >>> from pyms.Deconvolution.BillerBiemann.Function import \
    ... rel_threshold, num_ions_threshold
    >>> pl = rel_threshold(peak_list, percent=2)
    >>> new_peak_list = num_ions_threshold(pl, n=3, cutoff=10000)
    >>> len(new_peak_list)
    146

The number of detected peaks is now more realistic of what would be expected in
the test sample.

Noise analysis for peak filtering
==================================

.. note:: This example is in :ref:`pyms-demo/54 <demo-54>`

In the previous section the cutoff parameter for peak filtering was set by the user.
This can work well for individual data files, but can cause problems when applied to large experiments with many individual data files.
Where experimental conditions have changed slightly between experimental runs, the ion intensity over the GC-MS run may also change.
This means that an inflexible cutoff value can work for some data files, while excluding too many, or including too many peaks in other files.

An alternative to manually setting the value for cutoff is to use the :meth:`window_analyzer() <pyms.Noise.Analysis.window_analyzer>` function.
This function examines a Total Ion Chromatogram (TIC) and computes a value for the median absolute deviation in troughs between peaks.
This gives an approximate threshold value above which false peaks from noise should be filtered out.

To compute this noise value:

    >>> from pyms.Noise.Analysis import window_analyzer
    >>> # data is a GCMS data object
    >>> tic = data.tic
    >>> noise_level = window_analyzer(tic)


Now the usual peak deconvolution steps are performed, and the peak list is filtered using this noise value as the cutoff:

    >>> peak_list = num_ions_threshold(pl, n, noise_level)
    >>> # pl is a peak list, n is number of ions above threshold

Peak area estimation
========================

.. note:: This example is in :ref:`pyms-demo/55 <demo-55>`

The :class:`~pyms.Peak.Class.Peak` object does not contain any information about the width or area of the peak when it is created.
This information can be added after the instantiation of a Peak object.
The area of the peak can be set with the attribute :attr:`~pyms.Peak.Class.Peak.area`, or with the method :meth:`set_ion_areas() <pyms.Peak.Class.Peak.set_ion_areas>`.

The total peak area can by obtained by the :meth:`peak_sum_area() <pyms.Peak.Function.peak_sum_area>` function in :mod:`pyms.Peak.Function`.
The function determines the total area as the sum of the ion intensities for all masses that apex at the given peak.
To calculate the peak area of a single mass, the intensities are added from the apex of the mass peak outwards.

Edge values are added until the following conditions are met:
* the added intensity adds less than 0.5\% to the accumulated area; or
* the added intensity starts increasing (i.e. when the ion is common to co-eluting compounds).

To avoid noise effects, the edge value is taken at the midpoint of three consecutive edge values.

Given a list of peaks, areas can be determined and added as follows:

    >>> from pyms.Peak.Function import peak_sum_area
    >>> for peak in peak_list:
    ...     area = peak_sum_area(intensity_matrix, peak)
    ...     peak.area = area
    ...


.. individual_ion_areas:

Individual Ion Areas
------------------------

.. note:: This example is in :ref:`pyms-demo/56 <demo-56>`

While the previous approach uses the sum of all areas in the peak to estimate the peak area, the user may also choose to record the area of each individual ion in each peak.

This can be useful when the intention is to later perform quantitation based on the area of a single characteristic ion for a particular compound.
It is also essential if using the Common Ion Algorithm for quantitation, outlined in the section :ref:`common-ion`.

To set the area of each ion for each peak, the following code is used:

    >>> from pyms.Peak.Function import peak_top_ion_areas
    >>> for peak in peak_list:
    ...     area_dict = peak_top_ions_areas(intensity_matrix, peak)
    ...     peak.set_ion_areas(area_dict)
    ...

This will set the areas of the 5 most abundant ions in each peak.
If it is desired to record more than the top five ions, the argument ``num_ions=x`` should be supplied, where ``x`` is the number of most abundant ions to be recorded.
For example:

.. code-block:: python

    ...     area_dict = peak_top_ions_areas(intensity_matrix, peak, num_ions=10)

will record the 10 most abundant ions for each peak.

The individual ion areas can be set instead of, or in addition to the total area for each peak.

Reading the area of a single ion in a peak
-------------------------------------------

If the individual ion areas have been set for a peak, it is possible to read the area of an individual ion for the peak.
For example:

>>> peak.get_ion_area(101)

will return the area of the :math:`m/z` value 101 for the peak.
If the area of that ion has not been set (i.e. it was not one of the most abundant ions), the function will return ``None``.

References
============

.. [1] Biller JE and Biemann K. Reconstructed mass spectra, a novel approach for the utilization of gas chromatograph–mass spectrometer data. `Anal. Lett.`, 7:515–528, 1974
