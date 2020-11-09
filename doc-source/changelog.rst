===============
Changelog
===============

Changes in v2.3.0
--------------------------

* All functions, classes and methods now have :pep:`484` type hints. Contributed by Chris Davis-Foster in :pull:`4`

* All modules now implement ``__all__`` to limit the objects imported when using ``*`` imports.

* Removed the following deprecated functions:

  .. list-table::
    :header-rows: 1

    * - Removed object
      - Suggested replacement

    * - :meth:`pyms.Experiment.Experiment.get_expr_code`
      - :meth:`pyms.Experiment.Experiment.expr_code`

    * - :meth:`pyms.Experiment.Experiment.get_peak_list`
      - :meth:`pyms.Experiment.Experiment.peak_list`

    * - :meth:`pyms.Experiment.Experiment.store`
      - :meth:`pyms.Experiment.Experiment.dump`

    * - :func:`pyms.Experiment.store_expr`
      - :meth:`pyms.Experiment.Experiment.dump`

    * - :meth:`pyms.GCMS.Class.GCMS_data.get_scan_list`
      - :meth:`pyms.GCMS.Class.GCMS_data.scan_list`

    * - :meth:`pyms.GCMS.Class.GCMS_data.get_tic`
      - :meth:`pyms.GCMS.Class.GCMS_data.tic`

    * - :meth:`pyms.Gapfill.Class.MissingPeak.get_common_ion`
      - :meth:`pyms.Gapfill.Class.MissingPeak.common_ion`

    * - :meth:`pyms.Gapfill.Class.MissingPeak.get_common_ion_area`
      - :attr:`pyms.Gapfill.Class.MissingPeak.common_ion_area`

    * - :meth:`pyms.Gapfill.Class.MissingPeak.get_exact_rt`
      - :meth:`pyms.Gapfill.Class.MissingPeak.exact_rt`

    * - :meth:`pyms.Gapfill.Class.MissingPeak.get_qual_ion1`
      - :meth:`pyms.Gapfill.Class.MissingPeak.qual_ion1`

    * - :meth:`pyms.Gapfill.Class.MissingPeak.get_qual_ion2`
      - :meth:`pyms.Gapfill.Class.MissingPeak.qual_ion2`

    * - :meth:`pyms.Gapfill.Class.MissingPeak.get_rt`
      - :meth:`pyms.Gapfill.Class.MissingPeak.rt`

    * - :meth:`pyms.Gapfill.Class.MissingPeak.set_common_ion_area`
      - :attr:`pyms.Gapfill.Class.MissingPeak.common_ion_area`

    * - :meth:`pyms.Gapfill.Class.MissingPeak.set_exact_rt`
      - :attr:`pyms.Gapfill.Class.MissingPeak.exact_rt`

    * - :meth:`pyms.Gapfill.Class.Sample.get_missing_peaks`
      - :meth:`pyms.Gapfill.Class.Sample.missing_peaks`

    * - :meth:`pyms.Gapfill.Class.Sample.get_mp_rt_area_dict`
      - :meth:`pyms.Gapfill.Class.Sample.rt_areas`

    * - :meth:`pyms.Gapfill.Class.Sample.get_name`
      - :meth:`pyms.Gapfill.Class.Sample.name`

    * - :func:`pyms.Gapfill.Function.transposed`
      -

    * - :meth:`pyms.IonChromatogram.IonChromatogram.get_mass`
      - :meth:`pyms.IonChromatogram.IonChromatogram.mass`

    * - :meth:`pyms.IonChromatogram.IonChromatogram.get_time_step`
      - :meth:`pyms.IonChromatogram.IonChromatogram.time_step`

    * - :meth:`pyms.IonChromatogram.IonChromatogram.set_intensity_array`
      - :meth:`pyms.IonChromatogram.IonChromatogram.intensity_array`

    * - :meth:`pyms.Mixins.MaxMinMassMixin.get_max_mass`
      - :meth:`pyms.Mixins.MaxMinMassMixin.max_mass`

    * - :meth:`pyms.Mixins.MaxMinMassMixin.get_min_mass`
      - :meth:`pyms.Mixins.MaxMinMassMixin.min_mass`

    * - :meth:`pyms.Mixins.MassListMixin.get_mass_list`
      - :meth:`pyms.Mixins.MassListMixin.mass_list`

    * - :meth:`pyms.Mixins.TimeListMixin.get_time_list`
      - :meth:`pyms.Mixins.TimeListMixin.time_list`

    * - :meth:`pyms.Mixins.IntensityArrayMixin.get_intensity_array`
      - :meth:`pyms.Mixins.IntensityArrayMixin.intensity_array`

    * - :meth:`pyms.Mixins.IntensityArrayMixin.get_matrix_list`
      - :meth:`pyms.Mixins.IntensityArrayMixin.intensity_array_list`

    * - :meth:`pyms.Peak.Class.Peak.get_area`
      - :meth:`pyms.Peak.Class.Peak.area`

    * - :meth:`pyms.Peak.Class.Peak.get_ic_mass`
      - :meth:`pyms.Peak.Class.ICPeak.ic_mass`

    * - :meth:`pyms.Peak.Class.Peak.get_ion_areas`
      - :meth:`pyms.Peak.Class.Peak.ion_areas`

    * - :meth:`pyms.Peak.Class.Peak.get_mass_spectrum`
      - :meth:`pyms.Peak.Class.Peak.mass_spectrum`

    * - :meth:`pyms.Peak.Class.Peak.get_pt_bounds`
      - :meth:`pyms.Peak.Class.Peak.bounds`

    * - :meth:`pyms.Peak.Class.Peak.get_rt`
      - :meth:`pyms.Peak.Class.Peak.rt`

    * - :meth:`pyms.Peak.Class.Peak.get_UID`
      - :meth:`pyms.Peak.Class.Peak.UID`

    * - :meth:`pyms.Peak.Class.Peak.set_area`
      - :meth:`pyms.Peak.Class.Peak.area`

    * - :meth:`pyms.Peak.Class.Peak.set_ic_mass`
      - :meth:`pyms.Peak.Class.ICPeak.ic_mass`

    * - :meth:`pyms.Peak.Class.Peak.set_ion_areas`
      - :meth:`pyms.Peak.Class.Peak.ion_areas`

    * - :meth:`pyms.Peak.Class.Peak.set_mass_spectrum`
      - :meth:`pyms.Peak.Class.Peak.mass_spectrum`

    * - :meth:`pyms.Peak.Class.Peak.set_pt_bounds`
      - :attr:`pyms.Peak.Class.Peak.pt_bounds`

    * - :func:`pyms.Utils.Utils.is_positive_int`
      -

    * - :func:`pyms.Utils.Utils.is_list_of_dec_nums`
      -


* Renamed :func:`pyms.Gapfill.Function.file2matrix` to :func:`pyms.Gapfill.Function.file2dataframe`.
  The function now returns a Pandas DataFrame.

* Split :class:`pyms.IntensityMatrix.IntensityMatrix` into two classes:
  :class:`pyms.IntensityMatrix.BaseIntensityMatrix` and :class:`pyms.IntensityMatrix.IntensityMatrix`.
  This makes subclassing easier.

* Split :class:`pyms.Peak.Class.Peak` into three classes:
  :class:`pyms.Peak.Class.AbstractPeak`, :class:`pyms.Peak.Class.Peak`, :class:`pyms.Peak.Class.ICPeak`.
  :class:`~pyms.Peak.Class.ICPeak` is returned when a mass is passed to the :class:`~pyms.Peak.Class.Peak` constructor
  instead of a mass spectrum.

* Added the following functions and classes:

  .. autosummary::

    pyms.Gapfill.Function.MissingPeakFiletype
    pyms.IntensityMatrix.AsciiFiletypes
    pyms.IntensityMatrix.IntensityMatrix.bpc
    pyms.IonChromatogram.IonChromatogram.is_eic
    pyms.IonChromatogram.IonChromatogram.is_bpc
    pyms.IonChromatogram.ExtractedIonChromatogram
    pyms.IonChromatogram.BasePeakChromatogram
    pyms.Spectrum.array_as_numeric
    pyms.Utils.Utils.is_path
    pyms.Utils.Utils.is_sequence
    pyms.Utils.Utils.is_sequence_of
    pyms.Utils.Utils.is_number
    pyms.eic

* The ``ia`` parameter of :class:`pyms.IonChromatogram.IonChromatogram`` was renamed to ``intensity_list``.


Changes in v2.2.22-beta2
--------------------------

* :class:`pyms.Spectrum.Scan` and :class:`pyms.Spectrum.MassSpectrum` can now accept any values for ``mass`` and ``intensity`` that that can be converted to a :class:`float` or :class:`int`. This includes strings representing numbers. Previously only :class:`int` and :class:`float` values were permitted.

* If the mass and intensity values supplied to a :class:`pyms.Spectrum.Scan` or a :class:`pyms.Spectrum.MassSpectrum` are :class:`float`, :class:`int`, or a data type derived from :class:`numpy.number`, the data is stored in that type. For other data types, such as strings, :class:`decimal.Decimal` etc., the data is stored as :class:`float`.

  If the data contains values in mixed types then, in most cases, all values will be converted to :class:`float`. If you wish to control this behaviour you should construct a :class:`numpy.ndarray` with the desired type. See https://numpy.org/devdocs/user/basics.types.html for a list of types.

* A :exc:`TypeError` is no longer raised when creating a :class:`pyms.Spectrum.Scan` or a :class:`pyms.Spectrum.MassSpectrum` with a :class:`float`, :class:`int` etc. rather than a sequence. Instead, value is treated as being the sole element in a list.

* Passing a non-numeric string or a list of non-numeric strings to :class:`pyms.Spectrum.Scan` or :class:`pyms.Spectrum.MassSpectrum` now raises a :exc:`ValueError` and not a :exc:`TypeError` as in previous versions.

* :meth:`pyms.Peak.Class.Peak.ion_areas` now accepts dictionary keys as :class:`float` as well as :class:`int`.


Changes in v2.2.22-beta1
--------------------------

* :func:`~pyms.GCMS.IO.ANDI.ANDI_reader` and :class:`pyms.Spectrum.Scan` were modified to allow ANDI-MS files to be read if the data either:

    - had the *m/z* data stored from highest *m/z* to lowest; or
    - contained 0-length scans.
