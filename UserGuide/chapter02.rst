.. _chapter02:

******************************
PyMS tutorial and examples
******************************

.. contents:: Table of Contents

A tutorial illustrating various |pkgname| features in detail is provided
in subsequent chapters of this User Guide. The commands executed
interactively are grouped together by example, and can be found
:ref:`here <pyms-demo>`.

.. If you are viewing this source, the examples can be found in the pyms-demo directory, and the data files in pyms-data

The data used in the |pkgname| documentation and examples is available
:ref:`here <pyms-demo/data-files>`.

In the ":ref:`Demos and Examples <pyms-demo>`" section there
is a page corresponding to each example, coded with the chapter number
(ie. ``pyms-demo/20a/`` corresponds to the Example 21a, from Chapter 2).

Each example has a script named 'proc.py' which contains the commands given in the example.
These scripts can be run with the following command:

.. code-block:: bash

    $ python3 proc.py

Example processing GC-MS data
==============================

Download the file ``gc01\_0812\_066.jdx`` and save it in the folder ``Data``.
This file contains GC-MS data in the the JCAMP-DX format.

First the raw data is loaded:

    >>> from pyms.GCMS.IO.JCAMP.Function import JCAMP_reader
    >>> jcamp_file = "Data/gc01_0812_066.jdx"
    >>> data = JCAMP_reader(jcamp_file)
    -> Reading JCAMP file 'Data/gc01_0812_066.jdx'

The intensity matrix object is then built by binning:

    >>> from pyms.IntensityMatrix import build_intensity_matrix_i
    >>> im = build_intensity_matrix_i(data)

In this example, we show how to obtain the dimensions of the
newly created intensity matrix, then loop over all ion chromatograms,
and for each ion chromatogram apply Savitzky-Golay noise filter
and tophat baseline correction:

    >>> n_scan, n_mz = im.get_size()
    >>> from pyms.Noise.SavitzkyGolay import savitzky_golay
    >>> from pyms.Baseline.TopHat import tophat
    >>> for ii in range(n_mz):
    ...     print("working on IC", ii)
    ...     ic = im.get_ic_at_index(ii)
    ...     ic1 = savitzky_golay(ic)
    ...     ic_smooth = savitzky_golay(ic1)
    ...     ic_base = tophat(ic_smooth, struct="1.5m")
    ...     im.set_ic_at_index(ii, ic_base)

The resulting noise and baseline corrected ion chromatogram is saved
back into the intensity matrix.
