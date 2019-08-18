************
PyMassSpec
************

A Python toolkit for processing of chromatography--mass spectrometry data
===========================================================================

Originally by Andrew Isaac, Sean O'Callaghan and Vladimir Likić
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Forked from the original PyMS Repository: https://github.com/ma-bio21/pyms.

The project seems to have been abandoned as there has been no activity in 18 months.

|

Introduction
==============

PyMassSpec provides a framework and a set of components for rapid development
and testing of methods for processing of chromatography--mass spectrometry data.
PyMassSpec can be used interactively through the Python shell, or the functions
can be collected into scripts and run non-interactively when it is preferable
to perform data processing in the batch mode.

PyMassSpec consists of modules which are loaded when needed,
and different functions are completely decoupled from one another.
If desired, new functions (such as a test or prototype of a new algorithm)
can be implemented efficiently and ensuring that this will not break any
existing functionality.

The original publication can be found here: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-115


The PyMS project
=================

The directory structure of PyMassSpec is as follows:

.. code-block:: text

    /
    ├── pyms:      The PyMassSpec code
    │
    ├── docs: The PyMassSpec User Guide and Documentation
    │
    ├── pyms-test: Examples of how to use PyMassSpec
    │
    └── UserGuide: Sphinx source for pyms-docs

Features
=========

Installation
==============

PyMassSpec can be installed with the following command:

.. code-block:: bash

    $ python3 -m pip --user install PyMassSpec

This will also install the following dependencies:

.. literalinclude:: ../requirements.txt


Example
=======

Download the file ``gc01\_0812\_066.jdx`` and save it in the folder ``Data``.
This file contains GC-MS data in the the JCAMP-DX format.

First the raw data is loaded:

    >>> from pyms.GCMS.IO.JCAMP.Function import JCAMP_reader
    >>> jcamp_file = "Data/gc01_0812_066.jdx"
    >>> data = JCAMP_reader(jcamp_file)
    -> Reading JCAMP file 'Data/gc01_0812_066.jdx'

The intensity matrix object is then built by binning:

    >>> from pyms.GCMS.Function import build_intensity_matrix_i
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

Further examples can be found in the `documentation`_

Contributing
==============

Contributions are very welcome. Tests can be run with `pytest`_. Please
ensure the coverage stays at least the same before you submit a pull
request.

License
=========
PyMassSpec is Free and Open Source software released under the `GNU General Public License version 2 <GPL_>`__.


Issues
========

If you encounter any problems, please `file an issue`_ along with a
detailed description.


.. _`documentation`: https://pymassspec.readthedocs.io
.. _`pytest`: https://pytest.org/latest/contents.html
.. _`file an issue`: https://github.com/domdfcoding/pymassspec/issues
.. _Python: https://www.python.org/
.. _GPL: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html