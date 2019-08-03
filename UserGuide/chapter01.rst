.. _chapter01:

************
Installation
************

.. contents:: Table of Contents


PyMassSpec is a Python_ package for processing chromatography-mass spectrometry data.
Some of its functionality depends on other Python libraries,
such as 'NumPy_' (the Python library for numerical computing),
or 'Matplotlib_' (the Python library for plotting).
In general, PyMassSpec is available for all operating systems.

.. Attention:: The 'pycdf' library works only under Linux/Unix.
    Functionality that depends on this library, namely reading
    ``ANDI-MS`` format data files, works only under Unix operating systems.

    'pycdf_' - a python interface to Unidata netCDF library -
    was written by Andre Gosselin of the Institut Maurice-Lamontagne, Canada.

.. _Python: https://www.python.org/
.. _NumPy: http://www.numpy.org/
.. _Matplotlib: https://matplotlib.org/
.. _pycdf: http://pysclint.sourceforge.net/pycdf/

There are several ways to install PyMassSpec depending your computer
configuration and personal preferences. These installation
instructions assume that Python is already installed and can be
invoked with the ``python3`` command. Modify the instructions
given if your Python installation is invoked with a different
command, such as ``py`` on Windows.

Installation with pip
======================

.. The latest source package can be downloaded :download:`here <../dist/PyMassSpec-2.1.0.tar.gz>`

The recommended installation method is to use ``pip``.
PyMassSpec can be installed with the following command:

.. code-block:: bash

    $ python3 -m pip --user install PyMassSpec

This will also install the following dependencies:

.. literalinclude:: ../requirements.txt

* The package NumPy provides numerical capabilities to Python.

    This package is used throughout PyMassSpec, and is also required for some
    external packages used in PyMS.

    Further information can be found on the NumPy website: http://www.numpy.org/

* The package 'scipy.ndimage' provides the TopHat baseline corrector

    Further information can be found on the SciPy website: https://www.scipy.org/

* The package matplotlib is used to display Ion Chromatograms and detected peaks.

    Further information can be found on the matplotlib website: https://matplotlib.org/

* The package BioPython provides an alternative to Pycluster, which is used for peak alignment.

    Further information can be found on the BioPython website: https://biopython.org/

Additional Libraries
====================

For the full PyMassSpec functionality several other libraries are required.

Package 'netCDF4-python' (required for reading ANDI-MS files)
-----------------------------------------------------------------

.. caution:: This functionality has not been implemented in PyMassSpec fully

For further information see http://unidata.github.io/netcdf4-python/netCDF4/index.html


Package 'Pycluster' (required for peak alignment by dynamic programming)
-------------------------------------------------------------------------

.. Attention:: Pycluster is only required if BioPython is not installed

The peak alignment by dynamic programming is located in the subpackage
:mod:`pyms.Peak.List.DPA`. This subpackage uses the Python package 'Pycluster_'
as the clustering engine. Pycluster and its installation instructions
can be found here: http://bonsai.hgc.jp/~mdehoon/software/cluster/index.html

.. _Pycluster: http://bonsai.hgc.jp/~mdehoon/software/cluster/index.html

PyCluster can be installed with the following commands:

.. code-block:: bash

    $ tar xvf Pycluster-1.57.tar.gz
    $ cd Pycluster-1.57
    $ python3 setup.py install

Package 'mpi4py' (required for parallel processing)
-------------------------------------------------------

This package is required for parallel processing with PyMS.
Installation instructions and download links can be found at:
https://mpi4py.readthedocs.io/en/stable/


Since 'mpi4py' provides only Python bindings, it requires an MPI implementation.
We recommend using 'mpich': https://www.mpich.org/

* On Ubuntu or Debian, mpich can be installed with the following command:

    .. code-block:: bash

        sudo apt install mpich

* On Fedora, CentOS, or RHEL, mpich can be installed with the following command:

    .. code-block:: bash

        sudo yum install mpich

* Downloads for other Linux distributions, Windows, and macOS can be found here:

    https://www.mpich.org/downloads/

With this completed, 'mpi4py' can now be installed.

The recommended installation method is to use ``pip``, with the following command:

.. code-block:: bash

    $ python3 -m pip install --user mpi4py

To check that the installation of 'mpi4py' was successful:

    .. code-block:: bash

        $ python3
        Python 3.6.7 (default, Oct 22 2018, 11:32:17)
        [GCC 8.2.0] on linux
        Type "help", "copyright", "credits" or "license" for more information.
        >>> import mpi4py
        >>>

If the above command import produced no output, mpi4py is installed
properly and ready to use.

You may need to install the package ``python3-dev`` if the installation of 'mpi4py' fails.
