===========
PyMassSpec
===========

.. start short_desc

**Python Toolkit for Mass Spectrometry**

.. end short_desc

.. start shields

.. list-table::
	:stub-columns: 1
	:widths: 10 90

	* - Docs
	  - |docs| |docs_check|
	* - Tests
	  - |travis| |actions_windows| |actions_macos| |coveralls| |codefactor|
	* - PyPI
	  - |pypi-version| |supported-versions| |supported-implementations| |wheel|
	* - Anaconda
	  - |conda-version| |conda-platform|
	* - Activity
	  - |commits-latest| |commits-since| |maintained|
	* - Other
	  - |license| |language| |requires| |pre_commit|

.. |docs| rtfd-shield::
	:project: pymassspec
	:alt: Documentation Build Status

.. |docs_check| actions-shield::
	:workflow: Docs Check
	:alt: Docs Check Status

.. |travis| travis-shield::
	:travis-site: org
	:alt: Travis Build Status

.. |actions_windows| actions-shield::
	:workflow: Windows Tests
	:alt: Windows Tests Status

.. |actions_macos| actions-shield::
	:workflow: macOS Tests
	:alt: macOS Tests Status

.. |requires| requires-io-shield::
	:alt: Requirements Status

.. |coveralls| coveralls-shield::
	:alt: Coverage

.. |codefactor| codefactor-shield::
	:alt: CodeFactor Grade

.. |pypi-version| pypi-shield::
	:project: PyMassSpec
	:version:
	:alt: PyPI - Package Version

.. |supported-versions| pypi-shield::
	:project: PyMassSpec
	:py-versions:
	:alt: PyPI - Supported Python Versions

.. |supported-implementations| pypi-shield::
	:project: PyMassSpec
	:implementations:
	:alt: PyPI - Supported Implementations

.. |wheel| pypi-shield::
	:project: PyMassSpec
	:wheel:
	:alt: PyPI - Wheel

.. |conda-version| image:: https://img.shields.io/conda/v/domdfcoding/PyMassSpec?logo=anaconda
	:target: https://anaconda.org/domdfcoding/PyMassSpec
	:alt: Conda - Package Version

.. |conda-platform| image:: https://img.shields.io/conda/pn/domdfcoding/PyMassSpec?label=conda%7Cplatform
	:target: https://anaconda.org/domdfcoding/PyMassSpec
	:alt: Conda - Platform

.. |license| github-shield::
	:license:
	:alt: License

.. |language| github-shield::
	:top-language:
	:alt: GitHub top language

.. |commits-since| github-shield::
	:commits-since: v2.2.21
	:alt: GitHub commits since tagged version

.. |commits-latest| github-shield::
	:last-commit:
	:alt: GitHub last commit

.. |maintained| maintained-shield:: 2020
	:alt: Maintenance

.. |pre_commit| pre-commit-shield::
	:alt: pre-commit

.. end shields

Installation
-------------

.. start installation

.. installation:: PyMassSpec
	:pypi:
	:github:
	:anaconda:
	:conda-channels: bioconda, conda-forge, domdfcoding

.. end installation


.. toctree::
	:hidden:

	Home<self>

.. toctree::
	:numbered:
	:caption: User Guide
	:maxdepth: 2

	10_gcms_raw_data_model
	20_gcms_data_derived_objects
	30_data_filtering
	40_peak_detection_and_representation
	50_peak_alignment_by_dynamic_programming
	60_display

.. toctree::
	:maxdepth: 2
	:caption: Documentation

	pyms/documentation
	pyms/Base
	pyms/BillerBiemann
	pyms/Display
	pyms/DPA
	pyms/Experiment
	pyms/Gapfill
	pyms/GCMS
	pyms/IntensityMatrix
	pyms/IonChromatogram
	pyms/json
	pyms/Mixins
	pyms/Spectrum
	pyms/Noise
	pyms/Peak
	pyms/Simulator
	pyms/TopHat
	pyms/Utils
	changelog


.. toctree::
	:maxdepth: 2
	:caption: Demos and Examples

	../../pyms-demo/introduction
	pyms-demo/data-files

.. toctree::
	:hidden:

	pyms-demo/20a
	pyms-demo/20b
	pyms-demo/20c
	pyms-demo/20d
	pyms-demo/30a
	pyms-demo/30b
	pyms-demo/30c
	pyms-demo/31
	pyms-demo/32
	pyms-demo/40a
	pyms-demo/40b
	pyms-demo/41a
	pyms-demo/41b
	pyms-demo/42a
	pyms-demo/42b
	pyms-demo/43
	pyms-demo/50
	pyms-demo/51
	pyms-demo/52
	pyms-demo/53
	pyms-demo/54
	pyms-demo/60
	pyms-demo/61a
	pyms-demo/61b
	pyms-demo/62
	pyms-demo/63
	pyms-demo/64
	pyms-demo/70a
	pyms-demo/70b
	pyms-demo/71
	pyms-demo/90
	pyms-demo/91
	pyms-demo/92
	pyms-demo/93
	pyms-demo/94
	pyms-demo/95
	pyms-demo/A1
	pyms-demo/A2
	pyms-demo/x10


.. toctree::
	:maxdepth: 3
	:caption: Contributing

	contributing
	Source
	StyleGuide


.. start links

View the :ref:`Function Index <genindex>` or browse the `Source Code <_modules/index.html>`__.

`Browse the GitHub Repository <https://github.com/domdfcoding/PyMassSpec>`__

.. end links
