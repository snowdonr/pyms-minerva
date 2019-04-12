from setuptools import setup, find_packages

setup(
	name="PyMassSpec",
	version='2.1.0',
    	author='Sean O\'Callaghan',
	author_email='spoc@unimelb.edu.au',	
	packages=find_packages(),
	license="GNU General Public License v2.0",
	url="http://github.com/domdfcoding/pyms/",
	description='Python Toolkit for Mass Spectrometry',
	long_description=('PyMassSpec is a collection of PyMS libraries '
                      'which can be used to perform peak picking ',
                      'alignment by dynamic programming, and '
                      'gap-filling in Gas Chromatography Mass '
                      'Spectrometry experiments'),
	install_requires=[
	        "numpy >= 1.16.2",
	        "scipy >= 1.2.1",
	        "pymzml >= 2.2.1",
		"matplotlib >= 3.0.2
        	],
)
