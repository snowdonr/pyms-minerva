#  This file is managed by 'repo_helper'. Don't edit it directly.
#  Copyright © 2020 Dominic Davis-Foster <dominic@davis-foster.co.uk>
#
#  This file is distributed under the same license terms as the program it came with.
#  There will probably be a file called LICEN[S/C]E in the same directory as this file.
#
#  In any case, this program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This script based on https://github.com/rocky/python-uncompyle6/blob/master/__pkginfo__.py
#

# stdlib
import pathlib

__all__ = [
		"__copyright__",
		"__version__",
		"modname",
		"pypi_name",
		"__license__",
		"__author__",
		"short_desc",
		"author",
		"author_email",
		"github_username",
		"web",
		"github_url",
		"repo_root",
		"install_requires",
		"extras_require",
		"project_urls",

		"import_name",
		]

__copyright__ = """
2019-2020 Dominic Davis-Foster <dominic@davis-foster.co.uk>
"""

__version__ = "2.2.21"
modname = "PyMassSpec"
pypi_name = "PyMassSpec"
import_name = "pyms"
__license__ = "GNU General Public License v2 (GPLv2)"
short_desc = 'Python Toolkit for Mass Spectrometry'
__author__ = author = 'Dominic Davis-Foster'
author_email = 'dominic@davis-foster.co.uk'
github_username = "domdfcoding"
web = github_url = "https://github.com/domdfcoding/PyMassSpec"
repo_root = pathlib.Path(__file__).parent
install_requires = (repo_root / "requirements.txt").read_text(encoding="utf-8").split('\n')
extras_require = {'all': []}

long_description = long_description.replace(':ref:`here <pyms-demo>`.', '`here <pyms_demo_>`__.')
long_description = long_description.replace(':ref:`here <pyms-demo/data-files>`.', '`here <datafiles_>`__.')
long_description = long_description.replace(':ref:`Demos and Examples <pyms-demo>`', '`Demos and Examples <pyms_demo_>`__')
long_description = long_description.replace(':ref:`pyms-demo/20a/`', '`pyms-demo/20a/`')
long_description += '''

.. _pyms_demo: https://pymassspec.readthedocs.io/en/master/pyms-demo/introduction.html#pyms-demo
.. _datafiles: https://pymassspec.readthedocs.io/en/master/pyms-demo/data-files.html
'''

conda_description = """Python Toolkit for Mass Spectrometry


Before installing please ensure you have added the following channels: bioconda, conda-forge, domdfcoding"""
__all__.append("conda_description")


project_urls = {
		"Documentation": "https://PyMassSpec.readthedocs.io",
		"Issue Tracker": f"{github_url}/issues",
		"Source Code": github_url,
		}
