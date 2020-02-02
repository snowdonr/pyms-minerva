#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Built rst from jupyter notebooks
import sys
sys.path.append(".")
import ipynb2rst

#### No need to change anything in this file ####

import os
import re
import sys

sys.path.append(os.path.abspath('.'))
sys.path.append(os.path.abspath('..'))

from sphinx.locale import _

from pyms import name, __author__, __version__, __copyright__
from __pkginfo__ import github_username, modname

github_url = f"https://github.com/{github_username}/{name}"

rst_prolog = f""".. |pkgname| replace:: {name}
.. |pkgname2| replace:: ``{name}``
.. |browse_github| replace:: `Browse the GitHub Repository <{github_url}>`__
.. |ghurl| replace:: {github_url}
"""

project = name
slug = re.sub(r'\W+', '-', modname.lower())
version = __version__
release = __version__
#author = __author__
author = u'PyMassSpec Developers'
copyright = __copyright__
language = 'en'
nitpicky = True

extensions = [
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinxcontrib.httpdomain',
    'autodocsumm',
    'nbsphinx',
]

templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []

master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'

intersphinx_mapping = { # Is this where those mystery links are specified?
    'rtd': ('https://docs.readthedocs.io/en/latest/', None),
    'sphinx': ('http://www.sphinx-doc.org/en/stable/', None),
    'Python 3': ('https://docs.python.org/3', None),
    'NumPy [latest]': ('http://docs.scipy.org/doc/numpy/', None),
    'SciPy [latest]': ('http://docs.scipy.org/doc/scipy/reference', None),
    'matplotlib [latest]': ('http://matplotlib.org', None),
    'h5py [latest]': ('http://docs.h5py.org/en/latest/', None),
    'Sphinx [stable]': ('http://www.sphinx-doc.org/en/stable/', None),
    'Django [latest?]': ('http://docs.djangoproject.com/en/dev/', 'https://docs.djangoproject.com/en/dev/_objects/'),
    'sarge [latest]': ('http://sarge.readthedocs.io/en/latest/', None),
    'attrs [stable]': ('http://www.attrs.org/en/stable/', None),
}

autodoc_default_options = {
    'autosummary': True,
}

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': False,  # True will show just the logo
    'includehidden': False
}
html_theme_path = ["../.."]
#html_logo = "logo/pyms.png"
html_show_sourcelink = False    # True will show link to source

html_context = {
    # Github Settings
    "display_github": True, # Integrate GitHub
    "github_user": github_username, # Username
    "github_repo": name, # Repo name
    "github_version": "master", # Version
    "conf_py_path": "/", # Path in the checkout to the docs root
}

htmlhelp_basename = slug

latex_documents = [
  ('index', '{0}.tex'.format(slug), project, author, 'manual'),
]

man_pages = [
    ('index', slug, project, [author], 1)
]

texinfo_documents = [
  ('index', slug, project, author, slug, project, 'Miscellaneous'),
]


# Extensions to theme docs
def setup(app):
    from sphinx.domains.python import PyField
    from sphinx.util.docfields import Field

    app.add_object_type(
        'confval',
        'confval',
        objname='configuration value',
        indextemplate='pair: %s; configuration value',
        doc_field_types=[
            PyField(
                'type',
                label=_('Type'),
                has_arg=False,
                names=('type',),
                bodyrolename='class'
            ),
            Field(
                'default',
                label=_('Default'),
                has_arg=False,
                names=('default',),
            ),
        ]
    )
