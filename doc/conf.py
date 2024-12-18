# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import sys
sys.path.insert(0, os.path.abspath('../python/'))

github_doc_root = 'https://github.com/ICB_DCM/parPE/tree/master/doc/'

autodoc_mock_imports = ['amici', 'libsbml']

# RTD

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme

    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# -- Project information -----------------------------------------------------

project = 'parPE'
copyright = '2021, The parPE developers'
author = 'Daniel Weindl'

# The full version, including alpha/beta/rc tags
release = '-'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # 'sphinx.ext.pngmath'
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx.ext.napoleon',
    'breathe',
    'exhale',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
    # 'recommonmark',
    'nbsphinx',
    'IPython.sphinxext.ipython_console_highlighting',
    'myst_parser',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    'logo/LICENSE.md',
    'README.md',
    '**/collapsible-lists/**',
]
source_suffix = ['.rst', '.md']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

# html_static_path = ['_static']


# breathe settings
breathe_projects = {
    "parPE": "_exhale_doxygen_output/xml/",
}

breathe_default_project = "parPE"

# exhale settings

exhale_args = {
    # These arguments are required
    "containmentFolder": "./cpp_api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "parPE C++ API",
    "doxygenStripFromPath": "..",
    # Suggested optional arguments
    "createTreeView": True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": "\n".join([
        "INPUT = ../include",
        "BUILTIN_STL_SUPPORT    = YES",
        "PREDEFINED            += EXHALE_DOXYGEN_SHOULD_SKIP_THIS"
    ]),
    "afterTitleDescription":
        "parPE C++ library functions",
    "verboseBuild": True,
}

# Tell sphinx what the primary language being documented is.
# primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

html_logo = 'logo/parPE.png'

autosummary_generate = True

