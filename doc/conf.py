# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import subprocess

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# ?? To fix:
# WARNING: Pygments lexer name 'markdown' is not known
from recommonmark.transform import AutoStructify
github_doc_root = 'https://github.com/ICB_DCM/parPE/tree/master/doc/'


def setup(app):
    app.add_config_value('recommonmark_config', {
            'url_resolver': lambda url: github_doc_root + url,
            'auto_toc_tree_section': 'Contents',
            }, True)
    app.add_transform(AutoStructify)


# RTD

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# -- Project information -----------------------------------------------------

project = 'parPE'
copyright = '2019, Daniel Weindl'
author = 'Daniel Weindl'

# The full version, including alpha/beta/rc tags
release = '-'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    #'sphinx.ext.pngmath'
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'breathe',
    'exhale',
    #'recommonmark',
    'IPython.sphinxext.ipython_console_highlighting',
    'm2r',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

# html_static_path = ['_static']


# breathe settings
breathe_projects = {
    "parPE":"doxy/xml/",
    }

breathe_default_project = "parPE"


# exhale settings

exhale_args = {
    # These arguments are required
    "containmentFolder":     "./exhale_cpp_api",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "parPE API",
    "doxygenStripFromPath":  "..",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    #"exhaleExecutesDoxygen": True,
    #    "exhaleDoxygenStdin":    "INPUT = ../include"
    "verboseBuild": False,
}

# Tell sphinx what the primary language being documented is.
# primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

html_logo = 'logo/parPE.png'


# -- RTD custom build --------------------------------------------------------

# only execute those commands when running from RTD
if on_rtd:
    parpe_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    doc_dir = os.path.dirname(os.path.abspath(__file__))

    print("Running in", os.getcwd())
    print("Generating Doxyfile")

    # # need cmake to update doxyfile
    # subprocess.run(['cmake', '-B' 'build', '-DBUILD_EXAMPLES=OFF',
    #                 '-DPARPE_ENABLE_CERES=OFF', '-DPARPE_ENABLE_DLIB=OFF',
    #                 '-DPARPE_ENABLE_FSQP=OFF', '-DPARPE_ENABLE_IPOPT=OFF',
    #                 '-DPARPE_ENABLE_MPI=OFF', '-DPARPE_ENABLE_TOMS611=OFF',
    #                 #f'-DAmici_DIR={parpe_dir}/deps/AMICI'
    #                 ],
    #                cwd=parpe_dir)

    # FIXME: workaround until we have cmake on rtd:
    #  https://github.com/readthedocs/readthedocs-docker-images/issues/127
    replacements = {
        "@GIT_VERSION@": "",
        "@CMAKE_CURRENT_LIST_DIR@": parpe_dir
    }
    with open(os.path.join(doc_dir, "Doxyfile.in"), "rt") as fin:
        with open(os.path.join(doc_dir, "Doxyfile"), "wt") as fout:
            for line in fin:
                for needle, replacement in replacements.items():
                    fout.write(line.replace(needle, replacement))

    print("Generating Doxygen docs")
    subprocess.run(['doxygen'], cwd=doc_dir, check=True)
