# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath("..")) 

project = 'ler'
copyright = '2023, Phurailatpam Hemantakumar'
author = 'Phurailatpam Hemantakumar'
release = '0.1.3'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "nbsphinx",
    "sphinx.ext.mathjax",
    "numpydoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    "sphinx.ext.inheritance_diagram",
    "sphinx_tabs.tabs",
    "autoapi.extension",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '.ipynb_checkpoints','.ipynb']
numpydoc_show_class_members = False



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_theme = "sphinx_rtd_theme"

# -- Configure autoapi -------------------------------------------------------
autoapi_type = "python"
autoapi_dirs = ["../ler/"]
autoapi_add_toctree_entry = False
autoapi_options = ["members", "show-inheritance", "show-module-summary"]
