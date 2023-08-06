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
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath("../ler")) 


# -- Project information -----------------------------------------------------

project = 'ler'
copyright = '2023, hemanta'
author = 'hemanta'

# The full version, including alpha/beta/rc tags
release = '0.2.2'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "nbsphinx",
    "sphinx.ext.mathjax",
    "numpydoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.inheritance_diagram",
    "sphinx_tabs.tabs",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "myst_parser",
    "sphinx_copybutton",
    "autoapi.extension",
]
autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = ['.rst', '.md', '.txt', '.ipynb']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '.ipynb_checkpoints','.ipynb', "venv", ".*"]
autodoc_member_order = 'bysource'
numpydoc_show_class_members = False
autoapi_add_toctree_entry = False
# The master toctree document.
master_doc = 'index'
# -- Napoleon options
napoleon_include_special_with_doc = True
pygments_style = 'sphinx'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Plausible support
ENABLE_PLAUSIBLE = os.environ.get("READTHEDOCS_VERSION_TYPE", "") in ["branch", "tag"]
html_context = {"enable_plausible": ENABLE_PLAUSIBLE}

# -- Configure autoapi -------------------------------------------------------
autodoc_typehints = "signature"  # autoapi respects this

autoapi_type = "python"
autoapi_dirs = ["../ler"]
autoapi_template_dir = "_templates/autoapi"
autoapi_options = [
    "members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
    "special-members",
]
# autoapi_python_use_implicit_namespaces = True
autoapi_keep_files = True
# autoapi_generate_api_docs = False


# -- Napoleon options
napoleon_include_special_with_doc = True

def skip_member(app, what, name, obj, skip, options):
    if what == "ler.ler.C":
        skip = True
        
    if "ler.ler" in name:
        if obj.name in [
            "C",
        ]:
            skip = True
    return skip

def setup(sphinx):
    sphinx.connect("autoapi-skip-member", skip_member)
































