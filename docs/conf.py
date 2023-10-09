# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath("../ler")) 

project = 'ler'
copyright = '2023, Phurailatpam Hemantakumar'
author = 'Phurailatpam Hemantakumar'
release = '0.2.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

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

# nbsphinx_custom_formats = {
#     ".md": ["jupytext.reads", {"fmt": "mystnb"}],
# }

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '.ipynb_checkpoints','.ipynb', "venv", ".*"]
autodoc_member_order = 'bysource'
numpydoc_show_class_members = False
autoapi_add_toctree_entry = False
# -- Napoleon options
napoleon_include_special_with_doc = True
pygments_style = 'sphinx'

# Don't mess with double-dash used in CLI options
smartquotes_action = "qe"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = 'sphinx_rtd_theme'
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
