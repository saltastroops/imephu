import os
import sys
sys.path.insert(0, os.path.abspath("../src"))

project = 'imephu'
copyright = '2022, Southern African Large Telescope (SALT)'
author = 'Southern African Large Telescope (SALT)'
extensions = [
    "sphinx.ext.autodoc",
    "myst_parser",
    "numpydoc"
]
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
html_static_path = ['_static']
myst_enable_extensions = [
    "colon_fence"
]
html_theme = "sphinx_book_theme"
html_sidebars = {
    "**": ["sbt-sidebar-nav.html", "sbt-sidebar-footer.html"]
}
numpydoc_validation_checks = {"all"}
intersphinx_mapping = {
    'python': ('http://docs.python.org/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('http://docs.scipy.org/doc/scipy/reference/', None),
    'matplotlib': ('http://matplotlib.org/', None),
    'astropy': ('http://docs.astropy.org/en/stable/', None)
}
