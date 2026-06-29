"""Sphinx configuration for the geokmeans documentation."""
import importlib.metadata

project = "geokmeans"
author = "Parichit Sharma"
copyright = "2026, Parichit Sharma"

try:
    release = importlib.metadata.version("geokmeans")
except importlib.metadata.PackageNotFoundError:
    release = "0.1.0"
version = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "numpydoc",
]

autosummary_generate = True
numpydoc_show_class_members = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "sklearn": ("https://scikit-learn.org/stable/", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_title = f"geokmeans {version}"
html_static_path = []
