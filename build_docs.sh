pip install sphinx sphinx_rtd_theme sphinx-argparse
cd docs
sphinx-build -b html -d _build/doctrees  . _build/html
cd ../
touch docs/_build/html/.nojekyll
