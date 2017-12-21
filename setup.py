import os
import sys

import re

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

# Get the version number from _version.py, and exe_path
verstrline = open(os.path.join('tombo', '_version.py'), 'r').read()
vsre = r"^TOMBO_VERSION = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "tombo/_version.py".'.format(__pkg_name__))

def readme():
    with open('README.rst') as f:
        return f.read()

try:
    import numpy as np
    include_dirs = [np.get_include()]
except:
    sys.stderr.write(
        '*' * 60 + '\nINSTALLATION ERROR:\n'
        '\tNeed to install numpy before tombo installation.\n' +
        '\tThis is required in order to get maximum efficincy from ' +
        'cython code optimizations.\nTo install run:\n$ pip install numpy\n' +
        '*' * 60 + '\n')
    sys.exit(1)

if not sys.version_info[0] == 2:
    sys.exit("Sorry, Python 3 is not supported (yet)")

ext_modules = [
    Extension("tombo.dynamic_programming",
              ["tombo/dynamic_programming.pyx"],
              include_dirs=include_dirs,
              language="c++"),
    Extension("tombo.c_helper",
              ["tombo/c_helper.pyx"],
              include_dirs=include_dirs,
              language="c++")
]

for e in ext_modules:
    e.cython_directives = {"embedsignature": True}

setup(
    name = "ont-tombo",
    version = __version__,
    packages = ["tombo"],
    python_requires = '<3',
    install_requires = ['h5py', 'numpy', 'scipy', 'cython', 'setuptools >= 18.0'],
    extras_require={'plot':['rpy2<=2.8.6'], 'alt_est':['scikit-learn'],
                    'full':['rpy2<=2.8.6', 'scikit-learn']},

    author = "Marcus Stoiber",
    author_email = "marcus.stoiber@nanoporetech.com",
    description='Analysis of raw nanopore sequencing data.',
    long_description = readme(),
    license = "mpl-2.0",
    keywords = "nanopore high-throughput sequencing correction genome",
    url = "https://github.com/nanoporetech/tombo",

    entry_points={
        'console_scripts': [
            'tombo = tombo.__main__:main'
        ]
    },
    include_package_data=True,
    ext_modules=ext_modules,
    test_suite='nose2.collector.collector',
    tests_require=['nose2'],
)
