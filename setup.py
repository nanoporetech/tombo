import sys

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

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

setup(
    name = "ont-tombo",
    version = "1.0.1",
    packages = ["tombo"],
    python_requires = '<3',
    install_requires = ['h5py', 'numpy', 'scipy', 'cython', 'setuptools >= 18.0'],
    extras_require={'plot':['rpy2'], 'alt_est':['scikit-learn'],
                    'full':['rpy2', 'scikit-learn']},

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
    ext_modules=[
        Extension("tombo.dynamic_programming",
                  ["tombo/dynamic_programming.pyx"],
                  include_dirs=include_dirs,
                  language="c++"),
        Extension("tombo.c_helper",
                  ["tombo/c_helper.pyx"],
                  include_dirs=include_dirs,
                  language="c++")
    ],
    test_suite='nose2.collector.collector',
    tests_require=['nose2'],
)
