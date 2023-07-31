#!/usr/bin/env python

import ast
import re
from setuptools import setup, find_packages

_version_re = re.compile(r'__version__\s*=\s*(.*)')

with open('cell_line_mapper/__init__.py', 'rt') as f:
    version = str(ast.literal_eval(_version_re.search(
        f.read()).group(1)))

setup(name='cell_line_mapper',
      version=version,
      description='Functions for mapping between cell line identifiers',
      author='Phoebe Moh',
      author_email='pmoh@broadinstitute.org',
      install_requires=['pandas', 'requests'],
      packages=find_packages()
     )

