#!/usr/bin/env python
from distutils.core import setup

import inspect, os
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory

setup(name='Inchman',
      version='2.0',
      description='Inchman Utilities',
      author='Matthias Vigelius',
      author_email='matthias.vigelius@monash.edu',
      url='http://www.csse.monash.edu.au/~berndm/inchman/',
      packages=['gpgmp', 'gpgmp.test', 'gpgmp.models', 'gpgmp.inchman', 'gpgmp.common'],
      package_dir = {'': currentdir}
     )
