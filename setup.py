# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
setup(name = 'osmosisbackend', version = '1.0',  \
   ext_modules = [Extension('osmosisbackend', ['osmosisbackend.c'])])