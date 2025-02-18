from setuptools import setup, Extension

module = Extension('sho', sources=['sho.cpp'])

setup(
	name='sho',
	ext_modules=[module]
)
