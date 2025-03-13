from setuptools import setup, Extension
import glob

object_files = glob.glob("../cpl/build/*.o")
module = Extension(
		'pendulum', 
		sources=['pendulum.cpp'],
		include_dirs=['../cpl'],
		extra_objects=object_files
)

setup(
	name='pendulum',
	ext_modules=[module]
)
