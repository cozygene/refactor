from setuptools import setup
#from distutils.core import setup

setup(name='refactor',
	version='1.0',
	description='',
	url='https://github.com/cozygene/refactor',
	package=['python'],
	scripts=['python/refactor.py'], 
	# dependencies on PyPI
	install_requires=[
#          'numpy',
#          'scipy',
#          'sklearn',
#          'matplotlib',
#          'statsmodels'
	],
	# dependencies not on PyPI
	dependency_links=['https://pypi.python.org/pypi/scipy'],
	include_package_data=True,
	#zip_safe=True,
)
