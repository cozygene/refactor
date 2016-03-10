from setuptools import setup
setup(name='refactor',
	version='1.0',
	description='',
	url='https://github.com/cozygene/refactor',
	packages=['refactor_lib'],
	scripts=['refactor.py'], 
	#install_requires=[], # leave this commented out.declaring dependencies that are already install can cause problems
	include_package_data=True,
	zip_safe=False,
)
