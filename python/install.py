#!/usr/bin/env python

import os
import sys
from setuptools import setup

REFACTOR_OBLIGATORY_DEPENDENCIES = ['numpy', 'scipy', 'sklearn']
REFACTOR_OPTIONAL_DEPENDENCIES = ['matplotlib','statsmodels'] # for demo only

def setup():
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

def run_function_without_prints(func):
    """
    sets stdout and stderr to null before running function func
    restores original stdout and stderr after function func is done
    returns founction func output
    """
    def inner(*args):

        stdout = sys.stdout
        stderr = sys.stderr
        null = open(os.devnull, 'w')
        sys.stderr = null
        sys.stdout = null

        #run function
        output = func(*args)
        
        # restore stdout stderr
        sys.stderr = stderr
        sys.stdout = stdout
        return output
        
    return inner

@run_function_without_prints
def pip_install(module_name):
    """
    installs module_name using pip
    pip prints are disable
    """
    return pip.main(['install', module_name])

def load_pip():
    """
    imports pip if possible
    returns True if import was successfull
    returns False otherwise
    """
    try:
        import pip
        global pip
        return True
    except Exception, e:
        return False

def install(module_name):
    """
    trying to install module_name and import it after successfull installation
    returns True if installation was successfull (and module can be imported)
    returns False otherwise
    """
    exit_code = pip_install(module_name)

    if exit_code == 0: # success
        try:
            __import__(module_name)
            print "Package %s installed successfully!" % module_name
            return True
            
        except Exception, e:
            pass

    print "Package %s installation failed" % module_name
    return False


def already_installed(module_name):
    """
    checks if module_name ia already installed by trying to import it
    """
    try:
        __import__(module_name)
        return True

    # if __import__ excepts - module is not installled
    except Exception, e:
        return False

def user_installation_confirmation(modules):
    """
    ask user to install modules
    returns True if user confirmed installation
    returns False otherwise
    """
    request = raw_input ("Packages %s aren't installed. Do you want to start installation now? (y/n) " % str(modules))
    while(request.lower() != 'y' and request.lower != 'n'):
        request = raw_input ("Packages %s aren't installed. Do you want to start installation now? (y/n) " % str(modules))

    return request.lower() == 'y'


def check_dependencies(dependencies_list):
    """
    if a module in dependencies_list is not installed, tries to install it using pip (if user confirms installation)
    returns the list of dependencies left to be installed
    """
    dependencies_to_install = []

    # check if deoendencies are installed
    print "Checking if dependencies are installed..."
    for module in dependencies_list:
        if not already_installed(module):
            dependencies_to_install.append(module)

    # ask user to install dependencies. load pip and try to install if user confirms
    if len(dependencies_to_install) > 0:
        if user_installation_confirmation(dependencies_to_install):
            if not load_pip():
                print "pip (installation package) wasn't found. Can't install dependencies"
            else:
                dependencies_to_install = [module for module in dependencies_to_install if install(module) == False]

    return dependencies_to_install


def install_refactor():
    """
    runs setup.py script to put refactor in PATH
    installes dependencies if they arn't already installed
    """

    print "Installing ReFACTOR..."
    setup()
    
    # install dependencies
    dependencies_to_install = check_dependencies(REFACTOR_OBLIGATORY_DEPENDENCIES + REFACTOR_OPTIONAL_DEPENDENCIES)

    if len(dependencies_to_install) == 0:
        print "All dependencies are installed! Done."
    else:
        obligatory_dependencies = []
        optional_dependencied = []
        for module in dependencies_to_install:
            if module in REFACTOR_OBLIGATORY_DEPENDENCIES:
                obligatory_dependencies.append(module)
            else:
                optional_dependencied.append(module)
        
        if obligatory_dependencies:
            print "To run ReFACTOR you must install the following packages: %s" % str(obligatory_dependencies)
        if optional_dependencied:
            print "To run ReFACTOR's demo you need to install: %s " % str(optional_dependencied)

        raw_input("Press any key to continue")
        

if __name__ == '__main__':

    # run as root (root privilege is required to run python setup.py script and pip install)
    euid = os.geteuid()
    if euid != 0:
        print "Installation must be started as root. Running sudo..."
        args = ['sudo', sys.executable] + sys.argv + [os.environ]
        # the next line replaces the currently-running process with the sudo
        os.execlpe('sudo', *args)

    install_refactor()
	
