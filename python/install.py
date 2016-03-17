#!/usr/bin/env python

import os
import sys
from setuptools import setup

if (sys.version_info > (3, 0)):
    VERSION = 3
else:
    VERSION = 2

if os.name == 'nt':
        
    class FOREGROUND:
        BLACK     = 0x0000
        BLUE      = 0x0001
        GREEN     = 0x0002
        CYAN      = 0x0003
        RED       = 0x0004
        PURPLE    = 0x0005
        YELLOW    = 0x0006
        WHITE     = 0x0007
        GREY      = 0x0008

    class BACKGROUND:
        BLACK     = 0x0000
        BLUE      = 0x0010
        GREEN     = 0x0020
        CYAN      = 0x0030
        RED       = 0x0040
        PURPLE    = 0x0050
        YELLOW    = 0x0060
        GREY      = 0x0080

    import ctypes
    STD_OUTPUT_HANDLE = ctypes.windll.kernel32.GetStdHandle(-11) #windows stdout device 

    def get_csbi_attributes(handle):
        import struct
        csbi = ctypes.create_string_buffer(22)
        res = ctypes.windll.kernel32.GetConsoleScreenBufferInfo(handle, csbi)
        assert res
        (_, _, _, _, wattr, _, _, _, _, _, _) = struct.unpack("hhhhHhhhhhh", csbi.raw)
        return wattr

    RESET = get_csbi_attributes(STD_OUTPUT_HANDLE)


elif os.name == 'posix':
    class FOREGROUND:
        BLACK     = '\033[30m'
        BLUE      = '\033[94m'
        GREEN     = '\033[92m'
        CYAN      = '\033[96m'
        RED       = '\033[91m'
        PURPLE    = '\033[35m'
        YELLOW    = '\033[33m'
        WHITE     = '\033[37m'
        GREY      = '\033[90m'
               
    class BACKGROUND:
        BLACK     = '\033[40m'
        BLUE      = '\033[44m'
        GREEN     = '\033[42m'
        CYAN      = '\033[46m'
        RED       = '\033[41m'
        PURPLE    = '\033[45m'
        YELLOW    = '\033[43m'
        GREY      = '\033[47m'

    RESET ='\033[0m'
    
    # # unix styles :
    #     bold='\033[01m'
    #     disable='\033[02m'
    #     underline='\033[04m'
    #     reverse='\033[07m'
    #     strikethrough='\033[09m'
    #     invisible='\033[08m'



    
REFACTOR_OBLIGATORY_DEPENDENCIES = ['numpy', 'scipy', 'sklearn']
REFACTOR_OPTIONAL_DEPENDENCIES = ['matplotlib','statsmodels'] # for demo only

def setup_refactor():
    """
    Adds refactor.py to the PATH so it can be executed from anywhere.
    NOTE: on windows refactor.py is added to C:\PythonXX\Scripts, but this folder is not in the PATH by default.
    """
    # add "install" command for the setup script. "mimic" user choice: setup install 
    sys.argv.append("install")
    

    # TODO: add  C:\PythonXX\Scripts to registry or environment in a way it will stay there after script terminates
    # could be something like that:
    #
    # aKey = _winreg.OpenKey(aReg, r"SYSTEM\CurrentControlSet\Control\Session Manager\Environment", 0, _winreg.KEY_SET_VALUE)
    # try:   
    #     _winreg.SetValueEx(key = aKey,value_name="Path",type=_winreg.REG_EXPAND_SZ, value= os.path.join(os.path.dirname(sys.executable), "Scripts") ) 
    # except EnvironmentError:                                          
    #     print("Encountered problems writing into the Registry...")
    # _winreg.CloseKey(aKey)
    # _winreg.CloseKey(aReg)        
    
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
    except Exception as e:
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
            color_print("Package %s was installed successfully!" % module_name, FOREGROUND.GREEN)
            return True
            
        except Exception as e:
            pass

    color_print("Package %s auto-installation FAILED" % module_name, FOREGROUND.RED)
    return False


def already_installed(module_name):
    """
    checks if module_name ia already installed by trying to import it
    """
    try:
        __import__(module_name)
        return True

    # if __import__ excepts - module is not installled
    except Exception as e:
        return False

def user_installation_confirmation(modules):
    """
    ask user to install modules
    returns True if user confirmed installation
    returns False otherwise
    """
    print("Packages %s aren't installed" % str(modules))
    request = _input("Do you want to start installation now? (y/n) ")
    while(request.lower() not in ['y', 'n']):
        request = _input("Please select 'y' or 'n': (y/n) ")

    return request.lower() == 'y'


def check_dependencies(dependencies_list):
    """
    if a module in dependencies_list is not installed, tries to install it using pip (if user confirms installation)
    returns the list of dependencies left to be installed
    """
    dependencies_to_install = []

    # check if deoendencies are installed
    print("Checking if required dependencies are installed...")
    for module in dependencies_list:
        if not already_installed(module):
            dependencies_to_install.append(module)

    # ask user to install dependencies. load pip and try to install if user confirms
    if len(dependencies_to_install) > 0:
        if user_installation_confirmation(dependencies_to_install):
            if not load_pip():
                print("pip (installation package) wasn't found. Can't install dependencies.")
            else:
                dependencies_to_install = [module for module in dependencies_to_install if install(module) == False]

    return dependencies_to_install


def install_refactor():
    """
    runs setup.py script to put refactor in PATH
    installes dependencies if they arn't already installed
    """

    print("Installing ReFACTor...")
    setup_refactor()
    
    # install dependencies
    dependencies_to_install = check_dependencies(REFACTOR_OBLIGATORY_DEPENDENCIES + REFACTOR_OPTIONAL_DEPENDENCIES)

    if len(dependencies_to_install) == 0:
        color_print("All dependencies are installed! Done.", FOREGROUND.GREEN)
    else:
        obligatory_dependencies = []
        optional_dependencied = []
        for module in dependencies_to_install:
            if module in REFACTOR_OBLIGATORY_DEPENDENCIES:
                obligatory_dependencies.append(module)
            else:
                optional_dependencied.append(module)
        
        if obligatory_dependencies:
            color_print("To run ReFACTor you must install the following packages: %s" % str(obligatory_dependencies), BACKGROUND.BLACK + FOREGROUND.YELLOW)
        if optional_dependencied:
            color_print("To run ReFACTor's demo.py you must install: %s " % str(optional_dependencied), BACKGROUND.BLACK + FOREGROUND.YELLOW)

def isUserAdmin():

    # check for root on Windows
    if os.name == 'nt':
        try:
            return ctypes.windll.shell32.IsUserAnAdmin()
        except:
            print("Admin check failed, assuming not an admin.")
            return False
    # check for root on Posix
    elif os.name == 'posix':
        return os.getuid() == 0
    else:
        _input("Make sure you're running this as Administrator user. Press any key to continue...")
        return True
        

def color_print(string, color):
    """
    prints 'string' with color 'color'
    NOTE: both on Unix and Windows if you want to add FOREGROUND and BACKGROUND just add them by "+"
          e.g: FOREGROUND.GREEN + BACKGROUND.YELLOW
    """
    if os.name == 'nt':
        ctypes.windll.kernel32.SetConsoleTextAttribute(STD_OUTPUT_HANDLE, color)
        print(string)
        ctypes.windll.kernel32.SetConsoleTextAttribute(STD_OUTPUT_HANDLE, RESET)
    elif os.name == 'posix':
        print(color + string + RESET)
    # color_print on other OS is not supported
    else:
        print(string)

def _input(msg):
    if VERSION == 3:
        return input(msg)
    else:
        return raw_input(msg)

if __name__ == '__main__':

    # run as root (root privilege is required to run python setup.py script and pip install)
    if not isUserAdmin():
        # on Windows
        if os.name == 'nt':
            # TODO: make sure if windows require Admin, if so find a way to change to admin and print a message
            # I think those lines are not working. but should do something like that
            # import subprocess
            # subprocess.call(['runas', '/noprofile', '/user:Administrator' ," ".join(sys.argv)])
            pass
        # on unix    
        elif os.name == 'posix':
            print("Installation must be started as root. Running sudo...")
            args = ['sudo', sys.executable] + sys.argv + [os.environ]
            # the next line replaces the currently-running process with the sudo, won't work the same way on windows
            os.execlpe('sudo', *args)
    
    install_refactor()
    _input("\nPress any key to continue")
	
