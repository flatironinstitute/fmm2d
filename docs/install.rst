Installation
============

Obtaining fmm2d
***************

The source code can be downloaded from https://github.com/flatironinstitute/fmm2d 


Dependencies
************

This library is supported for unix/linux, Mac OSX, and Windows.

For the basic libraries

* Fortran compiler, such as ``gfortran`` packaged with GCC
* GNU make


Optional:

* for building Python wrappers you will need ``python3``, ``pip3`` and ``numpy``.
* for building standard MATLAB wrappers: MATLAB
* for modifying MATLAB wrappers (experts only): ``mwrap``


Quick install instructions
*********************************************

Make sure you have dependencies installed, and `cd` into your fmm2d
directory. 

-  For linux, run ``make install``.
-  For Mac OSX, run ``cp make.inc.macos.gnu make.inc`` followed by ``make install``.
-  For Windows, run ``cp make.inc.windows.mingw make.inc`` followed by ``make install`` 

This should compile the static library
in ``lib-static/``, the dynamic library in ``lib/`` and copy the dynamic 
library to ``$(HOME)/lib`` on Linux,  to
``/usr/local/lib`` on Mac OSX, and to ``C:\lib`` on Windows.
The location of the default installation directory can be changed by
running::

    make install PREFIX=(INSTALL_DIR)


In order to link against the dynamic library, you will have to update
the ``PATH`` environment variable on Windows, ``LD_LIBRARY_PATH`` environment
variable on Linux and ``DYLD_LIBRARY_PATH`` environment variable on Mac OSX
to the installation directory.
You may then link to the FMM library using the ``-lfmm2d`` option.

.. note :: 
   On MacOSX, /usr/local/lib is included by default in the
   DYLD_LIBRARY_PATH.


To verify successful compilation of the program, run ``make test``
which compiles some fortran test drivers in ``test/`` linked against
the static library, after which it
runs the test programs. The last 17 lines of the terminal output should be::

   cat print_testreshelm.txt
   Successfully completed 27 out of 27 tests in hfmm2d testing suite
   Successfully completed 27 out of 27 tests in hf-hfmm2d testing suite
   Successfully completed 27 out of 27 tests in hfmm2d vec testing suite
   Successfully completed 1 out of 1 tests in helm2d_mps testing suite
   cat print_testreslap.txt
   Successfully completed  9 out of  9 tests in cfmm2d testing suite
   Successfully completed  9 out of  9 tests in cfmm2d vec testing suite
   Successfully completed  9 out of  9 tests in lfmm2d testing suite
   Successfully completed  9 out of  9 tests in lfmm2d vec testing suite
   Successfully completed  9 out of  9 tests in rfmm2d testing suite
   Successfully completed  9 out of  9 tests in rfmm2d vec testing suite
   cat print_testresbh.txt
   Successfully completed  1 out of  1 tests in bhfmm2d testing suite
   rm print_testreshelm.txt
   rm print_testreslap.txt
   rm print_testresbh.txt


To verify successful installation of the program, and the correct
setting for environment variables, run ``make test-dyn`` which compiles
some fortran test drivers in ``test/`` linked against the dynamic
library, after which it runs teh test prgram. The output ofshould be the
same as above.

If ``make test`` fails, see more detailed instructions below.

If ``make test-dyn`` fails with an error about not finding
``-llibfmm2d_dll`` or ``-lfmm2d`` or ``libfmm2d`` make sure that the
appropriate environment variables have been set. If it fails with other
issues, see more detailed instructions below.

Type ``make`` to see a list of other build options (language
interfaces, etc). Please see `Fortran interfaces <fortran-c.html>`__ and look in
``examples/`` for sample drivers.

If there is an error in testing on a standard set-up,
please file a bug report as a New Issue at https://github.com/flatironinstitute/fmm2d/issues

.. _custom-install:

Custom library compilation options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the (default) easy-to-install version,
the library is compiled  without using the optimized direct evaluation kernels.

In order to disable multi-threading, append ``OMP=OFF`` to the make task.


All of these different libraries are
built with the same name, so you will have to move them to other
locations, or build a 2nd copy of the repo, if you want to keep both
versions.

You *must* do at least ``make objclean`` before changing to the openmp
/fast direct kernel evaluation options.


Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``make examples`` to compile and run the examples for calling from Fortran.

The ``examples`` directory is a good place to see usage 
examples for Fortran.
There are two sample Fortran drivers  
for Helmholtz FMMs,
one which demonstrates the use of FMMs, and one which demonstrates
the use of vectorized FMMs. 

The sample Helmholtz drivers are ``hfmm2d_example.f``,
``hfmm2d_vec_example.f``, and ``hfmm2d_legacy_example.f``.
The corresponding makefiles are ``hfmm2d_example.make``, 
``hfmm2d_vec_example.make``, and ``hfmm2d_legacy_example.make``.


Building Python wrappers
****************************

First make sure you have python (version 3 or higher), pip and numpy installed. 

You may then execute ``make python`` (after copying over the
operating system specific make.inc.* file to make.inc) which calls
pip for the install and then runs some tests.

To rerun the tests, you may run ``pytest`` in ``python/`` 
or alternatively run ``python python/test_hfmm.py`` and 
``python python/test_lfmm.py``.

See ``python/hfmmexample.py`` and ``python/lfmmexample.py`` to see
usage examples for the Python wrappers.

.. note::
   On windows, you will need to update ``distutils.cfg`` located in 
   ``(PYTHON_INSTALL_DIR)\Lib\distutils`` and set it to::
 
       [build]
       compiler=mingw32

       [build_ext]
       compiler=mingw32

   which forces python to use the mingw compiler for building its
   modules. In case you wish to revert to using VC/C++ for building python
   modules, make sure to update distutils.cfg appropriately.


A few words about Python environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There can be confusion and conflicts between various versions of Python and installed packages. It is therefore a very good idea to use virtual environments. Here's a simple way to do it (after installing python-virtualenv)::

  Open a terminal
  virtualenv -p /usr/bin/python3 env1
  . env1/bin/activate

Now you are in a virtual environment that starts from scratch. All pip installed packages will go inside the env1 directory. (You can get out of the environment by typing ``deactivate``)

It's advisable to install numpy into a virtual environment using pip as.

  virtualenv -p /usr/bin/python3 env1
  . env1/bin/activate
  pip install numpy


Building the MATLAB wrappers
****************************

First make sure you have MATLAB installed. 

Then run ``make matlab`` (after copying over the operating
system specific make.inc.* file to make.inc) which links the .m files to
the .c file in the matlab folder. 

To run tests, you can run ``matlab test_hfmm3d.m`` and 
``matlab test_lfmm3d.m`` and it should return with $0$ crashes.

Example codes for demonstrating the Helmholtz and Laplace
interfaces are ``hfmm3d_example.m`` and ``lfmm3d_example.m``.

In order to build the MATLAB routines with the optimized direct kernel 
evaluation routines on a Mac, we recommend building mex with gcc
instead of clang. The relevant xml files for configuring mex to use
gcc can be found at https://github.com/danfortunato/matlab-gcc.
Follow the instructions there to configure mex with gcc, and set
CC = ``gcc-<version number>`` in your make.inc file. 


Tips for installing dependencies
**********************************

On Ubuntu linux
~~~~~~~~~~~~~~~~

On Ubuntu linux (assuming python3 as opposed to python)::

  sudo apt-get install make build-essential gfortran  


On Fedora/CentOS linux
~~~~~~~~~~~~~~~~~~~~~~~~

On a Fedora/CentOS linux system, these dependencies can be installed as 
follows::

  sudo yum install make gcc gcc-c++ gcc-gfortran libgomp 

.. _mac-inst:

On Mac OSX
~~~~~~~~~~~~~~~~~~~~~~~~

First setup Homebrew as follows. If you don't have Xcode, install
Command Line Tools by opening a terminal (from /Applications/Utilities/)
and typing::

  xcode-select --install

Then install Homebrew by pasting the installation command from
https://brew.sh

Then do::
  
  brew install gcc 


On Windows
~~~~~~~~~~~~~~~

Download 64 bit mingw (available `here <http://mingw-w64.org/doku.php>`_). 
Follow the install instructions and append to the environment variable ``PATH`` the
location of the bin directory of your mingw installation.

Download  and install ``make`` for windows 
(Available `here <http://gnuwin32.sourceforge.net/packages/make.htm>`_).

Download and install ``git`` for windows
(Available `here <https://git-scm.com/download/win>`_).


Tips for installing optional dependencies
******************************************

Installing python and pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On Ubuntu linux
##################

::

  sudo apt-get install python3 python3-pip


On Mac OSX
############

Make sure you have homebrew installed. See `Tips for installing dependencies -> On Mac OSX <install.html#mac-inst>`__ 

::
  
  brew install python3


On Windows
###########

Download and install python3.7 from python.org.

Configuring MATLAB
~~~~~~~~~~~~~~~~~~~

On Windows
############

Update ``MINGW_LPATH`` in ``make.inc.windows.mingw`` to point to the
appropriate installation directory (it should be the one within the
``gcc`` folder).

To setup mingw as the C compiler on MATLAB run ``configuremingw.p``
(which can be downloaded from 
`here <https://www.mathworks.com/matlabcentral/answers/uploaded_files/88639/configuremingw.p>`_)
and choose the mingw directory. To verify successful setup run ``mex
-setup`` from matlab and it should be configured to compile with mingw.


Installing MWrap
~~~~~~~~~~~~~~~~~~

If you make any changes to the 
fortran code, you will need to regenerate the .c files
from the .mw files for which mwrap is required.
This is not needed for most users.
`MWrap <http://www.cs.cornell.edu/~bindel/sw/mwrap>`_
is a very useful MEX interface generator by Dave Bindel.

Make sure you have ``flex`` and ``bison`` installed.
Download version 0.33.5 or later from https://github.com/zgimbutas/mwrap, un-tar the package, cd into it, then::
  
  make
  sudo cp mwrap /usr/local/bin/


