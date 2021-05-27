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
runs the test programs. The last 4 lines of the terminal output should be::

   cat print_testreshelm.txt
   Successfully completed 27 out of 27 tests in hfmm2d testing suite
   Successfully completed 27 out of 27 tests in hfmm2d vec testing suite
   rm print_testreshelm.txt


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


