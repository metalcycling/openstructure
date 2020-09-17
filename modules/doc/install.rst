Installing OpenStructure From Source
================================================================================

Brief Overview
--------------------------------------------------------------------------------

For a simple and portable way to use OpenStructure we recommend using a
container solution. We provide recipes to build images for
`Docker <https://www.docker.com/>`_ and
`Singularity <https://sylabs.io/singularity/>`_.
The latest recipes and instructions can be found on our
`GitLab site <https://git.scicore.unibas.ch/schwede/openstructure/>`_, including
a link to OpenStructure's own `GitLab Docker registry <https://git.scicore.unibas.ch/schwede/openstructure/container_registry>`_ (`Docker instructions`_ and
`Singularity instructions`_).

If you wish to compile OpenStructure outside of a container, you need to follow
the steps which we describe in detail below. In essence, these steps are:

* Installing the Dependencies
* Checking out the source code from GIT
* Configuring the build with cmake
* Compiling an Linking


Installing the Dependencies
--------------------------------------------------------------------------------

OpenStructure requires a C++11 enabled compiler (e.g. recent gcc/clang) and uses 
a bunch of open-source libraries. If you haven't already installed them, please 
install them now! Where appropriate, the minimally required version is given in 
parentheses.

* `CMake <http://cmake.org>`_ (3.10.2)
* `Python3 <http://python.org>`_ (3.6)
* `Boost <http://boost.org>`_ (1.65)
* `zlib <https://zlib.net/>`_ (usually comes with Boost or system)
* `Eigen3 <http://eigen.tuxfamily.org>`_ (3.2.0)
* `SQLite3 <https://www3.sqlite.org>`_ (3.7.13)

When you enable support for image processing, you will need:

* `FFTW3 <http://fftw.org>`_. By default, OpenStructure is compiled with single
  precision and thus also requires FFTW to be compiled with single precision.
  Most platforms offer this as a second package. If you are compiling manually,
  use the `--enable-single` option.
* `libtiff <http://www.libtiff.org>`_
* `libpng <http://www.libpng.org>`_ (also needed for GUI)

If you would like to use the info module, also install:

* `Qt5 <http://qt-project.org/>`_ 

If you would like to use the graphical user interface (GUI), also install:

* `SIP <http://www.riverbankcomputing.co.uk/software/sip/download>`_
* `PyQt5 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`_

If you would like to use the :mod:`molecular mechanics <ost.mol.mm>` module:

* `OpenMM <https://simtk.org/home/openmm>`_ (6.1)

We do not provide backwards compatibility to Python 2.7. The last
release supporting Python 2.7 is 1.11.0.


Getting the Source Code
--------------------------------------------------------------------------------

OpenStructure uses `git` as the revision control system. The main repository can
be browsed `here <https://git.scicore.unibas.ch/schwede/openstructure.git>`_. To
get the source code, use git clone:

.. code-block:: bash

  git clone https://git.scicore.unibas.ch/schwede/openstructure.git <directory-name>
  
The above command will clone OpenStructure into the directory specified by
`<directory-name>`. If omitted, the directory will be called openstructure. 

.. note::

  Some versions of curl have trouble with the certificate of the OpenStructure
  git server and fail to clone the repository. To work around this, disable the
  SSL certificate verification by setting the following environment variable:
  
  .. code-block:: bash

    export GIT_SSL_NO_VERIFY=1


Picking the right branch
--------------------------------------------------------------------------------

By default you are checking out the master branch. Master is by definition a
stable branch. It always points to the latest release. However, there are
several other branches at your disposal. The main development is happening in
the develop branch. It contains the newest features and bug fixes. However, we
don't make any guarantees that the develop branch is bug free and doesn't
contain major bugs. After all, it's in constant flux. If you are developing new
features, start your feature branch off develop. Besides that, there are several
smaller features branches that are used to group together commits for one
specific features. To change to a specific branch, use

.. code-block:: bash

  git checkout <branch-name>


Configuring
--------------------------------------------------------------------------------

OpenStructure uses `CMake <http://cmake.org>`_ for compiling and building the
project. The next required step is to configure the build environment using
cmake. You can do that by invoking `cmake` in the project directory.

.. code-block:: bash

  cmake . <options>

There are two kinds of options: Options that let you control the building
behaviour, enabling and disabling the compilation of certain modules and options
that let you tell CMake where to find the dependencies. All of them are passed
to CMake via `-D<opt>=<value>`.


Flag to choose build generator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CMake supports different build generators. On UNIX, i.e. macOS and Linux, the
default build generator is Makefiles, but it is also possible to use other
programs. For a list of supported build generators on your platform, run
`cmake` without parameters.


.. _cmake-flags:

Flags to Control the Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, `CMake <http://cmake.org>`_ searches the standard directories for
dependencies. However, on some systems, this might not be enough. Here is a
short description of how CMake figures out what dependencies to take and how you
can influence it.

* Boost is mainly controlled via the `BOOST_ROOT` option. If boost wasn't
  found, it should be set to the prefix of the boost installation. If for some
  reason, it is desirable to use the non-multithreaded boost libraries, you can
  switch `Boost_USE_MULTITHREADED` off (it is on by default).

* `PYTHON_ROOT` is the Python equivalent of BOOST_ROOT. It should be set to 
  the prefix path containing the python binary, headers and libraries.

* `SYS_ROOT` controls the general prefix for searching libraries and headers.
  By default, it is set to `/`.
  
* `COMPOUND_LIB` specifies the location of the compound library and
  activates the rule-based-builder. The compound library is based on 
  the component dictionary released by the PDB, and it specifies atoms
  of a certain residue or connectivities between atoms etc. The 
  :doc:`compound library <conop/compoundlib>` itself is created from the 
  component dictionary by calling the OpenStructure chemdict_tool. 
  By default this is switched off but it is highly recommended to provide a
  compound library to use all features of OpenStructure.

* `ENABLE_GUI` controls whether to build the graphical user interface module.
  By default, this is switched on.

* `ENABLE_GFX` controls whether to build the graphics module. By default, this
  is switched on. If it is switched off, it also switches `ENABLE_GUI` off.

* `ENABLE_INFO` controls whether to build the info module. By default, this is
  switched on. If it is switched off, it also switches `ENABLE_GFX` off and
  removes all dependencies to Qt.

* `QT_QMAKE_EXECUTABLE` defines the exact Qt installation to take. It should 
  be set to the full path to `qmake`. This is only needed if `ENABLE_INFO` is
  switched on.

* `COMPILE_TMTOOLS` will activate bindings for TMAlign and TMScore, which are 
  then available at python level. This option requires a Fortran compiler. 
  By default, this option is switched off.

* `USE_NUMPY` allows OpenStructure to pass back data in NumPy format. By 
  default, this is switched off.

* `ENABLE_IMG` controls whether to build the image processing module. This will
  enable support for density maps, and general image processing in 1, 2 an 3
  dimensions. By default, this is switched on.

* `ENABLE_MM` controls whether the molecular mechanics module is enabled. By
  default, this is switched off. If it is turned on, you should also set the
  paths to your local OpenMM installation:

  * `OPEN_MM_INCLUDE_DIR`: the include path
  * `OPEN_MM_LIBRARY`: the libOpenMM library
  * `OPEN_MM_PLUGIN_DIR`: the path for OpenMM plugins
  * see example below for commonly used paths

* Several paths to other libraries can be set if they are not in the expected
  locations:

  * `PYTHON_LIBRARIES` defines the location of the Python library (file name
    starting with `libpython`). This must be set if it is not in
    `$PYTHON_ROOT/lib`.
  * `EIGEN3_INCLUDE_DIR` defines the include folder of Eigen3 (contains `Eigen`
    folder with include files).
  * `FFTW_LIBRARY` defines the location of the FFTW3 library (file name starting
    with `libfftw3f` (or `libfftw3` if `USE_DOUBLE_PRECISION` is switched on))
  * `FFTW_INCLUDE_DIR` defines the include folder of FFTW3 (contains include
    files directly)
  * `PNG_LIBRARY` defines the location of the libpng library (file name starting
    with `libpng`)
  * `PNG_PNG_INCLUDE_DIR` defines the include folder of libpng (contains include
    files directly)
  * `ZLIB_LIBRARY` defines the location of the zlib library (file name starting
    with `libz`)
  * `ZLIB_INCLUDE_DIR` defines the include folder of zlib (contains include
    files directly)
  * `TIFF_LIBRARY` defines the location of the libtiff library (file name
    starting with `libtiff`)
  * `TIFF_INCLUDE_DIR` defines the include folder of libtiff (contains include
    files directly)
  * `SQLITE3_LIBRARY` defines the location of the SQLite3 library (file name starting
    with `libsqlite3`)
  * `SQLITE3_INCLUDE_DIR` defines the include folder of SQLite3 (contains include
    files directly)
  * Usually, you will receive errors for those variables when executing `cmake`
    and set them accordingly as needed.

* `OPENGLPREFERENCE_LEGACY` switches the GL implementation to be used by OpenGL.
  The default is what should be used on modern systems. But since there are some
  reports on the internet claiming that the default does not work everywhere,
  this switch enables the usage of the legacy implementation of GL.
  
Build Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `OPTIMIZE` can be switched on to build an optimised (-O3 -DNDEBUG) version of
  OpenStructure. By default, this is switched off.

* `USE_DOUBLE_PRECISION` will switch on double precision within OpenStructure. 
  By default, this is switched off.

* `ENABLE_STATIC` allows some parts of OpenStructure to be statically linked 
  and thus can be used more easily across a heterogeneous setup, e.g. older 
  systems and newer systems. Note that enabling this flag will not compile the
  full OpenStructure package and it is not guaranteed to lead to fully portable
  binaries. By default, this is switched off.

* For deployment of OpenStructure with `make install` there are two relevant
  settings to consider:

  * `PREFIX` or `CMAKE_INSTALL_PREFIX` are used to define the path where the
    OpenStructure `stage` folder will be installed to.
  * `USE_RPATH` can be switched on to embed rpath upon make install. By default,
    this option is switched off.

* Experimental settings (only change if you know what you are doing):

  * `USE_SHADER` controls whether to compile with shader support. By default,
    this is turned off.
  * `ENABLE_SPNAV` controls whether 3DConnexion devices should be supported. By
    default, this is turned off.
  * `PROFILE` can be switched on to enable a (very verbose) code profiler. By
    default, this is turned off.
  * `UBUNTU_LAYOUT` can be turned on to switch the directory layout of the
    `stage` folder to be more ubuntu-like. By default, this is switched off.
  * `HIDDEN_VISIBILITY` can be turned on to add "-fvisibility=hidden" to gcc's
    compile flags (only if GNU compiler used). By default, this is switched off.


Example Configurations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Generic Linux without GUI**

The simplest way to compile OpenStructure is to disable the GUI and any
dependency to Qt5. You can build an optimised OpenStructure without GUI as
follows:

.. code-block:: bash

  cmake . -DOPTIMIZE=ON -DENABLE_INFO=OFF

The molecular mechanics module can be enabled by installing OpenMM and adding
the appropriate flags as follows (replace `<OPENMM>` with the actual path to
OpenMM):

.. code-block:: bash

  cmake . -DOPTIMIZE=ON -DENABLE_INFO=OFF -DENABLE_MM=ON \
          -DOPEN_MM_LIBRARY=<OPENMM>/lib/libOpenMM.so \
          -DOPEN_MM_INCLUDE_DIR=<OPENMM>/include/ \
          -DOPEN_MM_PLUGIN_DIR=<OPENMM>/lib/plugins

Note that the OpenMM binaries available online may be incompatible with files
compiled using your gcc compiler (known as "Dual ABI" issue). This has been
observed for OpenMM versions 6.1 until 7.1.1 when compiling with gcc versions >=
5.1. In those cases, you cannot use the binaries and will have to install OpenMM
from source.


**Ubuntu 20.04 LTS / Debian 10 with GUI**

All the dependencies can be installed from the package manager as follows:

.. code-block:: bash

  sudo apt-get install cmake g++ libtiff-dev libfftw3-dev libeigen3-dev \
               libpng-dev python3-all python3-pyqt5 libboost-all-dev \
               qt5-qmake qtbase5-dev libpng-dev libsqlite3-dev

Now, all dependencies are located in standard locations and cmake will
automatically find them without the need to pass any additional parameters. The
only exception is the Python library which is put in a different path than
expected. Also, we add -DOPTIMIZE, which will tell cmake to build an optimised
version of OpenStructure.

.. code-block:: bash

  cmake . -DPYTHON_LIBRARIES=/usr/lib/x86_64-linux-gnu/libpython3.8.so \
          -DOPTIMIZE=ON

Be careful at -DPYTHON_LIBRARIES, Debian 10 comes with Python 3.7 so that needs
to be substituted (libpython3.8.so -> libpython3.7m.so).


**macOS (Catalina) with Homebrew**

.. note::

  When switching the Qt version used to compile OST with support for the
  graphical user interface, dng may start behaving weird. Symptoms are that the
  user interface starts being unresponsive to mouse clicks. An easy solution
  may be to close dng and remove
  $HOME/Library/Preferences/org.openstructure.dng.plist and start dng again.

`Homebrew <https://brew.sh/>`_ can be used to conveniently install all
dependencies. The current Python version, as of writing these instructions, is
3.8.5 but works so far. Boost comes as 1.72.0 which seems to be OK. Do not
forget to also install `boost-python3` (your system may have a lower version of
Python than 3.8.5 but it seems like `boost-python` was compiled for 3.8.5).
Eigen and SQLite also seem to be unproblematic concerning higher version numbers.

If you want to build the info module or the graphical user interface, make sure
you have the Xcode app installed. Just the Xcode command line tools which are
sufficient for Homebrew, will not work with Qt5.

Before running CMake, some environment variables need to be set on the command
line. If omitted, the linker will throw a bunch of warnings later:

.. code-block:: bash

  export SDKROOT=/Applications/Xcode.app/Contents/Developer/Platforms/\
  MacOSX.platform/Developer/SDKs/MacOSX.sdk

If building the info module or with graphical user interface, get the Qt
binaries in your Path for CMake to determine its configuration:

.. code-block:: bash

  export PATH="/usr/local/opt/qt/bin:$PATH"

Homebrew installs all the software under /usr/local. Thus we have to tell cmake
where to find Boost and Python. Also the Python headers and libraries are not
located as they are on Linux and hence they must be specified too. To get rid of
a ton of compilation warnings from third party software, we add some dedicated
C flags:

.. code-block:: bash

  cmake . -DPYTHON_INCLUDE_PATH=/usr/local/opt/python@3.8/Frameworks/\
  Python.framework/Versions/Current/include/python3.8/ \
          -DPYTHON_LIBRARIES=/usr/local/opt/python@3.8/Frameworks/\
  Python.framework/Versions/Current/lib/libpython3.8.dylib \
          -DPYTHON_ROOT=/usr/local/opt/python@3.8/ \
          -DBOOST_ROOT=/usr/local \
          -DSYS_ROOT=/usr/local \
          -DOPTIMIZE=ON \
          -DCMAKE_C_FLAGS="-isystem /Applications/Xcode.app/Contents/\
  Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/\
  Library/Frameworks/OpenGL.framework/Headers/ -isystem /usr/local/opt/\
  qt/lib/QtCore.framework/Headers/ -isystem /usr/local/opt/qt/lib/\
  QtWidgets.framework/Headers/ -isystem /Applications/Xcode.app/\
  Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/\
  MacOSX.sdk/System/Library/Frameworks/Security.framework/ \
  -isystem /usr/local/opt/qt/lib/QtGui.framework/Headers/" \
          -DCMAKE_CXX_FLAGS="-isystem /Applications/Xcode.app/Contents/\
  Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/\
  Library/Frameworks/OpenGL.framework/Headers/ -isystem /usr/local/opt/\
  qt/lib/QtCore.framework/Headers/ -isystem /usr/local/opt/qt/lib/\
  QtWidgets.framework/Headers/ -isystem /Applications/Xcode.app/\
  Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/\
  MacOSX.sdk/System/Library/Frameworks/Security.framework/ \
  -isystem /usr/local/opt/qt/lib/QtGui.framework/Headers/"

Building the Project
--------------------------------------------------------------------------------

Type ``make``. If you are using a multi-core machine, you can use the `-j` flag
to run multiple jobs at once.


What's next?
--------------------------------------------------------------------------------

On Linux and macOS, you can start dng from the command-line. The binaries are
all located in stage/bin:

.. code-block:: bash

  stage/bin/dng
  
or, to start the command-line interpreter:

.. code-block:: bash

  stage/bin/ost
  
If you repeatedly use OpenStructure, it is recommended to add
/path/to/ost/stage/bin to your path.

Getting the newest changes
--------------------------------------------------------------------------------

To get the newest changes from the central git repository, enter

.. code-block:: bash

  git pull

in your terminal. This will fetch the newest changes.


..  LocalWords:  Homebrew cmake CMake zlib SQLite FFTW libtiff libpng PyQt
..  LocalWords:  SSL macOS Makefiles PDB qmake PNG libz libsqlite OPTIMIZE
..  LocalWords:  DNDEBUG RPATH rpath SHADER shader SPNAV DConnexion profiler
..  LocalWords:  DOPTIMIZE DENABLE DOPEN DPYTHON DBOOST DSYS Xcode Eigen
