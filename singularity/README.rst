OST Singularity
===============

Building Singularity image
--------------------------

In order to build OST Singularity image:

.. code-block:: bash

  cd <OST ROOT>/singularity
  sudo singularity build ost.img Singularity

.. note::

  Running singularity build command requires root permissions (sudo).

One can chose any name for an image. For the purose of this file we will assume
that the image name is ``ost.img``.

Here we only keep the recipe for the most recent version of OpenStructure. To
build an image for a different version, edit the source line (``From:``) in the
recipe or look in the git history for an older recipe.

Available apps
--------------

This container includes the following apps:
 * **OST** - OpenStructure binary
 * **IPython** - OST-powered iPython shell
 * **Notebook** - A Jupyter notebook playground with OST and nglview
 * **lDDT** - The Local Distance Difference Test
 * **Molck** - Molecular checker
 * **ChemdictTool** - Creating or update a compound library

To see the help for each individual app run:

.. code-block:: bash

    singularity run-help --app <APP NAME> <PATH TO OST IMAGE>

Eg.:

.. code-block:: bash

    singularity run-help --app OST ost.img


Facilitating the usage
----------------------

For each of these apps it is useful to create an alias if they will be
frequently used. Eg. to create an alias for IPython app one can run:

.. code-block::

  alias ost_ipython="singularity run --app IPython <PATH TO OST IMAGE>"

Then (in the same terminal window) to invoke IPython app one can just type:

.. code-block::

  ost_ipython

To make the alias permanent put it into your ``.bashrc`` file or whatever file
you use to store the aliases.


The Compound Library
--------------------

You'll have the exact same problem with outdated compound libraries as in the
raw Docker image. You can find more information on that matter in the Docker
section of the documentation: :ref:`docker_compound_lib`.

The same trick of mounting an up to date compound library from the local host into
the container applies. The two relevant commands for Singularity are building
a new library and mount it.

Build a new library:

.. code-block:: bash

  singularity run --app ChemdictTool <IMAGE> create components.cif.gz \
  compounds.chemlib

Run some script with an updated compound library from localhost:

.. code-block:: bash

  singularity run \
  -B <COMPLIB_DIR_LOCALHOST>/compounds.chemlib:/compounds.chemlib \
  --env OST_COMPOUNDS_CHEMLIB=/compounds.chemlib \
  --app PM <IMAGE> my_script.py

<COMPLIB_DIR_LOCALHOST> is the directory that contains the compound lib with the
name compounds.chemlib that you created before. Make sure that everything works
as expected by executing the exact same lines of Python code as described
in the Docker documentation: :ref:`docker_compound_lib`.
