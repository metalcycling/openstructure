OST Docker
==========

.. note::

  For many docker installations it is required to run docker commands as root. As
  this depends on set up, we skip the ``sudo`` in all commands.

Obtain Docker image from the OST registry
-----------------------------------------

OST has its own container registry inside GitLab. There we try to keep one
image for the latest stable version of OST. You can import it by

.. code-block:: bash

  docker pull registry.scicore.unibas.ch/schwede/openstructure:<TAG>

and just start using it without the overhead to build it yourself.


Build Docker image
------------------

In order to build OST image:

.. code-block:: bash

  cd <PATH TO OST>/docker
  docker build --tag <IMAGE NAME> -f Dockerfile .

or if you downloaded the Dockerfile directly:

.. code-block:: bash

  docker build --tag <IMAGE NAME> -f <DOCKERFILE NAME> <PATH TO DOCKERFILE DIR>

You can chose any image name (tag) eg. ost.
Here we only keep the recipe for the most recent version of OpenStructure. To
build an image for a different version, you can either adapt the
``OPENSTRUCTURE_VERSION`` variable in the recipe or look in the git history for
an older recipe.

Testing the image
-----------------

One can find a exemplary script (``test_docker.py``) in the downloaded directory.
To run it do:

.. code-block::

  cd <PATH TO OST>/docker
  docker run --rm -v $(pwd):/home <IMAGE NAME> test_docker.py

As the last line you should see ``OST is working!``.

Run script and action with OST
------------------------------

.. note::

  If script or action requires some external files eg. PDBs, they have to be located in the
  path accessible via mounted volume and should be accessed via docker (NOT LOCAL)
  path. Eg. assuming that we have a struc.pdb file in /home/user/pdbs directory and
  a script.py in /home/user we could mount the /home/user to /home in docker as
  above by specifying -v /home/user:/home. To run the script we thus need to
  provide the (relative) path to the script and (relative) path to the file eg:

  .. code-block:: bash

    docker run --rm -v /home/user:/home <IMAGE NAME> script.py pdbs/struct.pdb

  or with absolute paths:

  .. code-block:: bash

    docker run --rm -v /home/user:/home <IMAGE NAME> /home/script.py /home/pdbs/struct.pdb
  
  An easy solution to mount a CWD is to use $(pwd) command in the -v option
  of the Docker. For an example see the action exemplary run.
  The same reasoning is valid for the output files.

Actions
#######

To see the list of available actions do:

  .. code-block::

    docker run --rm <IMAGE NAME> -h

To run chosen action do:

  .. code-block::

    docker run --rm <IMAGE NAME> <ACTION NAME>

 
Here is an example run of the compare-structures action:

.. code-block::

  docker run --rm -v $(pwd):/home <IMAGE NAME> compare-structures \
      --model model.pdb \
      --reference reference.pdb \
      --output output.json \
      --qs-score \
      --residue-number-alignment \
      --lddt \
      --structural-checks \
      --consistency-checks \
      --inclusion-radius 15.0 \
      --bond-tolerance 15.0 \
      --angle-tolerance 15.0 \
      --molck \
      --remove oxt hyd unk \
      --clean-element-column \
      --map-nonstandard-residues


In order to see all available options for this action run:

.. code-block::

  docker run --rm <IMAGE NAME> compare-structures -h

Scripts
#######

In order to run OST script do:

.. code-block:: bash

  docker run [DOCKER OPTIONS] --rm -v <PATH TO SCRIPT DIR>:/home <IMAGE NAME> /home/<SCRIPT NAME> [SCRIPT OPTIONS]

Run ost with utility command
###############################

One can also use provided utility bash script ``run_docker_ost`` to run basic
scripts and actions:

.. code-block:: bash

  <PATH TO OST>/docker/run_docker_ost <IMAGE_NAME> [<SCRIPT_PATH>] [SCRIPT OPTIONS]

One just needs to provide image name and optionally a script/action and its
options. It is useful to link the command to the binary directory eg. in linux:

.. code-block:: bash

  ln -s <PATH TO OST>/docker/run_docker_ost /usr/bin/run_docker_ost

In order to run an exemplary script (``test_docker.py``) do:

.. code-block::

  cd <PATH TO OST>/docker
  ./run_docker_ost <IMAGE NAME> test_docker.py

To see the help for compare-structures action run:

.. code-block::

  cd <PATH TO OST>/docker
  ./run_docker_ost <IMAGE NAME> compare-structures


Running other commands
----------------------

The default entrypoint of the Docker image is "ost" thus in order to run other
available commands (or other commands in general) one need to override
the entrypoint:

.. code-block::

  docker run --rm -ti --entrypoint <COMMAND> <IMAGE NAME> [COMMAND OPTIONS]

Eg. to run molck type:

.. code-block::

  docker run --rm -ti --entrypoint molck <IMAGE NAME> --help

.. note::

  Note how the options to the command are specified after the image name.
