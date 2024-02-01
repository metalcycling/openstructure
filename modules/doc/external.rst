Using External Programs within OpenStructure
================================================================================

Introduction
--------------------------------------------------------------------------------

It is often very useful to use external programs to do a specific task. In principle, this can be done by writing out files from OpenStructure and manually running an external program, however, for convenience, this can also be done directly from within OpenStructure using Python commands. 

This tutorial will give you some hints how to do this for a new external program. The process basically consists of four steps:

  * locate the executable of the external program
  * prepare all necessary files
  * execute the external program from python
  * read in generated output


Locating the Executable
--------------------------------------------------------------------------------

There is a helper function available to locate files, and especially executables: :func:`~ost.settings.Locate`. Using this, you can obtain the full path of an executable.

As an example, we would like to obtain the full path of the msms executable (a program to calculate molecular surfaces):

.. code-block:: python

  from ost import settings
  exe_path = settings.Locate('msms', search_paths=['/opt/app','/home/app'],
              env_name='MSMS', search_system_paths=True)
  print(exe_path)
  
The :func:`~ost.settings.Locate` command looks for the program with the name 
`msms`. If env_name is set, it first looks if an environment variable with the 
name `MSMS` is set. If not, all paths in *search_paths* are searched. If the 
executable could still not be found and *search_system_paths* is set to True, 
the binary search paths are searched. If the executable could not be found, a 
:exc:`~ost.FileNotFound` exception is raised with a detailed description where 
Locate was searching for the executable.
    
Prepare All Files
--------------------------------------------------------------------------------

The preparation of the necessary files is very dependent on the external 
program. Often it is useful to generate a temporary directory or file. For 
this, the Python module tempfile is very handy.

An example how to generate a temporary directory, open a file in this directory and write the position and radius of all atoms into this file is shown here:

.. code-block:: python

  import tempfile
  import os
  
  # generate a temporary directory
  tmp_dir_name = tempfile.mkdtemp()
  print('temporary directory:',tmp_dir_name)
  
  # generate and open a file in the temp directory
  tmp_file_name = os.path.join(tmp_dir_name,"entity")
  tmp_file_handle = open(tmp_file_name, 'w')
  print('temporary file:',tmp_file_handle)
  
  # write position and radius of all atoms to file
  for a in entity.GetAtomList():
    position = a.GetPos()
    tmp_file_handle.write('%8.3f %8.3f %8.3f %4.2f\n' % (position[0],
                          position[1], position[2], a.GetProp().radius))
                          
  # close the file
  tmp_file_handle.close()

Execute the External Program
--------------------------------------------------------------------------------

The external program can be executed from python using the python module subprocess.

To run the external program msms from the above example, with the temporary file generated before, we can use the following:

.. code-block:: python

  import subprocess

  # set the command to execute
  command = "%s -if %s -of %s" % (exe_path,
            tmp_file_name, tmp_file_name)
  print('command:',command)

  # run the executable with the command
  proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  stdout_value, stderr_value = proc.communicate()

  # check for successful completion of msms
  if proc.returncode != 0:
    print("WARNING: msms error\n", stdout_value)
    raise subprocess.CalledProcessError(proc.returncode, command)

  # print everything written to the command line (stdout)
  print(stdout_value)
    
Read Generated Output
--------------------------------------------------------------------------------

The last step includes reading of generated files (like in the case of msms) and/or processing of the generated command line output.

Here we first print the command line output and then load the generated msms surface and print the number of vertex points:

.. code-block:: python

  # print everything written to the command line (stdout)
  print(stdout_value)
  
  # read msms surface from file
  surface = io.LoadSurface(tmp_file_name, "msms")
  print('number of vertices:',len(surface.GetVertexIDList()))
