Observers
================================================================================

.. currentmodule:: ost.mol.mm

Observers can be registered to a :class:`Simulation` and get called at a 
defined interval.


Trajectory observers
--------------------------------------------------------------------------------

.. class:: TrajObserver(rhythm) 

  Stores and updates the md trajectory produced by a :class:`Simulation`
  in memory

  :param rhythm:        Frequency, at which the observer gets
                        notified by the simulation
  :type rhythm:         :class:`int`

  .. method:: GetTraj()

    Get current Trajectory

    :returns:           :class:`ost.mol.CoordGroupHandle`

  .. method:: Rhythm()

    :returns:           :class:`int` calling rythm


.. class:: TrajWriter(rythm, pdb_filename, dcd_filename)

  Writes the MD trajectory produced by a :class:`Simulation` directly to
  disk in dcd format 

  :param rythm:         Frequency, at which the observer gets
                        notified by the simulation
  :param pdb_filename:  Name of produced pdb file
  :param dcd_filename:  Name of produced dcd file

  :type rythm:          :class:`int`
  :type pdb_filename:   :class:`str`
  :type dcd_filename:   :class:`str`

  .. method:: Rhythm()

    :returns:           :class:`int` calling rhythm

  .. method:: Finalize()

    Obligatory function to be called after the last simulation step to
    write missing information into the dcd header

