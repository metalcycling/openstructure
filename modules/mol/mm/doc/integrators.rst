Integrators
================================================================================

.. currentmodule:: ost.mol.mm

.. class:: Integrator

  The :class:`Integrator` are a wrapper around the integrators provided by OpenMM and
  have to be attached to the :class:`Settings` when setting up a simulation.
  They define the algorithm used to integrate the equations of motion of the system.



Verlet Dynamics
--------------------------------------------------------------------------------

.. class:: VerletIntegrator(step_size)

  Implementation of Verlet dynamics

  :param step_size:     Defining timestep in ps
  :type step_size:      :class:`float`


.. class:: VariableVerletIntegrator(error_tol)

  Implementation of Verlet dynamics with variable time steps

  :param error_tol:     Error tolerance
  :type error_tol:      :class:`float`

  .. method:: GetErrorTolerance()

    :returns:           :class:`float`

  .. method:: SetErrorTolerance(tol)

    :param tol:         Error tolerance
    :type tol:          :class:`float`


Langevin Dynamics
--------------------------------------------------------------------------------

.. class:: LangevinIntegrator(temperature, friction_coeff, step_size)

  Implementation of the Langevin dynamics

  :param temperature:   Temperature of heat bath in K
  :param friction_coeff: Friction coefficient coupling the
                         system to the heat bath in 1/ps
  :param step_size:     Defining timestep in ps

  :type temperature:    :class:`float`
  :type friction_coeff: :class:`float`
  :type step_size:      :class:`float`

  .. method:: GetTemperature()

    :returns:         :class:`float` - the temperature in K

  .. method:: SetTemperature(temperature)

    :param temperature: Temperature in K
    :type temperature:  :class:`float` 

  .. method:: GetFriction()

    :returns:           :class:`float` - the friction in 1/ps

  .. method:: SetFriction(friction)

    :param friction:    Friction in 1/ps
    :type friction:     :class:`float`

  .. method:: GetRandomNumberSeed()

    :returns:          :class:`int`

  .. method:: SetRandomNumberSeed(seed)

    :param seed:        Random number seed
    :type seed:         :class:`int`




.. class:: VariableLangevinIntegrator(temperature, friction_coeff, error_tol)

  Implementation of Langevin dynamics with variable time steps

  :param temperature:   Temperature of heat bath in K
  :param friction_coeff: Friction coefficient coupling the
                         system to the heat bath in 1/ps
  :param error_tol:     Error tolerance

  :type temperature:    :class:`float`
  :type friction_coeff: :class:`float`
  :type error_tolerance: :class:`float`

  .. method:: GetTemperature()

    :returns:         :class:`float` - the temperature in K

  .. method:: SetTemperature(temperature)

    :param temperature: Temperature in K
    :type temperature:  :class:`float`

  .. method:: GetFriction()

    :returns:           :class:`float` - the friction in 1/ps

  .. method:: SetFriction(friction)

    :param friction:    Friction in 1/ps
    :type friction:     :class:`float`

  .. method:: GetRandomNumberSeed()

    :returns:          :class:`int`

  .. method:: SetRandomNumberSeed(seed)

    :param seed:        Random number seed
    :type seed:         :class:`int`

  .. method:: GetErrorTolerance()

    :returns:           :class:`float`

  .. method:: SetErrorTolerance(tol)

    :param tol:         Error tolerance
    :type tol:          :class:`float`


Brownian Dynamics
--------------------------------------------------------------------------------

.. class:: BrownianIntegrator(temperature, friction_coeff, step_size)

  Implementation of Brownian dynamics

  :param temperature:   Temperature of heat bath in K
  :param friction_coeff: Friction coefficient coupling the
                         system to the heat bath in 1/ps
  :param step_size:     Defining timestep in ps.

  :type temperature:    :class:`float`
  :type friction_coeff:    :class:`float`
  :type step_size:    :class:`float`

  .. method:: GetTemperature()

    :returns:           :class:`float` - the temperature in K

  .. method:: SetTemperature(temperature)

    :param temperature: Temperature in K
    :type temperature:  :class:`float`

  .. method:: GetFriction()

    :returns:           :class:`float` - the friction in 1/ps

  .. method:: SetFriction(friction)

    :param friction:    Friction in 1/ps
    :type friction:     :class:`float`

  .. method:: GetRandomNumberSeed()

    :returns:          :class:`int`

  .. method:: SetRandomNumberSeed(seed)

    :param seed:        Random number seed
    :type seed:         :class:`int`
