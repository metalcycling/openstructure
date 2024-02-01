//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------

#include <boost/python.hpp>
#include <OpenMM.h>

using namespace boost::python;



void export_OpenMM()
{

  //we have to tell boost, that the Integrator class is around...
  class_<OpenMM::Integrator, boost::noncopyable>("Integrator",no_init)
    .def("GetConstraintTolerance",&OpenMM::Integrator::getConstraintTolerance)
    .def("SetConstraintTolerance",&OpenMM::Integrator::setConstraintTolerance,(arg("tolerance")))
  ;

  class_<OpenMM::VerletIntegrator, bases<OpenMM::Integrator> >("VerletIntegrator", init<double>());

  class_<OpenMM::BrownianIntegrator, bases<OpenMM::Integrator> >("BrownianIntegrator", init<double,double,double>())
    .def("GetTemperature",&OpenMM::BrownianIntegrator::getTemperature)
    .def("SetTemperature",&OpenMM::BrownianIntegrator::setTemperature,(arg("temperature")))
    .def("GetFriction",&OpenMM::BrownianIntegrator::getFriction)
    .def("SetFriction",&OpenMM::BrownianIntegrator::setFriction,(arg("friction")))
    .def("GetRandomNumberSeed",&OpenMM::BrownianIntegrator::getRandomNumberSeed)
    .def("SetRandomNumberSeed",&OpenMM::BrownianIntegrator::setRandomNumberSeed,(arg("seed")))
  ;

  class_<OpenMM::LangevinIntegrator, bases<OpenMM::Integrator> >("LangevinIntegrator", init<double,double,double>())
    .def("GetTemperature",&OpenMM::LangevinIntegrator::getTemperature)
    .def("SetTemperature",&OpenMM::LangevinIntegrator::setTemperature,(arg("temperature")))
    .def("GetFriction",&OpenMM::LangevinIntegrator::getFriction)
    .def("SetFriction",&OpenMM::LangevinIntegrator::setFriction,(arg("friction")))
    .def("GetRandomNumberSeed",&OpenMM::LangevinIntegrator::getRandomNumberSeed)
    .def("SetRandomNumberSeed",&OpenMM::LangevinIntegrator::setRandomNumberSeed,(arg("seed")))
  ;

  class_<OpenMM::VariableVerletIntegrator, bases<OpenMM::Integrator> >("VariableVerletIntegrator",init<double>())
    .def("GetErrorTolerance", &OpenMM::VariableVerletIntegrator::getErrorTolerance)
    .def("SetErrorTolerance", &OpenMM::VariableVerletIntegrator::setErrorTolerance,(arg("tolerance")))
  ;

  class_<OpenMM::VariableLangevinIntegrator, bases<OpenMM::Integrator> >("VaribaleLangevinIntegrator", init<double,double,double>())
    .def("GetTemperature",&OpenMM::VariableLangevinIntegrator::getTemperature)
    .def("SetTemperature",&OpenMM::VariableLangevinIntegrator::setTemperature,(arg("temperature")))
    .def("GetFriction",&OpenMM::VariableLangevinIntegrator::getFriction)
    .def("SetFriction",&OpenMM::VariableLangevinIntegrator::setFriction,(arg("friction")))
    .def("GetRandomNumberSeed",&OpenMM::VariableLangevinIntegrator::getRandomNumberSeed)
    .def("SetRandomNumberSeed",&OpenMM::VariableLangevinIntegrator::setRandomNumberSeed,(arg("seed")))
    .def("GetErrorTolerance",&OpenMM::VariableLangevinIntegrator::getErrorTolerance)
    .def("SetErrorTolerance",&OpenMM::VariableLangevinIntegrator::setErrorTolerance,(arg("tolerance")))
  ;
}
