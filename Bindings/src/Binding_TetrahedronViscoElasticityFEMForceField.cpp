/******************************************************************************
*  THE SOFA VISCOELASTIC PLUGIN.                                              *
*                                                                             * 
* DESCRIPTION:                                                                *
* This plugin is made for the Simulation Open-Framework Architecture (SOFA)   *
* (c) 2006 INRIA, USTL, UJF, CNRS, MGH.                                       *
* The plugin consist in a Visco-Elastic force field for tetrahedral meshes.   * 
* Several rheological models are implemented.                                 *
*                                                                             *
* CONTRIBUTORS:                                                               *         
* The plugin is made by the collaboration beween the Robotics and Multibody   * 
* Mechanics Department (R&MM) Vrije Universiteit Brussel (VUB), Bruxelles     *
* Belgium, and the DEFROST Team of the INRIA - Lille, France.                 *
*                                                                             *
*                                                                             *
* LICENSE:                                                                    *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Author: Pasquale Ferrentino                                                 *
*                                                                             *
* Contact information: pasquale.ferrentino@vub.be                             *
******************************************************************************/

#include <SofaPython3/Sofa/Core/Binding_Base.h>
#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_BaseObject.h>
#include <Binding_TetrahedronViscoElasticityFEMForceField.h>
#include <SofaViscoElastic/TetrahedronViscoelasticityFEMForceField.h>

using namespace sofa::SofaViscoElastic;

/// Makes an alias for the pybind11 namespace to increase readability.
namespace py { using namespace pybind11; }
typedef TetrahedronViscoelasticityFEMForceField<sofa::defaulttype::Vec3Types> Myclass;
namespace sofapython3 {

// Coord
py::object getShapeVector(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_shapeVector;

    return pybind11::cast(tmp);
}

// Coord
py::object getFiberDirection(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_fiberDirection;

    return pybind11::cast(tmp);
}


// Real
py::object getVolume(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_volume;

    return pybind11::cast(tmp);
}

// Real
py::object getRestVolume(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_restVolume;

    return pybind11::cast(tmp);
}

// Real
py::object getVolScale(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_volScale;

    return pybind11::cast(tmp);
}


//Matrix 3
py::object getF(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_deformationGradient;

    return pybind11::cast(tmp);
}


//Vec 6
py::object getCauchyStress(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_CauchyStress;

    return pybind11::cast(tmp);
}

// Real
py::object getVonMisesStress(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_VonMisesStress;

    return pybind11::cast(tmp);
}


void moduleAddTetrahedronViscoelasticityFEMForceField(py::module &m)
{
    PythonFactory::registerType<Myclass>([](sofa::core::objectmodel::Base* object){
        return py::cast(dynamic_cast<Myclass*>(object));
    });

    py::class_<Myclass, sofa::core::objectmodel::BaseObject, py_shared_ptr<Myclass>> p(m, "TetrahedronViscoelasticityFEMForceField");
    
    p.def("getShapeVector", getShapeVector, "get the shape vector of a tetrahedron");
    p.def("getFiberDirection", getFiberDirection, "get the fiber direction of a tetrahedron");
    p.def("getVolume", getVolume, "get the volume of a tetrahedron");
    p.def("getRestVolume", getRestVolume, "get the rest volume of a tetrahedron");
    p.def("getVolScale", getVolScale, "get the volume Scale of a tetrahedron");
    p.def("getF", getF, "get the deformation gradient of a tetrahedron");
    p.def("getCauchyStress", getCauchyStress, "get the Cauchy Stresses of a tetrahedron");
    p.def("getVonMisesStress", getVonMisesStress, "get the Von Mises Stresses of a tetrahedron");

}
 
}
