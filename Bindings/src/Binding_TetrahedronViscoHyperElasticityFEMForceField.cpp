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
#include <Binding_TetrahedronViscoHyperElasticityFEMForceField.h>
#include <SofaViscoElastic/TetrahedronViscoHyperelasticityFEMForceField.h>

using namespace sofa::SofaViscoElastic;

/// Makes an alias for the pybind11 namespace to increase readability.
namespace py { using namespace pybind11; }
typedef TetrahedronViscoHyperelasticityFEMForceField<sofa::defaulttype::Vec3Types> Myclass1;
namespace sofapython3 {

// Coord
py::object getShapeVector(Myclass1& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_shapeVector;

    return pybind11::cast(tmp);
}

// Coord
py::object getFiberDirection(Myclass1& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_fiberDirection;

    return pybind11::cast(tmp);
}


// Real
py::object getVolume(Myclass1& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_volume;

    return pybind11::cast(tmp);
}

// Real
py::object getRestVolume(Myclass1& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_restVolume;

    return pybind11::cast(tmp);
}

// Real
py::object getVolScale(Myclass1& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_volScale;

    return pybind11::cast(tmp);
}


//Matrix 3
py::object getF(Myclass1& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_deformationGradient;

    return pybind11::cast(tmp);
}

//Vec 6
py::object getSPKStress(Myclass1& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_SPKStress;

    return pybind11::cast(tmp);
}

//Vec 6
py::object getCauchyStress(Myclass1& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_CauchyStress;

    return pybind11::cast(tmp);
}

// Real
py::object getVonMisesStress(Myclass1& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_VonMisesStress;

    return pybind11::cast(tmp);
}


void moduleAddTetrahedronViscoHyperelasticityFEMForceField(py::module &s)
{
    PythonFactory::registerType<Myclass1>([](sofa::core::objectmodel::Base* object){
        return py::cast(dynamic_cast<Myclass1*>(object));
    });

    py::class_<Myclass1, sofa::core::objectmodel::BaseObject, py_shared_ptr<Myclass1>> c(s, "TetrahedronViscoHyperelasticityFEMForceField");
    
    c.def("getShapeVector", getShapeVector, "get the shape vector of a tetrahedron");
    c.def("getFiberDirection", getFiberDirection, "get the fiber direction of a tetrahedron");
    c.def("getVolume", getVolume, "get the volume of a tetrahedron");
    c.def("getRestVolume", getRestVolume, "get the rest volume of a tetrahedron");
    c.def("getVolScale", getVolScale, "get the volume Scale of a tetrahedron");
    c.def("getF", getF, "get the deformation gradient of a tetrahedron");
    c.def("getSPKStress", getSPKStress, "get the Second Piola Kirchhoff Stresses of a tetrahedron");
    c.def("getCauchyStress", getCauchyStress, "get the Cauchy Stresses of a tetrahedron");
    c.def("getVonMisesStress", getVonMisesStress, "get the Von Mises Stresses of a tetrahedron");

}
 
}
