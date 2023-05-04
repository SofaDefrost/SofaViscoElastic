/******************************************************************************
*                              SofaPython3 plugin                             *
*                  (c) 2021 CNRS, University of Lille, INRIA                  *
*                                                                             *
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
* Contact information: contact@sofa-framework.org                             *
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

py::object getVolume(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_volume;

    return pybind11::cast(tmp);
}

py::object getF(Myclass& self, int i)
{
    auto tmp = self.m_tetrahedronInfo.getValue()[i].m_deformationGradient;

    return pybind11::cast(tmp);
}

void moduleAddTetrahedronViscoelasticityFEMForceField(py::module &m)
{
    PythonFactory::registerType<Myclass>(
                [](sofa::core::objectmodel::Base* object)
    {
        return py::cast(dynamic_cast<Myclass*>(object));
    });

    py::class_<Myclass, sofa::core::objectmodel::BaseObject, py_shared_ptr<Myclass>> p(m, "TetrahedronViscoelasticityFEMForceField");

    p.def("getVolume", getVolume, "get the volume of a tetrahedron");
    p.def("getF", getF, "get the SPK stress tensor of a tetrahedron");

}
 
}
    //
    //  for tetra_info in obj.tetra_infos:
    //     print(tetra_info.stress)

    //  for i in range(0, obj.getTetraCount()):
    //     print(obj.getStressSymAt(i))
    //     print(obj.getMatrixSymAt(i))
