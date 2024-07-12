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

#include <pybind11/eval.h>
namespace py = pybind11;
#include <Binding_TetrahedronViscoElasticityFEMForceField.h>
#include <Binding_TetrahedronViscoHyperElasticityFEMForceField.h>

namespace sofapython3
{

PYBIND11_MODULE(SofaViscoElastic, m) {
    m.doc() = R"doc(
              Binding for the SofaViscoElastic plugin
              -----------------------------------

              Provides python bindings for the SofaViscoElastic module

              Example of use:

              .. code-block:: python

                import SofaViscoElastic

              .. autosummary::
                  :toctree: _autosummary/_autosummary

                  SofaViscoElastic.STLExporter
                  SofaViscoElastic.VisualModelOBJExporter

              )doc";

    py::module::import("Sofa.Core");
    py::module::import("SofaTypes");

     moduleAddTetrahedronViscoelasticityFEMForceField(m);
     moduleAddTetrahedronViscoHyperelasticityFEMForceField(m);
}

}  // namespace sofapython3

