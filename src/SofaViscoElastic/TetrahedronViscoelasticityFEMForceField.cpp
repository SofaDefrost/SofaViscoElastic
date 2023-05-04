/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
* Authors: Pasquale Ferrentino The SOFA Team(see Authors.txt)                 *
*                                                                             *
* Contact information: contact@sofa-framework.org & pasquale.ferrentino@vub.be*
******************************************************************************/

#define SOFA_COMPONENT_FORCEFIELD_TETRAHEDRONVISCOELASTICITYFEMFORCEFIELD_CPP

#include <SofaViscoElastic/TetrahedronViscoelasticityFEMForceField.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa::SofaViscoElastic
{

using namespace sofa::defaulttype;

//////////****************To register in the factory******************

// Register in the Factory
int TetrahedronViscoelasticityFEMForceFieldClass = core::RegisterObject("Generic Tetrahedral finite elements")
.add< TetrahedronViscoelasticityFEMForceField<Vec3Types> >();

template class SOFAVISCOELASTIC_API TetrahedronViscoelasticityFEMForceField<Vec3Types>;


} // namespace sofa::component::solidmechanics::fem::hyperelastic
