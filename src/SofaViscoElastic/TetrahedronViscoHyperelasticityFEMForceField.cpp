/******************************************************************************
*  THE SOFA VISCOELASTIC PLUGIN.                                              *
*                                               							  * 
* DESCRIPTION:                                                                *
* This plugin is made within the SOFA, Simulation Open-Framework Architecture *
* (c) 2006 INRIA, USTL, UJF, CNRS, MGH.                                       *
* The plugin consist in a Visco-Elastic force field for tetrahedral meshes.	  * 
* Several rheological models are implemented.                                 *
* 	                                                                          *
* CONTRIBUTORS:																  *			
* The plugin is made by the collaboration beween the Robotics and Multibody   * 
* Mechanics Department (R&MM) Vrije Universiteit Brussel (VUB), Bruxelles     *
* Belgium, and the DEFROST Team of the INRIA - Lille, France.    		      *
* 																			  *
* 																		  	  *
* LICENSE:    							  									  *
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


#define SOFA_COMPONENT_FORCEFIELD_TETRAHEDRONVISCOHYPERELASTICITYFEMFORCEFIELD_CPP

#include <SofaViscoElastic/TetrahedronViscoHyperelasticityFEMForceField.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa::SofaViscoElastic
{

using namespace sofa::defaulttype;

//////////****************To register in the factory******************

// Register in the Factory
int TetrahedronViscoHyperelasticityFEMForceFieldClass = core::RegisterObject("Generic Tetrahedral finite elements")
.add< TetrahedronViscoHyperelasticityFEMForceField<Vec3Types> >();

template class SOFAVISCOELASTIC_API TetrahedronViscoHyperelasticityFEMForceField<Vec3Types>;


} // namespace 
