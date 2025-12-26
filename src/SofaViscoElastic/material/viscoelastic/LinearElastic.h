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
#pragma once

#include <SofaViscoElastic/config.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/TopologyData.inl>

#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>
#include <string>

namespace sofa::SofaViscoElastic::material
{

/* a Class that describe a generic Viscoelastic material : example of  Standard Linear Solid Maxwell Representaion.
The material is described based on continuum mechanics and the description is independent
to any discretization method like the finite element method.

For the explanation of the algorithm, the user can see the documentation on this two links:

FEniCS : https://comet-fenics.readthedocs.io/en/latest/demo/viscoelasticity/linear_viscoelasticity.html

COMSOL : https://doc.comsol.com/5.5/doc/com.comsol.help.sme/sme_ug_theory.06.26.html 

*/  






template<class DataTypes>
class LinearElastic : public BaseViscoelasticMaterial<DataTypes>{
public:

  static constexpr std::string_view Name = "LinearElastic";
    typedef typename DataTypes::Coord::value_type Real;
    typedef type::Mat<3,3,Real> Matrix3;
    typedef type::MatSym<3,Real> MatrixSym;


    void deriveCauchyGreenStressTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param, MatrixSym &CauchyStressTensor, SReal& dt) override
    {
        Real mu = param.parameterArray[0];
        Real lambda = param.parameterArray[1];

        MatrixSym ID;
        ID.identity();

       
        // the equation in 3D uses  the lame' parameters, in the equation is considered the transformation from (E, nu) to (mu, lambda) which are the Lame' parameters
        // The equation is always  sigma =  2 * mu * E + lambda * trE * ID

        CauchyStressTensor = 2 * mu * sinfo->Edev + lambda* sinfo->trE * ID;


    }

    void applyElasticityTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,const MatrixSym& inputTensor, MatrixSym &outputTensor, SReal& t) override

    {
        Real mu = param.parameterArray[0];
        Real lambda = param.parameterArray[1];

        MatrixSym ID;
        ID.identity();


        const Real trH = sofa::type::trace(inputTensor);

        outputTensor = ID * (trH * lambda / 2.0) + inputTensor * mu;

    }

};


} // namespace 
