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

/* a Class that describe a generic Viscoelastic material : example of  Standard Linear Solid Maxwell Representaion with a Neo-Hookean 
hyperelastic spring working in parallel.
The material is described based on continuum mechanics and the description is independent
to any discretization method like the finite element method.

For the explanation of the algorithm, the user can see the documentation on this two links:

FEniCS : https://comet-fenics.readthedocs.io/en/latest/demo/viscoelasticity/linear_viscoelasticity.html

COMSOL : https://doc.comsol.com/5.5/doc/com.comsol.help.sme/sme_ug_theory.06.26.html 

*/  






template<class DataTypes>
class SLSNeoHookeanFirstOrder : public BaseViscoHyperelasticMaterial<DataTypes>{

public:
    static constexpr std::string_view Name = "SLSNeoHookeanFirstOrder";

    typedef typename DataTypes::Coord::value_type Real;
    typedef type::Mat<3,3,Real> Matrix3;
    typedef type::MatSym<3,Real> MatrixSym;


    void deriveSPKTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,MatrixSym &SPKTensorGeneral, SReal& dt) override
    {
        Real mu = param.parameterArray[0];
        Real G1 = param.parameterArray[1];
        Real tau = param.parameterArray[2];
        Real lambda = param.parameterArray[3];


        MatrixSym ID;
        ID.identity();
        
        // inverse of the right Cauchy-Green deformation tensor
        MatrixSym inverse_C;
        sofa::type::invertMatrix(inverse_C, sinfo->C);

        //det(F) = J
        const Real J = sinfo->J;

        // viscous strain
        sinfo->Evisc1 = (1 / (1 + ( tau / dt ))) * (( tau / dt ) * sinfo->Evisc_prev1 + sinfo->E);

        //second Piola-Kirchhoff stress tensor
        SPKTensorGeneral = mu * ID + (lambda * std::log(J) - mu) * inverse_C + 2 * G1 * (sinfo->E - sinfo->Evisc1);

        /// Store the viscous strain every time step.
        sinfo->Evisc_prev1 = sinfo->Evisc1;


    }

    void applyElasticityTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,const MatrixSym& inputTensor, MatrixSym &outputTensor, SReal& t) override

    {
        Real mu = param.parameterArray[0];
        Real G1 = param.parameterArray[1];
        Real tau = param.parameterArray[2];
        Real lambda = param.parameterArray[3];


 // inverse of the right Cauchy-Green deformation tensor
        MatrixSym inverse_C;
        sofa::type::invertMatrix(inverse_C, sinfo->C);

        // trace(C^-1 * H)
        Real trHC = inputTensor[0] * inverse_C[0] + inputTensor[2] * inverse_C[2] + inputTensor[5] * inverse_C[5]
                + 2 * inputTensor[1] * inverse_C[1] + 2 * inputTensor[3] * inverse_C[3] + 2 *
                inputTensor[4] * inverse_C[4];

        // C^-1 * H * C^-1
        MatrixSym Firstmatrix;
        MatrixSym::Mat2Sym(inverse_C * (inputTensor * inverse_C), Firstmatrix);

        outputTensor = Firstmatrix * (mu - lambda * log(sinfo->J)) + inverse_C * (lambda * trHC / 2);



    }

};


} // namespace 
