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

/* a Class that describe a generic Viscoelastic material : example of Standard Linear solid Maxwell representation of Second Order.
The material is described based on continuum mechanics and the description is independent
to any discretization method like the finite element method.

For the explanation of the algorithm, the user can see the documentation on this two links:

FEniCS : https://comet-fenics.readthedocs.io/en/latest/demo/viscoelasticity/linear_viscoelasticity.html

COMSOL: https://doc.comsol.com/5.5/doc/com.comsol.help.sme/sme_ug_theory.06.26.html 

*/  

  



template<class DataTypes>
class SLSMaxwellSecondOrder : public BaseViscoelasticMaterial<DataTypes>{
public:

  static constexpr std::string_view Name = "SLSMaxwellSecondOrder";
    typedef typename DataTypes::Coord::value_type Real;
    typedef type::Mat<3,3,Real> Matrix3;
    typedef type::MatSym<3,Real> MatrixSym;

    virtual void deriveCauchyGreenStressTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,MatrixSym &CauchyStressTensor, SReal& dt) override
    {

        Real G0 = param.parameterArray[0];
        Real G1 = param.parameterArray[1];
        Real tau1 = param.parameterArray[2];
        Real G2 = param.parameterArray[3];
        Real tau2 = param.parameterArray[4];
        Real lambda = param.parameterArray[5];        

        MatrixSym ID;
        ID.identity();

        sinfo->Evisc1 = (1 / (1 + ( tau1 / dt ))) * (( tau1 / dt ) * sinfo->Evisc_prev1 + sinfo->E);
        sinfo->Evisc2 = (1 / (1 + ( tau2 / dt ))) * (( tau2 / dt ) * sinfo->Evisc_prev2 + sinfo->E);        

        /// The equation of the Cauchy Stress tensor for the Maxwell Model.
        CauchyStressTensor = 2 * G0 * sinfo->E + 2 * G1 *(sinfo->E-sinfo->Evisc1) + 2 * G2 * (sinfo->E - sinfo->Evisc2) + lambda * sinfo->trE * ID;

        /// Store the viscous strain every time step.
        sinfo->Evisc_prev1 = sinfo->Evisc1;
        sinfo->Evisc_prev2 = sinfo->Evisc2;
        
    }

    virtual void applyElasticityTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,const MatrixSym& inputTensor, MatrixSym &outputTensor, SReal& dt) override
      {

        Real G0 = param.parameterArray[0];
        Real G1 = param.parameterArray[1];
        Real tau1 = param.parameterArray[2];
        Real G2 = param.parameterArray[3];
        Real tau2 = param.parameterArray[4];
        Real lambda = param.parameterArray[5];

        MatrixSym ID;
        ID.identity();

        const Real trH = sofa::type::trace(inputTensor);

 

        // The 4th order tensor of elasticity is always approximated to the one in case of pure Linear elasticity (Long-term elasticity tensor) 

        outputTensor = ID * (trH * lambda / 2.0) + inputTensor * G0;


    }

};



} // namespace sofa::component::solidmechanics::fem::hyperelastic::material
