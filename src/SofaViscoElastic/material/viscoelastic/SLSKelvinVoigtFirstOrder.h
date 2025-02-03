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

/* a Class that describe a generic Viscoelastic material : example of  Standard Linear Solid Kelvin Representaion.
The material is described based on continuum mechanics and the description is independent
to any discretization method like the finite element method.

For the explanation of the algorithm, the user can see the documentation on this two links:

FEniCS : https://comet-fenics.readthedocs.io/en/latest/demo/viscoelasticity/linear_viscoelasticity.html

COMSOL: https://doc.comsol.com/5.5/doc/com.comsol.help.sme/sme_ug_theory.06.26.html 

*/  

template<class DataTypes>
class SLSKelvinVoigtFirstOrder : public BaseViscoelasticMaterial<DataTypes>{

public:

  static constexpr std::string_view Name = "SLSKelvinVoigtFirstOrder";
    typedef typename DataTypes::Coord::value_type Real;
    typedef type::Mat<3,3,Real> Matrix3;
    typedef type::MatSym<3,Real> MatrixSym;


    void deriveSPKTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,MatrixSym &SPKTensorGeneral,MatrixSym &CauchyStressTensor, SReal& dt) override
    {

        Real E0=param.parameterArray[0];
        Real E1=param.parameterArray[1];
        Real tau=param.parameterArray[2];
        Real nu=param.parameterArray[3];        

        MatrixSym inversematrix;
        invertMatrix(inversematrix,sinfo->C);
        MatrixSym ID;
        ID.identity();
        Real trE = sinfo->E(0,0) + sinfo->E(1,1) +sinfo->E(2,2);

        /// Calculation Viscous strain
        sinfo->Evisc1 = (1/(1+(((E0*dt)/(E1*tau))+(dt/tau))))*(sinfo->Evisc_prev1+((E0*dt)/(E1*tau))*sinfo->E);


        /// The equation of the Cauchy Stress tensor for the SLS Kelvin Model.
        CauchyStressTensor = E0*(sinfo->E-sinfo->Evisc1)+((E0)/(3*(1-2*nu)))*trE*ID; 

        /// store the viscous strain every time step
        sinfo->Evisc_prev1 = sinfo->Evisc1;
        

        /// Do the Multiplication for C^-1 to obtain the Second Piola Kirchhoff stress tensor
        SPKTensorGeneral.Mat2Sym(inversematrix.SymSymMultiply(CauchyStressTensor), SPKTensorGeneral);
    }
 

    virtual void applyElasticityTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,const MatrixSym& inputTensor, MatrixSym &outputTensor, SReal& t) override
    {
        Real E0=param.parameterArray[0];
        Real E1=param.parameterArray[1];
        Real tau=param.parameterArray[2];
        Real nu=param.parameterArray[3];
        MatrixSym inversematrix;
        invertMatrix(inversematrix,sinfo->C);
        MatrixSym ID;
        ID.identity();


        Real trHC=inputTensor[0]*inversematrix[0]+inputTensor[2]*inversematrix[2]+inputTensor[5]*inversematrix[5]
                    +2*inputTensor[1]*inversematrix[1]+2*inputTensor[3]*inversematrix[3]+2*inputTensor[4]*inversematrix[4];




        MatrixSym Thirdmatrix;
        Thirdmatrix.Mat2Sym(inversematrix.SymMatMultiply(inputTensor.SymSymMultiply(inversematrix)),Thirdmatrix);

        Real alpha = E0+E1;

        outputTensor = Thirdmatrix*(0.5*alpha-((E0+E1)/(3*(1-2*nu)))*log(sinfo->J)*0.5)+ inversematrix*((E0+E1)/(3*(1-2*nu)))*trHC*0.5;

    }
};


} // namespace sofa::component::solidmechanics::fem::hyperelastic::material
