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
#pragma once

#include <SofaViscoElastic/config.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/TopologyData.inl>

#include <SofaViscoElastic/material/ViscoelasticMaterial.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>
#include <string>

namespace sofa::SofaViscoElastic::material
{

/** a Class that describe a generic Viscoelastic material : exemple of Maxwell First Order
The material is described based on continuum mechanics and the description is independent
to any discretization method like the finite element method.

**/  

template<class DataTypes>
class SLSKelvinVoigtFirstOrder : public BaseViscoelasticMaterial<DataTypes>{
    typedef typename DataTypes::Coord::value_type Real;
    typedef type::Mat<3,3,Real> Matrix3;
    typedef type::Mat<6,6,Real> Matrix6;
    typedef type::MatSym<3,Real> MatrixSym;


    //Comment(dmarchal): dt should be SReal as it is sofa defined...
    void deriveSPKTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,MatrixSym &SPKTensorGeneral, Real& dt) override
    {

        Real E0=param.parameterArray[0];
        Real E1=param.parameterArray[1];
        Real tau=param.parameterArray[2];
        MatrixSym inversematrix;
        invertMatrix(inversematrix,sinfo->C);
        MatrixSym ID;
        MatrixSym E = sinfo->E;
        MatrixSym Edot = sinfo->Edot;

        ID.identity();

        // 2) In this code we will apply the Newmark scheme integration to discretize the differential equations of the Linear ViscoElastic Materials.
        // In particular the Newmark Scheme is applied to the Stress rate tensor calculation.
        int k = 0, l = 0;

        // Calculation of Stress rate tensor according to Newmark    SPKdot(t+dt) = SPKdot(t) + 0.5*dt*(a(t+dt)+a(t)) Where a is the Stress acceleration
        for (k = 0; k < 3; ++k)
        {
            for (int l = 0; l < 3; ++l)
            {
                sinfo->SPKdot(k,l) = sinfo->SPKdotprev(k,l) + 0.5*dt*(sinfo->acc_SPK(k,l)+sinfo->prevacc_SPK(k,l));
            }
        }


        // Calculation of Stress acceleration.
        for (k = 0; k < 3; ++k)
        {
            for (l = 0; l < 3; ++l)
            {
                sinfo->acc_SPK(k,l) = (sinfo->SPKdot(k,l)-sinfo->SPKdotprev(k,l))/dt;
            }
        }

        // Differential equation: SPK = E0 * E + (E0+E1)* tau * Edot - tau * SPKdot

        Real beta = 1+(E0/E1);

        SPKTensorGeneral = (1/beta)*(E0*E+E0*tau*Edot-tau*sinfo->SPKdot);

        // Store the value of Stress rate every time step

        for (k = 0; k < 3; ++k)
        {
            for (l = 0; l < 3; ++l)
            {
                sinfo->SPKdotprev(k,l) = sinfo->SPKdot(k,l);
            }
        }




        // Store the  value of the Stress Acceleration every Time step.

        for (k = 0; k < 3; ++k)
        {
            for (l = 0; l < 3; ++l)
            {
                sinfo->prevacc_SPK(k,l) = sinfo->acc_SPK(k,l);
            }
        }

        // Do the Multiplication C^-1 * SPK

        SPKTensorGeneral.Mat2Sym(inversematrix.SymSymMultiply(SPKTensorGeneral), SPKTensorGeneral);

    }



    virtual void applyElasticityTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,const MatrixSym& inputTensor, MatrixSym &outputTensor, Real& dt)
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
        for(int k = 0; k<3; ++k){
            for(int l=0; l<3; ++l){
                if(sinfo->Edot(k,l) >= (1/tau)){
                    Real alpha = E0+E1;

                    outputTensor = Thirdmatrix*(0.5*alpha-((E0+E1)/(3*(1-2*nu)))*log(sinfo->J)*0.5)+ inversematrix*((E0+E1)/(3*(1-2*nu)))*trHC*0.5;
                }
                else{
                    Real alpha = E0+E1/(1-exp(-dt/tau));

                    outputTensor = Thirdmatrix*(0.5*alpha-((E0+E1)/(3*(1-2*nu)))*log(sinfo->J)*0.5)+ inversematrix*((E0+E1)/(3*(1-2*nu)))*trHC*0.5;
                }
            }
        }
    }
};


} // namespace sofa::component::solidmechanics::fem::hyperelastic::material
