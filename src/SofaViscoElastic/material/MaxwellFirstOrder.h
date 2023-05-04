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
#include "ViscoelasticMaterial.h"
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>
#include <string>



namespace sofa::SofaViscoElastic::material
{

/** a Class that describe a generic Viscoelastic material : example of Maxwell First Order
The material is described based on continuum mechanics.

**/  



template<class DataTypes>
class MaxwellFirstOrder : public ViscoelasticMaterial<DataTypes>{

  typedef typename DataTypes::Coord::value_type Real;
  typedef type::Mat<3,3,Real> Matrix3;
  typedef type::MatSym<3,Real> MatrixSym;





// 2) We discretize directly the equation with an Euler Scheme.
 
    virtual void deriveSPKTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param, MatrixSym &SPKTensorGeneral, Real& dt){

    Real E1=param.parameterArray[0];
    Real tau=param.parameterArray[1];


        MatrixSym inversematrix;
        invertMatrix(inversematrix,sinfo->C);
        MatrixSym Edot = sinfo-> Edot;
        MatrixSym ID;
        ID.identity();


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


//Differential Equation For Maxwell Model: SPK =  E1*tau*Edot - tau*SPKdot 

        SPKTensorGeneral = E1*tau*Edot -tau*sinfo->SPKdot;

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
  
    virtual void applyElasticityTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,const MatrixSym& inputTensor, MatrixSym &outputTensor, Real& dt)  {
		Real E1=param.parameterArray[0];
		Real tau=param.parameterArray[1];
		Real nu=param.parameterArray[2];
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
              outputTensor = Thirdmatrix*(E1-(E1/(3*(1-2*nu)))*log(sinfo->J))+ inversematrix*(E1/(3*(1-2*nu)))*trHC;
          }
          else{

          outputTensor = Thirdmatrix*(E1*exp(-dt/tau)-(E1/(3*(1-2*nu)))*log(sinfo->J))+ inversematrix*(E1/(3*(1-2*nu)))*trHC;

        }  
      }
    }
    
	}

};



} // namespace sofa::component::solidmechanics::fem::hyperelastic::material
