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

/* a Class that describe a generic Viscohyperelastic material : example of  Standard Linear Solid Maxwell Representaion with Mooney-Rivlin hyperelastic 
spring working in parallel.
The material is described based on continuum mechanics and the description is independent
to any discretization method like the finite element method.

For the explanation of the algorithm, the user can see the documentation on this two links:

FEniCS : https://comet-fenics.readthedocs.io/en/latest/demo/viscoelasticity/linear_viscoelasticity.html

COMSOL : https://doc.comsol.com/5.5/doc/com.comsol.help.sme/sme_ug_theory.06.26.html 

*/  

template <class DataTypes>
class SLSMooneyRivlinSecondOrder : public BaseViscoHyperelasticMaterial<DataTypes>
{
public:
  static constexpr std::string_view Name = "SLSMooneyRivlinSecondOrder";

    typedef typename DataTypes::Coord::value_type Real;
    typedef type::Mat<3, 3, Real> Matrix3;
    typedef type::Mat<6, 6, Real> Matrix6;
    typedef type::MatSym<3, Real> MatrixSym;


    void deriveSPKTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,MatrixSym &SPKTensorGeneral, SReal& dt) override
    {
        MatrixSym inversematrix;
        MatrixSym C = sinfo->C;
        MatrixSym E = sinfo->E;
        invertMatrix(inversematrix, C);
        Real I1 = sinfo->trC;
        Real I1square = (Real)(C[0]*C[0]+C[2]*C[2]+C[5]*C[5]+2*(C[1]*C[1]+C[3]*C[3]+C[4]*C[4]));
        Real I2 = (Real)((pow(I1, (Real)2)-I1square)/2);
        Real c1 = param.parameterArray[0];
        Real c2 = param.parameterArray[1];
        Real G1 = param.parameterArray[2];
        Real tau1 = param.parameterArray[3];
        Real G2 = param.parameterArray[4];
        Real tau2 = param.parameterArray[5];        
        Real k0 = param.parameterArray[6];

        MatrixSym ID;
        ID.identity();

        /// Calculation Viscous strain 
        sinfo->Evisc1 = (1 / (1 + ( tau1 / dt ))) * (( tau1 / dt ) * sinfo->Evisc_prev1 + sinfo->E);
        sinfo->Evisc2 = (1 / (1 + ( tau2 / dt ))) * (( tau2 / dt ) * sinfo->Evisc_prev2 + sinfo->E);


        SPKTensorGeneral = (-1*inversematrix*I1/3+ID)*(2*c1*pow(sinfo->J, (Real)(-2.0/3.0)))+(
                               -1*inversematrix*2*I2/3+ID*I1-C)*(
                               2*c2*pow(sinfo->J, (Real)(-4.0/3.0)))+inversematrix*(k0*log(sinfo->J)) + 2 * G1 * (sinfo->E - sinfo->Evisc1) + 2 * G2 * (sinfo->E - sinfo->Evisc2);




        /// store the viscous strain every time step
        sinfo->Evisc_prev1 = sinfo->Evisc1;
        sinfo->Evisc_prev2 = sinfo->Evisc2;


    }


    void applyElasticityTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,const MatrixSym& inputTensor, MatrixSym &outputTensor, SReal& dt) override
    {
        MatrixSym inversematrix;
        MatrixSym C = sinfo->C;
        invertMatrix(inversematrix, C);
        Real I1 = sinfo->trC;
        Real I1square = (Real)(C[0]*C[0]+C[2]*C[2]+C[5]*C[5]+2*(C[1]*C[1]+C[3]*C[3]+C[4]*C[4]));
        Real I2 = (Real)((pow(I1, (Real)2)-I1square)/2);
        Real c1 = param.parameterArray[0];
        Real c2 = param.parameterArray[1];
        Real G1 = param.parameterArray[2];
        Real tau1 = param.parameterArray[3];
        Real G2 = param.parameterArray[4];
        Real tau2 = param.parameterArray[5];        
        Real k0 = param.parameterArray[6];

        MatrixSym ID;
        ID.identity();
        // C-1:H
        Real _trHC = inputTensor[0]*inversematrix[0]+inputTensor[2]*inversematrix[2]+inputTensor[5]*
                     inversematrix[5]
                     +2*inputTensor[1]*inversematrix[1]+2*inputTensor[3]*inversematrix[3]+2*
                     inputTensor[4]*inversematrix[4];
        MatrixSym Firstmatrix;
        //C-1HC-1 convert to sym matrix
        Firstmatrix.Mat2Sym(inversematrix.SymMatMultiply(inputTensor.SymSymMultiply(inversematrix)),
                            Firstmatrix);
        //C:H
        Real trHC = inputTensor[0]*C[0]+inputTensor[2]*C[2]+inputTensor[5]*C[5]
                    +2*inputTensor[1]*C[1]+2*inputTensor[3]*C[3]+2*inputTensor[4]*C[4];

        //trH
        Real trH = inputTensor[0]+inputTensor[2]+inputTensor[5];
   
        outputTensor = ((ID-inversematrix*I1/(Real)3.0)*(-_trHC)/(Real)3.0+Firstmatrix*I1/(Real)3.0-
                        inversematrix*trH/(Real)3.0)*(Real)2.0*c1*pow(sinfo->J, (Real)(-2.0/3.0))
                       +((inversematrix*(Real)(-2.0)*I2/(Real)3.0+ID*I1-C)*(Real)(-2.0)*_trHC/(Real)
                         3.0+Firstmatrix*(Real)2.0*I2/(Real)3.0-inversematrix*(Real)2.0*(
                             I1*trH-trHC)/(Real)3.0+ID*trH-inputTensor)*(Real)2.0*c2*pow(
                           sinfo->J, (Real)(-4.0/3.0))
                       +inversematrix*_trHC*k0/(Real)2.0-Firstmatrix*k0*log(sinfo->J);          

    }

};
} // namespace 


