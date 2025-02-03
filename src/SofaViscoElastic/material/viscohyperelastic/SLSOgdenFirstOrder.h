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
#include <sofa/type/MatSym.h>

#include <string>

#include <Eigen/QR>
#include <Eigen/Eigenvalues>

namespace sofa::SofaViscoElastic::material
{

/* a Class that describe a generic Viscoelastic material : example of  Standard Linear Solid Maxwell Representaion with a Stable Neo-Hookean 
hyperelastic spring working in parallel.
The material is described based on continuum mechanics and the description is independent
to any discretization method like the finite element method.

For the explanation of the algorithm, the user can see the documentation on this two links:

FEniCS : https://comet-fenics.readthedocs.io/en/latest/demo/viscoelasticity/linear_viscoelasticity.html

COMSOL : https://doc.comsol.com/5.5/doc/com.comsol.help.sme/sme_ug_theory.06.26.html 

*/  






template<class DataTypes>
class SLSOgdenFirstOrder : public BaseViscoHyperelasticMaterial<DataTypes>{

public:
    static constexpr std::string_view Name = "SLSOgdenFirstOrder";

    typedef typename DataTypes::Coord::value_type Real;
    typedef type::Mat<3,3,Real> Matrix3;
    typedef type::MatSym<3,Real> MatrixSym;
    typedef type::Vec<3,Real> Vect;
    typedef typename Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real,3,3> >::MatrixType EigenMatrix;
    typedef typename Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real,3,3> >::RealVectorType CoordEigen;

    void deriveSPKTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,MatrixSym &SPKTensorGeneral, MatrixSym &CauchyStressTensor, SReal& dt) override
    {

        Real k0 = param.parameterArray[0];
        Real mu1 = param.parameterArray[1];
        Real alpha1 = param.parameterArray[2];
        Real E1 = param.parameterArray[3];
        Real tau = param.parameterArray[4];

        MatrixSym C=sinfo->C;
        EigenMatrix CEigen;
        CEigen(0,0)=C[0]; CEigen(0,1)=C[1]; CEigen(1,0)=C[1]; CEigen(1,1)=C[2]; CEigen(1,2)=C[4]; CEigen(2,1)=C[4];
        CEigen(2,0)=C[3]; CEigen(0,2)=C[3]; CEigen(2,2)=C[5];

        Eigen::SelfAdjointEigenSolver<EigenMatrix> Vect(CEigen,true);
        EigenMatrix Evect=Vect.eigenvectors();
        CoordEigen Evalue=Vect.eigenvalues();

        Real trCalpha=pow(Evalue[0],alpha1/(Real)2)+pow(Evalue[1],alpha1/(Real)2)+pow(Evalue[2],alpha1/(Real)2);
        Matrix3 Pinverse;
        Pinverse(0,0)=Evect(0,0); Pinverse(1,1)=Evect(1,1); Pinverse(2,2)=Evect(2,2); Pinverse(0,1)=Evect(1,0); Pinverse(1,0)=Evect(0,1); Pinverse(2,0)=Evect(0,2);
        Pinverse(0,2)=Evect(2,0); Pinverse(2,1)=Evect(1,2); Pinverse(1,2)=Evect(2,1);
        MatrixSym Dalpha_1=MatrixSym(pow(Evalue[0],alpha1/(Real)2.0-(Real)1.0),0,pow(Evalue[1],alpha1/(Real)2.0-(Real)1.0),0,0,pow(Evalue[2],alpha1/(Real)2.0-(Real)1.0));
        MatrixSym Calpha_1; Matrix3 Ca;
        Ca=Pinverse.transposed()*Dalpha_1.SymMatMultiply(Pinverse);
        Calpha_1.Mat2Sym(Ca,Calpha_1);

        MatrixSym product;
        product.Mat2Sym(C.SymSymMultiply(Calpha_1), product);

        MatrixSym inversematrix;
        invertMatrix(inversematrix,sinfo->C);
        MatrixSym ID;
        ID.identity();
        
        Real J = sinfo->J;

        /// Calculation Viscous strain 
        sinfo->Evisc1 = (1/(1+(dt/tau)))*(sinfo->Evisc_prev1+ (dt/tau)*sinfo->E);

        // C-1:Evisc_previous_step
        MatrixSym P; 
        P.Mat2Sym(inversematrix.SymSymMultiply(sinfo->Evisc_prev1),P);


        SPKTensorGeneral=(-(Real)1.0/(Real)3.0*trCalpha*inversematrix+Calpha_1)*(mu1/alpha1*pow(sinfo->J,-alpha1/(Real)3.0))+inversematrix*(k0*log(sinfo->J)
        -0.5*(E1/(1+(dt/tau))))+0.5*(E1/(1+(dt/tau)))*ID -(E1/(1+(dt/tau)))*P;



        /// store the viscous strain every time step
        sinfo->Evisc_prev1 = sinfo->Evisc1;



    }

    void applyElasticityTensor(StrainInformation<DataTypes> *sinfo, const MaterialParameters<DataTypes> &param,const MatrixSym& inputTensor, MatrixSym &outputTensor, SReal& t) override

    {
        Real k0 = param.parameterArray[0];
        Real mu1 = param.parameterArray[1];
        Real alpha1 = param.parameterArray[2];
        Real E1 = param.parameterArray[3];
        Real tau = param.parameterArray[4];

        MatrixSym C=sinfo->C;
        EigenMatrix CEigen;
        CEigen(0,0)=C[0]; CEigen(0,1)=C[1]; CEigen(1,0)=C[1]; CEigen(1,1)=C[2]; CEigen(1,2)=C[4]; CEigen(2,1)=C[4];
        CEigen(2,0)=C[3]; CEigen(0,2)=C[3]; CEigen(2,2)=C[5];
        Eigen::SelfAdjointEigenSolver<EigenMatrix> Vect(CEigen,true);
        EigenMatrix Evect=Vect.eigenvectors();
        CoordEigen Evalue=Vect.eigenvalues();

        Real trCalpha=pow(Evalue[0],alpha1/(Real)2)+pow(Evalue[1],alpha1/(Real)2)+pow(Evalue[2],alpha1/(Real)2);
        Matrix3 Pinverse;
        Pinverse(0,0)=Evect(0,0); Pinverse(1,1)=Evect(1,1); Pinverse(2,2)=Evect(2,2); Pinverse(0,1)=Evect(1,0); Pinverse(1,0)=Evect(0,1); Pinverse(2,0)=Evect(0,2);
        Pinverse(0,2)=Evect(2,0); Pinverse(2,1)=Evect(1,2); Pinverse(1,2)=Evect(2,1);
        MatrixSym Dalpha_1=MatrixSym(pow(Evalue[0],alpha1/(Real)2.0-(Real)1.0),0,pow(Evalue[1],alpha1/(Real)2.0-(Real)1.0),0,0,pow(Evalue[2],alpha1/(Real)2.0-(Real)1.0));
        MatrixSym Calpha_1; Matrix3 Ca;
        Ca=Pinverse.transposed()*Dalpha_1.SymMatMultiply(Pinverse);
        Calpha_1.Mat2Sym(Ca,Calpha_1);
        MatrixSym Dalpha_2=MatrixSym(pow(Evalue[0],alpha1/(Real)4.0-(Real)1.0),0,pow(Evalue[1],alpha1/(Real)4.0-(Real)1.0),0,0,pow(Evalue[2],alpha1/(Real)4.0-(Real)1.0));
        MatrixSym Calpha_2;
        Calpha_2.Mat2Sym(Pinverse.transposed()*Dalpha_2.SymMatMultiply(Pinverse),Calpha_2);
        MatrixSym inversematrix;
        invertMatrix(inversematrix,sinfo->C);
        MatrixSym ID;
        ID.identity();

        Real J = sinfo->J;

        Real _trHCalpha_1=inputTensor[0]*Calpha_1[0]+inputTensor[2]*Calpha_1[2]+inputTensor[5]*Calpha_1[5]
                +2*inputTensor[1]*Calpha_1[1]+2*inputTensor[3]*Calpha_1[3]+2*inputTensor[4]*Calpha_1[4];
        Real _trHC=inputTensor[0]*inversematrix[0]+inputTensor[2]*inversematrix[2]+inputTensor[5]*inversematrix[5]
                +2*inputTensor[1]*inversematrix[1]+2*inputTensor[3]*inversematrix[3]+2*inputTensor[4]*inversematrix[4];

        MatrixSym Firstmatrix;
        Firstmatrix.Mat2Sym(inversematrix.SymMatMultiply(inputTensor.SymSymMultiply(inversematrix)),Firstmatrix);
        MatrixSym Secondmatrix;
        Secondmatrix.Mat2Sym(Calpha_2.SymMatMultiply(inputTensor.SymSymMultiply(Calpha_2)),Secondmatrix);

        outputTensor =  (_trHC*(-alpha1/(Real)6.0)*(-(Real)1.0/(Real)3.0*inversematrix*trCalpha+Calpha_1)+(Real)1.0/(Real)3.0*Firstmatrix*trCalpha-(Real)1.0/(Real)3.0*inversematrix*_trHCalpha_1*alpha1/(Real)2.0
                +(alpha1/(Real)2.0-(Real)1)*Secondmatrix) * (mu1/alpha1*pow(sinfo->J,-alpha1/(Real)3.0))
                +k0/(Real)2.0*_trHC*inversematrix-(Real)(k0*log(sinfo->J))*Firstmatrix + Firstmatrix*E1*exp(-t/tau);


    }

};


} // namespace 
