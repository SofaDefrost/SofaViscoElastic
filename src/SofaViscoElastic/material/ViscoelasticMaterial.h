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
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once
#include <SofaViscoElastic/config.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/TopologyData.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>
#include <sofa/type/MatSym.h>
#include <string>

//#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
//#include <opencv2/imgproc.hpp>
//#include <opencv2/imgcodecs.hpp>
//#include <opencv2/highgui.hpp>

namespace sofa::SofaViscoElastic::material
{

template<typename Real>
class StrainInformation;

template<typename DataTypes>
struct MaterialParameters;

/** a Class that describe a generic hyperelastic material .
The material is described based on continuum mechanics and the description is independent
to any discretization method like the finite element method. 
A material is generically described by a strain energy function and its first and second derivatives.
*/
template<class DataTypes>
class ViscoelasticMaterial
{
public:

  typedef typename DataTypes::Coord Coord;
  typedef typename Coord::value_type Real;
  typedef type::MatSym<3,Real> MatrixSym;
  typedef type::Mat<3,3,Real> Matrix3;
 


   virtual ~ViscoelasticMaterial(){}

  



  /** computes the second Piola Kirchhoff stress tensor of the current configuration */
    virtual void deriveSPKTensor(StrainInformation<DataTypes> *, const  MaterialParameters<DataTypes> &,MatrixSym &, Real& )  {}

  /** computes the Elasticity Tensor of the current configuration */

    virtual void applyElasticityTensor(StrainInformation<DataTypes> *, const  MaterialParameters<DataTypes> &,const MatrixSym& , MatrixSym &, Real&)  {}




};

/** structure that store the parameters required to that are necessary to compute the strain energy
The material parameters might be constant in space (homogeneous material) or not */
template<typename DataTypes>
struct MaterialParameters {
  typedef typename DataTypes::Coord Coord;
  typedef typename Coord::value_type Real;

  /** an array of Real values that correspond to the material parameters : the size depends on the material,
  e.g. 2 Lame coefficients for St-Venant Kirchhoff materials */
  std::vector<Real> parameterArray;
  /** the direction of anisotropy in the rest configuration  : the size of the array is 0 if the material is
  isotropic, 1 if it is transversely isotropic and 2 for orthotropic materials (assumed to be orthogonal to each other)*/
  std::vector<Coord> anisotropyDirection;

};

template<typename DataTypes>
class StrainInformation
{
public:

  typedef typename DataTypes::Coord Coord;
  typedef typename Coord::value_type Real;
  typedef type::MatSym<3,Real> MatrixSym;

  /// Trace of C = I1
  Real trC;
  Real J;
  Real lambda;

  /// boolean indicating whether the invariants have been computed
  bool hasBeenInitialized;
  /// right Cauchy-Green deformation tensor C (gradPhi^T gradPhi)
  MatrixSym C;
  MatrixSym SPKdot;//stress rate
  MatrixSym SPKdotprev;
  MatrixSym acc_SPK;
  MatrixSym prevacc_SPK;

  Real logJ;
  MatrixSym E; //strain tensor


  MatrixSym Edot;// strain rate
  MatrixSym Edotprev;
  MatrixSym acc_E;
  MatrixSym prevacc_E;

  StrainInformation() : trC(0), J(0), lambda(0), hasBeenInitialized(false), C(), SPKdot(), SPKdotprev(), logJ(0),E(), Edot(), Edotprev(),acc_E(), prevacc_E() {}
  virtual ~StrainInformation() {}
};

} // namespace 
