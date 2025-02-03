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
#include <SofaViscoElastic/config.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/TopologyData.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>
#include <sofa/type/MatSym.h>
#include <string>




namespace sofa::SofaViscoElastic::material
{

template<typename DataTypes>
struct MaterialParameters;

template<typename DataTypes>
class StrainInformation;

/** a Class that describe a generic Viscoelastic material .
The material is described based on continuum mechanics and the description is independent
to any discretization method like the finite element method. 
*/
template<class DataTypes>
class BaseViscoelasticMaterial
{
public:

  typedef typename DataTypes::Coord Coord;
  typedef typename Coord::value_type Real;
  typedef type::MatSym<3,Real> MatrixSym;
  typedef type::Mat<3,3,Real> Matrix3;

  virtual ~BaseViscoelasticMaterial(){}

  void myCommonFunction(){
      ///real content
  };

  /** computes the second Piola Kirchhoff stress tensor of the current configuration */
   virtual void deriveSPKTensor(StrainInformation<DataTypes> *, const  MaterialParameters<DataTypes> &,MatrixSym &,MatrixSym &, SReal&)  = 0;

  /** computes the Elasticity Tensor of the current configuration */
   virtual void applyElasticityTensor(StrainInformation<DataTypes> *, const  MaterialParameters<DataTypes> &,const MatrixSym& , MatrixSym &, SReal&)  = 0;
};

/** structure that store the parameters required to compute the constitutive law.
The material parameters might be constant in space (homogeneous material) or not */
template<typename DataTypes>
struct MaterialParameters {
  typedef typename DataTypes::Coord Coord;
  typedef typename Coord::value_type Real;

  /** an array of Real values that correspond to the material parameters : the size depends on the material */
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

  Real logJ;
  MatrixSym E; //strain tensor
  MatrixSym Eprev; // strain tensor previous time step. 

  MatrixSym Evisc1; // viscous strain
  MatrixSym Evisc_prev1; // viscous strain tensor at previous time step

  MatrixSym Evisc2; // viscous strain (second dashpot)
  MatrixSym Evisc_prev2; // viscous strain tensor at previous time step (second dashpot)


  StrainInformation() : trC(0), J(0), lambda(0), hasBeenInitialized(false), C(), logJ(0),E(), Eprev(), Evisc1(), Evisc_prev1(), Evisc2(), Evisc_prev2() {}
  virtual ~StrainInformation() {}
};

} // namespace 
