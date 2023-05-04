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
#include <SofaViscoElastic/material/ViscoelasticMaterial.h>

#include <sofa/type/Mat.h>
#include <sofa/type/MatSym.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/TopologyData.h>

namespace sofa::SofaViscoElastic
{

namespace
{
    using namespace sofa::type;
    using namespace sofa::defaulttype;
    using namespace sofa::core::topology;
}

//***************** Tetrahedron FEM code for several elastic models: TotalLagrangianForceField************************//

/** Compute Finite Element forces based on tetrahedral elements.
*/
template<class DataTypes>
class TetrahedronViscoelasticityFEMForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TetrahedronViscoelasticityFEMForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;

    typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv;
    typedef core::objectmodel::Data<VecCoord>    DataVecCoord;

    typedef Mat<3,3,Real> Matrix3;
    typedef MatSym<3,Real> MatrixSym;
    typedef type::Mat<6,6,Real> Matrix6;
    typedef type::Mat<6,3,Real> Matrix63;

    typedef std::pair<MatrixSym,MatrixSym> MatrixPair;
    typedef std::pair<Real,MatrixSym> MatrixCoeffPair;

    typedef type::vector<Real> SetParameterArray;
    typedef type::vector<Coord> SetAnisotropyDirectionArray;

    typedef core::topology::BaseMeshTopology::Index Index;
    typedef core::topology::BaseMeshTopology::Tetra Element;
    typedef core::topology::BaseMeshTopology::SeqTetrahedra VecElement;
    typedef sofa::core::topology::Topology::Tetrahedron Tetrahedron;
    typedef sofa::core::topology::Topology::TetraID TetraID;
    typedef sofa::core::topology::Topology::Tetra Tetra;
    typedef sofa::core::topology::Topology::Edge Edge;
    typedef sofa::core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
    typedef sofa::core::topology::BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;
    typedef sofa::core::topology::BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;

    material::MaterialParameters<DataTypes> globalParameters;

    /// data structure stored for each tetrahedron
    class TetrahedronRestInformation : public material::StrainInformation<DataTypes>
    {
    public:
        /// shape vector at the rest configuration
        Coord m_shapeVector[4];
        /// fiber direction in rest configuration
        Coord m_fiberDirection;
        /// rest volume
        Real m_restVolume{};
        /// current tetrahedron volume
        Real m_volScale{};
        Real m_volume{};
        /// Second Piola Kirchhoff stress tensor 
        MatrixSym m_SPKTensorGeneral;
        /// deformation gradient = gradPhi
        Matrix3 m_deformationGradient;
        Real m_strainEnergy{};
        
        /// Output stream
        inline friend std::ostream& operator<< ( std::ostream& os, const TetrahedronRestInformation& eri ) {  os << eri.m_SPKTensorGeneral; return os;  }
        /// Input stream
        inline friend std::istream& operator>> ( std::istream& in, TetrahedronRestInformation& /*eri*/ ) { return in; }

        TetrahedronRestInformation() = default;
    };

    /// data structure stored for each edge
    class EdgeInformation
    {
    public:
        /// store the stiffness edge matrix
        Matrix3 DfDx;

        /// Output stream
        inline friend std::ostream& operator<< (std::ostream& os, const EdgeInformation& /*eri*/ ) {  return os;  }
        /// Input stream
        inline friend std::istream& operator>> (std::istream& in, EdgeInformation& /*eri*/ ) { return in; }

        EdgeInformation() = default;
    };

 public :
    Data<bool> d_stiffnessMatrixRegularizationWeight; ///< Regularization of the Stiffness Matrix (between true or false)
    Data<std::string> d_materialName; ///< the name of the material
    Data<SetParameterArray> d_parameterSet; ///< The global parameters specifying the material
    Data<SetAnisotropyDirectionArray> d_anisotropySet; ///< The global directions of anisotropy of the material

    TetrahedronData<sofa::type::vector<TetrahedronRestInformation> > m_tetrahedronInfo; ///< Internal tetrahedron data
    EdgeData<sofa::type::vector<EdgeInformation> > m_edgeInfo; ///< Internal edge data
   
    /// Link to be set to the topology container in the component graph.
    SingleLink<TetrahedronViscoelasticityFEMForceField<DataTypes>, sofa::core::topology::BaseMeshTopology, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_topology;

    //Comment(dmarchal): the material shouldn't be using string but sofa::helper::OptionsGroup, see examples of use
    void setMaterialName(std::string materialName);
    void setparameter(const SetParameterArray& param);
    void setdirection(const SetAnisotropyDirectionArray& direction);

    /**
     * Method to initialize @sa TetrahedronRestInformation when a new Tetrahedron is created.
     * Will be set as creation callback in the TetrahedronData @sa m_tetrahedronInfo
     */
    void createTetrahedronRestInformation(Index, TetrahedronRestInformation& t, const Tetrahedron&,
        const sofa::type::vector<Index>&, const sofa::type::vector<SReal>&);

public:
    void init() override;    
    void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v) override;
    void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx) override;
    SReal getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const override;
    void addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal k, unsigned int &offset) override;

    void draw(const core::visual::VisualParams* vparams) override;
    void computeBBox(const core::ExecParams* params, bool onlyVisible) override;

    Mat<3,3, SReal> getPhi( int tetrahedronIndex);

protected:
    TetrahedronViscoelasticityFEMForceField();
    ~TetrahedronViscoelasticityFEMForceField() override;


    /// the array that describes the complete material energy and its derivatives
    std::unique_ptr<material::BaseViscoelasticMaterial<DataTypes> > m_myMaterial;

    void testDerivatives();
    void updateTangentMatrix();
    void instantiateMaterial();

private:
    //Comment(dmarchal): it is better to use SinglLinks to store reference to instance of object.
    core::topology::BaseMeshTopology* m_topology;
    VecCoord m_initialPoints;   /// the intial positions of the points
    bool m_updateMatrix;

};

#if  !defined(SOFA_COMPONENT_FORCEFIELD_TetrahedronViscoelasticityFEMForceField_CPP)
extern template class SOFAVISCOELASTIC_API TetrahedronViscoelasticityFEMForceField<defaulttype::Vec3Types>;
#endif //  !defined(SOFA_COMPONENT_FORCEFIELD_TetrahedronViscoelasticityFEMForceField_CPP)

} // namespace sofa::component::solidmechanics::fem::hyperelastic
