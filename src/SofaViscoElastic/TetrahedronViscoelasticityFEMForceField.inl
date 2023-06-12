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

#include <SofaViscoElastic/TetrahedronViscoelasticityFEMForceField.h>
#include <SofaViscoElastic/TetrahedronViscoelasticityFEMDrawing.h>

#include <SofaViscoElastic/material/MaxwellFirstOrder.h>
#include <SofaViscoElastic/material/SLSMaxwellFirstOrder.h>
#include <SofaViscoElastic/material/KelvinVoigtFirstOrder.h>
#include <SofaViscoElastic/material/SLSKelvinVoigtFirstOrder.h>
#include <SofaViscoElastic/material/Burgers.h>




#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/TopologyData.inl>

namespace sofa::SofaViscoElastic
{

using namespace sofa::defaulttype;
using namespace core::topology;
using namespace sofa::SofaViscoElastic::material;

template <class DataTypes> TetrahedronViscoelasticityFEMForceField<DataTypes>::TetrahedronViscoelasticityFEMForceField()
    : m_topology(nullptr)
    , m_initialPoints(0)
    , m_updateMatrix(true)
    , d_stiffnessMatrixRegularizationWeight(initData(&d_stiffnessMatrixRegularizationWeight, (bool)false,"matrixRegularization","Regularization of the Stiffness Matrix (between true or false)"))
    , d_materialName(initData(&d_materialName,std::string("ArrudaBoyce"),"materialName","the name of the material to be used"))
    , d_parameterSet(initData(&d_parameterSet,"ParameterSet","The global parameters specifying the material"))
    , d_anisotropySet(initData(&d_anisotropySet,"AnisotropyDirections","The global directions of anisotropy of the material"))
    , m_tetrahedronInfo(initData(&m_tetrahedronInfo, "tetrahedronInfo", "Internal tetrahedron data"))
    , m_edgeInfo(initData(&m_edgeInfo, "edgeInfo", "Internal edge data"))
    , d_stressSPK(initData(&d_stressSPK, "stressSPK","The stress of the Second Piola Kirchhoff Stress Tensor per Element"))   
    , d_Cauchystress(initData(&d_Cauchystress, "CauchyStress","The stress of the Cauchy Stress Tensor per Element"))
    , d_stressVonMisesElement(initData(&d_stressVonMisesElement, "stressVonMisesElement","The stress of the Von Mises Stress per Element"))
    , l_topology(initLink("topology", "link to the topology container"))
{
}

template <class DataTypes> TetrahedronViscoelasticityFEMForceField<DataTypes>::~TetrahedronViscoelasticityFEMForceField()
    = default;

template <class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::instantiateMaterial()
{
    const std::string& material = d_materialName.getValue();

    if (material == "MaxwellFirstOrder")
    {        

        m_myMaterial = std::make_unique<MaxwellFirstOrder<DataTypes>>();
    }
    else if (material == "SLSMaxwellFirstOrder")
    {        

        m_myMaterial = std::make_unique<SLSMaxwellFirstOrder<DataTypes>>();
    }
    else if (material == "KelvinVoigtFirstOrder")
    {        

        m_myMaterial = std::make_unique<KelvinVoigtFirstOrder<DataTypes>>();
    } 
    else if (material == "SLSKelvinVoigtFirstOrder")
    {        

        m_myMaterial = std::make_unique<SLSKelvinVoigtFirstOrder<DataTypes>>();
    } 
    else if (material == "Burgers")
    {        

        m_myMaterial = std::make_unique<Burgers<DataTypes>>();
    }
                 
    else
    {
        msg_error() << "material name " << material <<
            " is not valid (should be MaxwellFirstOrder, SLSMaxwellFirstOrder, KelvinVoigtFirstOrder, SLSKelvinVoigtFirstOrder, Burgers)";
    }

    if (m_myMaterial)
    {
        msg_info() << "The model is " << material;
    }
}

template <class DataTypes> void TetrahedronViscoelasticityFEMForceField<DataTypes>::init()
{
    msg_info() << "initializing TetrahedronViscoelasticityFEMForceField";

    this->Inherited::init();

    /// parse the parameter set 
    const SetParameterArray& paramSet = d_parameterSet.getValue();
    if (!paramSet.empty())
    {
        globalParameters.parameterArray.resize(paramSet.size());
        std::copy(paramSet.begin(), paramSet.end(), globalParameters.parameterArray.begin());
    }

    /** parse the anisotropy Direction set */
    const SetAnisotropyDirectionArray& anisotropySet = d_anisotropySet.getValue();
    if (!anisotropySet.empty())
    {
        globalParameters.anisotropyDirection.resize(anisotropySet.size());
        std::copy(anisotropySet.begin(), anisotropySet.end(),
                  globalParameters.anisotropyDirection.begin());
    }

    /// We check if we need to search from the context as there is no topology explicitely given
    if (l_topology.empty())
    {
        msg_info() << "link to Topology container should be set to ensure right behavior. First Topology found in current context will be used.";
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    m_topology = l_topology.get();

    msg_info() << "Topology path used: '" << l_topology.getLinkedPath() << "'";

    /// If there is still no topology we consider it as invalid case.
    if (m_topology == nullptr)
    {
        msg_error() << "No topology component found at path: " << l_topology.getLinkedPath() << ", nor in current context: " << this->getContext()->name;
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    if (!m_topology->getNbTetrahedra())
    {
        msg_error() << "object must have a Tetrahedral Set Topology.";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    /** parse the input material name */
    instantiateMaterial();

    auto tetrahedronInf = sofa::helper::getWriteAccessor(m_tetrahedronInfo);

    /// prepare to store info in the triangle array
    tetrahedronInf.resize(m_topology->getNbTetrahedra());

    sofa::helper::getWriteAccessor(m_edgeInfo).resize(m_topology->getNbEdges());
    m_edgeInfo.createTopologyHandler(m_topology);

    // get restPosition
    if (m_initialPoints.empty())
    {
        m_initialPoints = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
    }

    /// initialize the data structure associated with each tetrahedron
    for (Topology::TetrahedronID i = 0; i < m_topology->getNbTetrahedra(); ++i)
    {
        createTetrahedronRestInformation(i, tetrahedronInf[i],
                                         m_topology->getTetrahedron(i),
                                         (const type::vector<Index>)0,
                                         (const type::vector<SReal>)0);
    }

    /// set the call back function upon creation of a tetrahedron
    m_tetrahedronInfo.createTopologyHandler(m_topology);
    m_tetrahedronInfo.setCreationCallback([this](Index tetrahedronIndex, TetrahedronRestInformation& tetraInfo,
                                                 const core::topology::BaseMeshTopology::Tetrahedron& tetra,
                                                 const sofa::type::vector< Index >& ancestors,
                                                 const sofa::type::vector< SReal >& coefs)
                                          {
                                              createTetrahedronRestInformation(tetrahedronIndex, tetraInfo, tetra, ancestors, coefs);
                                          });

    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}


template <class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::setMaterialName(
    const std::string materialName)
{
    d_materialName.setValue(materialName);
}

template <class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::setparameter(const SetParameterArray& param)
{
    d_parameterSet.setValue(param);
}

template <class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::setdirection(
    const SetAnisotropyDirectionArray& direction)
{
    d_anisotropySet.setValue(direction);
}

template< class DataTypes >
void TetrahedronViscoelasticityFEMForceField<DataTypes>::createTetrahedronRestInformation(Index tetrahedronIndex,
                                                                                          TetrahedronRestInformation& tinfo,
                                                                                          const Tetrahedron&,
                                                                                          const sofa::type::vector<Index>&,
                                                                                          const sofa::type::vector<SReal>&)
{
    const sofa::type::vector< Tetrahedron >& tetrahedronArray = m_topology->getTetrahedra();
    const sofa::type::vector< Edge>& edgeArray = m_topology->getEdges();
    unsigned int j;

    typename DataTypes::Real volume;
    typename DataTypes::Coord point[4];
    const VecCoord& restPosition = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();

    ///describe the indices of the 4 tetrahedron vertices
    const Tetrahedron& t = tetrahedronArray[tetrahedronIndex];
    BaseMeshTopology::EdgesInTetrahedron te = m_topology->getEdgesInTetrahedron(tetrahedronIndex);

    /// store the point position
    for (j = 0; j < 4; ++j)
        point[j] = restPosition[t[j]];

    /// compute 6 times the rest volume
    volume = dot(cross(point[2] - point[0], point[3] - point[0]), point[1] - point[0]);

    /// store the rest volume
    tinfo.m_volScale = (Real)(1.0 / volume);
    tinfo.m_restVolume = fabs(volume / 6);

    /// store shape vectors at the rest configuration
    for (j = 0; j < 4; ++j)
    {
        if (!(j % 2))
            tinfo.m_shapeVector[j] = -cross(point[(j + 2) % 4] - point[(j + 1) % 4],
                                            point[(j + 3) % 4] - point[(j + 1) % 4]) / volume;
        else
            tinfo.m_shapeVector[j] = cross(point[(j + 2) % 4] - point[(j + 1) % 4],
                                           point[(j + 3) % 4] - point[(j + 1) % 4]) / volume;;
    }

    for (j = 0; j < 6; ++j)
    {
        Edge e = m_topology->getLocalEdgesInTetrahedron(j);
        int k = e[0];
        if (edgeArray[te[j]][0] != t[k])
        {
            k = e[1];
        }
    }
}

template <class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::addForce(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& /* d_v */)
{
    if(!this->mstate)
    {
        msg_error() << "object must have a Tetrahedral Set Topology.";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }


    auto f = sofa::helper::getWriteAccessor(d_f);
    const VecCoord& x = d_x.getValue();

    /// get the time step of simulation
    Real dt = this->getContext()->getDt();

    unsigned int j = 0, k = 0, l = 0;
    const unsigned int nbTetrahedra = m_topology->getNbTetrahedra();

    auto tetrahedronInf = sofa::helper::getWriteAccessor(m_tetrahedronInfo);

    Coord dp[3], x0, sv;

    vector<Vec6d> vecStressSPK;
    vector<Vec6d> vecStressCauchy;

    vector<Real> vecStressVonMisesElement;


    for (unsigned int i = 0; i < nbTetrahedra; i++)
    {
        TetrahedronRestInformation* tetInfo = &tetrahedronInf[i];
        const Tetrahedron& ta = m_topology->getTetrahedron(i);

        x0 = x[ta[0]];

        /// compute the deformation gradient
        /// deformation gradient = sum of tensor product between vertex position and shape vector
        /// optimize by using displacement with first vertex
        dp[0] = x[ta[1]] - x0;
        sv = tetInfo->m_shapeVector[1];
        for (k = 0; k < 3; ++k)
        {
            for (l = 0; l < 3; ++l)
            {
                tetInfo->m_deformationGradient[k][l] = dp[0][k] * sv[l];
            }
        }
        for (j = 1; j < 3; ++j)
        {
            dp[j] = x[ta[j + 1]] - x0;
            sv = tetInfo->m_shapeVector[j + 1];

            for (k = 0; k < 3; ++k)
            {
                for (l = 0; l < 3; ++l)
                {
                    tetInfo->m_deformationGradient[k][l] += dp[j][k] * sv[l];
                }
            }
        }

        /// compute the right Cauchy-Green deformation matrix
        for (k = 0; k < 3; ++k)
        {
            for (l = k; l < 3; ++l)
            {
                tetInfo->C(k, l) =
                    tetInfo->m_deformationGradient(0, k) * tetInfo->m_deformationGradient(0, l) +
                    tetInfo->m_deformationGradient(1, k) * tetInfo->m_deformationGradient(1, l) +
                    tetInfo->m_deformationGradient(2, k) * tetInfo->m_deformationGradient(2, l);

            }
        }

        MatrixSym ID;
        ID.identity();


        /// Definition the Cauchy-Green strain Tensor E

        tetInfo->E = 0.5*(tetInfo->C-ID);


        if (globalParameters.anisotropyDirection.size() > 0)
        {
            tetInfo->m_fiberDirection = globalParameters.anisotropyDirection[0];
            Coord vectCa = tetInfo->C * tetInfo->m_fiberDirection;
            Real aDotCDota = dot(tetInfo->m_fiberDirection, vectCa);
            tetInfo->lambda = (Real)sqrt(aDotCDota);
        }
        const Coord areaVec = cross( dp[1], dp[2] );

        tetInfo->J = dot(areaVec, dp[0]) * tetInfo->m_volScale;
        tetInfo->trC = (Real)(tetInfo->C(0, 0) + tetInfo->C(1, 1) +
                               tetInfo->C(2, 2));



        tetInfo->m_SPKTensorGeneral.clear();
        tetInfo->m_CauchyStressTensor.clear();

        tetInfo->m_SPKStress.clear();
        tetInfo->m_CauchyStress.clear();

        m_myMaterial->deriveSPKTensor(tetInfo, globalParameters, tetInfo->m_SPKTensorGeneral,tetInfo->m_CauchyStressTensor,dt);

        /// Stress map in simulation (Voigt notation of the Stress tensor)
        vecStressSPK.push_back(tetInfo->m_SPKTensorGeneral.getVoigt());
        vecStressCauchy.push_back(tetInfo->m_CauchyStressTensor.getVoigt());

        /// Python Bindings
        tetInfo->m_SPKStress = tetInfo->m_SPKTensorGeneral.getVoigt();
        tetInfo->m_CauchyStress = tetInfo->m_CauchyStressTensor.getVoigt();

        /// Calculation of the General Von mises Stress per Element
        tetInfo->m_VonMisesStress = sqrt(0.5*((tetInfo->m_CauchyStress(0)-tetInfo->m_CauchyStress(1))*(tetInfo->m_CauchyStress(0)-tetInfo->m_CauchyStress(1)) +
            (tetInfo->m_CauchyStress(1)-tetInfo->m_CauchyStress(2))*(tetInfo->m_CauchyStress(1)-tetInfo->m_CauchyStress(2)) +
            (tetInfo->m_CauchyStress(2)-tetInfo->m_CauchyStress(0))*(tetInfo->m_CauchyStress(2)-tetInfo->m_CauchyStress(0)) +
            3*(tetInfo->m_CauchyStress(3)*tetInfo->m_CauchyStress(3)+tetInfo->m_CauchyStress(4)*tetInfo->m_CauchyStress(4) + 
            tetInfo->m_CauchyStress(5)*tetInfo->m_CauchyStress(5))
            ));

        vecStressVonMisesElement.push_back(tetInfo->m_VonMisesStress);


        for (l = 0; l < 4; ++l)
        {
            f[ta[l]] -= tetInfo->m_deformationGradient* (
                            tetInfo->m_SPKTensorGeneral * tetInfo->m_shapeVector[l]) * tetInfo->m_restVolume;
        }
        
    }



    d_stressSPK.beginEdit();
    d_stressSPK.setValue(vecStressSPK);
    d_stressSPK.endEdit();

    d_Cauchystress.beginEdit();
    d_Cauchystress.setValue(vecStressCauchy);
    d_Cauchystress.endEdit();

    d_stressVonMisesElement.beginEdit();
    d_stressVonMisesElement.setValue(vecStressVonMisesElement);
    d_stressVonMisesElement.endEdit();



    /// indicates that the next call to addDForce will need to update the stiffness matrix
    m_updateMatrix = true;
}

template <class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::updateTangentMatrix()
{

    SReal dt = this->getContext()->getDt();

    unsigned int k = 0, l;
    const unsigned int nbEdges = m_topology->getNbEdges();
    const type::vector<Edge>& edgeArray = m_topology->getEdges();

    auto edgeInf = sofa::helper::getWriteAccessor(m_edgeInfo);
    auto tetrahedronInf = sofa::helper::getWriteAccessor(m_tetrahedronInfo);

    const unsigned int nbTetrahedra = m_topology->getNbTetrahedra();
    const type::vector<Tetrahedron>& tetrahedronArray = m_topology->getTetrahedra();

    for (l = 0; l < nbEdges; l++)
    {
        edgeInf[l].DfDx.clear();
    }
    for (unsigned int i = 0; i < nbTetrahedra; i++)
    {
        TetrahedronRestInformation* tetInfo = &tetrahedronInf[i];
        Matrix3& df = tetInfo->m_deformationGradient;
        const BaseMeshTopology::EdgesInTetrahedron& te = m_topology->getEdgesInTetrahedron(i);

        /// describe the jth vertex index of triangle no i
        const Tetrahedron& ta = tetrahedronArray[i];
        for (unsigned int j = 0; j < 6; j++)
        {
            EdgeInformation* einfo = &edgeInf[te[j]];
            Edge e = m_topology->getLocalEdgesInTetrahedron(j);

            k = e[0];
            l = e[1];
            if (edgeArray[te[j]][0] != ta[k])
            {
                k = e[1];
                l = e[0];
            }
            Matrix3 &edgeDfDx = einfo->DfDx;

            const Coord& svl = tetInfo->m_shapeVector[l];
            const Coord& svk = tetInfo->m_shapeVector[k];

            Matrix3  M, N;
            MatrixSym outputTensor;
            N.clear();
            type::vector<MatrixSym> inputTensor;
            inputTensor.resize(3);

            for (int m = 0; m < 3; m++)
            {
                for (int n = m; n < 3; n++)
                {
                    inputTensor[0](m, n) = svl[m] * df[0][n] + df[0][m] * svl[n];
                    inputTensor[1](m, n) = svl[m] * df[1][n] + df[1][m] * svl[n];
                    inputTensor[2](m, n) = svl[m] * df[2][n] + df[2][m] * svl[n];
                }
            }

            for (int m = 0; m < 3; m++)
            {
                m_myMaterial->applyElasticityTensor(tetInfo, globalParameters, inputTensor[m],
                                                    outputTensor, dt);
                Coord vectortemp = df * (outputTensor * svk);
                Matrix3 Nv;
                for (int u = 0; u < 3; u++)
                {
                    Nv[u][m] = vectortemp[u];
                }
                N += Nv.transposed();
            }

            /// Now M
            Real productSD = 0;
            const Coord vectSD = tetInfo->m_SPKTensorGeneral * svk;
            productSD = dot(vectSD, svl);
            M[0][1] = M[0][2] = M[1][0] = M[1][2] = M[2][0] = M[2][1] = 0;
            M[0][0] = M[1][1] = M[2][2] = (Real)productSD;

            edgeDfDx += (M+N)*tetInfo->m_restVolume;
        }
    }
    m_updateMatrix=false;
}

template <class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
{
    auto df = sofa::helper::getWriteAccessor(d_df);
    const VecDeriv& dx = d_dx.getValue();
    const Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams, this->rayleighStiffness.getValue());

    const unsigned int nbEdges=m_topology->getNbEdges();
    const type::vector< Edge> &edgeArray=m_topology->getEdges() ;

    auto edgeInf = sofa::helper::getWriteAccessor(m_edgeInfo);

    /// if the  matrix needs to be updated
    if (m_updateMatrix)
    {
        this->updateTangentMatrix();
    }

    Deriv deltax;
    Deriv dv0,dv1;

    for (unsigned int l = 0; l < nbEdges; l++)
    {
        EdgeInformation* einfo = &edgeInf[l];
        unsigned int v0 = edgeArray[l][0];
        unsigned int v1 = edgeArray[l][1];

        deltax = dx[v0] - dx[v1];
        dv0 = einfo->DfDx * deltax;

        /// do the transpose multiply:
        dv1[0] = (Real)(deltax[0] * einfo->DfDx[0][0] + deltax[1] * einfo->DfDx[1][0] + deltax[2] * einfo->DfDx[2][0]);
        dv1[1] = (Real)(deltax[0] * einfo->DfDx[0][1] + deltax[1] * einfo->DfDx[1][1] + deltax[2] * einfo->DfDx[2][1]);
        dv1[2] = (Real)(deltax[0] * einfo->DfDx[0][2] + deltax[1] * einfo->DfDx[1][2] + deltax[2] * einfo->DfDx[2][2]);

        /// add forces
        df[v0] += dv1 * kFactor;
        df[v1] -= dv0 * kFactor;
    }
}

template<class DataTypes>
SReal TetrahedronViscoelasticityFEMForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const
{
    msg_warning() << "Method getPotentialEnergy not implemented yet.";
    return 0.0;
}

template <class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal k, unsigned int &offset)
{

    if (m_updateMatrix)
    {
        this->updateTangentMatrix();
    }

    const unsigned int nbEdges=m_topology->getNbEdges();
    const type::vector< Edge> &edgeArray=m_topology->getEdges() ;
    auto edgeInf = sofa::helper::getWriteAccessor(m_edgeInfo);

    for (unsigned int l = 0; l < nbEdges; l++)
    {
        EdgeInformation* einfo = &edgeInf[l];
        const Index noeud0 = edgeArray[l][0];
        const Index noeud1 = edgeArray[l][1];
        const unsigned int N0 = offset + 3 * noeud0;
        const unsigned int N1 = offset + 3 * noeud1;

        for (unsigned int i = 0; i < 3; i++)
        {
            for (unsigned int j = 0; j < 3; j++)
            {
                mat->add(N0 + i, N0 + j, + einfo->DfDx[j][i] * k);
                mat->add(N0 + i, N1 + j, - einfo->DfDx[j][i] * k);
                mat->add(N1 + i, N0 + j, - einfo->DfDx[i][j] * k);
                mat->add(N1 + i, N1 + j, + einfo->DfDx[i][j] * k);
            }
        }
    }
}

template<class DataTypes>
Mat<3,3,SReal> TetrahedronViscoelasticityFEMForceField<DataTypes>::getPhi(int tetrahedronIndex)
{
    auto tetrahedronInf = sofa::helper::getWriteAccessor(m_tetrahedronInfo);
    TetrahedronRestInformation* tetInfo = &tetrahedronInf[tetrahedronIndex];



    return tetInfo->m_deformationGradient;
}

template<class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::computeBBox(const core::ExecParams*, bool onlyVisible)
{
    if( !onlyVisible ) return;
    if (!this->mstate) return;

    this->f_bbox.setValue(this->mstate->computeBBox());
}

template<class DataTypes>
void TetrahedronViscoelasticityFEMForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields()) return;
    if (!this->mstate) return;

    const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);

    drawViscoelasticTets(vparams, x, m_topology, d_materialName.getValue());

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);
}

} // namespace 
