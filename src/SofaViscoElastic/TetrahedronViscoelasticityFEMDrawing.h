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
#include <sofa/core/visual/VisualParams.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/Topology.h>


namespace sofa::SofaViscoElastic
{

template <class VecCoord>
void drawViscoelasticTets(const core::visual::VisualParams* vparams, const VecCoord& x, core::topology::BaseMeshTopology* topology, const std::string& materialName)
{
    std::vector<type::Vec3 > points[4];
    for(auto& p : points)
    {
        p.reserve(3 * topology->getNbTetrahedra());
    }

    for(core::topology::Topology::TetrahedronID i = 0 ; i < topology->getNbTetrahedra();++i)
    {
        const auto t = topology->getTetrahedron(i);

        Index a = t[0];
        Index b = t[1];
        Index c = t[2];
        Index d = t[3];
        const auto center = (x[a] + x[b] + x[c] + x[d]) * 0.125;
        const auto pa = (x[a] + center) * 0.666667;
        const auto pb = (x[b] + center) * 0.666667;
        const auto pc = (x[c] + center) * 0.666667;
        const auto pd = (x[d] + center) * 0.666667;

        points[0].push_back(pa);
        points[0].push_back(pb);
        points[0].push_back(pc);

        points[1].push_back(pb);
        points[1].push_back(pc);
        points[1].push_back(pd);

        points[2].push_back(pc);
        points[2].push_back(pd);
        points[2].push_back(pa);

        points[3].push_back(pd);
        points[3].push_back(pa);
        points[3].push_back(pb);
    }

    sofa::type::RGBAColor color1;
    sofa::type::RGBAColor color2;
    sofa::type::RGBAColor color3;
    sofa::type::RGBAColor color4;

    if (materialName=="MaxwellFirstOrder") {
        color1 = sofa::type::RGBAColor(0.0,1.0,0.0,1.0);
        color2 = sofa::type::RGBAColor(0.5,1.0,0.0,1.0);
        color3 = sofa::type::RGBAColor(1.0,1.0,0.0,1.0);
        color4 = sofa::type::RGBAColor(1.0,1.0,0.5,1.0);
    }
    else if (materialName=="KelvinVoigtFirstOrder"){
        color1 = sofa::type::RGBAColor(1.0,0.0,0.0,1.0);
        color2 = sofa::type::RGBAColor(1.0,0.0,0.5,1.0);
        color3 = sofa::type::RGBAColor(1.0,1.0,0.0,1.0);
        color4 = sofa::type::RGBAColor(1.0,0.5,1.0,1.0);
    }
    else if (materialName=="SLSMaxwellFirstOrder"){
        color1 = sofa::type::RGBAColor(0.0,1.0,1.0,1.0);
        color2 = sofa::type::RGBAColor(0.5,0.0,1.0,1.0);
        color3 = sofa::type::RGBAColor(1.0,0.0,1.0,1.0);
        color4 = sofa::type::RGBAColor(1.0,0.5,1.0,1.0);
    }
    else if (materialName=="SLSKelvinVoigtFirstOrder"){
        color1 = sofa::type::RGBAColor(0.0,1.0,0.0,1.0);
        color2 = sofa::type::RGBAColor(0.0,1.0,0.5,1.0);
        color3 = sofa::type::RGBAColor(0.0,1.0,1.0,1.0);
        color4 = sofa::type::RGBAColor(0.5,1.0,1.0,1.0);
    }
    else if (materialName=="Burgers"){
        color1 = sofa::type::RGBAColor(0.0,1.0,0.0,1.0);
        color2 = sofa::type::RGBAColor(0.5,1.0,0.0,1.0);
        color3 = sofa::type::RGBAColor(1.0,1.0,0.0,1.0);
        color4 = sofa::type::RGBAColor(1.0,1.0,0.5,1.0);
    }
    else {
        color1 = sofa::type::RGBAColor(0.0,1.0,0.0,1.0);
        color2 = sofa::type::RGBAColor(0.5,1.0,0.0,1.0);
        color3 = sofa::type::RGBAColor(1.0,1.0,0.0,1.0);
        color4 = sofa::type::RGBAColor(1.0,1.0,0.5,1.0);
    }

    vparams->drawTool()->drawTriangles(points[0], color1);
    vparams->drawTool()->drawTriangles(points[1], color2);
    vparams->drawTool()->drawTriangles(points[2], color3);
    vparams->drawTool()->drawTriangles(points[3], color4);
}

}
