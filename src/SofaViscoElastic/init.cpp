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
#include <SofaViscoElastic/init.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/PluginManager.h>
using sofa::helper::system::PluginManager;

#include <fstream>

namespace sofa::SofaViscoElastic
{
    
extern "C" {
    SOFAVISCOELASTIC_API void initExternalModule();
    SOFAVISCOELASTIC_API const char* getModuleName();
    SOFAVISCOELASTIC_API const char* getModuleVersion();
    SOFAVISCOELASTIC_API const char* getModuleComponentList();
}

void initExternalModule()
{
    init();
}

const char* getModuleName()
{
    return "SofaViscoElastic";
}

const char* getModuleVersion()
{
    return "1.0";
}

void init()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

const char* getModuleComponentList()
{
    /// string containing the names of the classes provided by the plugin
    static std::string classes = core::ObjectFactory::getInstance()->listClassesFromTarget(sofa_tostring(SOFA_TARGET));
    return classes.c_str();
}
} // namespace sofa::component::solidmechanics::fem::viscoelastic
