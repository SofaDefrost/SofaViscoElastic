cmake_minimum_required(VERSION 3.12)
project(BindingsSofaViscoElastic)



sofa_find_package(SofaPython3 REQUIRED COMPONENTS Bindings.Sofa)
sofa_find_package(SofaViscoElastic REQUIRED)


find_file(SofaPython3Tools NAMES "SofaPython3/lib/cmake/SofaPython3/SofaPython3Tools.cmake")

if(SofaPython3_FOUND)
set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Binding_TetrahedronViscoElasticityFEMForceField.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Module_TetrahedronViscoElasticityFEMForceField.cpp

)

set(HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Binding_TetrahedronViscoElasticityFEMForceField.h

)
endif()



SP3_add_python_module(
    TARGET       ${PROJECT_NAME}
    PACKAGE      Bindings
    MODULE       SofaViscoElastic
    DESTINATION  .
    SOURCES      ${SOURCE_FILES}
    HEADERS      ${HEADER_FILES}
    DEPENDS      Sofa.Core Sofa.GUI.Common SofaPython3::Plugin SofaPython3::Bindings.Sofa.Core SofaViscoElastic ${SUPPORTED_GUIS} 
)




sofa_create_component_in_package_with_targets(
    COMPONENT_NAME ${PROJECT_NAME}
    COMPONENT_VERSION ${SofaViscoElastic_VERSION}
    PACKAGE_NAME SofaViscoElastic
    TARGETS ${PROJECT_NAME}
    )