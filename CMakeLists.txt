cmake_minimum_required(VERSION 3.12)
project(SofaViscoElastic VERSION 1.0)


set(SOFAVISCOELASTIC_SOURCE_DIR src/SofaViscoElastic)

set(HEADER_FILES
    ${SOFAVISCOELASTIC_SOURCE_DIR}/config.h.in
    ${SOFAVISCOELASTIC_SOURCE_DIR}/init.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoelasticityFEMDrawing.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoelasticityFEMForceField.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoelasticityFEMForceField.inl
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoHyperelasticityFEMForceField.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoHyperelasticityFEMForceField.inl

    ## VISCOELASTIC MATERIALS
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/BaseViscoelasticMaterial.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/MaxwellFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/MaxwellSecondOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/SLSMaxwellSecondOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/KelvinVoigtFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/KelvinVoigtSecondOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/SLSKelvinVoigtSecondOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/SLSKelvinVoigtFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/SLSMaxwellFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/Burgers.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscoelastic/LinearElastic.h

    ## VISCOHYPERELASTIC MATERIALS 
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscohyperelastic/BaseViscoHyperelasticMaterial.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscohyperelastic/SLSNeoHookeanFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscohyperelastic/SLSNeoHookeanSecondOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscohyperelastic/SLSStableNeoHookeanFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscohyperelastic/SLSStableNeoHookeanSecondOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscohyperelastic/SLSMooneyRivlinFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscohyperelastic/SLSMooneyRivlinSecondOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscohyperelastic/SLSOgdenFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/viscohyperelastic/SLSOgdenSecondOrder.h


)

set(SOURCE_FILES
    ${SOFAVISCOELASTIC_SOURCE_DIR}/init.cpp
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoelasticityFEMForceField.cpp
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoHyperelasticityFEMForceField.cpp

    )

find_package(Sofa.Simulation.Core REQUIRED)

sofa_find_package(SofaPython3 QUIET)


add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Simulation.Core )


sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${Sofa_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR "${PROJECT_NAME}"
     RELOCATABLE "plugins"
)

if(SofaPython3_FOUND)
    add_subdirectory(Bindings)
endif()
