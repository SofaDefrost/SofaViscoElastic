cmake_minimum_required(VERSION 3.12)
project(SofaViscoElastic VERSION 1.0)


set(SOFAVISCOELASTIC_SOURCE_DIR src/SofaViscoElastic)

set(HEADER_FILES
    ${SOFAVISCOELASTIC_SOURCE_DIR}/config.h.in
    ${SOFAVISCOELASTIC_SOURCE_DIR}/init.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoelasticityFEMDrawing.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoelasticityFEMForceField.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoelasticityFEMForceField.inl


    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/BaseViscoelasticMaterial.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/MaxwellFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/KelvinVoigtFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/SLSKelvinVoigtFirstOrder.h
    ${SOFAVISCOELASTIC_SOURCE_DIR}/material/SLSMaxwellFirstOrder.h
        ${SOFAVISCOELASTIC_SOURCE_DIR}/material/Burgers.h

)

set(SOURCE_FILES
    ${SOFAVISCOELASTIC_SOURCE_DIR}/init.cpp
    ${SOFAVISCOELASTIC_SOURCE_DIR}/TetrahedronViscoelasticityFEMForceField.cpp
    )

sofa_find_package(Sofa.Simulation.Core REQUIRED)

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
