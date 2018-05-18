APPEND_SOURCE_FILES(SOURCES)

add_library(SmallDeformationNonlocalHydroMechanics ${SOURCES})
target_link_libraries(SmallDeformationNonlocalHydroMechanics PUBLIC ProcessLib)

include(Tests.cmake)
