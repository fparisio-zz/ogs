AddTest(
    NAME SDN_HydroMechanics_HML_square_1e2_quad8_confined_compression
    PATH SmallDeformationNonlocalHydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB square_1e2_pcs_0_ts_*.vtu displacement displacement 1e-15 1e-15
    GLOB square_1e2_pcs_0_ts_*.vtu pressure pressure 1e-15 1e-15
    GLOB square_1e2_pcs_0_ts_*.vtu velocity velocity 1e-15 1e-15
    GLOB square_1e2_pcs_0_ts_*.vtu HydraulicFlow HydraulicFlow 1e-15 1e-15
    GLOB square_1e2_pcs_0_ts_*.vtu NodalForces NodalForces 1e-15 1e-15
)
