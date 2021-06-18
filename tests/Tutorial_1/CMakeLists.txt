

# if (MSVC AND BUILD_SHARED_LIBS)
#   # copy dll
#   add_custom_command(TARGET test_c POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:IPhreeqc> $<TARGET_FILE_DIR:test_c>
#   )
# endif()

configure_file(Case1.cgn Case1.cgn COPYONLY)

configure_file(iRICZone.1.edge iRICZone.1.edge COPYONLY)

configure_file(iRICZone.1.v.edge iRICZone.1.v.edge COPYONLY)

# test compile and run
message(STATUS $<TARGET_FILE:SToRM>)
##add_test(Test-Tutorial_1 $<TARGET_FILE:SToRM> Case1.cgn)
