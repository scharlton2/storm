# copy files into build directory for testing
configure_file(Case1.cgn Case1.cgn COPYONLY)
configure_file(iRICZone.1.edge iRICZone.1.edge COPYONLY)
configure_file(iRICZone.1.v.edge iRICZone.1.v.edge COPYONLY)

# test run
add_test(NAME Test-Tutorial_1
  COMMAND $<TARGET_FILE:SToRM> Case1.cgn)
