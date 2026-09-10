include(TestingMacros)

add_regression_test(reks-sec66_02 "psi;reks;noc1;long")

# Heavy test: skipped by default. To enable, drop this line.
set_tests_properties(reks-sec66_02 PROPERTIES DISABLED TRUE)
