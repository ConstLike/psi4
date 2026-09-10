include(TestingMacros)

add_regression_test(reks-h4_g2_4sa_9si "psi;reks;noc1;long")

# Heavy test: skipped by default. To enable, drop this line.
set_tests_properties(reks-h4_g2_4sa_9si PROPERTIES DISABLED TRUE)
