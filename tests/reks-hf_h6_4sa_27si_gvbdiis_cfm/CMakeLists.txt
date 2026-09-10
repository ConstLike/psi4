include(TestingMacros)

add_regression_test(reks-hf_h6_4sa_27si_gvbdiis_cfm "psi;reks;noc1;long")

# Heavy test: skipped by default. To enable, drop this line.
set_tests_properties(reks-hf_h6_4sa_27si_gvbdiis_cfm PROPERTIES DISABLED TRUE)
