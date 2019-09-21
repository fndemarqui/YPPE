

devtools::load_all(".")
devtools::document()
devtools::install()
#devtools::install(quick=TRUE)


devtools::build()
#devtools::build_vignettes()
devtools::check()
devtools::missing_s3()
devtools::release_checks()
# devtools::check_win_devel()
# devtools::spell_check()
devtools::test()
usethis::use_build_ignore(c("script_gastric.R", "script_ovarian.R", "script_ipass.R"))
# devtools::reload()
# devtools::run_examples()


#install.packages("crossSurv_0.0.0.9000.tar.gz", repos = FALSE, dependencies=TRUE)
