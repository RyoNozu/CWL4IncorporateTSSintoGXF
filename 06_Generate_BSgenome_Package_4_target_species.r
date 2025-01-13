# generate BSgenome Package for target species
BSgenomeForge::forgeBSgenomeDataPkg("./BSgenome.Htrimaculatus.Htriv1.1-seed.txt", replace=TRUE)
devtools::build("./BSgenome.Htrimaculatus.inhouse.Htriv1")
devtools::check_built("BSgenome.Htrimaculatus.inhouse.Htriv1_1.0.0.tar.gz")
devtools::install_local("BSgenome.Htrimaculatus.inhouse.Htriv1_1.0.0.tar.gz")
