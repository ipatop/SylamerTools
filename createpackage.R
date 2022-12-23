#create th package
setwd("~/Documents/shRNA_design_circRNA/SylamerTools/")
usethis::create_package("SylamerTools")
#this creates a new project with the requieremnts for R package
#NAMESpace store the functinos that will be imported etc
#create an r script, this will be the name of the R script with
usethis::use_r("Sylamer_Tools")


#edit rprofile to add name etc
usethis::edit_r_profile()

#add needed package for IMPORT
# usethis::use_package("assertthat")
# usethis::use_package("Seurat")
# usethis::use_package("SeuratDisk")
# usethis::use_package("grDevices")
# usethis::use_package("ggplot2")
# usethis::use_package("hdf5r")
# usethis::use_package("loomR")
# usethis::use_package("scater")
# usethis::use_package("patchwork")
# usethis::use_package("SeuratObject")
# usethis::use_package("dplyr")
# usethis::use_package("tidyr")
# usethis::use_package("cowplot")

usethis::use_package("dplyr")

usethis::use_package("tidyr")

usethis::use_package("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
usethis::use_package("BSgenome.Dmelanogaster.UCSC.dm6")
usethis::use_package("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
usethis::use_package("BSgenome.Dmelanogaster.UCSC.dm3")
usethis::use_package("TxDb.Hsapiens.UCSC.hg19.knownGene")
usethis::use_package("BSgenome.Hsapiens.UCSC.hg19")
usethis::use_package("TxDb.Hsapiens.UCSC.hg38.knownGene")
usethis::use_package("BSgenome.Hsapiens.UCSC.hg38")
usethis::use_package("TxDb.Mmusculus.UCSC.mm10.knownGene")
usethis::use_package("BSgenome.Mmusculus.UCSC.mm10")
usethis::use_package("TxDb.Rnorvegicus.UCSC.rn4.ensGene")
usethis::use_package("BSgenome.Rnorvegicus.UCSC.rn4")
usethis::use_package("Rsubread")
usethis::use_package("GenomicFeatures")
usethis::use_package("Biostrings")

#put on github
usethis::git_sitrep()
usethis::use_git()
usethis::use_github()
#usethis::create_from_github("ipatop/circRNAshRNAdesign", fork = FALSE)


###document###
"Add roxygen comments to your .R files.
Run devtools::document() to convert roxygen comments to .Rd files.
Load the current version fo the package with devtools::load_all()
Preview documentation with ?.
Rinse and repeat until the documentation looks the way you want."
devtools::load_all()
devtools::document()

#create the vignet
usethis::use_vignette("Sylamer_tools",title = "Sylamer tools")
#buiold the vignete
devtools::build_vignettes()

#createreadme
usethis::use_readme_md()

#add external file
usethis::use_data_raw()

#add data https://r-pkgs.org/data.html
usethis::use_data(my_pkg_data)


#modify
#load it in the memory
devtools::load_all()
#test it
devtools::check()
#license
usethis::use_mit_license("Ines Patop")


