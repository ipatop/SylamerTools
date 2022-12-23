## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
#install library
#devtools::install_github("ipatop/SylamerTools")
library(SylamerTools)

## -----------------------------------------------------------------------------
#This is how the input looks like
head(read.table(file = paste0(system.file("extdata", package = "SylamerTools"),"/res.txt"),header = T))

## -----------------------------------------------------------------------------
sort_forSylamer(DEres = paste0(system.file("extdata", package = "SylamerTools"),"/res.txt"),sortBy = "log2FoldChange",outName =paste0(system.file("extdata", package = "SylamerTools") ,"/ForSylamer.txt"))

head(read.table(file = paste0(system.file("extdata", package = "SylamerTools"),"/ForSylamer.txt") ,header = T))

## -----------------------------------------------------------------------------
create_sylamer_script(utrFile =  paste0(system.file("extdata", package = "SylamerTools"),"/dm3_flybase_3utrs.fa"), overwrite = T, DEsorted=paste0(system.file("extdata", package = "SylamerTools"),c("/ForSylamer.txt","/ForSylamer_2.txt")))

## -----------------------------------------------------------------------------
read.delim("sylamer_command.sh",header = F)

## -----------------------------------------------------------------------------
read.delim(paste0(system.file("extdata", package = "SylamerTools"),"/L2FcSh_sylamer.output.tab"))[1:3,1:3]

## -----------------------------------------------------------------------------
sylamerplots(oligosChosen=c("TGTAAA","GATGCT"),sylamer.out.tab=paste0(system.file("extdata", package = "SylamerTools"),"/L2FcSh_sylamer.output.tab"),outName ="shRNA_fake")

## -----------------------------------------------------------------------------
SylamerTools::sylamerplots(oligosChosen=c("TGTAAA","GATGCT"),sylamer.out.tab=paste0(system.file("extdata", package = "SylamerTools"),"/L2FcSh_sylamer.output.tab"),greylines = T,outName = "shRNA_fake_markingGrey",free_ylim = T)

