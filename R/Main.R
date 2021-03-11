#Main
#library(tidyverse)

#' Main function to run PARIS pipeline
#' @param imp.score importance score methods; NB raw/Gini for Boruta; Ginicorrected RF
#' @param omic Features data to use: Expression or Mutation
#' @param outdir Write Boruta/Random Forest Feature selection models and selected gene paris in the outdir
#' @param genelist_dep Gene names vector to use as response variables (dependency scores), ALL if NULL
#' @param genelist_feat Gene names vector to use as features (Expression or Mutation), ALL if NULL
#' @return Write Boruta/Random Forest Feature selection models and selected gene paris in the outdir
#' @export
#' @example
#' \dontrun{
#' run_PARIS(imp.score = "raw", omic = "expression", outdir = "/PARIS_test", genelist_dep = c("WRN"))
#' }
run_PARIS <- function(imp.score = c("raw", "Gini", "Ginicorrected"), omic = c("mutation", "expression"), outdir, genelist_dep, genelist_feat){

  imp.score <- match.arg(imp.score, c("raw", "Gini", "Ginicorrected"))
  omic <- match.arg(omic, c("mutation", "expression"))

  dir.create(outdir, showWarnings = T)

  if (missing(genelist_dep) & missing(genelist_feat)) { #ALL dependecies vs ALL features

    print("No genelists selected")

    myselectedmerge <- mynewmerge

    ###remove cell lines not screened
    mynewexp <- semi_join(mynewexp, myselectedmerge, by = "DepMap_ID")

    ###remove cell lines not screened
    mynewmut <- semi_join(mynewmut, myselectedmerge, by = "DepMap_ID")
  }

  if (missing(genelist_feat) & !missing(genelist_dep)) { #subset dependencies vs ALL features

    print("Only dependencies genelist selected")

    #subset
    myselectedmerge <- filter(mynewmerge, gene %in% genelist_dep)

    ###remove cell lines not screened
    mynewexp <- semi_join(mynewexp, myselectedmerge, by = "DepMap_ID")


    ###remove cell lines not screened
    mynewmut <- semi_join(mynewmut, myselectedmerge, by = "DepMap_ID")
  }

  if (!missing(genelist_feat) & missing(genelist_dep)) { #ALL dependencies vs subset features

    print("Only features genelist selected")

    myselectedmerge <- mynewmerge

    ###remove cell lines not screened
    mynewexp <- semi_join(mynewexp, myselectedmerge, by = "DepMap_ID")

    mynewexp <- filter(mynewexp, gene %in% genelist_feat)

    ###remove cell lines not screened
    mynewmut <- semi_join(mynewmut, myselectedmerge, by = "DepMap_ID")

    mynewmut <- filter(mynewmut, gene %in% genelist_feat)
  }

  if (!missing(genelist_feat) & !missing(genelist_dep)) { #subset dependencies vs subset features

    print("Dependencies and Features genelists selected")

    #subset
    myselectedmerge <- filter(mynewmerge, gene %in% genelist_dep)

    ###remove cell lines not screened
    mynewexp <- semi_join(mynewexp, myselectedmerge, by = "DepMap_ID")

    mynewexp <- filter(mynewexp, gene %in% genelist_feat)

    ###remove cell lines not screened
    mynewmut <- semi_join(mynewmut, myselectedmerge, by = "DepMap_ID")

    mynewmut <- filter(mynewmut, gene %in% genelist_feat)
  }

  if (omic == "mutation") {
    if (imp.score == "raw") {
      RF_mut(data = myselectedmerge, mut = mynewmut, imp.score = "raw", outdir = outdir)
    }
    if (imp.score == "Gini") {
      RF_mut(data = myselectedmerge, mut = mynewmut, imp.score = "Gini", outdir = outdir)
    }
    if (imp.score == "Ginicorrected") {
      RF_mut(data = myselectedmerge, mut = mynewmut, imp.score = "Ginicorrected", outdir = outdir)
    }
  }
  if (omic == "expression") {
    if (imp.score == "raw") {
      RF_exp(data = myselectedmerge, exp = mynewexp, imp.score = "raw", outdir = outdir)
    }
    if (imp.score == "Gini") {
      RF_exp(data = myselectedmerge, exp = mynewexp, imp.score = "Gini", outdir = outdir)
    }
    if (imp.score == "Ginicorrected") {
      RF_exp(data = myselectedmerge, exp = mynewexp, imp.score = "Ginicorrected", outdir = outdir)
    }
  }
}

