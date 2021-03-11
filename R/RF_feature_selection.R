#RF feature selection
# library(Boruta)
# library(ranger)

#' Random Forest Feature Selection - Mutated genes -> Dependent gene
#' @param data Dataframe with dependecy scores
#' @param mut Dataframe with mutational data
#' @param imp.score importance score methods; NB raw/Gini for Boruta; Ginicorrected RF
#' @param outdir Directory where to write Models
#' @return Write Boruta/Random Forest Feature selection models and selected gene paris in the outdir
#' @export
#' @example
#' \dontrun{
#' RF_mut(data = myselectedmerge, mut = mynewmut, imp.score = "raw", outdir = outdir)
#' }
RF_mut <- function(data, mut, imp.score = c("raw", "Gini", "Ginicorrected"), outdir){
  imp.score <- match.arg(imp.score)

  myselectedmerge <- data
  mynewmut <- mut
  myoutdir <- outdir

  if (imp.score == "raw") {

    dir.create(paste(myoutdir, "/Boruta_models_mut_ImpRfRaw", sep = ""),showWarnings = F)
    ###run Boruta mut ImpRfRaw
    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      ###merge dependency and mutation
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct() %>% inner_join(mynewmut, by = "DepMap_ID")
      ####some tranformation
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, gene.y)
      mybigmergeuniq$mut <- 1

      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, mut, fill = 0) %>% select(-DepMap_ID, -gene.x)
      results <- Boruta(y = mybigene$value, x = mybigene[2:ncol(mybigene)], maxRuns = 500, getImp = getImpRfRaw)
      save(results, file = paste(paste(myoutdir, "/Boruta_models_mut_ImpRfRaw/model", sep = ""), i, "rds", sep = "."))
      pdf(paste(paste(myoutdir, "/Boruta_models_mut_ImpRfRaw/model", sep = ""), i, "pdf", sep = "."), width = 20)
      plotboruta(results, las=2, cex.axis=0.5, main = i, xlab = "Mutations")
      dev.off()
    }

    ####calculate correlation with important mutations
    res <- data.frame()

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, mut, fill = 0) %>% select(-DepMap_ID, -gene.x)
      load(file = paste(paste(myoutdir, "/Boruta_models_mut_ImpRfRaw/model", sep = ""), i, "rds", sep = "."))
      muts <- getSelectedAttributes(results, withTentative = F)
      for (c in muts) {
        coeff <- cor.test(mybigene$value, mybigene[[c]])$estimate
        pval <- cor.test(mybigene$value, mybigene[[c]])$p.value
        #if (coeff > 0.3 | coeff < -0.3) {
        res <- rbind(res, data.frame(i, c, coeff, pval))
        #}
      }
    }
    colnames(res)[1] <- "gene"
    res$group <- "mr"
    imp <- load_imp_scores2(paste(myoutdir, "/Boruta_models_mut_ImpRfRaw/", sep = ""))
    res <- inner_join(res, imp)
    write_csv(res, path = paste(myoutdir, "/PairsMut_ImpRfRaw.csv", sep = ""))
  }

  if (imp.score == "Gini") {

    dir.create(paste(myoutdir, "/Boruta_models_mut_ImpRfGini", sep = ""), showWarnings = F)
    ###run Boruta mut ImpRfGini
    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      ###merge dependency and mutation
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct() %>% inner_join(mynewmut, by = "DepMap_ID")
      ####some tranformation
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, gene.y)
      mybigmergeuniq$mut <- 1

      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, mut, fill = 0) %>% select(-DepMap_ID, -gene.x)
      results <- Boruta(y = mybigene$value, x = mybigene[2:ncol(mybigene)], maxRuns = 500, getImp = getImpRfGini)
      save(results, file = paste(paste(myoutdir, "/Boruta_models_mut_ImpRfGini/model", sep = ""), i, "rds", sep = "."))
      pdf(paste(paste(myoutdir, "/Boruta_models_mut_ImpRfGini/model", sep = ""), i, "pdf", sep = "."), width = 20)
      plotboruta(results, las=2, cex.axis=0.5, main = i, xlab = "Mutations")
      dev.off()
      #mylist[[i]] <- getSelectedAttributes(results, withTentative = F)

    }

    ####calculate correlation with important mutations

    res <- data.frame()

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, mut, fill = 0) %>% select(-DepMap_ID, -gene.x)
      load(file = paste(paste(myoutdir, "/Boruta_models_mut_ImpRfGini/model", sep = ""), i, "rds", sep = "."))
      muts <- getSelectedAttributes(results, withTentative = F)
      for (c in muts) {
        coeff <- cor.test(mybigene$value, mybigene[[c]])$estimate
        pval <- cor.test(mybigene$value, mybigene[[c]])$p.value
        #if (coeff > 0.3 | coeff < -0.3) {
        res <- rbind(res, data.frame(i, c, coeff, pval))
        #}
      }
    }
    colnames(res)[1] <- "gene"
    res$group <- "mg"
    imp <- load_imp_scores2(paste(myoutdir, "/Boruta_models_mut_ImpRfGini/", sep = ""))
    res <- inner_join(res, imp)
    write_csv(res, path = paste(myoutdir, "/PairsMut_ImpRfGini.csv", sep = ""))
  }

  if (imp.score == "Ginicorrected") {

    dir.create(paste(myoutdir, "/ranger_models_mut_Ginicorrected", sep = ""),showWarnings = F)
    ###run ranger impurity_corrected mut

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      ###merge dependency and mutation
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct() %>% inner_join(mynewmut, by = "DepMap_ID")
      ####some tranformation
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, gene.y)
      mybigmergeuniq$mut <- 1
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, mut, fill = 0) %>% select(-DepMap_ID, -gene.x)
      results <- ranger(dependent.variable.name = "value", data = mybigene, importance = 'impurity_corrected')
      save(results, file = paste(paste(myoutdir, "/ranger_models_mut_Ginicorrected/model", sep = ""), i, "rds", sep = "."))
      #mylist[[i]] <- getSelectedAttributes(results, withTentative = F)

    }

    ####calculate correlation with important mutations

    res <- data.frame()

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      ###merge dependency and mutation
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct() %>% inner_join(mynewmut, by = "DepMap_ID")
      ####some tranformation
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, gene.y)
      mybigmergeuniq$mut <- 1
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, mut, fill = 0) %>% select(-DepMap_ID, -gene.x)
      load(file = paste(paste(myoutdir, "/ranger_models_mut_Ginicorrected/model", sep = ""), i, "rds", sep = "."))
      inv <- try(importance_pvalues(results), silent=TRUE)
      if ('try-error' %in% class(inv)) next
      muts <- importance_pvalues(results) %>% as.data.frame() %>% subset(pvalue == 0) %>% rownames()
      for (c in muts) {
        coeff <- cor.test(mybigene$value, mybigene[[c]])$estimate
        pval <- cor.test(mybigene$value, mybigene[[c]])$p.value
        #if (coeff > 0.3 | coeff < -0.3) {
        res <- rbind(res, data.frame(i, c, coeff, pval, results$variable.importance[c]))
        #}
      }
    }
    colnames(res)[1] <- "gene"
    res$group <- "mc"
    colnames(res)[5] <- "importance"
    write_csv(res, path = paste(myoutdir, "/PairsMut_imp_corr.csv", sep = ""))

  }
}

#' Random Forest Feature Selection - Genes expression -> Dependent gene
#' @param data Dataframe with dependecy scores
#' @param exp Dataframe with expression data
#' @param imp.score importance score methods; NB raw/Gini for Boruta; Ginicorrected RF
#' @param outdir Directory where to write Models
#' @return Write Boruta/Random Forest Feature selection models and selected gene paris in the outdir
#' @export
#' @example
#' \dontrun{
#' RF_exp(data = myselectedmerge, mut = mynewexp, imp.score = "raw", outdir = outdir)
#' }
RF_exp <- function(data, exp, imp.score = c("raw", "Gini", "Ginicorrected"), outdir){
  imp.score <- match.arg(imp.score)

  myselectedmerge <- data
  mynewexp <- exp
  myoutdir <- outdir

  if (imp.score == "raw") {

    dir.create(paste(myoutdir, "/Boruta_models_exp_ImpRfRaw", sep = ""), showWarnings = F)
    ###run Boruta exp ImpRfRaw

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct()  %>% inner_join(mynewexp, by = "DepMap_ID")
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, TPM, gene.y)
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, TPM, fill = 0) %>% select(-DepMap_ID, -gene.x)
      results <- Boruta(y = mybigene$value, x = mybigene[2:ncol(mybigene)], maxRuns = 500, getImp = getImpRfRaw)
      save(results, file = paste(paste(myoutdir, "/Boruta_models_exp_ImpRfRaw/model", sep = ""), i, "rds", sep = "."))
      pdf(paste(paste(myoutdir, "/Boruta_models_exp_ImpRfRaw/model", sep = ""), i, "pdf", sep = "."), width = 20)
      plotboruta(results, las=2, cex.axis=0.5, main = i, xlab = "Gene")
      dev.off()
      #mylist[[i]] <- getSelectedAttributes(results, withTentative = F)

    }

    ####calculate correlation with important exp genes

    res <- data.frame()

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct()  %>% inner_join(mynewexp, by = "DepMap_ID")
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, TPM, gene.y)
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, TPM, fill = 0) %>% select(-DepMap_ID, -gene.x)
      load(file = paste(paste(myoutdir, "/Boruta_models_exp_ImpRfRaw/model", sep = ""), i, "rds", sep = "."))
      muts <- getSelectedAttributes(results, withTentative = F)
      for (c in muts) {
        coeff <- cor.test(mybigene$value, mybigene[[c]])$estimate
        pval <- cor.test(mybigene$value, mybigene[[c]])$p.value
        #if (coeff > 0.3 | coeff < -0.3) {
        res <- rbind(res, data.frame(i, c, coeff, pval))
        #}
      }
    }
    colnames(res)[1] <- "gene"
    res$group <- "er"
    imp <- load_imp_scores2(paste(myoutdir, "/Boruta_models_exp_ImpRfRaw/", sep = ""))
    res <- inner_join(res, imp)
    write_csv(res, path = paste(myoutdir, "/PairsExp_ImpRfRaw.csv", sep = ""))

  }

  if (imp.score == "Gini") {

    dir.create(paste(myoutdir, "/Boruta_models_exp_ImpRfGini", sep = ""), showWarnings = F)
    ###run Boruta exp ImpRfGini

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct()  %>% inner_join(mynewexp, by = "DepMap_ID")
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, TPM, gene.y)
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, TPM, fill = 0) %>% select(-DepMap_ID, -gene.x)
      results <- Boruta(y = mybigene$value, x = mybigene[2:ncol(mybigene)], maxRuns = 500, getImp = getImpRfGini)
      save(results, file = paste(paste(myoutdir, "/Boruta_models_exp_ImpRfGini/model", sep = ""), i, "rds", sep = "."))
      pdf(paste(paste(myoutdir, "/Boruta_models_exp_ImpRfGini/model", sep = ""), i, "pdf", sep = "."), width = 20)
      plotboruta(results, las=2, cex.axis=0.5, main = i, xlab = "Gene")
      dev.off()
      #mylist[[i]] <- getSelectedAttributes(results, withTentative = F)

    }

    ####calculate correlation with important exp genes

    res <- data.frame()

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct()  %>% inner_join(mynewexp, by = "DepMap_ID")
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, TPM, gene.y)
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, TPM, fill = 0) %>% select(-DepMap_ID, -gene.x)
      load(file = paste(paste(myoutdir, "/Boruta_models_exp_ImpRfGini/model", sep = ""), i, "rds", sep = "."))
      muts <- getSelectedAttributes(results, withTentative = F)
      for (c in muts) {
        coeff <- cor.test(mybigene$value, mybigene[[c]])$estimate
        pval <- cor.test(mybigene$value, mybigene[[c]])$p.value
        #if (coeff > 0.3 | coeff < -0.3) {
        res <- rbind(res, data.frame(i, c, coeff, pval))
        #}
      }
    }
    colnames(res)[1] <- "gene"
    res$group <- "eg"
    imp <- load_imp_scores2(paste(myoutdir, "/Boruta_models_exp_ImpRfGini/", sep = ""))
    res <- inner_join(res, imp)
    write_csv(res, path = paste(myoutdir, "/PairsExp_ImpRfGini.csv", sep = ""))
  }

  if (imp.score == "Ginicorrected") {

    dir.create(paste(myoutdir, "/ranger_models_exp_Ginicorrected", sep = ""), showWarnings = F)
    ###run ranger impurity_corrected exp

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct()  %>% inner_join(mynewexp, by = "DepMap_ID")
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, TPM, gene.y)
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, TPM, fill = 0) %>% select(-DepMap_ID, -gene.x)
      results <- ranger(dependent.variable.name = "value", data = mybigene, importance = 'impurity_corrected')
      save(results, file = paste(paste(myoutdir, "/ranger_models_exp_Ginicorrected/model", sep = ""), i, "rds", sep = "."))
      #mylist[[i]] <- getSelectedAttributes(results, withTentative = F)

    }

    ####calculate correlation with important exp genes

    res <- data.frame()

    for (i in (unique(myselectedmerge$gene))) {
      set.seed(43)
      mybigmerge <- myselectedmerge %>% na.omit() %>% select(DepMap_ID, gene, value) %>% filter(gene == i) %>% distinct()  %>% inner_join(mynewexp, by = "DepMap_ID")
      mybigmergeuniq <- distinct(mybigmerge, DepMap_ID, gene.x, value, TPM, gene.y)
      mybigene <- mybigmergeuniq %>% filter(gene.x == i) %>%  spread( gene.y, TPM, fill = 0) %>% select(-DepMap_ID, -gene.x)
      load(file = paste(paste(myoutdir, "/ranger_models_exp_Ginicorrected/model", sep = ""), i, "rds", sep = "."))
      inv <- try(importance_pvalues(results), silent=TRUE)
      if ('try-error' %in% class(inv)) next
      muts <- importance_pvalues(results) %>% as.data.frame() %>% subset(pvalue == 0) %>% rownames()
      for (c in muts) {
        coeff <- cor.test(mybigene$value, mybigene[[c]])$estimate
        pval <- cor.test(mybigene$value, mybigene[[c]])$p.value
        #if (coeff > 0.3 | coeff < -0.3) {
        res <- rbind(res, data.frame(i, c, coeff, pval, results$variable.importance[c]))
        #}
      }
    }
    colnames(res)[1] <- "gene"
    res$group <- "ec"
    colnames(res)[5] <- "importance"
    write_csv(res, path = paste(myoutdir, "/PairsExp_imp_corr.csv", sep = ""))

  }

}

