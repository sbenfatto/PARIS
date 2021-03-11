#Utilis

# library(Boruta)
# library(tidyverse)
# library(STRINGdb)

#' Scaling function
#' @param x vector of importance scores
#' @return The function return scaled importance scores between 0 and 1
#' @export
#' @example
#' \dontrun{
#' scale01(importance)
#' }
scale01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}

#' Load importance scores
#' @param path directory to Boruta models
#' @param group group to be added
#' @return Return dataframe importance score gene pairs
#' @export
#' @example
#' \dontrun{
#' load_imp_scores("/Boruta_models_exp_ImpRfGini/", "mr")
#' }
load_imp_scores <- function(path, group = c("mr","mg","er","eg")){

  group <- match.arg(group)

  dir <- path
  res <- tibble()
  mylist <- list.files(dir, pattern = "*.rds")

  for (i in (mylist)) {
    set.seed(43)
    #mybigene <-  filter(mybigmergeuniqsp, gene.x == i) %>% select(-DepMap_ID, -gene.x)
    load(file = paste(dir, i, sep = ""))
    muts <- subset(attStats(results), decision == "Confirmed")
    muts$feature <- rownames(muts)
    rownames(muts) <- NULL
    if (length(complete.cases(muts)) > 0 ){
      muts$group <- group
      muts$KO <- gsub(".rds", "", gsub("model.", "", i))
      res <- bind_rows(res,muts)
    }}

  colnames(res)[7] <- "c"
  colnames(res)[9] <- "gene"

  res <- select(res, c, gene, medianImp , group)

  return(res)
}

#' Load importance scores
#' @param path directory to Boruta models
#' @return Return dataframe importance score gene pairs
#' @export
#' @example
#' \dontrun{
#' load_imp_scores2("/Boruta_models_exp_ImpRfGini/")
#' }
load_imp_scores2 <- function(path){

  dir <- path
  res <- tibble()
  mylist <- list.files(dir, pattern = "*.rds")

  for (i in (mylist)) {
    set.seed(43)
    #mybigene <-  filter(mybigmergeuniqsp, gene.x == i) %>% select(-DepMap_ID, -gene.x)
    load(file = paste(dir, i, sep = ""))
    muts <- subset(attStats(results), decision == "Confirmed")
    muts$feature <- rownames(muts)
    rownames(muts) <- NULL
    if (length(complete.cases(muts)) > 0 ){
      muts$KO <- gsub(".rds", "", gsub("model.", "", i))
      res <- bind_rows(res,muts)
    }}

  colnames(res)[7] <- "c"
  colnames(res)[8] <- "gene"

  res <- select(res, c, gene, medianImp)
  colnames(res)[3] <- "importance"

  return(res)
}

#' Add group to Pairs file
#' @param pairfile File with selected genes pairs
#' @param group group to be added
#' @return Overwrite the file adding group
#' @export
#' @example
#' \dontrun{
#' add_group("/PairsMut_ImpRfGini.csv", "mg")
#' }
add_group <- function(pairfile, g = c("mr","mg","er","eg")){

  g <- match.arg(g)
  pairs <- read_csv(pairfile)
  pairs$group <- g

  write_csv(pairs, pairfile)

}

#' Merge selected gene pairs with their importance scores
#' @param impfile file with importance scores
#' @param pairsfile file with selected gene pairs
#' @return dataframe with selected gene pairs along with their importance scores
#' @export
#' @example
#' \dontrun{
#' add_imp_scores("/Imp.csv", "/PairsExp_ImpRfGini.csv")
#' }
add_imp_scores <- function(impfile, pairsfile){

  imp <- read_csv(impfile)
  pairs <- read_csv(pairsfile) %>% select(-1)
  colnames(pairs)[1] <- "gene"
  colnames(imp)[3] <- "importance"

  res <- inner_join(pairs, imp)

  return(res)
}

#' Merge selected gene pairs with their importance scores
#' @param path directory containing Selected Pairs files and Importance file
#' @return Overwrite Selected gene pairs with their importance scores
#' @export
#' @example
#' \dontrun{
#' add_imp_scores_loop("/Models/")
#' }
add_imp_scores_loop <- function(path){

  myfiles <- list.files(path, "Pairs.*Imp.*", full.names = T)

  for (i in myfiles) {
    pairs_imp <- add_imp_scores(paste(path, "Imp.csv", sep = ""), i)

    write_csv(pairs_imp, i, append = F)
  }
}

#' Filter essential genes
#' @param data dataframe with selected gene pairs
#' @param cutoff value to use as a threshold for the dependecy score coefficient of variantion
#' @return Dataframe with essential genes removed
#' @export
#' @example
#' \dontrun{
#' filter_essential_genes(mydata_DDR, 0.3)
#' }
filter_essential_genes <- function(data, cutoff){

  data <- filter(data, !(gene %in% filter(mystat, cv < cutoff)$gene) & !(c %in% filter(mystat, cv < cutoff)$gene)) #filter out essential genes

  return(data)
}

#' Define GenerateCol as inside Boruta
#' @noRd
generateCol<-function(x,colCode,col,numShadow){
  #Checking arguments
  if(is.null(col) & length(colCode)!=4)
    stop('colCode should have 4 elements.')
  #Generating col
  if(is.null(col)){
    rep(colCode[4],length(x$finalDecision)+numShadow)->cc
    cc[c(x$finalDecision=='Confirmed',rep(FALSE,numShadow))]<-colCode[1]
    cc[c(x$finalDecision=='Tentative',rep(FALSE,numShadow))]<-colCode[2]
    cc[c(x$finalDecision=='Rejected',rep(FALSE,numShadow))]<-colCode[3]
    col=cc
  }
  return(col)
}

#' plot Boruta top 10 features
#' @param x Boruta object
#' @return Plot top 10 features of a Boruta object
#' @export
#' @example
#' \dontrun{
#' plotboruta(results, las=2, cex.axis=0.5, main = i, xlab = "Mutations")
#' }
plotboruta<-function(x,colCode=c('green','yellow','red','blue'),sort=TRUE,whichShadow=c(TRUE,TRUE,TRUE),
                     col=NULL,xlab='Attributes',ylab='Importance',...){
  #Checking arguments
  if(class(x)!='Boruta')
    stop('This function needs Boruta object as an argument.')
  if(is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.')

  #Removal of -Infs and conversion to a list
  lz<-lapply(1:ncol(x$ImpHistory),function(i) x$ImpHistory[is.finite(x$ImpHistory[,i]),i])
  colnames(x$ImpHistory)->names(lz)

  #Selection of shadow meta-attributes
  numShadow<-sum(whichShadow)
  lz[c(rep(TRUE,length(x$finalDecision)),whichShadow)]->lz

  #Generating color vector
  col<-generateCol(x,colCode,col,numShadow)

  #Ordering boxes due to attribute median importance
  if(sort){
    ii<-order(sapply(lz,stats::median))
    lz[ii]->lz; col<-col[ii]
  }

  #Final plotting
  graphics::boxplot(tail(lz,50),xlab=xlab,ylab=ylab,col=tail(col,50),...)
  invisible(x)
}

#' Read selected gene pairs from a directory
#' @param path directory with selected gene pairs files
#' @param coeff.filter Boolean, if to apply the Pearson coefficient-based filtering step
#' @return dataframe with all selected gene pairs
#' @export
#' @example
#' \dontrun{
#' select_pairs("/Models/", coeff.filter = T)
#' }
select_pairs <- function(path, coeff.filter=TRUE){


  files <- list.files(path, "^Pairs*", full.names = T)

  mydata <- tibble()

  for (i in files) {
    print(paste("Reading......",i))
    current_tibble <- read_csv(i)
    mydata <- bind_rows(mydata, current_tibble)
  }

  mydata <- filter(mydata, gene != c) #remove self-pairs
  if (coeff.filter == T) {
    mydata <- filter(mydata, (grepl("e", group) & coeff < 0) | (grepl("m", group) & coeff > 0))
  }
  mydata <- mydata %>% group_by(group) %>% mutate(scaled = scale01(importance))

  return(mydata)
}

#' Read selected gene pairs from a directory
#' @param path directory with selected gene pairs files
#' @param coeff.filter Boolean, traditional Pearson coefficient-based filtering step or the other one
#' @return dataframe with all selected gene pairs
#' @export
#' @example
#' \dontrun{
#' select_pairs2("/Models/", coeff.filter = T)
#' }
select_pairs2 <- function(path, coeff.filter=TRUE){


  files <- list.files(path, "^Pairs*", full.names = T)

  mydata <- tibble()

  for (i in files) {
    print(paste("Reading......",i))
    current_tibble <- read_csv(i)
    mydata <- bind_rows(mydata, current_tibble)
  }

  mydata <- filter(mydata, gene != c) #remove self-pairs
  if (coeff.filter == T) {
    mydata <- filter(mydata, (grepl("e", group) & coeff < 0) | (grepl("m", group) & coeff > 0))
  }else{
    mydata <- filter(mydata, (grepl("e", group) & coeff > 0) | (grepl("m", group) & coeff < 0))
  }
  mydata <- mydata %>% group_by(group) %>% mutate(scaled = scale01(importance))

  return(mydata)
}

#' Retrieve StringDB interactions for gene pairs
#' @param dataframe Dataframe with a subset of gene pairs from \link{subset_by_threshold}
#' @return Dataframe with StringDB interactions (if any) of the selected gene pairs
#' @export
#' @example
#' \dontrun{
#' get_stringdb_interactions(mydata_DDR_sub)
#' }
get_stringdb_interactions <- function(data){

  # create a new STRING_db object ##NB till version 10 is supported in R
  string_db <- STRINGdb$new(version="10",species=9606)

  res <- tibble()
  for (i in unique(data$group)) {

    subdata <- filter(data, group == i)

    genes_mapped <- string_db$map(subdata %>%
                                        ungroup() %>%
                                        gather(mygenes, name, 1:2) %>%
                                        select(name) %>%
                                        distinct() %>%
                                        as.data.frame(), "name", removeUnmappedRows = TRUE )

    myinteractions <- string_db$get_interactions(genes_mapped$STRING_id)

    #add gene name to interaction table
    colnames(myinteractions)[1] <- "STRING_id"
    myinteractions <- inner_join(genes_mapped, myinteractions, by = "STRING_id")
    colnames(myinteractions)[c(2,3)] <- c("from", "STRING_id")
    myinteractions <- inner_join(genes_mapped, myinteractions, by = "STRING_id")
    colnames(myinteractions)[c(1,2,3,4)] <- c("node1", "to", "node2", "from")

    #duplicate row to have also the interactions in the opposite direction
    myinteractions2 <- myinteractions
    myinteractions2$node1 <- myinteractions$node2
    myinteractions2$node2 <- myinteractions$node1
    myinteractions <- bind_rows(myinteractions, myinteractions2)

    #filter interactions not present in tha data pairs
    myinteractions_filt <- filter(myinteractions, interaction(node1, node2) %in% interaction(subdata$gene, subdata$c))

    myinteractions_filt$group <- i

    res <- bind_rows(res, myinteractions_filt)
  }

  return(res)
}

#' Subset a dataframe of gene pairs by threshold
#' @param mydata Dataframe of selected gene pairs from \link{select_pairs}
#' @param g group: expression or mutation
#' @param threshold cutoff for the scaled importance score
#' @return subset dataframe
#' @export
#' @example
#' \dontrun{
#' subset_by_threshold(mydata_DDR, "e", 0.4)
#' }
subset_by_threshold <- function(data, g = c("e", "m"), threshold){

  g <- match.arg(g)

  high <- data %>%
            separate(group, c("source", "imp"), sep = 1) %>%
            group_by(gene, c, source) %>%
            mutate(n=n()) %>% filter(n ==3) %>%
            unite("group", source:imp, sep = "") %>%
            filter(grepl(g, group)) %>%
            filter(scaled > threshold)

   low <- data %>%
            separate(group, c("source", "imp"), sep = 1) %>%
            group_by(gene, c, source) %>%
            mutate(n=n()) %>% filter(n ==3) %>%
            unite("group", source:imp, sep = "") %>%
            filter(grepl(g, group)) %>% group_by(gene,c) %>% filter(all(scaled <= threshold)) %>%
            mutate(group = "low")

   res <- bind_rows(high, low)

   return(res)
}

#' Calculate fraction of interacting gene pairs over total
#' @param data Dataframe with a subset of gene pairs from \link{subset_by_threshold}
#' @param interactions Dataframe with StringDB interactions from \link{get_stringdb_interactions}
#' @return Dataframe with fraction of interacting gene pairs over total by group
#' @export
#' @example
#' \dontrun{
#' compute_interaction_freq(mydata_DDR_sub, mydata_DDR_sub_string)
#' }
compute_interaction_freq <- function(data, interactions){

  submerge <- data
  submerge$subgroup <- "not interaction"
  submerge <- submerge %>% ungroup() %>% select(group, subgroup)
  submerge <- submerge %>% group_by(group) %>% mutate(tot=n()) %>% distinct()

  bins2 <- interactions %>% select(group)
  bins2$subgroup <- "interaction"
  bins2 <- bins2 %>% group_by(group) %>% mutate(count=n()) %>% distinct()

  dd <- inner_join(submerge, bins2, by="group") %>%
    mutate(interacting=count/tot, not_interacting=(tot-count)/tot) %>%
    select(group, interacting, not_interacting) %>% gather(subgroup, frequency, 2:3)

  return(dd)
}

human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#' Retrieve Paralogs from gene vector
#' @param x vector of gene names
#' @return Dataframe of paralogs in the gene vector
#' @export
#' @example
#' \dontrun{
#' get_paralog_biomart(unique(c(mydata_DDR$gene)))
#' }
get_paralog_biomart <- function(x){

  res <- biomaRt::getBM(attributes = c("external_gene_name",
                           "hsapiens_paralog_associated_gene_name"),
            filters = "external_gene_name",
            values = x,
            mart = human)
  return(res)

}
