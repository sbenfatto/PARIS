#plot utils

# library(tidyverse)
# library(ggplot2)
# library(ggpubr)
# library(wesanderson)
# library(viridis)
# library(UpSetR)

#' Plot scatterplot dependecy vs expression
#' @param depgene Dependecy score gene name
#' @param expgene Expression gene name
#' @return Plot scatterplot dependecy score vs expression
#' @export
#' @example
#' \dontrun{
#' plot.sca.exp.dep(TYMS, TYMP)
#' }
plot.sca.exp.dep <- function(depgene, expgene) {

  depgene <- deparse(substitute(depgene))
  expgene <- deparse(substitute(expgene))

  exp <- filter(mynewexp, gene == expgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- inner_join(dep, exp)

  ggplot(merge, aes(TPM, value))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "pearson")+
    theme_light()+
    theme(text = element_text(size = 20), axis.text = element_text(size = 15, color = "black"))+
    xlab(paste(expgene, "expression"))+
    ylab(paste(depgene, "dependecy"))
}

#' Plot scatterplot dependecy vs ratio expression
#' @param depgene Dependecy score gene name
#' @param expgenes Expression gene name
#' @return Plot scatterplot dependecy score of one gene vs the ratio of expression of two genes
#' @export
#' @example
#' \dontrun{
#' plot.sca.exp.dep.ratio(TYMS)
#' }
plot.sca.exp.dep.ratio <- function(depgene, expgenes = c("TYMS", "TYMP")) {

  depgene <- deparse(substitute(depgene))

  exp <- filter(mynewexp, gene %in% expgenes) %>% spread(gene, TPM) %>% mutate(ratio = TYMP/TYMS)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- inner_join(dep, exp)
  print(expgenes)
  ggplot(merge, aes(ratio, value))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "pearson")+
    theme_light()+
    theme(text = element_text(size = 20), axis.text = element_text(size = 15, color = "black"))+
    xlab("TYMP/TYMS TMP")+
    ylab(paste(depgene, "dependecy"))
}

#' Plot scatterplot dependecy vs mutation
#' @param depgene Dependecy score gene name
#' @param mutgene Expression gene name
#' @return Plot scatterplot dependecy score vs mutation
#' @export
#' @example
#' \dontrun{
#' plot.sca.mut.dep(TYMS, CDKN2A)
#' }
plot.sca.mut.dep <- function(depgene, mutgene) {

  depgene <- deparse(substitute(depgene))
  mutgene <- deparse(substitute(mutgene))

  mut <- filter(mynewmut, gene == mutgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- dep %>% mutate(mut=ifelse(DepMap_ID %in% mut$DepMap_ID, 1, 0))

  ggplot(merge, aes(mut, value))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "pearson")+
    theme_light()+
    theme(text = element_text(size = 20), axis.text = element_text(size = 15, color = "black"))+
    xlab(paste(mutgene, "mutated"))+
    ylab(paste(depgene, "dependecy"))
}

#' Plot rank dependecy vs expression
#' @param depgene Dependecy score gene name
#' @param expgene Expression gene name
#' @return Plot ranked barplot dependecy score, colored by expression
#' @export
#' @example
#' \dontrun{
#' plot.rank.exp.dep(TYMS, TYMP)
#' }
plot.rank.exp.dep <- function(depgene, expgene) {

  depgene <- deparse(substitute(depgene))
  expgene <- deparse(substitute(expgene))

  exp <- filter(mynewexp, gene == expgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- inner_join(dep, exp) %>% na.omit()

  ggplot(merge, aes(reorder(DepMap_ID, value), value-mean(value), fill=TPM))+
    geom_bar(stat = "identity", width = 1)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15))+
    ylab(paste(depgene, " dependency mean centered"))+
    xlab("")+
    scale_fill_viridis(name = paste(expgene, " TPM"), direction = -1)
  #scale_fill_gradientn(colors = wes_palette("Zissou1", 10, type = "continuous"), name = paste(expgene, " TPM"))

}

#' Plot rank dependecy vs expression
#' @param depgene Dependecy score gene name
#' @param expgene Expression gene name
#' @param cells Cell lines DepMap to subset
#' @return Plot ranked barplot dependecy score, colored by expression, cell lines subset
#' @export
#' @example
#' \dontrun{
#' plot.rank.exp.dep.sub(TYMS, TYMP, c())
#' }
plot.rank.exp.dep.sub <- function(depgene, expgene, cells) {

  depgene <- deparse(substitute(depgene))
  expgene <- deparse(substitute(expgene))

  exp <- filter(mynewexp, gene == expgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- inner_join(dep, exp)
  merge <- filter(merge, DepMap_ID %in% cells) %>% inner_join(myinfo)

  ggplot(merge, aes(reorder(DepMap_ID, value), value-mean(value), fill=TPM))+
    geom_bar(stat = "identity", width = 1)+
    theme_minimal()+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15))+
    ylab(paste(depgene, " dependency mean centered"))+
    xlab("")+
    scale_fill_viridis(name = paste(expgene, " TPM"), direction = -1)+
    #scale_fill_gradientn(colors = wes_palette("Zissou1", 10, type = "continuous"), name = paste(expgene, " TPM"))+
    geom_text(aes(y=0,label=stripped_cell_line_name), size=5)

}

#' Plot boxplot dependecy vs mutation
#' @param depgene Dependecy score gene name
#' @param mutgene Expression gene name
#' @return Plot boxplot dependecy score, by mutational status
#' @export
#' @example
#' \dontrun{
#' plot.mut.dep(TYMS, CDKN2A)
#' }
plot.mut.dep <- function(depgene, mutgene) {

  depgene <- deparse(substitute(depgene))
  mutgene <- deparse(substitute(mutgene))

  mut <- filter(mynewmut, gene == mutgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- dep %>% mutate(mut=ifelse(DepMap_ID %in% mut$DepMap_ID, 1, 0))

  ggplot(merge, aes(mutgene, value, fill=factor(mut)))+
    geom_boxplot()+
    ylab(paste(depgene, "Dependency"))+
    theme_light()+
    xlab("")+
    theme(text = element_text(size = 20))+
    guides(fill=guide_legend(title="Mutated"))+
    stat_compare_means(method = "wilcox", label.y = 1.3, size=5)+
    scale_fill_manual(values = wes_palette("Zissou1", 2, type = "continuous"))

}

#' Plot rank dependecy vs mutation
#' @param depgene Dependecy score gene name
#' @param mutgene Expression gene name
#' @return Plot ranked barplot dependecy score, colored by mutational status
#' @export
#' @example
#' \dontrun{
#' plot.rank.mut.dep(TYMS, CDKN2A)
#' }
plot.rank.mut.dep <- function(depgene, mutgene) {

  depgene <- deparse(substitute(depgene))
  mutgene <- deparse(substitute(mutgene))

  mut <- filter(mynewmut, gene == mutgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- dep %>% mutate(mut=ifelse(DepMap_ID %in% mut$DepMap_ID, 1, 0))

  ggplot(merge, aes(reorder(DepMap_ID, value), value-mean(value), fill=as.factor(mut)))+
    geom_bar(stat = "identity", width = 1)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15))+
    ylab(paste(depgene, " dependency mean centered"))+
    xlab("")+
    guides(fill=guide_legend(title=paste(mutgene, " Mutated")))+
    scale_fill_manual(values = wes_palette("Zissou1", 2, type = "continuous"))

}

#' Plot rank dependecy vs mutation
#' @param depgene Dependecy score gene name
#' @param mutgene Expression gene name
#' @param cells Cell lines DepMap to subset
#' @return Plot ranked barplot dependecy score, colored by mutational status
#' @export
#' @example
#' \dontrun{
#' plot.rank.mut.dep.sub(TYMS, CDKN2A, c())
#' }
plot.rank.mut.dep.sub <- function(depgene, mutgene, cells) {

  depgene <- deparse(substitute(depgene))
  mutgene <- deparse(substitute(mutgene))

  mut <- filter(mynewmut, gene == mutgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- dep %>% mutate(mut=ifelse(DepMap_ID %in% mut$DepMap_ID, 1, 0))
  merge <- filter(merge, DepMap_ID %in% cells) %>% inner_join(myinfo)

  ggplot(merge, aes(reorder(DepMap_ID, value), value-mean(value), fill=as.factor(mut)))+
    geom_bar(stat = "identity", width = 1)+
    theme_minimal()+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15))+
    ylab(paste(depgene, " dependency mean centered"))+
    xlab("")+
    guides(fill=guide_legend(title=paste(mutgene, " Mutated")))+
    scale_fill_manual(values = wes_palette("Zissou1", 2, type = "continuous"))+
    geom_text(aes(y=0,label=stripped_cell_line_name), size=5)

}

#' Plot boxplot dependecy vs mutation
#' @param depgene Dependecy score gene name
#' @param mutgene Expression gene name
#' @return Plot boxplot dependecy score, by mutational status plus jitter
#' @export
#' @example
#' \dontrun{
#' plot.box.mut.dep(TYMS, CDKN2A)
#' }
plot.box.mut.dep <- function(depgene, mutgene) {

  depgene <- deparse(substitute(depgene))
  mutgene <- deparse(substitute(mutgene))

  mut <- filter(mynewmut, gene == mutgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- dep %>% mutate(mut=ifelse(DepMap_ID %in% mut$DepMap_ID, 1, 0))

  ggplot(merge, aes(as.factor(mut), value, fill=as.factor(mut)))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.3)+
    theme_light()+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20),
          axis.text.y = element_text(size = 15, color = "black"))+
    ylab(paste(depgene, " dependency score"))+
    xlab(mutgene)+
    guides(fill=guide_legend(title="Mutated"))+
    scale_fill_manual(values = wes_palette("Zissou1", 5)[c(1,5)])

}

#' Plot violin dependecy vs mutation
#' @param depgene Dependecy score gene name
#' @param mutgene Expression gene name
#' @return Plot violin dependecy score, by mutational status plus jitter
#' @export
#' @example
#' \dontrun{
#' plot.vl.mut.dep(TYMS, CDKN2A)
#' }
plot.vl.mut.dep <- function(depgene, mutgene) {

  depgene <- deparse(substitute(depgene))
  mutgene <- deparse(substitute(mutgene))

  mut <- filter(mynewmut, gene == mutgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- dep %>% mutate(mut=ifelse(DepMap_ID %in% mut$DepMap_ID, 1, 0))

  ggplot(merge, aes(as.factor(mut), value, fill=as.factor(mut)))+
    geom_violin()+
    geom_jitter(width = 0.3)+
    theme_light()+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20),
          axis.text.y = element_text(size = 15, color = "black"))+
    ylab(paste(depgene, " dependency score"))+
    xlab(mutgene)+
    guides(fill=guide_legend(title="Mutated"))+
    scale_fill_manual(values = wes_palette("Zissou1", 5)[c(1,5)])

}

#' Plot boxplot dependecy vs deficient
#' @param depgene Dependecy score gene name
#' @param defgene Expression gene name
#' @return Plot boxplot dependecy score, by deficiency status (mutated or copy number loss)
#' @export
#' @example
#' \dontrun{
#' plot.box.def.dep(TYMS, CDKN2A)
#' }
plot.box.def.dep <- function(depgene, defgene) {

  depgene <- deparse(substitute(depgene))
  defgene <- deparse(substitute(defgene))

  mut <- filter(mynewmut, gene == defgene) %>% select(-gene)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  cn  <- filter(mynewcn, gene == defgene) %>% select(-gene)
  #select cells with CN loss
  cn <- filter(cn, CN <= 0.01)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- dep %>% mutate(def=ifelse((DepMap_ID %in% mut$DepMap_ID
                                      | DepMap_ID %in% cn$DepMap_ID), 1, 0))

  ggplot(merge, aes(as.factor(def), value, fill=as.factor(def)))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter()+
    theme_light()+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20),
          axis.text.y = element_text(size = 15, color = "black"))+
    ylab(paste(depgene, " dependency score"))+
    xlab(mutgene)+
    guides(fill=guide_legend(title="Mutated"))+
    scale_fill_manual(values = wes_palette("Zissou1", 5)[c(1,5)])

}

#' Plot top rank selected genes pairs
#' @param data dataframe with selected pairs from \link{select_pairs}
#' @param g group to look for
#' @param cutoff cutoff of scaled importance score
#' @return Plot top rank selected gene pairs for a group with scaled importance score above the cutoff
#' @export
#' @example
#' \dontrun{
#' plot.top.rank(mydata_DDR, "mr", 0.5)
#' }
plot.top.rank <- function(data, g = c("mr","mg","mc","er","eg","ec"), cutoff){

  g <- match.arg(g)

  data <- filter_essential_genes(data, 0.3) #remove essential genes based on the cv of the dep distribution

  subdata <- filter(data, group == g & scaled > cutoff)

  ggplot(subdata, aes(fct_reorder(interaction(gene, c, sep = "/"), scaled), scaled))+
    geom_bar(stat="identity", fill=wes_palette("Zissou1", 2)[2])+
    coord_flip()+
    theme_minimal()+
    theme(text = element_text(size = 20), axis.text = element_text(colour = "black"))+
    xlab("")+
    ylab("scaled importance score")+
    labs(fill="")
}

#' Plot top rank selected genes pairs colored
#' @param data dataframe with selected pairs from \link{select_pairs}
#' @param g group to look for
#' @param color second group
#' @param cutoff cutoff of scaled importance score
#' @return Plot top rank selected gene pairs for a group with scaled importance score above the cutoff,
#'  colored by another scaled importance score group
#' @export
#' @example
#' \dontrun{
#' plot.top.rank.color(mydata_DDR, "mr", "mc", 0.5)
#' }
plot.top.rank.color <- function(data, g = c("mr","mg","mc","er","eg","ec"), color = c("mr","mg","mc","er","eg","ec"), cutoff){

  g <- match.arg(g)
  color <- match.arg(color)

  data <- filter_essential_genes(data, 0.3) #remove essential genes based on the cv of the dep distribution

  subdata <- filter(data, group == g & scaled > cutoff)
  subdata2 <- filter(data, group == color)

  left_join(subdata, subdata2, by = c("gene", "c")) %>%
  ggplot(aes(fct_reorder(interaction(gene, c, sep = "/"), scaled.x), scaled.x, fill=scaled.y))+
    geom_bar(stat="identity")+
    coord_flip()+
    theme_minimal()+
    theme(text = element_text(size = 20), axis.text = element_text(colour = "black"))+
    xlab("")+
    ylab("scaled importance score")+
    labs(fill="")+
    scale_fill_gradientn(colours =  wes_palette("Zissou1", 5, type = "continuous"), limits=c(0,1))
}

#' Compute and plot HTBreaks
#' @param data dataframe with selected pairs from \link{select_pairs}
#' @param by.group Boolean, if group by
#' @return Compute Head Tail Breaks and plot density distribution and breaks as vertical lines
#' @export
#' @example
#' \dontrun{
#' HTBreaks(mydata_DDR)
#' }
HTBreaks <- function(data, by.group=TRUE ){

  if (by.group==T) {
    p <- list()
    count = 0
    for (g in unique(data$group)) {

      subdata <- filter(data, group == g)
      myvector <- subdata$scaled
      print(g)

      i=1
      class <- vector()

      m <- mean(myvector)
      h <- myvector[myvector > m]

      while (length(h)/length(myvector) < 0.4 & length(h)/length(myvector) != 0) {

        class[i] = m
        myvector <- h

        m <- mean(myvector)
        h <- myvector[myvector > m]
        t <- myvector[myvector < m]
        print(length(h)/length(myvector))

        i = i + 1

      }

      print((class))
      #plot breaks
      count = count +1

      if(length(class) == 0){
      p[[count]] <- (ggplot(subdata, aes(scaled))+geom_density()+
            theme_light()+theme(text = element_text(size = 20))+
            xlab("scaled importance scores")+
            geom_vline(xintercept = 0)+labs(title = g))
      }else{
        p[[count]] <- (ggplot(subdata, aes(scaled))+geom_density()+
                         theme_light()+theme(text = element_text(size = 20))+
                         xlab("scaled importance scores")+
                         geom_vline(xintercept = class)+labs(title = g))
      }
    }
    gridExtra::grid.arrange(grobs=p)
    }else{

        myvector <- data$scaled

        i=1
        class <- vector()

        m <- mean(myvector)
        h <- myvector[myvector > m]

        while (length(h)/length(myvector) < 0.4 & length(h)/length(myvector) != 0) {

          class[i] = m
          myvector <- h

          m <- mean(myvector)
          h <- myvector[myvector > m]
          t <- myvector[myvector < m]
          print(length(h)/length(myvector))

          i = i + 1

        }

        print(class)
        #plot breaks

        ggplot(data, aes(scaled))+geom_density()+
          theme_light()+theme(text = element_text(size = 20))+
          xlab("scaled importance scores")+
          geom_vline(xintercept = class)
    }
}

#' Compute and plot HTBreaks + histogram
#' @param data dataframe with selected pairs from \link{select_pairs}
#' @param by.group Boolean, if group by
#' @return Compute Head Tail Breaks and plot density distribution, histogram and breaks as vertical lines
#' @export
#' @example
#' \dontrun{
#' HTBreaks.hist(mydata_DDR)
#' }
HTBreaks.hist <- function(data, by.group=TRUE ){

  if (by.group==T) {
    p <- list()
    count = 0
    for (g in unique(data$group)) {

      subdata <- filter(data, group == g)
      myvector <- subdata$scaled
      print(g)

      i=1
      class <- vector()

      m <- mean(myvector)
      h <- myvector[myvector > m]

      while (length(h)/length(myvector) < 0.4 & length(h)/length(myvector) != 0) {

        class[i] = m
        myvector <- h

        m <- mean(myvector)
        h <- myvector[myvector > m]
        t <- myvector[myvector < m]
        print(length(h)/length(myvector))

        i = i + 1

      }

      print((class))
      #plot breaks
      count = count +1

      if(length(class) == 0){
        p[[count]] <- (ggplot(subdata, aes(scaled))+geom_histogram(aes(y=..density..))+geom_density()+
                         theme_light()+theme(text = element_text(size = 20))+
                         xlab("scaled importance scores")+
                         geom_vline(xintercept = 0)+labs(title = g))
      }else{
        p[[count]] <- (ggplot(subdata, aes(scaled))+geom_histogram(aes(y=..density..))+geom_density()+
                         theme_light()+theme(text = element_text(size = 20))+
                         xlab("scaled importance scores")+
                         geom_vline(xintercept = class)+labs(title = g))
      }
    }
    gridExtra::grid.arrange(grobs=p)
  }else{

    myvector <- data$scaled

    i=1
    class <- vector()

    m <- mean(myvector)
    h <- myvector[myvector > m]

    while (length(h)/length(myvector) < 0.4 & length(h)/length(myvector) != 0) {

      class[i] = m
      myvector <- h

      m <- mean(myvector)
      h <- myvector[myvector > m]
      t <- myvector[myvector < m]
      print(length(h)/length(myvector))

      i = i + 1

    }

    print(class)
    #plot breaks

    ggplot(data, aes(scaled))+geom_histogram(aes(y=..density..))+geom_density()+
      theme_light()+theme(text = element_text(size = 20))+
      xlab("scaled importance scores")+
      geom_vline(xintercept = class)
  }
}

#' Plot StringDB combined scores
#' @param interactions dataframe with interaction for each gene pairs from \link{get_stringdb_interactions}
#' @return Plot StringDB combined scores of interactiong gene pairs
#' @export
#' @example
#' \dontrun{
#' plot.stringdb.interactions.combined.score.by.group(mydata_DDR_sub_string)
#' }
plot.stringdb.interactions.combined.score.by.group <- function(interactions){

  ggplot(interactions, aes(group, combined_score/1000, fill=group))+
    geom_boxplot(outlier.shape = NA)+geom_jitter(width = 0.2)+
    scale_fill_manual(values = c(wes_palette("Zissou1", 5)[c(5,3,1)], "grey"))+
    theme_light()+
    theme(text = element_text(size = 25), axis.text = element_text(colour = "black"))+
    ylab("Combined score")
}

#' Plot percentage of interacting gene over total
#' @param interaction_freq dataframe from \link{compute_interaction_freq}
#' @return Plot percentage of interacting gene over total by group
#' @export
#' @example
#' \dontrun{
#' plot.stringdb.interactions.by.group(mydata_DDR_sub_string_freq)
#' }
plot.stringdb.interactions.by.group <- function(interaction_freq){

  ggplot(interaction_freq, aes(x = group, y = frequency, fill=subgroup))+
    geom_bar(stat = "identity", color="black")+
    scale_fill_manual(values = c("black", "grey95"))+
    theme_light()+
    theme(text = element_text(size = 25), axis.text = element_text(colour = "black"))+
    xlab("group")+
    labs(fill="subgroup")
}

#' Plot scatterplot Gini corrected vs raw permutation for a group
#' @param data dataframe with selected pairs from \link{select_pairs}
#' @param g group: expression or mutation
#' @param imp_gini gini or corrected gini
#' @return Plot scatterplot Gini corrected vs raw permutation for expression or mutation data
#' @export
#' @example
#' \dontrun{
#' plot.common.pairs(mydata_DDR, "m")
#' }
plot.common.pairs <- function(data, g = c("e", "m"), imp_gini = c("g", "c")){

  g <- match.arg(g, c("e", "m"))
  imp_gini <- match.arg(imp_gini)

  subdata <- data %>% filter(grepl(g, group)) %>%
              select(gene, c, scaled, group) %>%
              spread(group, scaled) %>%
              na.omit()


  if (g == "e") {
    if (imp_gini == "c") {

     ggplot(subdata, aes(er, ec))+
        geom_point(size=5)+
        theme_light()+
        theme(text = element_text(size = 20), axis.text = element_text(colour = "black"))+
        stat_cor(method = "pearson", size = 8)
    }else{

      ggplot(subdata, aes(er, eg))+
        geom_point(size=5)+
        theme_light()+
        theme(text = element_text(size = 20), axis.text = element_text(colour = "black"))+
        stat_cor(method = "pearson", size = 8)
    }

  }else{

    if (imp_gini == "c") {

      ggplot(subdata, aes(mr, mc))+
        geom_point(size=5)+
        theme_light()+
        theme(text = element_text(size = 20), axis.text = element_text(colour = "black"))+
        stat_cor(method = "pearson", size = 8)
    }else{

      ggplot(subdata, aes(mr, mg))+
        geom_point(size=5)+
        theme_light()+
        theme(text = element_text(size = 20), axis.text = element_text(colour = "black"))+
        stat_cor(method = "pearson", size = 8)
    }
  }
}

#' Plot ranked barplot dependecy score, colored by expression and red dots in deficient cells
#' @param depgene Gene dependecy score
#' @param expgene Gene expression
#' @param defgene Gene deficient (mutated or copy number loss)
#' @return Plot ranked barplot dependecy score, colored by expression and red dots in deficient cells
#' @export
#' @example
#' \dontrun{
#' plot.rank.def.exp.dep(TYMS, TYMP, CDKN2A)
#' }
plot.rank.def.exp.dep <- function(depgene, expgene, defgene) {

  depgene <- deparse(substitute(depgene))
  defgene <- deparse(substitute(defgene))
  expgene <- deparse(substitute(expgene))

  exp <- filter(mynewexp, gene == expgene) %>% select(-gene)
  mut <- filter(mynewmut, gene == defgene) %>% select(-gene)
  cn  <- filter(mynewcn, gene == defgene) %>% select(-gene)
  #select cells with CN loss
  cn <- filter(cn, CN <= 0.01)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- dep %>% mutate(def=ifelse((DepMap_ID %in% mut$DepMap_ID
                                      | DepMap_ID %in% cn$DepMap_ID), 1, 0))

  merge <- inner_join(merge, exp) %>% inner_join(myinfo)

  ggplot(merge, aes(reorder(DepMap_ID, value), value-mean(value), fill=TPM))+
    geom_bar(stat = "identity", width = 1)+
    geom_point(aes(color=as.factor(def)))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15))+
    ylab(paste(depgene, " dependency mean centered"))+
    xlab("")+
    scale_fill_viridis(name = paste(expgene, " TPM"), direction = -1)+
    scale_color_manual(name = paste(defgene, " deficient"), values = wes_palette("Zissou1", 2, type = "continuous"))+
    facet_wrap(~disease, scales = "free")

}

#' Plot ranked barplot dependecy score, colored by ratio expression of two genes and red dots in deficient cells by disease
#' @param depgene Gene dependecy score
#' @param expgenes Genes ratio expression
#' @param defgene Gene deficient (mutated or copy number loss)
#' @return Plot ranked barplot dependecy score, colored by expression and red dots in deficient cells bydisease
#' @export
#' @example
#' \dontrun{
#' plot.rank.def.ratio.dep.by.disease(TYMS, defgene=CDKN2A)
#' }
plot.rank.def.ratio.dep.by.disease <- function(depgene, expgenes = c("TYMS", "TYMP"), defgene) {

  depgene <- deparse(substitute(depgene))
  defgene <- deparse(substitute(defgene))

  exp <- filter(mynewexp, gene %in% expgenes) %>% spread(gene, TPM) %>% mutate(ratio = TYMP/TYMS)
  mut <- filter(mynewmut, gene == defgene) %>% select(-gene)
  cn  <- filter(mynewcn, gene == defgene) %>% select(-gene)
  #select cells with CN loss
  cn <- filter(cn, CN <= 0.01)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- dep %>% mutate(def=ifelse((DepMap_ID %in% mut$DepMap_ID
                                      | DepMap_ID %in% cn$DepMap_ID), 1, 0))

  merge <- inner_join(merge, exp) %>% inner_join(myinfo)

  ggplot(merge, aes(reorder(DepMap_ID, value), value-mean(value), fill=ratio))+
    geom_bar(stat = "identity", width = 1)+
    geom_point(aes(color=as.factor(def)))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15))+
    ylab(paste(depgene, " dependency mean centered"))+
    xlab("")+
    scale_fill_viridis(name = "TYMP/TYMS TMP", direction = -1)+
    scale_color_manual(name = paste(defgene, " deficient"), values = wes_palette("Zissou1", 2, type = "continuous"))+
    facet_wrap(~disease, scales = "free")

}

#' Plot ranked barplot dependecy score, colored by ratio expression of two genes
#' @param depgene Gene dependecy score
#' @param expgenes Genes ratio expression
#' @return Plot ranked barplot dependecy score, colored by expression
#' @export
#' @example
#' \dontrun{
#' plot.rank.def.ratio.dep.by.disease(TYMS)
#' }
plot.rank.def.ratio.dep <- function(depgene, expgenes = c("TYMS", "TYMP")) {

  depgene <- deparse(substitute(depgene))
  defgene <- deparse(substitute(defgene))

  exp <- filter(mynewexp, gene %in% expgenes) %>% spread(gene, TPM) %>% mutate(ratio = TYMP/TYMS)
  dep <- filter(mynewmerge, gene == depgene) %>% select(-gene)
  merge <- inner_join(dep, exp) %>% na.omit()

  ggplot(merge, aes(reorder(DepMap_ID, value), value-mean(value), fill=ratio))+
    geom_bar(stat = "identity", width = 1)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15))+
    ylab(paste(depgene, " dependency mean centered"))+
    xlab("")+
    scale_fill_viridis(name = "TYMP/TYMS TMP", direction = -1)

}

#' Upsetplot
#' @param data dataframe with selected pairs from \link{select_pairs}
#' @param g group: expression or mutationa
#' @param cutoff cutoff of scaled importance score
#' @return Subset the dataframe based on the cutoff and plot groups intersections as upsetplot
#' @export
#' @example
#' \dontrun{
#' plot.upset(mydata_DDR, "e", 0.4)
#' }
plot.upset <- function(data, g = c("e", "m"), cutoff){

  subdata <- filter(data, grepl(g, group)) %>% dplyr::group_by(gene, c) %>%  add_tally() %>% filter(n == 3) %>% mutate(group=ifelse(scaled >= cutoff, paste(group, "high"), paste(group, "low")))

  subdata %>% select(gene, c, group) %>% mutate(value=1) %>% spread(group, value, fill = 0)  %>% as.data.frame() %>%
    upset(nsets = 6, text.scale =  c(2, 2, 2, 1.5, 2, 1.5))
}
