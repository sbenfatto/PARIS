#data documentation

#' Expression data.
#'
#' A dataset containing the gene expression data.
#'
#'
#' @format A data frame with 23164240 rows and 3 variables:
#' \describe{
#'   \item{DepMap_ID}{DepMap_ID}
#'   \item{gene}{gene name}
#'   \item{TPM}{TPM, Transcripts Per Kilobase Million}
#'   ...
#' }
#' @source \url{https://depmap.org/portal/}
"mynewexp"

#' Mutation data.
#'
#' A dataset containing the gene mutation data.
#'
#'
#' @format A data frame with 716951 rows and 2 variables:
#' \describe{
#'   \item{DepMap_ID}{DepMap_ID}
#'   \item{gene}{mutated gene name}
#'   ...
#' }
#' @source \url{https://depmap.org/portal/}
"mynewmut"

#' Copy number data.
#'
#' A dataset containing the gene copy number data.
#'
#'
#' @format A data frame with 1657 rows and 3 variables:
#' \describe{
#'   \item{DepMap_ID}{DepMap_ID}
#'   \item{gene}{gene name}
#'   \item{CN}{CN, Copy number}
#'   ...
#' }
#' @source \url{https://depmap.org/portal/}
"mynewcn"

#' Dependency data.
#'
#' A dataset containing the gene dependency data.
#'
#'
#' @format A data frame with 11458125 rows and 3 variables:
#' \describe{
#'   \item{DepMap_ID}{DepMap_ID}
#'   \item{gene}{gene name}
#'   \item{value}{value, gene dependency}
#'   ...
#' }
#' @source \url{https://depmap.org/portal/}
"mynewmerge"

#' Stats data.
#'
#' A dataset containing the gene dependency statistics data.
#'
#'
#' @format A data frame with 18333 rows and 6 variables:
#' \describe{
#'   \item{gene}{gene name}
#'   \item{med}{gene dependency median}
#'   \item{cv}{gene dependency, coefficient of variation}
#'   \item{sd}{gene dependency, standard deviation}
#'   \item{cvmod}{gene dependency, coefficient of variation modified}
#'   \item{r}{gene dependency, range}
#'   ...
#' }
#' @source \url{https://depmap.org/portal/}
"mystat"

#' Info data.
#'
#' A dataset containing the cell lines info data.
#'
#'
#' @format A data frame with 1736 rows and 23 variables:
#' \describe{
#'   \item{DepMap_ID}{DepMap_ID}
#'   \item{stripped_cell_line_name}{stripped_cell_line_name}
#'   \item{CCLE Name}{CCLE Name}
#'   \item{alias}{alias}
#'   \item{COSMIC_ID}{COSMIC_ID}
#'   \item{lineage}{lineage}
#'   \item{lineage_subtype}{lineage_subtype}
#'   \item{lineage_sub_subtype}{lineage_sub_subtype}
#'   \item{sex}{sex}
#'   \item{source}{source}
#'   \item{Achilles_n_replicates}{Achilles_n_replicates}
#'   \item{cell_line_NNMD}{cell_line_NNMD}
#'   \item{culture_type}{culture_type}
#'   \item{culture_medium}{culture_medium}
#'   \item{cas9_activity}{cas9_activity}
#'   \item{RRID}{RRID}
#'   \item{sample_collection_site}{sample_collection_site}
#'   \item{primary_or_metastasis}{primary_or_metastasis}
#'   \item{disease}{disease}
#'   \item{disease_subtype}{disease_subtype}
#'   \item{age}{age}
#'   \item{Sanger_model_ID}{Sanger_model_ID}
#'   \item{additional_info}{additional_info}
#'   
#'   ...
#' }
#' @source \url{https://depmap.org/portal/}
"myinfo"

