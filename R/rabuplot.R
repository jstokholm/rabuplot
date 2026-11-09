#' Relative Abundance Plot
#'
#' For creating nice microbiome plots
#'
#' @param phylo_ob Phyloseq object with metadata in sample_data.
#' @param predictor Predictor of interest for statistics/plotting in sample_data.
#' @param covariates Character vector of covariate column names in sample_data to adjust for; default is NULL.
#' @param type Taxonomic rank from tax_table, case insensitive; default is "genus".
#' @param relative_abun Use relative abundances, else absolute; default is TRUE.
#' @param min_abund Minimum relative abundance threshold (e.g., 0.001 for 0.1\%); default is 0.
#' @param min_prev Minimum prevalence fraction threshold (fraction of samples where taxon exceeds min_abund, e.g., 0.1 for 10\% of samples); default is 0.
#' @param id Define id variable for mixed models.
#' @param xlabs X-axis label
#' @param ylabs Y-axis label
#' @param main Title of plot
#' @param violin Use geom_violin for plotting, else boxplot; default is TRUE.
#' @param violin_scale Scale option for geom_violin; default is "width".
#' @param legend_title Legend title; default is name of predictor.
#' @param N_taxa Number of taxa to be plotted; default is 15.
#' @param By_median Order plot by median abundances, else mean abundances; default is TRUE.
#' @param no_other_type Taxa in lower abundances than top N_taxa, are grouped as "other", this will remove this group from the plot; default is FALSE.
#' @param legend_names Define variable names for legend text.
#' @param Time Time variable name for longitudinal datasets.
#' @param Timepoint Value in variable \code{Time} to select.
#' @param Strata Name of variable for stratification; default is NULL.
#' @param Strata_val Value in variable \code{Strata} to keep; default is 1.
#' @param no_legends Removes legend; default is FALSE.
#' @param no_names Removes taxa names; default is FALSE.
#' @param italic_names Taxa names will be in italic e.g. usable for family, genus, species levels; default is TRUE
#' @param Only_sig Only keep significant taxa; default is FALSE.
#' @param log Present plot on a log scale; default is TRUE.
#' @param log_max Maximum value of log-axis options:1, 10, 100; default is 100.
#' @param stat_out Outputs a data.frame with statistics to Global environment; default is FALSE.
#' @param p_val Displays p-values on plot; default is TRUE.
#' @param p_stars Shows stars instead of p-values; default is FALSE.
#' @param stats Select type of statistical test; options: "non-parametric", "parametric", "mixed", "mgs_feature", "maaslin2", "maaslin3", "ancombc2", "deseq2"; default is "non-parametric".
#' @param p_adjust adjust p-values; default is FALSE.
#' @param p_adjust_method options: "holm", "hochberg", "hommel", "bonferroni", "bh", "by", "fdr"; default is "fdr".
#' @param p_adjust_full correction applied for all taxa in the dataset; default is FALSE.
#' @param colors define list of colors for plot. If not color brewer will be used; default is NULL.
#' @param color_by define taxonomic rank to color by; default is NULL.
#' @param order Order by abundance, else alphabetically; default is TRUE.
#' @param reverse Flip taxa order; default is FALSE.
#' @param list_taxa A list of specific taxa names to be analyzed; default is NULL.
#' @param select_taxa Choose all taxa from one or more taxonomic variables, eg. "Staphylococcus" or "Staph" or "coccus" or c("staph",bifido"); default is NULL.
#' @param select_type Taxonomic rank of the \code{select_taxa}; default is "genus".
#' @param bar_chart Choose to make bar chart; default is FALSE.
#' @param bar_chart_stacked Produce stacked bar chart; default is FALSE
#' @param percent Print percentages on bar chart; default is FALSE.
#' @param facet_wrap Facet wrap chart by variable; eg. Time; default is NULL.
#' @param facet_label Facet wrap labels; default is NULL.
#' @param facet_n Show n for each facet; default is TRUE.
#' @param order_by Choose variable to order the selected taxa by; eg. Time; default is Time.
#' @param order_val Choose value for \code{order_by}; default is NULL.
#' @param text_angle_x Choose value for rotation of axis-text; default is 0.
#' @param rank Taxonomic rank from tax_table, replaces "type", case insensitive; default is "genus".
#' @param select_rank Taxonomic rank of the \code{select_taxa}, replaces "select_type"; default is "genus".
#' @param volcano_plot Produce volcano plot (log2FC vs -log10 P-value); default is FALSE.
#' @param p_cutoff P-value significance threshold cutoff line for volcano plot; default is 0.05.
#' @param adj_p_cutoff q-value threshold cutoff line for volcano plot (used when p_adjust=TRUE); default is 0.05.
#' @param label_taxa Display bacterial names on significant points in volcano plot; default is TRUE.
#' @param top_n_labels Max number of significant bacterial labels to plot on volcano plot; default is 15.
#'
#' @import ggplot2 phyloseq dplyr tidyr RColorBrewer lmerTest grid
#' @importFrom ggrepel geom_text_repel
#' @return A ggplot
#' @export

rabuplot <- function(phylo_ob,
                     predictor="none",
                     covariates=NULL,
                     type="genus",
                     relative_abun=TRUE,
                     min_abund=0,
                     min_prev=0,
                     id=NULL,
                     xlabs = "Relative abundance (%)",
                     ylabs = "Average relative abundance",
                     main = "Relative abundance plot",
                     violin=TRUE,
                     violin_scale = "width",
                     legend_title=predictor,
                     N_taxa=NULL,
                     By_median=TRUE,
                     no_other_type=FALSE,
                     legend_names=NULL,
                     Time="Time",
                     Timepoint=NULL,
                     Strata=NULL,
                     Strata_val="1",
                     no_legends = FALSE,
                     no_names=FALSE,
                     italic_names=TRUE,
                     Only_sig=FALSE,
                     log=TRUE,
                     log_max=100,
                     stat_out=FALSE,
                     p_val = TRUE,
                     p_stars=FALSE,
                     stats="non-parametric",
                     p_adjust=FALSE,
                     p_adjust_method="fdr",
                     p_adjust_full=FALSE,
                     colors=NULL,
                     color_by=NULL,
                     order=TRUE,
                     reverse=FALSE,
                     list_taxa=NULL,
                     select_taxa=NULL,
                     select_type="genus",
                     bar_chart=FALSE,
                     bar_chart_stacked=FALSE,
                     facet_wrap=NULL,
                     facet_label=NULL,
                     facet_n=TRUE,
                     percent=FALSE,
                     order_by="Time",
                     order_val=NULL,
                     text_angle_x=0,
                     rank=NULL,
                     select_rank=NULL,
                     volcano_plot=FALSE,
                     p_cutoff=0.05,
                     adj_p_cutoff=0.05,
                     label_taxa=TRUE,
                     top_n_labels=30
)
{
  if (!is.numeric(min_abund) || min_abund < 0 || min_abund > 1) {
    stop("`min_abund` must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(min_prev) || min_prev < 0 || min_prev > 1) {
    stop("`min_prev` must be a numeric value between 0 and 1.")
  }

  if(!is.null(list_taxa) & is.null(N_taxa)) N_taxa = length(list_taxa)
  if(is.null(N_taxa) & is.null(list_taxa)) N_taxa=15
  if(!is.null(rank)) type <- rank
  if(!is.null(select_rank)) select_type <- select_rank
  options(dplyr.summarize.inform = FALSE)
  if(volcano_plot == TRUE) {
    p_val = TRUE
    bar_chart = FALSE
    bar_chart_stacked = FALSE
    if (main == "Relative abundance plot") main <- ""
    if (xlabs == "Relative abundance (%)") xlabs <- "log2FC"
  }
  if(bar_chart_stacked==TRUE) {
    bar_chart=TRUE
    p_val=FALSE
  }
  if(predictor=="none") {
    sample_data(phylo_ob)$none <- "All samples"
    p_val=FALSE
    if(bar_chart_stacked==FALSE & is.null(color_by)) no_legends = TRUE
  }
  if (!is.null(Strata)) {
    strata_vals <- phyloseq::get_variable(phylo_ob, Strata)
    keep_samples <- !is.na(strata_vals) & as.character(strata_vals) == as.character(Strata_val)

    if (sum(keep_samples) == 0) {
      stop(sprintf("No samples remaining after filtering by Strata '%s' == '%s'.", Strata, Strata_val))
    }
    phylo_ob <- phyloseq::prune_samples(keep_samples, phylo_ob)
  }

  otu_raw <- as(otu_table(phylo_ob), "matrix")
  if (taxa_are_rows(phylo_ob)) otu_raw <- t(otu_raw)

  has_data <- apply(otu_raw, 1, function(row) any(!is.na(row) & row != 0))
  if (all(!has_data)) {
    stop("[rabuplot Error]: All samples in the dataset are empty (all 0 or NA).", call. = FALSE)
  }
  if (any(!has_data)) {
    message(sprintf("%d empty samples removed from dataset.", sum(!has_data)))
    phylo_ob <- phyloseq::prune_samples(names(has_data)[has_data], phylo_ob)
  }

  if (!is.null(covariates)) {
    missing_covs <- setdiff(covariates, phyloseq::sample_variables(phylo_ob))
    if (length(missing_covs) > 0) {
      stop(sprintf("Covariates not found in sample_data: %s", paste(missing_covs, collapse = ", ")), call. = FALSE)
    }
  }

  otu_mat <- as(otu_table(phylo_ob), "matrix")
  if(taxa_are_rows(phylo_ob)) otu_mat <- t(otu_mat)

  if(!is.null(facet_wrap)) {
    index <- !is.na(get_variable(phylo_ob, predictor)) & !is.na(get_variable(phylo_ob, facet_wrap)) & as.character(get_variable(phylo_ob, facet_wrap)) != "NA"
  } else {
    index <- !is.na(get_variable(phylo_ob, predictor))
  }

  if (!is.null(covariates)) {
    cov_mat <- as.data.frame(phyloseq::get_variable(phylo_ob, covariates))
    cov_na  <- apply(cov_mat, 1, function(row) any(is.na(row) | as.character(row) == "NA"))
    index   <- index & !cov_na
  }

  if(length(unique(index)) != 1) message(paste(length(which(index==F)), "samples have been removed from full dataset (predictor/facet_wrap/covariate NAs)"))
  otu_mat <- otu_mat[index, , drop = FALSE]

  keep_feats <- apply(otu_mat, 2, function(col) any(!is.na(col) & col != 0))
  otu_mat    <- otu_mat[, keep_feats, drop = FALSE]

  OTU_index <- colnames(otu_mat)
  tax <- as(tax_table(phylo_ob), "matrix") %>% data.frame(stringsAsFactors=FALSE)
  tax <- tax[rownames(tax) %in% OTU_index,,drop=FALSE]
  tax[is.na(tax)] <- "unclassified"
  tax[tax==""] <- "unclassified"
  names(tax) <- tolower(names(tax))
  type <- tolower(type)
  if(!is.null(select_type)) select_type <- tolower(select_type)
  tax$OTU <- rownames(tax)
  samp <- data.frame(sample_data(phylo_ob), stringsAsFactors=TRUE)
  samp <- samp[index,]
  if(is.null(facet_wrap)) samp$wrap <- ""
  if(!is.null(facet_wrap)) samp$wrap <- samp[,facet_wrap]
  if(!is.null(Timepoint)){
    index <- rownames(samp[(samp[,Time] ==Timepoint),])
    otu_mat <- otu_mat[rownames(otu_mat) %in% index,]
    otu_mat  <- otu_mat[, colSums(otu_mat != 0, na.rm = TRUE) > 0, drop = FALSE]
    OTU_index <- colnames(otu_mat)
    tax <- tax[rownames(tax) %in% OTU_index,]
    samp <- samp[rownames(samp) %in% index,]
  }
  taxa_vec <- as.character(tax[,type])
  unique_tax <- unique(taxa_vec)

  abund <- as.data.frame(matrix(rep(0,(length(unique_tax)*nrow(otu_mat))),ncol=length(unique_tax)))
  row.names(abund) <- row.names(otu_mat)
  names(abund) <- unique_tax
  for(i in names(abund)){
    if(is.array(otu_mat[,taxa_vec==i]))  abund[,i] <- rowSums(otu_mat[,taxa_vec== i], na.rm = TRUE)
    else   abund[,i] <- otu_mat[,taxa_vec== i]
  }
  abund_org <- abund
  if(relative_abun==TRUE) abund <- apply(abund,1,function(x) x/sum(x)) %>% t %>% as.data.frame()

  if (is.null(facet_wrap) && (min_abund > 0 || min_prev > 0)) {
    n_taxa_before <- ncol(abund)
    rel_matrix <- if (relative_abun) abund else apply(abund_org, 1, function(x) x / sum(x)) %>% t %>% as.data.frame()

    n_samples <- nrow(rel_matrix)
    prevalence_frac <- colSums(rel_matrix > min_abund, na.rm = TRUE) / n_samples
    keep_taxa_filter <- prevalence_frac >= min_prev

    abund <- abund[, keep_taxa_filter, drop = FALSE]
    abund_org <- abund_org[, keep_taxa_filter, drop = FALSE]
    unique_tax <- names(abund)

    n_taxa_after <- ncol(abund)
    n_removed <- n_taxa_before - n_taxa_after

    message(sprintf(
      "Global filter applied (min_abund = %g, min_prev = %.1f%%). Removed %d / %d taxa (%d taxa remaining).",
      min_abund,
      min_prev * 100,
      n_removed,
      n_taxa_before,
      n_taxa_after
    ))

    if (n_taxa_after == 0) {
      stop("All taxa were removed by the specified `min_abund` and `min_prev` thresholds. Please relax your criteria.")
    }
  }

  abund_all <- abund
  if (is.null(list_taxa) & !is.null(select_taxa)) {
    list_taxa <- NULL
    for(i in 1:length(select_taxa)){
      list_taxa <- c(list_taxa,(as.character(unique(tax[grep(select_taxa[[i]],tax[,select_type],ignore.case=TRUE),type]))))
    }
  }
  if (!is.null(list_taxa)) {
    if (is.null(N_taxa)) N_taxa <- length(list_taxa)
    abund <- abund[,colnames(abund) %in% list_taxa, drop = FALSE]
    unique_tax <- names(abund)
  }

  if(length(abund)>1 && !volcano_plot){
    index <- !is.na(rownames(samp))
    if (!is.null(order_val))  index <- samp[,order_by] ==order_val
    abund <- abund[,order(-colSums(abund[index, , drop = FALSE], na.rm = TRUE))]
    if (By_median)  abund <- abund[,order(-apply(abund[index, , drop = FALSE], 2, stats::median, na.rm = TRUE))]
    if("unclassified" %in% unique_tax) abund <- abund[c(setdiff(names(abund), "unclassified"),"unclassified")] #Move unclassified to end
    if(N_taxa<length(unique_tax)) abund <- abund[-(length(unique_tax)-(length(unique_tax)-N_taxa)+1):-length(unique_tax)]
    if(no_other_type==FALSE) abund[, paste("Other",type)] <- rowSums(abund_all[,!names(abund_all) %in% names(abund), drop = FALSE], na.rm = TRUE)
  }
  index <- !is.na(rownames(samp))
  samp2 <- samp[index, , drop = FALSE]
  if(p_val==TRUE & (bar_chart==FALSE | (bar_chart==TRUE & bar_chart_stacked==FALSE))){
    if(p_adjust_full ==TRUE | stats=="mgs_feature" | stats=="ancombc2" | stats=="deseq2" | volcano_plot==TRUE){
      abund2 <- abund_org[index, , drop = FALSE]
      if(relative_abun==TRUE & stats!="mgs_feature" & stats!="ancombc2" & stats!="deseq2") abund2 <- apply(abund2,1,function(x) x/sum(x)) %>% t %>% as.data.frame()
    }
    else abund2 <- abund[index, , drop = FALSE]

    if(stats=="mgs_feature" & length(levels(factor(samp2[,predictor])))>2){
      stats="non-parametric"
      message("MGS not available for >2 predictors, switching to non-parametric")
    }
    if(stats=="mixed" & is.null(id)){
      stats="non-parametric"
      message("No id variable for mixed model, switching to non-parametric")
    }
    if(stats == "mixed"){
      if (!requireNamespace("lmerTest", quietly = TRUE)) {
        stop("Package 'lmerTest' is required for mixed models. Install via install.packages('lmerTest').", call. = FALSE)
      }
      if (!id %in% colnames(samp2)) {
        stop(sprintf("id variable '%s' not found in sample_data.", id), call. = FALSE)
      }
      message("Mixed model statistics")

      pred_vec <- samp2[, predictor]
      id_vec   <- samp2[, id]

      has_facet <- !is.null(facet_wrap)
      if(has_facet) {
        message(paste0("Using ", facet_wrap, " as mixed model group variable"))
        facet_vec <- samp2[, "wrap"]
      } else {
        message("No group variable defined for mixed model in facet_wrap")
        facet_vec <- NULL
      }

      taxa_names <- colnames(abund2)
      pvals <- rep(NA_real_, length(taxa_names))

      for(j in seq_along(taxa_names)) {
        y_val <- abund2[, j]
        if(all(is.na(y_val)) || stats::var(y_val, na.rm = TRUE) == 0) {
          next
        }

        df_fit <- data.frame(
          y    = y_val,
          pred = pred_vec,
          id   = id_vec,
          stringsAsFactors = FALSE
        )
        if(has_facet) df_fit$facet <- factor(facet_vec)
        if(!is.null(covariates)) df_fit <- cbind(df_fit, samp2[, covariates, drop = FALSE])

        fml_terms <- c("pred", if(has_facet) "facet", covariates)
        fml <- stats::as.formula(paste("y ~", paste(fml_terms, collapse = " + "), "+ (1 | id)"))

        pvals[j] <- tryCatch({
          fit <- lmerTest::lmer(fml, data = df_fit)
          anv <- as.data.frame(stats::anova(fit))
          row_idx <- which(rownames(anv) == "pred")
          if(length(row_idx) == 0) row_idx <- 1
          as.numeric(anv[row_idx, "Pr(>F)"])
        }, error = function(e) NA_real_)
      }

      pval <- data.frame(variable = taxa_names, pval = pvals, stringsAsFactors = FALSE)
      pval <- pval %>% dplyr::mutate(wrap = "Mixed", p_adjust = stats::p.adjust(pval, method = p_adjust_method))

      if(has_facet) {
        pval$wrap <- factor(pval$wrap, levels = c(levels(factor(samp2[, facet_wrap])), "Mixed"))
      }
    }
    else {
      pval <- data.frame()
      for (i in 1:length(unique(samp2$wrap))){
        index_wrap <- samp2$wrap == unique(samp2$wrap)[[i]]
        abund3 <- abund2[index_wrap, , drop = FALSE]
        pred   <- samp2[index_wrap, predictor]

        has_var <- apply(abund3, 2, function(col) {
          val <- col[!is.na(col)]
          length(unique(val)) > 1 && any(val != 0)
        })
        abund3 <- abund3[, has_var, drop = FALSE]

        if (ncol(abund3) == 0) {
          warning(sprintf("No non-zero taxa found in facet '%s'. Skipping stats.", unique(samp2$wrap)[[i]]))
          next
        }

        if (!is.null(facet_wrap) && (min_abund > 0 || min_prev > 0)) {
          n_facet_samples <- nrow(abund3)
          n_taxa_before_facet <- ncol(abund3)

          rel_facet_mat <- if(relative_abun) abund3 else apply(abund3, 1, function(x) x / sum(x)) %>% t %>% as.data.frame()
          facet_prev_frac <- colSums(rel_facet_mat > min_abund, na.rm = TRUE) / n_facet_samples
          keep_facet_taxa <- facet_prev_frac >= min_prev

          abund3 <- abund3[, keep_facet_taxa, drop = FALSE]

          message(sprintf(
            "Facet '%s' filter (min_abund = %g, min_prev = %.1f%%): Kept %d / %d taxa.",
            unique(samp2$wrap)[[i]], min_abund, min_prev * 100, ncol(abund3), n_taxa_before_facet
          ))
        }

        if (ncol(abund3) == 0) {
          warning(sprintf("No taxa passed min_abund/min_prev filters in facet '%s'. Skipping facet stats.", unique(samp2$wrap)[[i]]))
          next
        }

        # --- METAGENOMESEQ (mgs_feature) ---
        if(stats=="mgs_feature"){
          if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
            stop(
              "Package 'metagenomeSeq' is required for this function. ",
              "Install it with BiocManager::install('metagenomeSeq').",
              call. = FALSE
            )
          }
          counts_sub <- abund_org[index_wrap, colnames(abund3), drop = FALSE]
          mgs <- metagenomeSeq::newMRexperiment(counts = t(counts_sub))
          mgsp <- metagenomeSeq::cumNormStat(mgs)
          mgs <- metagenomeSeq::cumNorm(mgs, mgsp)
          mod <- stats::model.matrix(~as.numeric(pred == unique(pred)[1]))
          if(length(unique(samp2$wrap))>1) message(paste0("MGS FeatureModel for facet_wrap = ",unique(samp2$wrap)[[i]]))
          else message("MGS FeatureModel")
          if (!is.null(covariates)) message("metagenomeSeq fitFeatureModel does not support covariate adjustment.")
          mgsfit <- metagenomeSeq::fitFeatureModel(obj=mgs,mod=mod)
          if (isS4(mgsfit)) {
            mgs_taxa <- slot(mgsfit, "taxa")
            mgs_pvals <- slot(mgsfit, "pvalues")
          } else if (is.list(mgsfit)) {
            mgs_taxa <- mgsfit$taxa
            mgs_pvals <- mgsfit$pvalues
          } else {
            mgs_taxa <- rownames(abund3)
            mgs_pvals <- metagenomeSeq::MRpvalues(mgsfit)
          }
          pval_tmp <- data.frame(variable = mgs_taxa, pval = mgs_pvals, stringsAsFactors = FALSE)
        }

        # --- NON-PARAMETRIC ---
        if(stats=="non-parametric"){
          if(i==1) {
            message("Non-parametric statistics")
            if(!is.null(covariates)) {
              message("Non-parametric test does not adjust for covariates.")
            }
          }
          pval_tmp <- data.frame(
            variable = colnames(abund3),
            pval     = apply(abund3, 2, function(x) stats::kruskal.test(x ~ pred)$p.value),
            stringsAsFactors = FALSE
          )
        }

        # --- PARAMETRIC ---
        if(stats=="parametric"){
          if(i==1) message("Parametric statistics")
          if (is.null(covariates)) {
            pval_tmp <- data.frame(
              variable = colnames(abund3),
              pval     = apply(abund3, 2, function(x) stats::oneway.test(x ~ pred)$p.value),
              stringsAsFactors = FALSE
            )
          } else {
            cov_df <- samp2[index_wrap, covariates, drop = FALSE]
            fml <- stats::as.formula(paste("y ~ pred +", paste(covariates, collapse = " + ")))
            p_vec <- apply(abund3, 2, function(y) {
              fit <- stats::lm(fml, data = cbind(data.frame(y = y, pred = pred), cov_df))
              anv <- stats::anova(fit)
              anv["pred", "Pr(>F)"]
            })
            pval_tmp <- data.frame(variable = colnames(abund3), pval = as.numeric(p_vec), stringsAsFactors = FALSE)
          }
        }

        # --- MAASLIN 2 ---
        if(stats=="maaslin2"){
          if (!requireNamespace("Maaslin2", quietly = TRUE)) {
            stop(
              "Package 'Maaslin2' is required for this function. ",
              "Install it with BiocManager::install('Maaslin2').",
              call. = FALSE
            )
          }
          if(i==1) message("Maaslin2 statistics")
          test_set <- setNames(abund3, paste0('X', seq_along(abund3)))
          rownames(test_set) <- rownames(abund3)

          meta_df <- data.frame(
            pred = pred,
            samp2[index_wrap, covariates, drop = FALSE],
            row.names = rownames(samp2)[index_wrap],
            stringsAsFactors = FALSE
          )

          capture.output({fit_data <- Maaslin2::Maaslin2(
            input_data = test_set,
            input_metadata = meta_df,
            output = tempdir(),
            fixed_effects = c('pred', covariates),
            normalization = 'NONE',
            transform = "LOG",
            min_abundance = 0.0,
            min_prevalence = 0.0,
            min_variance = 0.0,
            plot_scatter = FALSE,
            plot_heatmap = FALSE,
            standardize = FALSE)$results
          }, file = nullfile())

          fit_data_pred <- fit_data[fit_data$metadata == "pred", , drop = FALSE]
          pval_tmp <- data.frame(variable=names(abund3), pval=fit_data_pred[match(names(test_set), fit_data_pred$feature), "pval"], stringsAsFactors = FALSE)
        }

        # --- MAASLIN 3 ---
        if(stats=="maaslin3"){
          if (!requireNamespace("maaslin3", quietly = TRUE)) {
            stop("Package 'maaslin3' is required. Install via BiocManager::install('maaslin3') or remotes::install_github('cBioCenter/maaslin3').", call. = FALSE)
          }
          if(i==1) message("MaAsLin 3 statistics")

          test_set <- setNames(abund3, paste0('X', seq_along(abund3)))
          rownames(test_set) <- rownames(abund3)

          meta_df <- data.frame(
            pred = pred,
            samp2[index_wrap, covariates, drop = FALSE],
            row.names = rownames(samp2)[index_wrap],
            stringsAsFactors = FALSE
          )

          formula_str <- if (!is.null(covariates)) paste("~ pred +", paste(covariates, collapse = " + ")) else "~ pred"

          capture.output({fit_data <- maaslin3::maaslin3(
            input_data = test_set,
            input_metadata = meta_df,
            output = tempdir(),
            formula = formula_str,
            normalization = "NONE",
            transform = "LOG",
            min_abundance = 0.0,
            min_prevalence = 0.0,
            min_variance = 0.0,
            warn_prevalence = FALSE,
            plot_summary_plot = FALSE,
            plot_associations = FALSE,
            standardize = FALSE)$fit_data_abundance$results
          }, file = nullfile())

          fit_data_pred <- fit_data[fit_data$metadata == "pred", , drop = FALSE]
          pval_tmp <- data.frame(variable=names(abund3), pval=fit_data_pred[match(names(test_set), fit_data_pred$feature), "pval_joint"], stringsAsFactors = FALSE)
          cat(sprintf("Maaslin3: Using combined abundance + presence/absence model\n"))
          cat(sprintf("Total input taxa: %d\n", ncol(test_set)))
          cat(sprintf("Excluded taxa (zero variance / zero present): %d\n", ncol(test_set) - sum(!is.na(pval_tmp$pval))))
          pval_tmp <- pval_tmp[!is.na(pval_tmp$pval),]
          # Target fix for MaAsLin 3 multcomp precision artifact:
          # Convert hard 0s resulting from catastrophic cancellation to machine epsilon
          pval_tmp$pval[pval_tmp$pval == 0] <- 2e-16
        }

        # --- ANCOM-BC2 ---
        if(stats=="ancombc2"){
          if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
            stop("Package 'ANCOMBC' is required. Install via BiocManager::install('ANCOMBC').", call. = FALSE)
          }
          if (!requireNamespace("microbiome", quietly = TRUE)) {
            stop("Package 'microbiome' is required by ANCOM-BC2. Install via BiocManager::install('microbiome').", call. = FALSE)
          }
          if(i==1) message("ANCOM-BC2 statistics")

          counts_sub <- abund_org[index_wrap, colnames(abund3), drop = FALSE]
          samp_sub   <- samp2[index_wrap, , drop = FALSE]

          ps_sub <- phyloseq::phyloseq(
            phyloseq::otu_table(as.matrix(t(counts_sub)), taxa_are_rows = TRUE),
            phyloseq::sample_data(samp_sub)
          )

          formula_str <- if(!is.null(covariates)) paste(c(predictor, covariates), collapse = " + ") else predictor

          capture.output({
            ancom_out <- ANCOMBC::ancombc2(
              data = ps_sub,
              fix_formula = formula_str,
              p_adj_method = p_adjust_method,
              prv_cut = 0,
              lib_cut = 0,
              verbose = FALSE
            )
          }, file = nullfile())

          res_df <- ancom_out$res
          p_col <- grep(paste0("^p_", predictor), colnames(res_df), value = TRUE)
          if(length(p_col) == 0) {
            p_cols <- grep("^p_", colnames(res_df), value = TRUE)
            p_cols <- setdiff(p_cols, "p_val")
            p_col  <- p_cols[1]
          }

          pval_tmp <- data.frame(
            variable = names(abund3),
            pval     = res_df[match(names(abund3), res_df$taxon), p_col],
            stringsAsFactors = FALSE
          )

          cat(sprintf("ANCOM-BC2: Using bias-corrected log-linear model\n"))
          cat(sprintf("Total input taxa: %d\n", ncol(counts_sub)))
          cat(sprintf("Excluded taxa (zero variance / fitting error): %d\n", sum(is.na(pval_tmp$pval))))

          pval_tmp <- pval_tmp[!is.na(pval_tmp$pval), ]
        }

        # --- DESEQ2 ---
        if(stats=="deseq2"){
          if (!requireNamespace("DESeq2", quietly = TRUE)) {
            stop("Package 'DESeq2' is required. Install via BiocManager::install('DESeq2').", call. = FALSE)
          }
          if(i==1) message("DESeq2 statistics")

          counts_sub <- abund_org[index_wrap, colnames(abund3), drop = FALSE]
          samp_sub   <- samp2[index_wrap, , drop = FALSE]
          cts        <- t(counts_sub)

          design_terms   <- if(!is.null(covariates)) c(covariates, predictor) else predictor
          design_formula <- stats::as.formula(paste("~", paste(design_terms, collapse = " + ")))

          dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = cts,
            colData   = samp_sub,
            design    = design_formula
          )

          capture.output({
            dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
            dds <- DESeq2::DESeq(dds, quiet = TRUE)
            res_name <- grep(paste0("^", predictor), DESeq2::resultsNames(dds), value = TRUE)[1]
            res <- DESeq2::results(dds, name = res_name)
          }, file = nullfile())

          res_df <- as.data.frame(res)
          res_df$taxon <- rownames(res_df)

          pval_tmp <- data.frame(
            variable = names(abund3),
            pval     = res_df[match(names(abund3), res_df$taxon), "pvalue"],
            stringsAsFactors = FALSE
          )

          cat(sprintf("DESeq2: Using Negative Binomial Wald test (poscounts size factors)\n"))
          cat(sprintf("Total input taxa: %d\n", ncol(counts_sub)))
          cat(sprintf("Excluded taxa (zero variance / fitting error): %d\n", sum(is.na(pval_tmp$pval))))

          pval_tmp <- pval_tmp[!is.na(pval_tmp$pval), ]
        }

        pval_tmp <- pval_tmp %>%
          dplyr::mutate(wrap=unique(samp2$wrap)[[i]],p_adjust=stats::p.adjust(pval, p_adjust_method))
        pval <- rbind(pval,pval_tmp)
      }
    }
    if(p_adjust) message(paste(p_adjust_method,"correction applied for",length(unique(pval$variable)),"taxa"))
  }

  if(volcano_plot == TRUE) {
    vol_df <- data.frame()
    if(!is.null(facet_wrap)) {
      samp2 <- samp2[!is.na(samp2[[facet_wrap]]) & samp2[[facet_wrap]] != "NA", , drop = FALSE]
      samp2[[facet_wrap]] <- droplevels(factor(samp2[[facet_wrap]]))
      v_wrap_vec <- as.character(samp2[[facet_wrap]])
      wrap_levels <- levels(samp2[[facet_wrap]])
    } else {
      v_wrap_vec <- rep("All", nrow(samp2))
      wrap_levels <- "All"
    }

    for(i in 1:length(wrap_levels)) {
      w <- wrap_levels[[i]]
      if(is.na(w) || w == "NA" || w == "Mixed") next

      index_wrap <- v_wrap_vec == w
      samp_sub <- samp2[index_wrap, , drop = FALSE]

      pval_sub <- if(stats == "mixed" && !is.null(facet_wrap)) pval else {
        if(!is.null(facet_wrap)) pval[pval$wrap == w, ] else pval
      }
      if (nrow(pval_sub) == 0) next

      facet_taxa_names <- as.character(pval_sub$variable)
      abund_sub <- abund2[index_wrap, facet_taxa_names, drop = FALSE]

      pred_factor <- factor(samp_sub[, predictor])
      pred_levels <- levels(pred_factor)
      if(length(pred_levels) != 2) {
        stop(paste0("volcano_plot requires a binary predictor variable (exactly 2 groups). Error in group/facet: ", w))
      }

      g1_idx <- pred_factor == pred_levels[1]
      g2_idx <- pred_factor == pred_levels[2]

      is_pre_logged <- any(abund_sub < 0, na.rm = TRUE)

      if (is_pre_logged) {
        # Already on log scale (e.g. metabolomics): simple difference of means
        m1 <- colMeans(abund_sub[g1_idx, , drop = FALSE], na.rm = TRUE)
        m2 <- colMeans(abund_sub[g2_idx, , drop = FALSE], na.rm = TRUE)
        log2fc <- m2 - m1
      } else {
        # Universal geometric mean log2FC: log-transform BEFORE averaging.
        # This properly aligns effect size with non-parametric rank tests and GLMs.
        pseudocount <- if (relative_abun) {
          nz <- abund_sub[abund_sub > 0]
          if (length(nz) > 0) min(1e-4, min(nz, na.rm = TRUE) / 2) else 1e-5
        } else {
          1
        }

        m1 <- colMeans(log2(abund_sub[g1_idx, , drop = FALSE] + pseudocount), na.rm = TRUE)
        m2 <- colMeans(log2(abund_sub[g2_idx, , drop = FALSE] + pseudocount), na.rm = TRUE)
        log2fc <- m2 - m1
      }

      p_vals_vec <- as.numeric(unlist(pval_sub$pval[match(colnames(abund_sub), pval_sub$variable)]))
      padj_vals_vec <- as.numeric(unlist(pval_sub$p_adjust[match(colnames(abund_sub), pval_sub$variable)]))

      sub_df <- data.frame(
        Taxa = colnames(abund_sub),
        log2FC = as.numeric(log2fc[colnames(abund_sub)]),
        pval = p_vals_vec,
        padj = padj_vals_vec,
        wrap = as.character(w),
        stringsAsFactors = FALSE
      )
      vol_df <- rbind(vol_df, sub_df)
    }

    vol_df <- vol_df[!is.na(vol_df$wrap) & vol_df$wrap != "NA" & vol_df$wrap != "Mixed", , drop = FALSE]

    if(!is.null(facet_wrap)) {
      vol_df$wrap <- factor(vol_df$wrap, levels = wrap_levels)
      vol_df$wrap <- droplevels(vol_df$wrap)
    }
    vol_df$Taxa_clean <- gsub("_", " ", vol_df$Taxa)
    vol_df$negLogP <- -log10(vol_df$pval)

    # Status classification
    vol_df$Status <- "NS"
    vol_df$Status[vol_df$pval < p_cutoff & vol_df$log2FC > 0] <- "Enriched"
    vol_df$Status[vol_df$pval < p_cutoff & vol_df$log2FC < 0] <- "Depleted"
    vol_df$Status <- factor(vol_df$Status, levels = c("Enriched", "Depleted", "NS"))

    # Alpha transparency scaling
    if (p_adjust) {
      vol_df$alpha_val <- ifelse(vol_df$padj < adj_p_cutoff, 1.00,
                                 ifelse(vol_df$pval < p_cutoff, 0.45, 0.25))
    } else {
      vol_df$alpha_val <- ifelse(vol_df$pval < p_cutoff, 1.00, 0.25)
    }

    # Color mapping: color_by vs. default Status
    use_color_by <- FALSE
    status_cols  <- c("Enriched" = "#C0392B", "Depleted" = "#2980B9", "NS" = "#BDC3C7")

    if (!is.null(color_by)) {
      color_by_rank <- tolower(color_by)
      if (color_by_rank %in% names(tax)) {
        use_color_by <- TRUE
        tax_match <- match(vol_df$Taxa, tax[, type])
        tax_group <- as.character(tax[tax_match, color_by_rank])
        tax_group[is.na(tax_group) | tax_group == ""] <- "unclassified"

        # Non-significant taxa are assigned "NS" (greyish)
        vol_df$colvar <- ifelse(vol_df$Status == "NS", "NS", tax_group)

        tax_cats <- sort(setdiff(unique(vol_df$colvar), "NS"))

        default_cols <- c(
          RColorBrewer::brewer.pal(8, "Set1"),
          RColorBrewer::brewer.pal(7, "Dark2"),
          RColorBrewer::brewer.pal(7, "Set2"),
          RColorBrewer::brewer.pal(12, "Set3"),
          RColorBrewer::brewer.pal(7, "Accent"),
          RColorBrewer::brewer.pal(12, "Paired")
        )

        if (!is.null(colors)) {
          assigned_cols <- rep_len(colors, length(tax_cats))
        } else {
          if (length(tax_cats) <= length(default_cols)) {
            assigned_cols <- default_cols[seq_len(length(tax_cats))]
          } else {
            assigned_cols <- grDevices::colorRampPalette(default_cols)(length(tax_cats))
          }
        }

        active_cols <- stats::setNames(assigned_cols, tax_cats)
        active_cols["NS"] <- "#BDC3C7"

        vol_df$colvar <- factor(vol_df$colvar, levels = c(tax_cats, "NS"))

        # Push NS points to background so colored points render on top
        vol_df <- vol_df[order(vol_df$colvar == "NS", decreasing = TRUE), ]
      } else {
        warning(sprintf("[rabuplot Warning]: color_by rank '%s' not found in taxonomy table. Defaulting to Status colors.", color_by))
        active_cols <- status_cols
      }
    } else {
      active_cols <- status_cols
    }

    color_var        <- if (use_color_by) "colvar" else "Status"
    legend_title_vol <- if (use_color_by) tools::toTitleCase(color_by) else NULL

    adj_p_val <- NULL
    filtering_active <- (min_abund > 0 || min_prev > 0)

    if (p_adjust && (is.null(facet_wrap) || !filtering_active)) {
      sig_adj_p_taxa <- vol_df[vol_df$padj <= adj_p_cutoff & !is.na(vol_df$padj), ]
      if (nrow(sig_adj_p_taxa) > 0) {
        adj_p_val <- max(sig_adj_p_taxa$pval, na.rm = TRUE)
      } else {
        message(paste0("No taxa met adjusted p <= ", adj_p_cutoff, ". Omitting 'adjusted p' line."))
      }
    }

    # 1. Expand y-ceiling to give top points and ggrepel labels breathing room
    y_max  <- max(vol_df$negLogP, na.rm = TRUE)
    y_ceil <- max(y_max * 1.07, -log10(p_cutoff) + 0.5)

    # 2. Generate standard evenly spaced breaks up to the tighter ceiling
    std_breaks <- pretty(c(0, y_ceil))
    std_breaks <- std_breaks[std_breaks >= 0 & std_breaks <= y_ceil]

    # 3. Assemble breaks & labels data frame to guarantee sorted alignment
    break_df <- data.frame(
      brk = std_breaks,
      lbl = as.character(std_breaks),
      stringsAsFactors = FALSE
    )

    # Add p_cutoff tick mark and drop any standard break that collides with it (< 0.5 distance)
    p_cut_val <- -log10(p_cutoff)
    break_df  <- break_df[abs(break_df$brk - p_cut_val) > 0.5, , drop = FALSE]
    break_df  <- rbind(break_df, data.frame(brk = p_cut_val, lbl = paste0("p=", p_cutoff)))

    # Add adjusted p-value tick mark if available
    if (p_adjust && is.null(facet_wrap) && !is.null(adj_p_val)) {
      adj_p_val_num <- -log10(adj_p_val)
      adj_label_prefix <- switch(
        tolower(p_adjust_method),
        "fdr" = "FDR",
        "bh"  = "BH",
        "by"  = "BY",
        tools::toTitleCase(tolower(p_adjust_method))
      )
      break_df <- break_df[abs(break_df$brk - adj_p_val_num) > 0.5, , drop = FALSE]
      break_df <- rbind(break_df, data.frame(brk = adj_p_val_num, lbl = paste0(adj_label_prefix, " p=", adj_p_cutoff)))
    }

    # 4. Sort breaks and labels monotonically
    break_df      <- break_df[order(break_df$brk), , drop = FALSE]
    custom_breaks <- break_df$brk
    custom_labels <- break_df$lbl

    p <- ggplot(vol_df, aes(x = log2FC, y = negLogP, color = .data[[color_var]])) +
      geom_point(aes(alpha = alpha_val), size = 2.8) +
      scale_color_manual(values = active_cols, name = legend_title_vol) +
      scale_alpha_identity() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
      geom_hline(yintercept = -log10(p_cutoff), linetype = "dotted", color = "black", linewidth = 0.5) +
      theme_bw() +
      theme(
        panel.grid.major = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)
      ) +
      scale_y_continuous(
        limits = c(0, y_ceil),
        breaks = custom_breaks,
        labels = custom_labels,
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      xlab(xlabs) + ylab(expression(-log[10](p))) + ggtitle(main) +
      guides(
        color = guide_legend(title = legend_title_vol, override.aes = list(size = 4, alpha = 1))
      )

    if (p_adjust && is.null(facet_wrap) && !is.null(adj_p_val)) {
      p <- p + geom_hline(yintercept = -log10(adj_p_val), linetype = "dashed", color = "black", linewidth = 0.5)
    }

    if(!is.null(facet_wrap)) {
      if(is.null(facet_label)) label_names <- levels(factor(samp[,facet_wrap]))
      if(!is.null(facet_label)) label_names <- facet_label
      if(facet_n==TRUE){
        label_names_df <- samp2 %>%
          dplyr::filter(!is.na(.data[[facet_wrap]]) & .data[[facet_wrap]] != "NA") %>%
          dplyr::group_by(.data[[facet_wrap]]) %>%
          dplyr::summarize(n = dplyr::n(), .groups = "drop") %>%
          dplyr::mutate(pasted_label = paste0(.data[[facet_wrap]], ", n = ", n))
        label_names <- as.character(label_names_df$pasted_label)
        names(label_names) <- as.character(label_names_df[[facet_wrap]])
      } else {
        names(label_names) <- levels(factor(samp2[[facet_wrap]]))
      }
      if(stats=="mixed") {
        label_names <- c(label_names, "Mixed" = "Mixed")
      }
      p <- p + facet_grid(~wrap, labeller = labeller(wrap=label_names),scales = "free") + theme(strip.background = element_blank())
    }
    if (label_taxa) {
      sig_df <- data.frame()
      for(w in wrap_levels) {
        if(is.na(w) || w == "NA" || w == "" || w == "Mixed") next
        sub_v <- vol_df[vol_df$wrap == w & !is.na(vol_df$wrap), ]
        target_sig_df <- if(p_adjust) sub_v[sub_v$padj < adj_p_cutoff | sub_v$pval < p_cutoff,] else sub_v[sub_v$pval < p_cutoff,]
        sub_sig <- target_sig_df %>%
          dplyr::arrange(if(p_adjust) padj else pval) %>%
          head(top_n_labels)
        sig_df <- rbind(sig_df, sub_sig)
      }

      sig_df <- sig_df[!is.na(sig_df$wrap) & sig_df$wrap != "NA" & sig_df$wrap != "" & sig_df$wrap != "Mixed", , drop = FALSE]
      if(!is.null(facet_wrap)) {
        sig_df$wrap <- factor(sig_df$wrap, levels = wrap_levels)
      }

      if (nrow(sig_df) > 0) {
        font_style <- if (italic_names && type %in% c("genus", "family", "species")) "italic" else "plain"

        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = sig_df,
            aes(label = Taxa_clean),
            color = "black",
            size = 3.5,
            fontface = font_style,
            show.legend = FALSE,
            segment.color = "grey30",
            segment.size = 0.4,
            max.overlaps = 30,
            box.padding = 0.4
          )
        } else {
          p <- p + geom_text(
            data = sig_df,
            aes(label = Taxa_clean),
            color = "black",
            size = 3.5,
            fontface = font_style,
            vjust = -0.5,
            show.legend = FALSE
          )
        }
      }
    }

    if (no_legends) p <- p + theme(legend.position = "none")
    if (stat_out) {
      pval_out <<- vol_df
    }
    return(p)
  }

  bacteria <- rev(names(abund))
  subset <- cbind(samp[!names(samp) %in% bacteria], abund) #fjerner evt eksisterende navne fra dataset og merger;
  subset$predictor2 <-  as.factor(subset[,predictor])
  subset$ID <- rownames(subset)
  if(!is.null(facet_wrap)){
    subset$wrap <-  as.factor(subset[,facet_wrap])
    molten <- subset[, c("ID", paste(bacteria), "predictor2", "wrap")] %>% tidyr::pivot_longer(cols = -c(ID, predictor2, wrap), names_to = "variable", values_to = "value")
  }
  if(is.null(facet_wrap)){
    molten <- subset[, c("ID", bacteria, "predictor2")] %>% tidyr::pivot_longer(cols = -c(ID, predictor2), names_to = "variable", values_to = "value")
  }

  if (!is.null(color_by)) {
    other_label <- paste("Other", type)

    non_other <- molten$variable != other_label

    molten$colvar <- NA_character_

    tax_match <- match(
      molten$variable[non_other],
      tax[, type]
    )

    molten$colvar[non_other] <- as.character(
      tax[tax_match, color_by]
    )

    molten$colvar[!non_other] <- paste("Other", color_by)
  }

  molten$variable <- gsub("_", " ", molten$variable)
  ordered <- unique(molten$variable)

  if (!order) {
    ordered <- sort(ordered)
  }

  molten$variable <- factor(
    molten$variable,
    levels = ordered
  )

  if (is.null(color_by)) {
    molten$colvar <- as.character(molten$variable)
  }

  default_cols <- c(
    RColorBrewer::brewer.pal(8, "Set1"),
    RColorBrewer::brewer.pal(7, "Dark2"),
    RColorBrewer::brewer.pal(7, "Set2"),
    RColorBrewer::brewer.pal(12, "Set3"),
    RColorBrewer::brewer.pal(7, "Accent"),
    RColorBrewer::brewer.pal(12, "Paired"),
    "gray"
  )

  predictor_levels <- levels(factor(molten$predictor2))

  if (is.null(legend_names)) {
    legend_names <- predictor_levels
  }

  if (is.null(colors)) {
    cols <- default_cols[seq_len(length(predictor_levels))]
  } else {
    cols <- colors
  }

  ordered2 <- rev(unique(as.character(molten$colvar)))

  if (reverse) {
    if (!bar_chart) {
      molten$predictor2 <- factor(
        molten$predictor2,
        levels = rev(predictor_levels)
      )

      legend_names <- rev(legend_names)
      cols <- rev(cols)
    } else {
      molten$colvar <- factor(
        molten$colvar,
        levels = rev(unique(as.character(molten$colvar)))
      )

      molten$variable <- factor(
        molten$variable,
        levels = rev(levels(molten$variable))
      )

      cols <- rev(cols)
      ordered2 <- rev(ordered2)
    }
  }

  if (bar_chart) {
    log <- FALSE

    colvar_levels <- levels(factor(molten$colvar))

    if (is.null(colors)) {
      if (is.null(color_by) && !bar_chart_stacked) {
        cols <- default_cols[seq_len(length(predictor_levels))]
      } else {
        n_colvar <- length(colvar_levels)

        if (n_colvar <= 1L) {
          cols <- tail(default_cols, 1L)
        } else {
          cols <- c(
            default_cols[seq_len(n_colvar - 1L)],
            tail(default_cols, 1L)
          )
        }
      }

      if (
        (is.null(color_by) && !reverse) ||
        (!is.null(color_by) && reverse)
      ) {
        cols <- rev(cols)
      }
    }

    if (is.null(facet_wrap)) {
      molten$wrap <- ""
    }

    molten_mean <- molten %>%
      dplyr::group_by(variable, predictor2, wrap, colvar) %>%
      dplyr::summarise(
        value = mean(value, na.rm = TRUE),
        .groups = "drop"
      )

    molten_mean$colvar <- factor(
      molten_mean$colvar,
      levels = ordered2
    )
  }

  # --- ORIGINAL RETAINED: Calculate pvalue for outcomes ---
  if(p_val==TRUE & ((bar_chart==TRUE & bar_chart_stacked==FALSE) | bar_chart==FALSE) & is.null(color_by)){
    if(is.null(facet_wrap)) molten$wrap <- ""
    if(!is.null(facet_wrap)) {
      pval <- data.frame(pval=pval[gsub('_',' ',pval$variable) %in% ordered,]$pval,p_adjust=pval[gsub('_',' ',pval$variable) %in% ordered,]$p_adjust, variable=gsub('_',' ',pval[gsub('_',' ',pval$variable) %in% ordered,]$variable),wrap=pval[gsub('_',' ',pval$variable) %in% ordered,]$wrap)
    }
    else {
      pval <- data.frame(pval=pval[gsub('_',' ',pval$variable) %in% ordered,]$pval,p_adjust=pval[gsub('_',' ',pval$variable) %in% ordered,]$p_adjust, variable=gsub('_',' ',pval[gsub('_',' ',pval$variable) %in% ordered,]$variable))
      if(length(pval$variable)-length(ordered)<0) pval <- pval[match(pval$variable,ordered[length(pval$variable)-length(ordered)]),]
    }
    pval$predictor2 <- molten$predictor2[1]
    pval$pval <- ifelse(is.na(pval$pval),1,pval$pval)
    pval$p_adjust <- ifelse(is.na(pval$p_adjust),1,pval$p_adjust)
    if(Only_sig){
      index <- pval[pval$pval<0.05,"variable"]
      molten <- molten[molten$variable %in% index,]
      pval <- pval[pval$pval<0.05,]
    }

    if(stat_out){
      median_iqr <<- molten %>% dplyr::group_by(variable, predictor2) %>% dplyr::summarize( N = length(value),median = stats::median(value)*100,Q1=stats::quantile(value, 1/4)*100,Q3=stats::quantile(value, 3/4)*100, IQR = stats::IQR(value)) %>% as.data.frame
      pval_out <<- pval
      mean_sd <<- molten %>% dplyr::group_by(variable, predictor2) %>% dplyr::summarize( N = length(value),mean = mean(value)*100,sd=stats::sd(value)*100) %>% as.data.frame
    }
  }

  if(bar_chart==FALSE){
    if(ncol(tax)>=6) molten$value <- molten$value+1e-6 #add pseudocount for log scale 0;
    else   molten$value <- molten$value+0.001 #add pseudocount for log scale 0;
    ordered <- levels(factor(molten$colvar))
    p <- ggplot(molten, aes(x=variable, y=value, fill=predictor2)) +
      {if(violin){geom_violin(scale = violin_scale,width = 0.65, position=position_dodge(width=0.9),size=1, color="#00000000")} else {geom_boxplot(width = 0.55, position=position_dodge(width=0.8),size=0.3,outlier.size = 0,outlier.color = "grey")}}+
      {if(violin){stat_summary(fun=stats::median, fun.min = min, fun.max = max, geom="point", size=0.8, color="black", position=position_dodge(width=0.9))} else {stat_summary(fun=stats::median, fun.min = min, fun.max = max, geom="point", size=0.8, color="#00000000", position=position_dodge(width=0.9))}}+ theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),legend.text=element_text(size=12),legend.key.size = grid::unit(0.5, "cm"))+ coord_flip() +xlab(NULL)+ylab(xlabs)+ggtitle(main)
    if(length(unique(molten$variable))>1) p <- p+ geom_vline(xintercept=seq(1.5, length(unique(molten$variable))-0.5, 1),lwd=0.2, colour="grey")

    p <- p +  scale_fill_manual(values =cols,labels=legend_names) + guides(fill = guide_legend(title=legend_title, reverse = TRUE,override.aes = list(linetype=0, shape=16,color=rev(cols),size=5, bg="white")))
  }

  if(bar_chart==TRUE){
    if(bar_chart_stacked==TRUE)
      if(!is.null(color_by)) cols <- rev(cols)
    p <-  ggplot(molten_mean,aes(x=factor(predictor2,levels=legend_names,labels=legend_names),y=value, fill=variable)) + theme_bw()+geom_bar(stat="identity")+ theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),axis.title=element_text(size=14),legend.text=element_text(size=12), axis.text = element_text(size = 12),strip.text = element_text(size = 12),legend.key.size = grid::unit(0.5, "cm"),text=element_text(size=12)) +xlab(NULL)+ylab(ylabs)+ggtitle(main) +  scale_fill_manual(values =cols,labels=ordered) + guides(fill = guide_legend(title=NULL))
    if(bar_chart_stacked==FALSE){
      if(!is.null(color_by)) p <-   ggplot(molten_mean,aes(x=variable,y=value, fill=colvar,group=wrap))+geom_bar(stat="identity", position = position_dodge(width = 0.95))+ scale_fill_manual(values =cols,labels=ordered2)+ guides(fill = guide_legend(title=color_by))
      else {
        p <-   ggplot(molten_mean,aes(x=variable,y=value, fill=predictor2))+geom_bar(stat="identity", position = position_dodge(width = 0.95))+ scale_fill_manual(values =cols,labels=legend_names)+ guides(fill = guide_legend(title=legend_title))
      }
      if(length(unique(molten_mean$variable))>1) p <- p+ geom_vline(xintercept=seq(1.5, length(unique(molten_mean$variable))-0.5, 1),lwd=0.2, colour="grey")
      p <-  p+ theme_bw()  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),axis.title=element_text(size=14),legend.text=element_text(size=12), axis.text = element_text(size = 12),strip.text = element_text(size = 12),legend.key.size = grid::unit(0.5, "cm"),text=element_text(size=12)) +xlab(NULL)+ylab(ylabs)+ggtitle(main)+ theme(strip.background = element_blank()) +coord_flip()
    }
  }

  if(!is.null(facet_wrap))   {
    if(is.null(facet_label)) label_names <- levels(factor(samp[,facet_wrap]))
    if(!is.null(facet_label)) label_names <- facet_label
    if(facet_n==TRUE){
      label_names <- samp2 %>%
        dplyr::group_by(get(facet_wrap)) %>%
        dplyr::summarize(n = dplyr::n()) %>%
        dplyr::mutate(pasted_label = paste0(levels(factor(samp2[,facet_wrap])), ", n = ", n))
      label_names <- as.character(label_names$pasted_label)
    }
    names(label_names) <- levels(factor(samp2[,facet_wrap]))
    if(stats=="mixed") {
      label_names <- c(label_names,"Mixed")
      names(label_names) <- c(levels(factor(samp2[,facet_wrap])),"Mixed")
    }
    p <- p+ facet_grid(~wrap,labeller = labeller(wrap=label_names),scales = "free", space = "free")+  theme(strip.background = element_blank())
    if(bar_chart==FALSE) p$layers[4:5] <- NULL
  }

  if(italic_names==TRUE &  (bar_chart==FALSE | (bar_chart==TRUE & bar_chart_stacked==FALSE)))  {
    if(type=="genus" | type=="family" | type=="species")
      if(no_other_type==F)
        p <- suppressWarnings(p+ theme(axis.text.y=element_text(face = c("plain", rep("italic",length(unique(molten$variable))-1)))))
    else p <- p+ theme(axis.text.y=element_text(face = "italic"))
  }

  if(!is.null(color_by)) {
    if(color_by=="genus" | color_by=="family" | color_by=="species") p <- p+ theme(legend.text=element_text(face = "italic"))
    if(color_by==type & bar_chart_stacked==FALSE ) p <- p+theme(legend.position="none")
  }

  if(p_val==TRUE){
    if(log==FALSE){
      if(bar_chart==TRUE) pval$y <- max(molten_mean$value)*1.10
      else pval$y <- max(molten$value)*1.15
    }
    else pval$y <-ifelse(log_max==100,10,ifelse(log_max==10,0.126,0.0126))
    if(p_adjust==TRUE){
      if(log==FALSE & bar_chart==FALSE) pval$y_adjust <- 1.22
      if(log==FALSE & bar_chart==TRUE) pval$y_adjust <- max(molten_mean$value)*1.25
      if(log==TRUE) pval$y_adjust <- ifelse(log_max==100,105,ifelse(log_max==10,1.26,0.126))
    }
  }

  if(log==TRUE){
    if(p_val==FALSE){
      if(log_max == 100)  p <- p+ scale_y_log10(breaks=c(.000001,.001,.01,.1,1),labels=c("0%","0.1%","1%","10%","100%"))
      if(log_max == 10)  p <- p+ scale_y_log10(limits=c(0.001,0.13),breaks=c(.001,.01,.05,.1),labels=c("0%","1%","5%","10%"))
      if(log_max == 1)  p <- p+ scale_y_log10(limits=c(0.001,0.013),breaks=c(.001,.01),labels=c("0%","1%"))
    }
    if(p_val==TRUE){
      if(p_adjust){
        if(log_max == 100)  p <- p+ scale_y_log10(breaks=c(.000001,.001,.01,.1,1,7,70),labels=c("0%","0.1%","1%","10%","100%", "P-value", "q-value"))
        if(log_max == 10)  p <- p+ scale_y_log10(breaks=c(.001,.01,.05,0.1,0.126,1.26),labels=c("0%","1%","5%","10%", "P-value", "q-value"))
        if(log_max == 1)  p <- p+ scale_y_log10(breaks=c(.001,.01,0.0126,0.126),labels=c("0%","1%", "P-value", "q-value"))
      }
      else{
        if(log_max == 100)  p <- p+ scale_y_log10(breaks=c(.000001,.001,.01,.1,1,7),labels=c("0%","0.1%","1%","10%","100%", "P-value"))
        if(log_max == 10)  p <- p+ scale_y_log10(breaks=c(.001,.01,.05,0.10,0.126),labels=c("0%","1%","5%","10%", "P-value"))
        if(log_max == 1)  p <- p+ scale_y_log10(breaks=c(.001,.01,0.0126),labels=c("0%","1%", "P-value"))
      }
    }
  }

  if(log==FALSE){
    if(p_val==FALSE) p <- p + scale_y_continuous(breaks=c(0,.25,.50,.75,1),labels=c("0%","25%","50%","75%","100%"))
    if(p_val==TRUE){
      if(p_adjust==TRUE) {
        if(max(molten_mean$value)<=1 & max(molten_mean$value)>=0.75) p <- p + scale_y_continuous(breaks=c(0,.25,.50,.75,1,max(molten_mean$value)*1.07,max(molten_mean$value)*1.25),labels=c("0%","25%","50%","75%","100%", "P-value", "q-value"))
        if(max(molten_mean$value)<0.75 & max(molten_mean$value)>=0.50) p <- p + scale_y_continuous(breaks=c(0,.25,.50,max(molten_mean$value)*1.07,max(molten_mean$value)*1.25),labels=c("0%","25%","50%", "P-value", "q-value"))
        if(max(molten_mean$value)<0.50 & max(molten_mean$value)>=0.25) p <- p + scale_y_continuous(breaks=c(0,.1,.2,.3,.4,max(molten_mean$value)*1.07,max(molten_mean$value)*1.25),labels=c("0%","10%","20%","30%","40%", "P-value", "q-value"))
        if(max(molten_mean$value)<0.25) p <- p + scale_y_continuous(breaks=c(0,.05,.1,.15,.2,max(molten_mean$value)*1.07,max(molten_mean$value)*1.25),labels=c("0%","5%","10%","15%","20%", "P-value", "q-value"))
      }
      if(p_adjust==FALSE) {
        if(max(molten_mean$value)<=1 & max(molten_mean$value)>=0.75)  p <- p + scale_y_continuous(breaks=c(0,.25,.50,.75,1,max(molten_mean$value)*1.07),labels=c("0%","25%","50%","75%","100%", "P-value"))
        if(max(molten_mean$value)<0.75 & max(molten_mean$value)>=0.50) p <- p + scale_y_continuous(breaks=c(0,.25,.50,max(molten_mean$value)*1.07),labels=c("0%","25%","50%", "P-value"))
        if(max(molten_mean$value)<0.50 & max(molten_mean$value)>=0.25) p <- p + scale_y_continuous(breaks=c(0,.1,.2,.3,.4,max(molten_mean$value)*1.07),labels=c("0%","10%","20%","30%","40%", "P-value"))
        if(max(molten_mean$value)<0.25) p <- p + scale_y_continuous(breaks=c(0,.05,.1,.15,.2,max(molten_mean$value)*1.07),labels=c("0%","5%","10%","15%","20%", "P-value"))
      }
    }
  }

  p <-  p + theme(plot.background = element_blank(),panel.background = element_blank(),plot.title = element_text(hjust = 0.5))
  p <-  p + theme(axis.text.x = element_text(angle = text_angle_x, vjust = ifelse(text_angle_x<0 & text_angle_x>-90 , 0, ifelse(text_angle_x>=0 & text_angle_x<90, 1, 0.5)), hjust=ifelse(text_angle_x==0 | text_angle_x==180 | text_angle_x==-180 | text_angle_x==180, 0.5, ifelse((text_angle_x<0 & text_angle_x>=-90) | text_angle_x>=270, 0, 1))))
  if (bar_chart==TRUE & bar_chart_stacked==FALSE & percent==TRUE)  p <- p+  geom_text(aes(label = paste0(sprintf("%.2f",value*100), "%")), hjust = -.12, position=position_dodge(width=0.95))+scale_y_continuous(limits=c(0,max(molten_mean$value)+0.2),labels = scales::percent)
  if(no_legends) p <- p + theme(legend.position="none")
  if(no_names)  p <- p + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  stars.pval <- function (p.value)
  {    unclass(stats::symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**",  "*", "NS")))
  }
  if(p_stars==TRUE & p_val==TRUE) p <- p + geom_text(data=pval,aes(x=variable,y=y,label=paste(stars.pval(pval))) ,size=3,hjust=1)

  if(p_stars==FALSE & p_val==TRUE & (bar_chart==FALSE | (bar_chart==TRUE & bar_chart_stacked==FALSE))){
    p <- p + geom_text(data=pval,aes(x=variable,y=y,label=ifelse(pval<0.05, paste(format.pval(pval,1,0.001,nsmall=3)),"")) ,size=3,hjust=1,fontface="bold")
    p <- p + geom_text(data=pval,aes(x=variable,y=y,label=ifelse(pval>=0.05, paste(format.pval(pval,1,0.001,nsmall=3)),"")) ,size=3,hjust=1)
    if(p_adjust){
      p <- p + geom_text(data=pval,aes(x=variable,y=y_adjust,label=ifelse(p_adjust<0.05, paste(format.pval(p_adjust,1,0.001,nsmall=3)),"")) ,size=3,hjust=1,fontface="bold")
      p <- p + geom_text(data=pval,aes(x=variable,y=y_adjust,label=ifelse(p_adjust>=0.05, paste(format.pval(p_adjust,1,0.001,nsmall=3)),"")) ,size=3,hjust=1)
    }
    if(stats=="mixed" & !is.null(facet_wrap)){
      if(bar_chart==FALSE) p <- p + expand_limits(y = 2)
      if(bar_chart==TRUE) p <- p + expand_limits(y = max(molten_mean2$value))
    }
  }
  p
}
