#' Relative Abundance Plot
#'
#' For creating nice microbiome plots
#'
#' @param phylo_ob Phyloseq object with metadata in sample_data.
#' @param predictor Predictor of interest for statistics/plotting in sample_data.
#' @param type Taxonomic rank from tax_table, case insensitive; default is "genus".
#' @param relative_abun Use relative abundances, else absolute; default is TRUE.
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
#' @param Timepoint Value in variable @Time to select.
#' @param Strata Name of variable for stratification;
#' @param Strata_val Value in variable @Strata to keep; default is 1.
#' @param no_legends Removes legend; default is FALSE.
#' @param no_names Removes taxa names; default is FALSE.
#' @param italic_names Taxa names will be in italic e.g. usable for family, genus, species levels; default is TRUE
#' @param Only_sig Only keep significant taxa; default is FALSE.
#' @param log Present plot on a log scale; default is TRUE.
#' @param log_max Maximum value of log-axis options:1, 10, 100; default is 100.
#' @param stat_out Outputs a data.frame with statistics to Global environment; default is FALSE.
#' @param p_val Displays p-values on plot; default is TRUE.
#' @param p_stars Shows stars instead of p-values; default is FALSE.
#' @param stats Select type of statistical test; options: "non-parametric", "parametric", "mixed", "mgs_feature"; default is "non-parametric".
#' @param p_adjust adjust p-values; default is "FALSE.
#' @param p_adjust_method options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"; default is "fdr".
#' @param p_adjust_full correction applied for all taxa in the dataset; default is FALSE.
#' @param colors define list of colors for plot. If not color brewer will be used; default is NULL.
#' @param color_by define taxonomic rank to color by; default is NULL.
#' @param order Order by abundance, else alphabetically; default is TRUE.
#' @param reverse Flip taxa order; default is FALSE.
#' @param list_taxa A list of specific taxa names to be analyzed; default is NULL.
#' @param select_taxa Choose all taxa from one or more taxonomic variables, eg. "Staphylococcus" or "Staph" or "coccus" or c("staph",bifido"); default is NULL.
#' @param select_type Taxonomic rank of the @select_taxa; default is "genus".
#' @param bar_chart Choose to make bar chart; default is FALSE.
#' @param bar_chart_stacked Produce stacked bar chart; default is FALSE
#' @param percent Print percentages on bar chart; default is FALSE.
#' @param facet_wrap Facet wrap chart by variable; eg. Time; default is NULL.
#' @param facet_label Facet wrap labels; default is NULL.
#' @param facet_n Show n for each facet; default is TRUE.
#' @param order_by Choose variable to order the selected taxa by; eg. Time; default is Time.
#' @param order_val Choose value for @order_by; default is NULL.
#'
#' @import ggplot2 phyloseq metagenomeSeq dplyr tidyr RColorBrewer lmerTest
#' @return A ggplot
#' @export

rabuplot <- function(phylo_ob,
                     predictor="none",
                     type="genus",
                     relative_abun=TRUE,
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
                     order_val=NULL)
{
  if(!is.null(list_taxa) & is.null(N_taxa)) N_taxa = length(list_taxa)
  if(is.null(N_taxa) & is.null(list_taxa)) N_taxa=15
  options(dplyr.summarise.inform = FALSE)
  if(bar_chart_stacked==TRUE & bar_chart==FALSE) {
    bar_chart=TRUE
    p_val=FALSE
  }
  if(predictor=="none") {
    sample_data(phylo_ob)$none <- "All samples"
    p_val=FALSE
    if(bar_chart_stacked==FALSE & is.null(color_by)) no_legends = TRUE
  }
  phylo_ob <- prune_samples(sample_sums(phylo_ob)>0,phylo_ob) #removes empty samples;
  otu_mat <- as(otu_table(phylo_ob), "matrix")
  if(taxa_are_rows(phylo_ob)) otu_mat <- t(otu_mat)
  if(!is.null(facet_wrap)) index <- !is.na(get_variable(phylo_ob, predictor)) & !is.na(get_variable(phylo_ob, facet_wrap))
  else   index <- !is.na(get_variable(phylo_ob, predictor))
  if(length(unique(index)) !=1) message(paste(length(which(index==F)), "samples have been removed from full dataset (predictor/facet_wrap NAs)"))
  otu_mat <- otu_mat[index,]
  otu_mat  <- otu_mat[,colSums(otu_mat)>0] #removes empty OTUs;
  OTU_index <- colnames(otu_mat)
  tax <- as(tax_table(phylo_ob), "matrix") %>% data.frame(stringsAsFactors=FALSE)
  tax <- tax[rownames(tax) %in% OTU_index,]
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
    otu_mat  <- otu_mat[,colSums(otu_mat)>0] #removes empty OTUs;
    OTU_index <- colnames(otu_mat)
    tax <- tax[rownames(tax) %in% OTU_index,]
    samp <- samp[rownames(samp) %in% index,]
  }

  list <-as.character(tax[,type])
  unique_tax <- unique(list)

  abund <- as.data.frame(matrix(rep(0,(length(unique_tax)*nrow(otu_mat))),ncol=length(unique_tax)))
  row.names(abund) <- row.names(otu_mat)
  names(abund) <- unique_tax
  for(i in names(abund)){
    if(is.array(otu_mat[,list==i]))  abund[,i] <- rowSums(otu_mat[,list== i])
    else   abund[,i] <- otu_mat[,list== i]
  }
  abund_org <- abund
  if(relative_abun==TRUE) abund <- apply(abund,1,function(x) x/sum(x)) %>% t %>% as.data.frame()

  if (is.null(list_taxa) & !is.null(select_taxa)) {
    list_taxa <- NULL
    for(i in 1:length(select_taxa)){
      list_taxa <- c(list_taxa,(as.character(unique(tax[grep(select_taxa[[i]],tax[,select_type],ignore.case=TRUE),type]))))
    }
  }
  if (!is.null(list_taxa)) {
    no_other_type <- TRUE
    if (is.null(N_taxa)) N_taxa <- length(list_taxa)
    abund <- abund[,colnames(abund) %in% list_taxa, drop = FALSE]
    unique_tax <- names(abund)
  }

  if(length(abund)>1){
    index <- !is.na(rownames(samp))
    if (!is.null(order_val))  index <- samp[,order_by] ==order_val
    abund <- abund[,order(-colSums(abund[index,]))]
    if (By_median)  abund <- abund[,order(-apply(abund[index,], 2, median))]
    if("unclassified" %in% unique_tax) abund <- abund[c(setdiff(names(abund), "unclassified"),"unclassified")] #Move unclassified to end
    if(N_taxa<length(unique_tax) ){
      abund[,paste("Other",type)] <- rowSums(abund[(length(unique_tax)-(length(unique_tax)-N_taxa)+1):length(unique_tax)])
      abund <- abund[-(length(unique_tax)-(length(unique_tax)-N_taxa)+1):-length(unique_tax)]
    }
    if(no_other_type)  abund[,paste("Other",type)] <- NULL
  }
  index <- !is.na(rownames(samp))
  if(!is.null(Strata)) index <- samp[,Strata]==Strata_val
  samp2 <- samp %>% filter(index)
  if(p_val==TRUE & (bar_chart==FALSE | (bar_chart==TRUE & bar_chart_stacked==FALSE))){
    if(p_adjust_full ==TRUE | stats=="mgs_feature"){
      abund2 <- abund_org %>% filter(index)
      if(relative_abun==TRUE & stats!="mgs_feature") abund2 <- apply(abund2,1,function(x) x/sum(x)) %>% t %>% as.data.frame()
    }
    else abund2 <- abund %>% filter(index)
    if(stats=="mgs_feature" & length(levels(factor(samp2[,predictor])))>2){
      stats="non-parametric"
      message("MGS not available for >2 predictors, switching to non-parametric")
    }
    if(stats=="mixed" & is.null(id)){
      stats="non-parametric"
      message("No id variable for mixed model, switching to non-parametric")
    }
    if(stats=="mixed"){
      pred <- samp2[,predictor]
      id <- samp2[,id]
      message("Mixed model statistics")
      pval <- cbind(abund2,pred,id) %>% as_tibble() %>%
        gather(variable, value,-"pred",-"id") %>%
        group_by(variable) %>%
        summarize(pval = lmer(value ~ pred + (1 | id)) %>% anova %>% .$'Pr(>F)', .groups = 'drop') %>%
        mutate(wrap="Mixed",p_adjust=p.adjust(pval, p_adjust_method))
    }
    else {
      pval <- data.frame()
    for (i in 1:length(unique(samp2$wrap))){
      index <- samp2$wrap==unique(samp2$wrap)[[i]]
      abund3 <- abund2 %>% filter(index)
      pred <- samp2[index,predictor]
      # test with featureModel
      if(stats=="mgs_feature"){
        mgs <- metagenomeSeq::newMRexperiment(counts = t(abund3))
        mgsp <- metagenomeSeq::cumNormStat(mgs)
        mgs <- metagenomeSeq::cumNorm(mgs, mgsp)
        mod <- model.matrix(~as.numeric(pred == unique(pred)[1]))
        if(length(unique(samp2$wrap))>1) message(paste0("MGS FeatureModel for facet_wrap = ",unique(samp2$wrap)[[i]]))
        else message("MGS FeatureModel")
        mgsfit <- metagenomeSeq::fitFeatureModel(obj=mgs,mod=mod)
        pval_tmp <- data.frame(variable=mgsfit$taxa,pval=mgsfit$pvalues)
      }
      if(stats=="non-parametric"){   #Kruskal-Wallis
        if(i==1) message("Non-parametric statistics")
        pval_tmp <- cbind(abund3,pred) %>% as_tibble() %>%
          gather(variable, value,-"pred") %>%
          group_by(variable) %>%
          summarize(pval = kruskal.test(value ~ pred)$p.value, .groups = 'drop')
      }
      if(stats=="parametric"){
        if(i==1) message("Parametric statistics")
        pval_tmp <- cbind(abund3,pred) %>% as_tibble() %>%
          gather(variable, value,-"pred") %>%
          group_by(variable) %>%
          summarize(pval = oneway.test(value ~ pred)$p.value, .groups = 'drop')
      }
      pval_tmp <- pval_tmp %>%
        mutate(wrap=unique(samp2$wrap)[[i]],p_adjust=p.adjust(pval, p_adjust_method))
      pval <- rbind(pval,pval_tmp)
    }
    }
    if(p_adjust) message(paste(p_adjust_method,"correction applied for",length(unique(pval$variable)),"taxa"))
  }

  bacteria <- rev(names(abund))
  subset <- cbind(samp[!names(samp) %in% bacteria], abund) #fjerner evt eksisterende navne fra dataset og merger;
  subset$predictor2 <-  as.factor(subset[,predictor])
  subset$ID <- rownames(subset)
  if(!is.null(Strata)) subset[,Strata] <- as.factor(subset[,Strata])
  if(!is.null(facet_wrap)){
    subset$wrap <-  as.factor(subset[,facet_wrap])
    if(!is.null(Strata))
      molten <- subset[,c("ID",paste(bacteria),"predictor2",Strata,"wrap")] %>% gather(variable, value,-"predictor2",-"ID",-all_of(Strata),-"wrap")
    else
      molten <- subset[,c("ID",paste(bacteria),"predictor2","wrap")] %>% gather(variable, value,-"predictor2",-"ID",-"wrap")
  }
  if(is.null(facet_wrap)){
    if(!is.null(Strata))
      molten <- subset[,c("ID",paste(bacteria),"predictor2",Strata)] %>% gather(variable, value,-"predictor2",-"ID",-all_of(Strata))
    else
      molten <- subset[,c("ID",paste(bacteria),"predictor2")] %>% gather(variable, value,-"predictor2",-"ID")
  }
  if(!is.null(color_by)){
    molten[molten$variable != paste("Other",type),"colvar"] <- molten %>% dplyr::filter(variable != paste("Other",type)) %>% .[,"variable"] %>% match(tax[,type]) %>% tax[.,color_by] %>% as.character
    molten[molten$variable == paste("Other",type),"colvar"] <- paste("Other",color_by) %>% as.character
  }

  molten$variable <- gsub('_',' ',molten$variable)

  if(order)   ordered <- unique(molten$variable) #level order
  if(!order)   ordered <-sort(unique(molten$variable))#level order alphabetically

  molten$variable <- factor(molten$variable, levels=ordered)
  if(is.null(color_by))  molten$colvar <- molten$variable
  if(!is.null(Strata))  molten <- molten[which(molten[,Strata]==Strata_val), ]

  if(is.null(colors)){
    cols  <- c(brewer.pal(8,"Set1"), brewer.pal(7,"Dark2"),brewer.pal(7,"Set2"),brewer.pal(12,"Set3"),brewer.pal(7,"Accent"),brewer.pal(12,"Paired"),"gray")
    cols <- cols[1:length(levels(factor(molten$predictor2)))]
  }
  if(!is.null(colors)) cols <- colors

  if(bar_chart==TRUE & bar_chart_stacked==FALSE & is.null(legend_names))  legend_names <- as.character(levels(factor(molten$predictor2)))
  if(is.null(legend_names))  legend_names <- as.character(levels(factor(molten$predictor2)))
  ordered2<- rev(unique(molten$colvar))
  if(reverse){
    if(bar_chart==FALSE) {
      molten$predictor2 <- factor(molten$predictor2, levels=rev(levels(molten$predictor2)))#manual faceting for levels;
      legend_names <- rev(legend_names)
      cols <- rev(cols)
    }
    if(bar_chart==TRUE) {
      molten$colvar <- factor(molten$colvar, levels=rev(levels(factor(molten$colvar))))#manual faceting for levels;
      molten$variable <- factor(molten$variable, levels=rev(levels(factor(molten$variable))))
      cols <- rev(cols)
      ordered2<- rev(ordered2)
    }
  }

  if(bar_chart){
    log=FALSE
    cols  <- c(brewer.pal(8,"Set1"), brewer.pal(7,"Dark2"),brewer.pal(7,"Set2"),brewer.pal(12,"Set3"),brewer.pal(7,"Accent"),brewer.pal(12,"Paired"),"gray")
    #  ordered <- levels(factor(molten$colvar))
    if(is.null(color_by) & bar_chart_stacked==FALSE)   cols <- cols[1:length(levels(factor(molten$predictor2)))]

    else cols <- cols[c(1:length(levels(factor(molten$colvar)))-1,length(cols))]
    if(!is.null(colors)) cols <- colors
    if(is.null(color_by) & reverse==FALSE) cols <- rev(cols)
    if(!is.null(color_by) & reverse==TRUE) cols <- rev(cols)
    if(is.null(facet_wrap))  molten$wrap <- ""
    molten_mean <- molten %>%
      dplyr::group_by(variable,predictor2,wrap,colvar) %>%
      dplyr::summarize(value = mean(value))
    molten_mean$colvar <- factor(molten_mean$colvar, levels=ordered2)
  }
  #Calculate pvalue for outcomes
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
      median_iqr <<- molten %>% dplyr::group_by(variable, predictor2) %>% dplyr::summarize( N = length(value),median = median(value)*100,Q1=quantile(value, 1/4)*100,Q3=quantile(value, 3/4)*100, IQR = IQR(value)) %>% as.data.frame
      pval_out <<- pval
      mean_sd <<- molten %>% dplyr::group_by(variable, predictor2) %>% dplyr::summarize( N = length(value),mean = mean(value)*100,sd=sd(value)*100) %>% as.data.frame
    }
  }
  if(bar_chart==FALSE){
    if(ncol(tax)>=6) molten$value <- molten$value+1e-6 #add pseudocount for log scale 0;
    else   molten$value <- molten$value+0.001 #add pseudocount for log scale 0;
    ordered <- levels(factor(molten$colvar))
    p <- ggplot(molten, aes(x=variable, y=value, fill=predictor2)) +
      {if(violin){geom_violin(scale = violin_scale,width = 0.65, position=position_dodge(width=0.9),size=1, color="#00000000")} else {geom_boxplot(width = 0.55, position=position_dodge(width=0.8),size=0.3,outlier.size = 0,outlier.color = "grey")}}+
      {if(violin){stat_summary(fun=median, fun.min = min, fun.max = max, geom="point", size=0.8, color="black", position=position_dodge(width=0.9))} else {stat_summary(fun=median, fun.min = min, fun.max = max, geom="point", size=0.8, color="#00000000", position=position_dodge(width=0.9))}}+ theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),legend.text=element_text(size=12),legend.key.size = unit(0.5, "cm"))+ coord_flip() +xlab(NULL)+ylab(xlabs)+ggtitle(main)
    if(length(unique(molten$variable))>1) p <- p+ geom_vline(xintercept=seq(1.5, length(unique(molten$variable))-0.5, 1),lwd=0.2, colour="grey")

    p <- p +  scale_fill_manual(values =cols,labels=legend_names) + guides(fill = guide_legend(title=legend_title, reverse = TRUE,override.aes = list(linetype=0, shape=16,color=rev(cols),size=5, bg="white")))

  }
  if(bar_chart==TRUE){
    if(bar_chart_stacked==TRUE)
      p <-  ggplot(molten_mean,aes(x=factor(predictor2,labels=legend_names),y=value, fill=variable)) + theme_bw()+geom_bar(stat="identity")+ theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),axis.title=element_text(size=14),legend.text=element_text(size=12), axis.text = element_text(size = 12),strip.text = element_text(size = 12),legend.key.size = unit(0.5, "cm"),text=element_text(size=12)) +xlab(NULL)+ylab(ylabs)+ggtitle(main) +  scale_fill_manual(values =cols,labels=ordered) + guides(fill = guide_legend(title=NULL))
    if(bar_chart_stacked==FALSE){
      if(!is.null(color_by)) p <-   ggplot(molten_mean,aes(x=variable,y=value, fill=colvar,group=wrap))+geom_bar(stat="identity", position = position_dodge(width = 0.95))+ scale_fill_manual(values =cols,labels=ordered2)+ guides(fill = guide_legend(title=color_by))
      else {
        p <-   ggplot(molten_mean,aes(x=variable,y=value, fill=predictor2))+geom_bar(stat="identity", position = position_dodge(width = 0.95))+ scale_fill_manual(values =cols,labels=legend_names)+ guides(fill = guide_legend(title=legend_title))
      }
      if(length(unique(molten_mean$variable))>1) p <- p+ geom_vline(xintercept=seq(1.5, length(unique(molten_mean$variable))-0.5, 1),lwd=0.2, colour="grey")
      p <-  p+ theme_bw()  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),axis.title=element_text(size=14),legend.text=element_text(size=12), axis.text = element_text(size = 12),strip.text = element_text(size = 12),legend.key.size = unit(0.5, "cm"),text=element_text(size=12)) +xlab(NULL)+ylab(ylabs)+ggtitle(main)+ theme(strip.background = element_blank()) +coord_flip()
    }
  }
  if(!is.null(facet_wrap))   {
    if(is.null(facet_label)) label_names <- levels(factor(samp[,facet_wrap]))
    if(!is.null(facet_label)) label_names <- facet_label
    if(facet_n==TRUE){
      label_names <- samp2 %>%
        dplyr::group_by(get(facet_wrap)) %>%
        dplyr::summarise(n = n()) %>%
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
  if(italic_names==TRUE &  (bar_chart==FALSE | (bar_chart==TRUE & bar_chart_stacked==FALSE)))   p <- p+ theme(axis.text.y=element_text(face = "italic"))
  if(!is.null(color_by)) {
    # p <- p + facet_grid(~predictor2, scales = "free", space = "free")
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

  p <-  p+ theme(plot.background = element_blank(),panel.background = element_blank(),plot.title = element_text(hjust = 0.5))
  if (bar_chart==TRUE & bar_chart_stacked==FALSE & percent==TRUE)  p <- p+  geom_text(aes(label = paste0(sprintf("%.2f",value*100), "%")), hjust = -.12, position=position_dodge(width=0.95))+scale_y_continuous(limits=c(0,max(molten_mean$value)+0.2),labels = scales::percent)
  if(no_legends) p <- p + theme(legend.position="none")
  if(no_names)  p <- p + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  stars.pval <- function (p.value)
  {    unclass(symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**",  "*", "NS")))
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
