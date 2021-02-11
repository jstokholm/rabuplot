#' Relative Abundance Plot
#'
#' For creating nice microbiome plots
#'
#' @param phylo_ob Phyloseq object with metadata in sample_data.
#' @param predictor Predictor of interestfor stastics/plotting in sample_data.
#' @param type Taxonomic rank from tax_table, case insensitive; default is "genus".
#' @param relative_abun Use relative abundances, else absolute; default is TRUE.
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
#' @param log_max Maximum value of log-axis; default is 100.
#' @param stat_out Outputs a data.frame with statistics to Global environment; default is FALSE.
#' @param p_val Displays p-values on plot; default is TRUE.
#' @param p_stars Shows stars instead of p-values; default is FALSE.
#' @param stats Select type of statistical test; options "non-parametric", "parametric", "mgs_feature"; defaule is "non-parametric".
#' @param p_adjust FDR adjust p-values; default is FALSE.
#' @param colors define list of colors for plot. If not color brewer will be used; default is NULL.
#' @param color_by define taxonomic rank to color by; default is NULL.
#' @param order Order by abundance, else alphabetically; default is TRUE.
#' @param reverse Flip taxa order; default is FALSE.
#' @param list_taxa A list of specific taxa names to be analyzed; default is NULL.
#' @param list_type Taxonomic rank of the @list_taxa; default is NULL.
#' @param select_taxa Choose all taxa from a Taxonomic variable, eg. "Staphylocuccus" or "Staph" or "cuccus"; default is NULL.
#' @param select_type Taxonomic rank of the @select_taxa; default is NULL.
#' @param bar_chart Choose to make bar chart; default is FALSE.
#' @param bar_chart_stacked Produce stacked bar chart, when @bar_chart is TRUE; default is TRUE.
#' @param percent Print percentages on bar chart; default is FALSE.
#' @param facet_wrap Facet wrap chart by variable; eg. Time; default is NULL.
#' @param order_by Choose variable to order the selected taxa by; eg. Time; default is Time.
#' @param order_val Choose value for @order_by; default is NULL.
#'
#' @import ggplot2 phyloseq metagenomeSeq dplyr tidyr RColorBrewer
#' @return A ggplot
#' @export

rabuplot <- function(phylo_ob,
                     predictor,
                     type="genus",
                     relative_abun=TRUE,
                     xlabs = "Relative abundance (%)",
                     ylabs = "Average relative abundance",
                     main = "Relative abundance plot",
                     violin=TRUE,
                     violin_scale = "width",
                     legend_title=predictor,
                     N_taxa=15,
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
                     p_adjust = FALSE,
                     colors=NULL,
                     color_by=NULL,
                     order=TRUE,
                     reverse=FALSE,
                     list_taxa=NULL,
                     list_type=NULL,
                     select_taxa=NULL,
                     select_type=NULL,
                     bar_chart=FALSE,
                     bar_chart_stacked=TRUE,
                     facet_wrap=NULL,
                     percent=FALSE,
                     order_by="Time",
                     order_val=NULL)
{

  phylo_ob <- prune_samples(sample_sums(phylo_ob)>0,phylo_ob) #removes empty samples;
  otu_mat <- t(as(otu_table(phylo_ob), "matrix"))
  if(!is.null(facet_wrap)) index <- !is.na(get_variable(phylo_ob, predictor)) & !is.na(get_variable(phylo_ob, facet_wrap))
  else   index <- !is.na(get_variable(phylo_ob, predictor))
  if(length(unique(index)) !=1) message("NAs have been removed for predictor/facet_wrap variable(s)")
  otu_mat <- otu_mat[index,]
  otu_mat  <- otu_mat[,colSums(otu_mat)>0] #removes empty OTUs;
  OTU_index <- colnames(otu_mat)
  tax <- as(tax_table(phylo_ob), "matrix") %>% data.frame(stringsAsFactors=FALSE)
  tax <- tax[rownames(tax) %in% OTU_index,]
  tax[is.na(tax)] <- "unclassified"
  names(tax) <- tolower(names(tax))
  type <- tolower(type)
  tax$OTU <- rownames(tax)
  samp <- data.frame(sample_data(phylo_ob), stringsAsFactors=TRUE)
  samp <- samp[index,]

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

  if(stats=="mgs_feature" & p_val==TRUE){
    if(!is.null(Strata))  {
      samp2 <- samp[which(samp[,Strata]==Strata_val),]
      abund2 <- abund[which(samp[,Strata]==Strata_val),]
      pred <- samp2[,predictor]
    }
    else  {
      samp2 <- samp
      abund2 <- abund
      pred <- samp2[,predictor]
    }
    if(length(levels(factor(pred)))>2){
      stats="non-parametric"
      message("MGS not available for >2 predictors, switching to non-parametric")
    }
    if(length(levels(factor(pred)))==2 & is.null(facet_wrap)){
      # test with featureModel
      mgs <- metagenomeSeq::newMRexperiment(counts = t(abund2))
      mgsp <- metagenomeSeq::cumNormStat(mgs)
      mgs <- metagenomeSeq::cumNorm(mgs, mgsp)
      mod <- model.matrix(~as.numeric(pred == unique(pred)[1]))
      message("MGS FeatureModel")
      mgsfit <- metagenomeSeq::fitFeatureModel(obj=mgs,mod=mod)
      mgs_pvalues <- data.frame(variable=mgsfit$taxa,pvalues=mgsfit$pvalues) %>% mutate(p_adjust=p.adjust(pvalues, "fdr"))
    }
    if(length(levels(factor(pred)))==2 & !is.null(facet_wrap)){
      mgs_pvalues <- data.frame()
      for (i in 1:length(unique(samp2[,facet_wrap]))){
        index <- samp2[,facet_wrap]==unique(samp2[,facet_wrap])[[i]]
        samp3 <- samp2[index,]
        abund3 <- abund2[index,]
        pred <- samp3[,predictor]
        # test with featureModel
        mgs <- metagenomeSeq::newMRexperiment(counts = t(abund3))
        mgsp <- metagenomeSeq::cumNormStat(mgs)
        mgs <- metagenomeSeq::cumNorm(mgs, mgsp)
        mod <- model.matrix(~as.numeric(pred == unique(pred)[1]))
        message(paste0("MGS FeatureModel for facet_wrap = ",unique(samp2[,facet_wrap])[[i]]))
        mgsfit <- metagenomeSeq::fitFeatureModel(obj=mgs,mod=mod)
        pval_tmp <- data.frame(variable=mgsfit$taxa,pvalues=mgsfit$pvalues,wrap=unique(samp2[,facet_wrap])[[i]]) %>% mutate(p_adjust=p.adjust(pvalues, "fdr"))
        mgs_pvalues <- rbind(mgs_pvalues,pval_tmp)
      }
    }
  }
  #Sort by abundance
  abund$"reads" <- rowSums(abund)
  if(relative_abun){
    for(i in names(abund)){
      abund[,i] <- abund[,i] / abund$"reads"
    }
  }
  abund$"reads" <- NULL

  if (is.null(list_taxa) & !is.null(select_taxa))  list_taxa <- as.character(unique(tax[grep(select_taxa,tax[,select_type],ignore.case=TRUE),type]))

  if (!is.null(list_taxa)) {
    abund <- abund[,colnames(abund) %in% list_taxa, drop = FALSE]
    unique_tax <- names(abund)
  }
  index <- !is.na(rownames(samp))
  if(length(abund)>1){
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
  bacteria <- rev(names(abund))
  subset <- cbind(samp[!names(samp) %in% bacteria], abund) #fjerner evt eksisterende navne fra dataset og merger;
  subset$predictor2 <-  as.factor(subset[,predictor])
  subset$ID <- rownames(subset)
  if(!is.null(Strata)) subset[,Strata] <- as.factor(subset[,Strata])
  if(!is.null(facet_wrap)){
    subset$wrap <-  as.factor(subset[,facet_wrap])
    if(!is.null(Strata))
      molten <- subset[,c("ID",paste(bacteria),"predictor2",Strata,"wrap")] %>% gather(variable, value,-"predictor2",-"ID",-Strata,-"wrap")
    else
      molten <- subset[,c("ID",paste(bacteria),"predictor2","wrap")] %>% gather(variable, value,-"predictor2",-"ID",-"wrap")
  }
  if(is.null(facet_wrap)){
    if(!is.null(Strata))
      molten <- subset[,c("ID",paste(bacteria),"predictor2",Strata)] %>% gather(variable, value,-"predictor2",-"ID",-Strata)
    else
      molten <- subset[,c("ID",paste(bacteria),"predictor2")] %>% gather(variable, value,-"predictor2",-"ID")
  }
  if(!is.null(color_by)){
    molten[molten$variable != paste("Other",type),"colvar"] <- molten %>% dplyr::filter(variable != paste("Other",type)) %>% .[,"variable"] %>% match(tax[,type]) %>% tax[.,"phylum"] %>% as.character
    molten[molten$variable == paste("Other",type),"colvar"] <- "Other" %>% as.character
  }

  molten$variable <- gsub('_',' ',molten$variable)

  if(order)   ordered<- unique(molten$variable) #level order
  if(!order)   ordered<-sort(unique(molten$variable))#level order alphabetically

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
    if(!is.null(facet_wrap)){
      molten_mean <- aggregate(molten$value,by=list(molten$variable, molten$predictor2, molten$wrap, molten$colvar),FUN=mean)
      names(molten_mean) <- c("type", "predictor2","wrap","colvar", "value")
    }
    if(is.null(facet_wrap)){
      molten_mean <- aggregate(molten$value,by=list(molten$variable, molten$predictor2,molten$colvar),FUN=mean)
      names(molten_mean) <- c("type", "predictor2","colvar", "value")
      molten_mean$wrap <- ""
    }
    molten_mean$colvar <- factor(molten_mean$colvar, levels=ordered2)

  }

  #Calculate pvalue for outcomes
  if(p_val==TRUE & ((bar_chart==TRUE & bar_chart_stacked==FALSE) | bar_chart==FALSE) & is.null(color_by)){
    if(is.null(facet_wrap)) molten$wrap <- ""
    if(stats=="mgs_feature"){
      if(!is.null(facet_wrap)) {
        pval <- data.frame(y=ifelse(log_max==100,26,ifelse(log_max==10,0.126,0.0126)), pval=mgs_pvalues[gsub('_',' ',mgs_pvalues$variable) %in% ordered,]$pvalues,p_adjust=mgs_pvalues[gsub('_',' ',mgs_pvalues$variable) %in% ordered,]$p_adjust, variable=gsub('_',' ',mgs_pvalues[gsub('_',' ',mgs_pvalues$variable) %in% ordered,]$variable),wrap=mgs_pvalues[gsub('_',' ',mgs_pvalues$variable) %in% ordered,]$wrap)
      }
      else {
        pval <- data.frame(y=ifelse(log_max==100,26,ifelse(log_max==10,0.126,0.0126)), pval=mgs_pvalues[gsub('_',' ',mgs_pvalues$variable) %in% ordered,]$pvalues,p_adjust=mgs_pvalues[gsub('_',' ',mgs_pvalues$variable) %in% ordered,]$p_adjust, variable=gsub('_',' ',mgs_pvalues[gsub('_',' ',mgs_pvalues$variable) %in% ordered,]$variable))
        if(length(pval$variable)-length(ordered)<0) pval <- pval[match(pval$variable,ordered[length(pval$variable)-length(ordered)]),]
      }
    }
    if(stats=="non-parametric"){
      if(length(levels(subset$predictor2))>2) {#Kruskal-Wallis for more than 2 groups
        pval <- data.frame()
        for (i in 1:length(unique(molten$wrap)))
        {
          test <- molten[molten$wrap==unique(molten$wrap)[[i]],]
          pval_tmp<- data.frame(y=ifelse(log_max==100,26,ifelse(log_max==10,0.126,0.0126)), pval=sapply(split(test, test$variable) , function(x) kruskal.test(value ~ predictor2, x)$p.value),variable=factor(paste(ordered)),wrap=unique(test$wrap))%>% mutate(p_adjust=p.adjust(pval, "fdr"))
          pval <- rbind(pval,pval_tmp)
        }
      }
      if(length(levels(subset$predictor2))==2) {  #Wilcoxon
        pval <- data.frame()
        for (i in 1:length(unique(molten$wrap)))
        {
          test <- molten[molten$wrap==unique(molten$wrap)[[i]],]
          pval_tmp<- data.frame(y=ifelse(log_max==100,26,ifelse(log_max==10,0.126,0.0126)), pval=sapply(split(test, test$variable) , function(x) wilcox.test(value ~ predictor2, x)$p.value),variable=factor(paste(ordered)),wrap=unique(test$wrap))%>% mutate(p_adjust=p.adjust(pval, "fdr"))
          pval <- rbind(pval,pval_tmp)
        }
      }
      message("Non-parametric")
    }
    if(stats=="parametric") {
      pval <- data.frame()
      for (i in 1:length(unique(molten$wrap)))
      {
        test <- molten[molten$wrap==unique(molten$wrap)[[i]],]
        pval_tmp<- data.frame(y=ifelse(log_max==100,26,ifelse(log_max==10,0.126,0.0126)), pval=sapply(split(test, test$variable) , function(x) oneway.test(value ~ predictor2, x)$p.value),variable=factor(paste(ordered)),wrap=unique(test$wrap)) %>% mutate(p_adjust=p.adjust(pval, "fdr"))
        pval <- rbind(pval,pval_tmp)
      }
      message("Parametric")
    }
    pval$predictor2 <- molten$predictor2[1]
    #  pval$pval = p.adjust(pval$pval, "fdr")
    pval$pval <- ifelse(is.na(pval$pval),1,pval$pval)
    pval$p_adjust <- ifelse(is.na(pval$p_adjust),1,pval$p_adjust)
    if(Only_sig){
      index <- rownames(pval[pval$pval<0.05,])
      molten <- molten[molten$variable %in% index,]
      pval <- pval[pval$pval<0.05,]
    }
    if (!is.null(list_type)){
      molten <- molten[molten$variable  %in% list_type,]
      pval <-    pval[pval$variable  %in% list_type,]
    }

    if(stat_out){
      median_iqr <<- molten %>% dplyr::group_by(variable, predictor2) %>% dplyr::summarize( N = length(value),median = median(value)*100,Q1=quantile(value, 1/4)*100,Q3=quantile(value, 3/4)*100, IQR = IQR(value)) %>% as.data.frame
      pval_out <<- pval
      mean_sd <<- molten %>% dplyr::group_by(variable, predictor2) %>% dplyr::summarize( N = length(value),mean = mean(value)*100,sd=sd(value)*100) %>% as.data.frame
    }
    if(log==FALSE & bar_chart==FALSE) pval$y <- max(molten$value)*1.10
    if(log==FALSE & bar_chart==TRUE) pval$y <- max(molten_mean$value)*1.20
    stars.pval <- function (p.value)
    {    unclass(symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**",  "*", "NS")))
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
    if(p_stars==TRUE & p_val==TRUE)
      p <- p + annotate("text",size=3,hjust=1,fontface="plain", x = pval$variable, y = pval$y,  label=paste(stars.pval(pval$pval)))

    p <- p +  scale_fill_manual(values =cols,labels=legend_names) + guides(fill = guide_legend(title=legend_title, reverse = TRUE,override.aes = list(linetype=0, shape=16,color=rev(cols),size=5, bg="white")))

  }
  if(bar_chart==TRUE){
    #   ordered  <- paste0(ordered, '   ')  #For legends in one line.
    #    ordered[length(ordered)]  <- paste0(ordered[length(ordered)], '             ')
    if(bar_chart_stacked==TRUE)
      p <-  ggplot(molten_mean,aes(x=factor(predictor2,labels=legend_names),y=value, fill=type)) + theme_bw()+geom_bar(stat="identity")+ theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),axis.title=element_text(size=14),legend.text=element_text(size=12), axis.text = element_text(size = 12),strip.text = element_text(size = 12),legend.key.size = unit(0.5, "cm"),text=element_text(size=12)) +xlab(NULL)+ylab(ylabs)+ggtitle(main) +  scale_fill_manual(values =cols,labels=ordered) + guides(fill = guide_legend(title=legend_title))+ scale_y_continuous(labels = scales::percent)
    if(bar_chart_stacked==FALSE){
      if(!is.null(color_by)) p <-   ggplot(molten_mean,aes(x=type,y=value, fill=colvar,group=wrap))+geom_bar(stat="identity", position = position_dodge(width = 0.95))+ scale_fill_manual(values =cols,labels=ordered2)
      else {
        p <-   ggplot(molten_mean,aes(x=type,y=value, fill=predictor2))+geom_bar(stat="identity", position = position_dodge(width = 0.95))+ scale_fill_manual(values =cols,labels=legend_names)
        if(p_stars==TRUE & p_val==TRUE)
          p <- p + annotate("text",size=3,hjust=1,fontface="plain", x = pval$variable, y = pval$y,  label=paste(stars.pval(pval$pval)))

      }
      p <-  p+ theme_bw()  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),axis.title=element_text(size=14),legend.text=element_text(size=12), axis.text = element_text(size = 12),strip.text = element_text(size = 12),legend.key.size = unit(0.5, "cm"),text=element_text(size=12)) +xlab(NULL)+ylab(ylabs)+ggtitle(main)+ guides(fill = guide_legend(title=legend_title)) + theme(strip.background = element_blank()) +coord_flip()

      if(is.null(color_by) & is.null(facet_wrap) & bar_chart_stacked==TRUE) p <- p+theme(legend.position="none")
    }
  }
  if(!is.null(facet_wrap))   { p <- p+ facet_grid(~wrap,scales = "free", space = "free")+ theme(strip.background = element_blank())
  if(bar_chart==FALSE) p$layers[4:5] <- NULL
  }
  if(italic_names==TRUE &  (bar_chart==FALSE | (bar_chart==TRUE & bar_chart_stacked==FALSE)))   p <- p+ theme(axis.text.y=element_text(face = "italic"))
  #  if(!is.null(color_by) & (color_by=="genus" | color_by=="family" | color_by=="species"))
  if(!is.null(color_by)) {
    p <- p + facet_grid(~predictor2, scales = "free", space = "free")
    if(color_by=="genus" | color_by=="family" | color_by=="species") p <- p+ theme(legend.text=element_text(face = "italic"))
  }
  if(log){
    if(log_max == 100)  p <- p+ scale_y_log10(breaks=c(.000001,.001,.01,.1,1),labels=c("0%","0.1%","1%","10%","100%"))
    else{
      if(log_max == 10)  p <- p+ scale_y_log10(limits=c(0.001,0.13),breaks=c(.001,.01,.05,.1),labels=c("0%","1%","5%","10%"))
      else  p <- p+ scale_y_log10(limits=c(0.001,0.013),breaks=c(.001,.01),labels=c("0%","1%"))
    }
  }
  p <-  p+ theme(plot.background = element_blank(),panel.background = element_blank(),plot.title = element_text(hjust = 0.5))
  if (bar_chart==TRUE & bar_chart_stacked==FALSE  & percent==FALSE) p <- p+ scale_y_continuous(labels = scales::percent)
  if (bar_chart==TRUE & bar_chart_stacked==FALSE & percent==TRUE) p <- p+  geom_text(aes(label = paste0(sprintf("%.2f",value*100), "%")), hjust = -.12, position=position_dodge(width=0.95))+scale_y_continuous(limits=c(0,max(molten_mean$value)+0.2),labels = scales::percent)
  if(no_legends) p <- p + theme(legend.position="none")
  if(no_names)  p <- p + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  # if(!is.null(facet_wrap)) p + guides(fill = guide_legend(title="legend_title", reverse = F))
  if(p_stars==FALSE & p_val==TRUE & (bar_chart==FALSE | (bar_chart==TRUE & bar_chart_stacked==FALSE)))
    p <- p + geom_text(data=pval,aes(x=variable,y=y,label=ifelse(pval<0.05, paste(format.pval(pval,1,0.001,nsmall=3)),"")) ,size=3,hjust=1,fontface="bold")
  p <- p + geom_text(data=pval,aes(x=variable,y=y,label=ifelse(pval>=0.05, paste(format.pval(pval,1,0.001,nsmall=3)),"")) ,size=3,hjust=1)
  if(p_adjust){
    if(stats=="mgs_feature") message(paste("FDR correction applied for",length(unique(mgs_pvalues$variable)),"taxa"))
    else  message(paste("FDR correction applied for",length(unique(pval$variable)),"taxa"))
    p <- p + scale_y_log10(breaks=c(.000001,.001,.01,.1,1,10,100),labels=c("0%","0.1%","1%","10%","100%", "P-value", "q-value"))
    p <- p + geom_text(data=pval,aes(x=variable,y=100,label=ifelse(p_adjust<0.05, paste(format.pval(p_adjust,1,0.001,nsmall=3)),"")) ,size=3,hjust=.5,fontface="bold")
    p <- p + geom_text(data=pval,aes(x=variable,y=100,label=ifelse(p_adjust>=0.05, paste(format.pval(p_adjust,1,0.001,nsmall=3)),"")) ,size=3,hjust=.5)
  }

  p
}
