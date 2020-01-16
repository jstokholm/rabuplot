#'Relative Abundance Plot
#'
#'For creating plots
#'
#'
#'@import ggplot2 phyloseq metagenomeSeq dplyr RColorBrewer
#'
#'@export
#'
#'@return Plot
#'

rabuplot <- function(phylo_ob,
                                predictor,
                                type="genus",
                                xlabs = "Proportional abundance (%)",
                                ylabs = "Average relative abundance",
                                main = "Proportional abundance plot",
                                violin=TRUE,
                                violin_scale = "width",
                                legend_title=predictor,
                                N_taxa=15,
                                No_other_type=FALSE,
                                legend_names=NULL,
                                Time="Time",
                                Timepoint=NULL,
                                Strata=NULL,
                                Strata_val="1",
                                No_legends = FALSE,
                                No_names=FALSE,
                                By_median=TRUE,
                                Only_sig=FALSE,
                                log=TRUE,
                                log_max=100,
                                stat_out=FALSE,
                                p_val = TRUE,
                                colors=NULL,
                                order=TRUE,
                                reverse=FALSE,
                                list_OTU=NULL,
                                list_type=NULL,
                                select_taxa=NULL,
                                select_type=NULL,
                                relative_abun=TRUE,
                                bar_chart=FALSE,
                                bar_wrap=NULL,
                                bar_chart_stacked=TRUE,
                                color_by=NULL,
                                percent=FALSE,
                                order_by="Time",
                                order_val=NULL,
                                group=NULL,
                                p_stars=FALSE,
                                remove_collapsed_taxa=FALSE,
                                stats="non-parametric",
                                p_adjust = FALSE)
{

  otu_mat <- t(as(otu_table(phylo_ob), "matrix"))
  otu_mat  <- otu_mat[,colSums(otu_mat)>1] #removes OTUs <1;
  index <- !is.na(get_variable(phylo_ob, predictor))
  otu_mat <- otu_mat[index,]
  OTU_index <- colnames(otu_mat)
  tax <- as(tax_table(phylo_ob), "matrix") %>% data.frame
  tax <- tax[rownames(tax) %in% OTU_index,]
  tax[is.na(tax)] <- as.factor("unclassified")
  org_tax <- names(tax)
  names(tax) <- tolower(names(tax))
  tax$OTU <- rownames(tax)
  samp <- data.frame(sample_data(phylo_ob), stringsAsFactors=TRUE)
  samp <- samp[index,]

  if(!is.null(Timepoint)){
    index <- rownames(samp[(samp[,Time] ==Timepoint),])
    otu_mat <- otu_mat[rownames(otu_mat) %in% index,]
    otu_mat  <- otu_mat[,colSums(otu_mat)>0] #fjerner tomme OTUs;
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
    if(length(levels(factor(pred)))==2){
      # test with featureModel

      mgs <- newMRexperiment(counts = t(abund2))
      mgsp <- cumNormStat(mgs)
      mgs <- cumNorm(mgs, mgsp)
      mod <- model.matrix(~as.numeric(pred == unique(pred)[1]))
      message("MGS FeatureModel")
      mgsfit <- metagenomeSeq::fitFeatureModel(obj=mgs,mod=mod)
      mgs_pvalues <- data.frame(variable=mgsfit$taxa,pvalues=mgsfit$pvalues)
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


  if (is.null(list_OTU) & !is.null(select_taxa))  list_OTU <- tax[tax[,select_type] %in% select_taxa,type]
  if (!is.null(list_OTU)) {
    abund <- abund[,colnames(abund) %in% list_OTU, drop = FALSE]
    unique_tax <- names(abund)
  }
  if(remove_collapsed_taxa){
    for(i in org_tax) {
      unique_tax <- unique_tax[!grepl(paste(i), unique_tax)]
    }
    unique_tax <- names(abund)
    abund1 <- abund[,(colnames(abund) %in% unique_tax), drop = FALSE]
    '%!in%' <- function(x,y)!('%in%'(x,y))
    unclassified <- data.frame(unclassified=rowSums(abund[,(colnames(abund) %!in% unique_tax), drop = FALSE]))
    abund <- cbind(abund1,unclassified)
    unique_tax <- c(unique_tax,"unclassified")
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

    if(No_other_type)  abund[,paste("Other",type)] <- NULL
  }
  bacteria <- rev(names(abund))
  subset <- cbind(samp[!names(samp) %in% bacteria], abund) #fjerner evt eksisterende navne fra dataset og merger;
  subset$predictor2 <-  as.factor(subset[,predictor])
  subset$ID <- rownames(subset)
  if(!is.null(Strata)) subset[,Strata] <- as.factor(subset[,Strata])
  if(!is.null(bar_wrap)){
    subset$wrap <-  as.factor(subset[,bar_wrap])
    molten <- melt(subset[,c("ID",paste(bacteria),"predictor2",Strata,"wrap")], )
  }
  if(is.null(bar_wrap)){
    #if (!is.null(group)) subset[,group] <-  as.factor(subset[,group])
    molten <- melt(subset[,c("ID",paste(bacteria),"predictor2",Strata)], )
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
    if(is.null(color_by) & reverse==FALSE) cols <- rev(cols)
    if(!is.null(color_by) & reverse==TRUE) cols <- rev(cols)
    if(!is.null(bar_wrap)){
      molten_mean <- aggregate(molten$value,by=list(molten$variable, molten$predictor2, molten$wrap, molten$colvar),FUN=mean)
      names(molten_mean) <- c("type", "predictor2","wrap","colvar", "value")
    }
    # if(!is.null(group)){
    #   molten_mean <- aggregate(molten$value,by=list(molten$variable, molten$predictor2,molten$colvar, molten[,group]),FUN=mean)
    #   names(molten_mean) <- c("type", "predictor2","colvar","group", "value")
    # }
    if(is.null(bar_wrap)){
      molten_mean <- aggregate(molten$value,by=list(molten$variable, molten$predictor2,molten$colvar),FUN=mean)
      names(molten_mean) <- c("type", "predictor2","colvar", "value")
      molten_mean$wrap <- ""
    }
    molten_mean$colvar <- factor(molten_mean$colvar, levels=ordered2)

  }


  #Calculate pvalue for outcomes
  if(p_val==TRUE & ((bar_chart==TRUE & bar_chart_stacked==FALSE) | bar_chart==FALSE) & is.null(color_by)){
    if(stats=="mgs_feature"){
      pval <- data.frame(y=ifelse(log_max==100,26,ifelse(log_max==10,0.126,0.0126)), pval=mgs_pvalues[gsub('_',' ',mgs_pvalues$variable) %in% ordered,]$pvalues, variable=gsub('_',' ',mgs_pvalues[gsub('_',' ',mgs_pvalues$variable) %in% ordered,]$variable))
      if(length(pval$variable)-length(ordered)<0) pval <- pval[match(pval$variable,ordered[length(pval$variable)-length(ordered)]),]
      else pval <- pval[match(pval$variable,ordered),]
    }
    if(stats=="non-parametric"){
      if(length(levels(subset$predictor2))>2) #Kruskal-Wallis for more than 2 groups
        pval <- data.frame(y=ifelse(log_max==100,26,ifelse(log_max==10,0.126,0.0126)), pval=sapply(split(molten, molten$variable), function(x) kruskal.test(value ~ predictor2, x)$p.value), variable=factor(paste(ordered)))
      if(length(levels(subset$predictor2))==2)   #Wilcoxon
        pval <- data.frame(y=ifelse(log_max==100,26,ifelse(log_max==10,0.126,0.0126)), pval=sapply(split(molten, molten$variable), function(x) wilcox.test(value ~ predictor2, x)$p.value), variable=factor(paste(ordered)))
      message("Non-parametric")
    }
    if(stats=="parametric") {
      pval <- data.frame(y=ifelse(log_max==100,26,ifelse(log_max==10,0.126,0.0126)), pval=sapply(split(molten, molten$variable), function(x) oneway.test(value ~ predictor2, x)$p.value), variable=factor(paste(ordered)))
      message("Parametric")
    }
    pval <<- pval
    #  pval$pval = p.adjust(pval$pval, "fdr")
    pval$pval <- ifelse(is.na(pval$pval),1,pval$pval)
    pval$pval_sig <- ifelse(pval$pval<0.05,pval$pval,NA)
    pval$pval_notsig <- ifelse(pval$pval>=0.05,pval$pval,NA)
    if(Only_sig){
      index <- rownames(pval[is.na(pval$pval_notsig),])
      molten <- molten[molten$variable %in% index,]
      pval <-   pval[is.na(pval$pval_notsig),]
    }
    #  if (!is.null(list_OTU)){
    #   list_tax <- as.character(tax[rownames(tax) %in% list_OTU,type])
    #  molten <- molten[molten$variable  %in% list_tax,]
    #    pval <-    pval[pval$variable  %in% list_tax,]
    #  }
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
    {if(violin){stat_summary(fun.y=median, fun.ymin = min, fun.ymax = max, geom="point", size=0.8, color="black", position=position_dodge(width=0.9))} else {stat_summary(fun.y=median, fun.ymin = min, fun.ymax = max, geom="point", size=0.8, color="#00000000", position=position_dodge(width=0.9))}}+ theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),legend.text=element_text(size=12),legend.key.size = unit(0.5, "cm"))+ coord_flip() +xlab(NULL)+ylab(xlabs)+ggtitle(main)
    if(length(unique(molten$variable))>1) p <- p+ geom_vline(xintercept=seq(1.5, length(unique(molten$variable))-0.5, 1),lwd=0.2, colour="grey")
    if(p_stars==TRUE & p_val==TRUE)
      p <- p + annotate("text",size=3,hjust=1,fontface="plain", x = pval$variable, y = pval$y,  label=paste(stars.pval(pval$pval)))
    if(p_stars==FALSE & p_val==TRUE)
      p <- p + annotate("text",size=3,hjust=1,fontface="plain", x = pval$variable, y = pval$y,  label=ifelse(is.na(pval$pval_notsig), "",paste(format.pval(pval$pval_notsig,1,0.001,nsmall=3)))) + annotate("text",size=3,hjust=1,fontface="bold", x = pval$variable, y = pval$y,  label=ifelse(is.na(pval$pval_sig), "",paste(format.pval(pval$pval_sig,1,0.001,nsmall=3))))
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
        if(p_stars==FALSE & p_val==TRUE)
          p <- p + annotate("text",size=3,hjust=1,fontface="plain", x = pval$variable, y = pval$y,  label=ifelse(is.na(pval$pval_notsig), "",paste(format.pval(pval$pval_notsig,1,0.001,nsmall=3)))) + annotate("text",size=3,hjust=1,fontface="bold", x = pval$variable, y = pval$y,  label=ifelse(is.na(pval$pval_sig), "",paste(format.pval(pval$pval_sig,1,0.001,nsmall=3))))
      }
      p <-  p+ theme_bw()  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.key = element_blank(),axis.title=element_text(size=14),legend.text=element_text(size=12), axis.text = element_text(size = 12),strip.text = element_text(size = 12),legend.key.size = unit(0.5, "cm"),text=element_text(size=12)) +xlab(NULL)+ylab(ylabs)+ggtitle(main)+ guides(fill = guide_legend(title=legend_title)) + theme(strip.background = element_blank()) +coord_flip()



      if(is.null(color_by) & is.null(bar_wrap) & bar_chart_stacked==TRUE) p <- p+theme(legend.position="none")
    }
  }
  if(!is.null(bar_wrap))   { p <- p+ facet_grid(~wrap)+ theme(strip.background = element_blank())
  if(bar_chart==FALSE) p$layers[4:5] <- NULL
  }
  if((type=="genus" | type=="family" | type=="species") &  bar_chart==FALSE | (bar_chart==TRUE & bar_chart_stacked==FALSE))   p <- p+ theme(axis.text.y=element_text(face = "italic"))
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
  if(No_legends) p <- p + theme(legend.position="none")
  if(No_names)  p <- p + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  # if(!is.null(bar_wrap)) p + guides(fill = guide_legend(title="legend_title", reverse = F))
  if(p_adjust){

    p <- p + scale_y_log10(breaks=c(.000001,.001,.01,.1,1,10,100),labels=c("0%","0.1%","1%","10%","100%", "P-value", "q-value"))
    p <- p + annotate("text",size=3,hjust=.5,fontface="bold", x = pval$variable, y = 100,label=ifelse(p.adjust(pval$pval, "fdr") < 0.05, paste(format.pval(p.adjust(pval$pval, "fdr"),1,0.001,nsmall=3)), ""))
    p <- p + annotate("text",size=3,hjust=.5,fontface="plain", x = pval$variable, y = 100,label=ifelse(p.adjust(pval$pval, "fdr") >= 0.05, paste(format.pval(p.adjust(pval$pval, "fdr"),1,0.001,nsmall=3)), ""))
  }

  p
}
