# results_dir <- "/project2/guiming/xsun/proteome_alzheimer/14.brain_schwarz_ukbb_joint_newscript/brain_output_preharmo_zt/"
# analysis_id <- "brain"

locus_plot <- function(region_tag, xlim=NULL, return_table=F, focus=NULL, label_panel="TWAS", label_genes=NULL, label_pos=NULL, plot_eqtl=NULL, rerun_ctwas=F, rerun_load_only=F, legend_side="right", legend_panel="cTWAS", twas_ymax=NULL,draw_gene_track=T){
  region_tag1 <- unlist(strsplit(region_tag, "_"))[1]
  region_tag2 <- unlist(strsplit(region_tag, "_"))[2]
  
  #ctwas_res <- ctwas_res[!duplicated(ctwas_res$id),]  ############
  
  a <- ctwas_res[ctwas_res$region_tag==region_tag,]
  regionlist <- readRDS(paste0(results_dir, "/", analysis_id, "_ctwas.regionlist.RDS"))
  region <- regionlist[[as.numeric(region_tag1)]][[region_tag2]]
  
  R_snp_info <- do.call(rbind, lapply(region$regRDS, function(x){data.table::fread(paste0(tools::file_path_sans_ext(x), ".Rvar"))}))
  
  if (isTRUE(rerun_ctwas)){
    ld_exprfs <- paste0(results_dir, "/", analysis_id, "_expr_chr", 1:22, ".expr.gz")
    temp_reg <- data.frame("chr" = paste0("chr",region_tag1), "start" = region$start, "stop" = region$stop)
    
    write.table(temp_reg, 
                #file= paste0(results_dir, "/", analysis_id, "_ctwas.temp.reg.txt") , 
                file= "temp_reg.txt",
                row.names=F, col.names=T, sep="\t", quote = F)
    
    load(paste0(results_dir, "/", analysis_id, "_expr_z_snp.Rd"))
    
    z_gene_temp <-  z_gene[z_gene$id %in% a$id[a$type=="gene"],]
    z_snp_temp <-  z_snp[z_snp$id %in% R_snp_info$id,]
    
    if (!rerun_load_only){
      ctwas::ctwas_rss(z_gene_temp, z_snp_temp, ld_exprfs, ld_pgenfs = NULL, 
                       ld_R_dir = dirname(region$regRDS)[1],
                       ld_regions_custom = "temp_reg.txt", thin = 1, 
                       outputdir = ".", outname = "temp", ncore = 1, ncore.rerun = 1, prob_single = 0,
                       group_prior = estimated_group_prior, group_prior_var = estimated_group_prior_var,
                       estimate_group_prior = F, estimate_group_prior_var = F)
    }
    
    a_bkup <- a         
    a <- as.data.frame(data.table::fread("temp.susieIrss.txt", header = T))
    
    rownames(z_snp_temp) <- z_snp_temp$id
    z_snp_temp <- z_snp_temp[a$id[a$type=="SNP"],]
    z_gene_temp <- z_gene_temp[a$id[a$type=="gene"],]
    
    a$genename <- NA
    a$gene_type <- NA
    
    a[a$type=="gene",c("genename", "gene_type")] <- a_bkup[match(a$id[a$type=="gene"], a_bkup$id),c("genename","gene_type")]
    
    a$z <- NA
    a$z[a$type=="SNP"] <- z_snp_temp$z
    a$z[a$type=="gene"] <- z_gene_temp$z
  }
  
  a_pos_bkup <- a$pos
  a$pos[a$type=="gene"] <- G_list$tss[match(sapply(a$genename[a$type=="gene"], function(x){unlist(strsplit(x, "[.]"))[1]}) ,G_list$hgnc_symbol)]
  a$pos[is.na(a$pos)] <- a_pos_bkup[is.na(a$pos)]
  a$pos <- a$pos/1000000
  
  if (!is.null(xlim)){
    
    if (is.na(xlim[1])){
      xlim[1] <- min(a$pos)
    }
    
    if (is.na(xlim[2])){
      xlim[2] <- max(a$pos)
    }
    
    a <- a[a$pos>=xlim[1] & a$pos<=xlim[2],,drop=F]
  }
  
  if (is.null(focus)){
    focus <- a$id[which.max(abs(a$z)[a$type=="gene"])]
  }
  
  if (is.null(label_genes)){
    label_genes <- focus
  }
  
  if (is.null(label_pos)){
    label_pos <- rep(3, length(label_genes))
  }
  
  # if (is.null(plot_eqtl)){
  #   plot_eqtl <- focus
  # }
  
  #### modified for qtl plot
  if (is.null(plot_eqtl)){
    gene_name_focus <- a$genename[a$id == focus]  
    plot_eqtl <- a$id[a$genename %in% gene_name_focus]
  }
  
  focus <- a$id[which(a$id==focus)]
  a$focus <- 0
  a$focus <- as.numeric(a$id==focus)
  
  a$PVALUE <- (-log(2) - pnorm(abs(a$z), lower.tail=F, log.p=T))/log(10)
  
  R_gene <- readRDS(region$R_g_file)
  R_snp_gene <- readRDS(region$R_sg_file)
  R_snp <- as.matrix(Matrix::bdiag(lapply(region$regRDS, readRDS)))
  
  rownames(R_gene) <- region$gid
  colnames(R_gene) <- region$gid
  rownames(R_snp_gene) <- R_snp_info$id
  colnames(R_snp_gene) <- region$gid
  rownames(R_snp) <- R_snp_info$id
  colnames(R_snp) <- R_snp_info$id
  
  a$r2max <- NA
  a$r2max[a$type=="gene"] <- R_gene[focus,a$id[a$type=="gene"]]
  a$r2max[a$type=="SNP"] <- R_snp_gene[a$id[a$type=="SNP"],focus]
  
  r2cut <- 0.4
  colorsall <- c("#7fc97f", "#beaed4", "#fdc086")
  
  start <- min(a$pos)
  end <- max(a$pos)
  
  #layout(matrix(1:3, ncol = 1), widths = 1, heights = c(1.5,1.75,0.75), respect = FALSE) without eqtl
  #### modified for qtl plot
  layout(matrix(1:(3+length(plot_eqtl)), ncol = 1), widths = 1, heights = c(1.5,rep(0.25,length(plot_eqtl)),1.75,0.75), respect = FALSE)
  
  par(mar = c(0, 4.1, 0, 2.1))
  
  if (is.null(twas_ymax)){
    twas_ymax <- max(a$PVALUE)*1.1
  }
  
  a <- a[!duplicated(a$id),]
  
  plot(a$pos[a$type=="SNP"], a$PVALUE[a$type == "SNP"], pch = 21, xlab=paste0("Chromosome ", region_tag1, " position (Mb)"), frame.plot=FALSE, bg = colorsall[1], ylab = "-log10(p value)", panel.first = grid(), ylim =c(0, twas_ymax), xaxt = 'n', xlim=c(start, end))
  
  abline(h=-log10(alpha/nrow(ctwas_gene_res)), col ="red", lty = 2)
  points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$PVALUE[a$type == "SNP"  & a$r2max > r2cut], pch = 21, bg = "purple")
  points(a$pos[a$type=="SNP" & a$focus == 1], a$PVALUE[a$type == "SNP" & a$focus == 1], pch = 21, bg = "salmon")
  points(a$pos[a$type=="gene" & a$group=="Expression"], a$PVALUE[a$type == "gene" & a$group=="Expression"], pch = 22, bg = colorsall[1], cex = 2)
  points(a$pos[a$type=="gene" & a$group=="Protein"], a$PVALUE[a$type == "gene" & a$group=="Protein"], pch = 23, bg = colorsall[1], cex = 2)
  
  points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Expression"], a$PVALUE[a$type == "gene"  & a$r2max > r2cut & a$group=="Expression"], pch = 22, bg = "purple", cex = 2)
  points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Protein"], a$PVALUE[a$type == "gene"  & a$r2max > r2cut & a$group=="Protein"], pch = 23, bg = "purple", cex = 2)
  
  points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Expression"], a$PVALUE[a$type == "gene" & a$focus == 1 & a$group=="Expression"], pch = 22, bg = "salmon", cex = 2)
  points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Protein"], a$PVALUE[a$type == "gene" & a$focus == 1 & a$group=="Protein"], pch = 23, bg = "salmon", cex = 2)
  
  if (legend_panel=="TWAS"){
    x_pos <- ifelse(legend_side=="right", max(a$pos)-0.2*(max(a$pos)-min(a$pos)), min(a$pos))
    legend(x_pos, y= twas_ymax*0.95, c("Expression","Protein","SNP","Lead TWAS Gene", "R2 > 0.4", "R2 <= 0.4"), pch = c(22,23,21,19,19,19), col = c("black","black","black", "salmon", "purple", colorsall[1]), cex=0.7, title.adj = 0)
  }
  
  
  label_genes <- a[a$id==label_genes,]$genename
  
  if (label_panel=="TWAS" | label_panel=="both"){
    for (i in 1:length(label_genes)){
      text(a$pos[a$genename %in% label_genes[i]], a$PVALUE[a$genename %in% label_genes[i]], labels=a$genename[a$genename %in% label_genes[i]], pos=label_pos[i], cex=0.7)
    }
  }
  
  #####for eqtl 
  
  #### modified for qtl plot
  type_plot <- a$group[a$id == plot_eqtl]
  
  eqtl_list<-list()
  for (i in 1:length(plot_eqtl)){
    par(mar = c(0.25, 4.1, 0.25, 2.1))
    
    # plot(NA, xlim = c(start, end), ylim = c(0, length(plot_eqtl)*length(all_cgenes)), 
    #      frame.plot = F, axes = F, xlab = NA, ylab = NA)
    
    plot(NA, xlim = c(start, end), ylim = c(0, 1), 
         frame.plot = F, axes = F, xlab = NA, ylab = NA)
    
    type_plot <- a$group[a$id == plot_eqtl[i]]
    
    cgene_list <- a[a$group ==  type_plot & a$id == plot_eqtl[i], ]$id
    load(paste0(results_dir, "/",analysis_id, "_expr_chr", region_tag1, ".exprqc.Rd"))
    
    for(k in 1:length(cgene_list)){
      cgene = cgene_list[k]
      #peak_id = a[a$id == cgene, ]$peak_id
      #print(peak_id)
      #peak_gr = strsplit(strsplit(peak_id, "[:]")[[1]][2], "-")[[1]]
      #cgene_start = as.numeric(peak_gr)[1]/1000000
      #cgene_end = as.numeric(peak_gr)[2]/1000000
      #col="grey"
      col="#c6e8f0"
      
      
      rect(xleft=start, ybottom =length(cgene_list)+1-k-0.8, 
           xright=end, ytop=length(cgene_list)+1-k-0.2, 
           col = col, border = T, lwd = 1)
      text(start, length(cgene_list)-k+0.5,  
           labels = paste0( type_plot,"QTL_", k), srt = 0, pos = 2, xpd = TRUE, cex=0.7)
      
      names(wgtlist)<-basename(names(wgtlist))
      eqtls <- rownames(wgtlist[[cgene]])
      eqtl_pos<-(bigSNP$map[bigSNP$map$marker.ID %in% eqtls, "physical.pos"])/1000000
      snps_to_plot<-wgtlist[[cgene]][eqtl_pos >= start & eqtl_pos <= end]
      eqtl_pos<-eqtl_pos[eqtl_pos >= start & eqtl_pos <= end]
      pos<-eqtl_pos*1000000
      eqtl_list[[cgene]]<-data.frame(snps_to_plot, 
                                     pos)
      #eqtl_pos <- a$pos[a$id %in% eqtls]
      
      if (length(eqtl_pos)>0){
        for (j in 1:length(eqtl_pos)){
          segments(x0=eqtl_pos[j], x1=eqtl_pos[j], y0=length(cgene_list)+1-k-0.2, length(cgene_list)+1-k-0.8, lwd=1.5)  
        }
      }
    }
  }
  
  
  
  ##### end eqtl
  
  
  
  
  par(mar = c(4.1, 4.1, 0, 2.1))
  
  plot(a$pos[a$type=="SNP"], a$susie_pip[a$type == "SNP"], pch = 19, xlab=paste0("Chromosome ", region_tag1, " position (Mb)"),frame.plot=FALSE, col = "white", ylim= c(0,1.1), ylab = "cTWAS PIP", xlim = c(start, end))
  
  grid()
  points(a$pos[a$type=="SNP"], a$susie_pip[a$type == "SNP"], pch = 21, xlab="Genomic position", bg = colorsall[1])
  points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$susie_pip[a$type == "SNP"  & a$r2max >r2cut], pch = 21, bg = "purple")
  points(a$pos[a$type=="SNP" & a$focus == 1], a$susie_pip[a$type == "SNP" & a$focus == 1], pch = 21, bg = "salmon")
  points(a$pos[a$type=="gene" & a$group=="Expression"], a$susie_pip[a$type == "gene" & a$group=="Expression"], pch = 22, bg = colorsall[1], cex = 2)
  points(a$pos[a$type=="gene" & a$group=="Protein"], a$susie_pip[a$type == "gene" & a$group=="Protein"], pch = 23, bg = colorsall[1], cex = 2)
  points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Expression"], a$susie_pip[a$type == "gene"  & a$r2max > r2cut & a$group=="Expression"], pch = 22, bg = "purple", cex = 2)
  points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Protein"], a$susie_pip[a$type == "gene"  & a$r2max > r2cut & a$group=="Protein"], pch = 23, bg = "purple", cex = 2)
  
  points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Expression"], a$susie_pip[a$type == "gene" & a$focus == 1 & a$group=="Expression"], pch = 22, bg = "salmon", cex = 2)
  points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Protein"], a$susie_pip[a$type == "gene" & a$focus == 1 & a$group=="Protein"], pch = 23, bg = "salmon", cex = 2)
  
  if (legend_panel=="cTWAS"){
    x_pos <- ifelse(legend_side=="right", max(a$pos)-0.2*(max(a$pos)-min(a$pos)), min(a$pos))
    legend(x_pos, y= 1 ,c("Expression","Protein","SNP","Lead TWAS Gene", "R2 > 0.4", "R2 <= 0.4"), pch = c(22,23,21,19,19,19), col = c("black","black","black", "salmon", "purple", colorsall[1]), cex=0.7, title.adj = 0)
  }
  
  
  
  if (label_panel=="cTWAS" | label_panel=="both"){
    for (i in 1:length(label_genes)){
      text(a$pos[a$genename==label_genes[i]], a$susie_pip[a$genename==label_genes[i]], labels=a$genename[a$genename==label_genes[i]], pos=label_pos[i], cex=0.7)
    }
  }
  
  if (return_table){
    return(a)
  }
  
  
  if (draw_gene_track){
    source("/project2/guiming/xsun/proteome_alzheimer/codes/trackplot.R")
    
    query_ucsc = TRUE
    build = "hg19"
    col = "gray70"
    txname = NULL
    genename = NULL
    collapse_txs = TRUE
    gene_model = "/project2/guiming/xsun/proteome_alzheimer/data_others/hg19.ensGene.gtf.gz"          #######
    isGTF = T
    
    ##########
    start <- min(a$pos)*1000000
    end <- max(a$pos)*1000000
    chr <- paste0("chr",as.character(unique(a$chrom)))
    
    #collect gene models
    if(is.null(gene_model)){
      if(query_ucsc){
        message("Missing gene model. Trying to query UCSC genome browser..")
        etbl = .extract_geneModel_ucsc(chr, start = start, end = end, refBuild = build, txname = txname, genename = genename)
      } else{
        etbl = NULL
      }
    } else{
      if(isGTF){
        etbl = .parse_gtf(gtf = gene_model, chr = chr, start = start, end = end, txname = txname, genename = genename)  
      } else{
        etbl = .extract_geneModel(ucsc_tbl = gene_model, chr = chr, start = start, end = end, txname = txname, genename = genename)  
      }
    }
    
    #draw gene models
    if(!is.null(etbl)){
      if(collapse_txs){
        etbl = .collapse_tx(etbl)
      }
      
      #subset to protein coding genes in ensembl and lincRNAs included in the analysis
      #etbl <- etbl[names(etbl) %in% G_list$ensembl_gene_id]
      etbl <- etbl[names(etbl) %in% c(G_list$ensembl_gene_id[G_list$gene_biotype=="protein_coding"], sapply(a$id[a$type=="gene"], function(x){unlist(strsplit(x, split="[.]"))[1]}))]
      
      for (i in 1:length(etbl)){
        ensembl_name <- attr(etbl[[i]], "gene")
        gene_name <- G_list$hgnc_symbol[match(ensembl_name, G_list$ensembl_gene_id)]
        if (gene_name==""){
          gene_name <- a$genename[sapply(a$id, function(x){unlist(strsplit(x,split="[.]"))[1]})==ensembl_name]
        }
        if (length(gene_name)>0){
          attr(etbl[[i]], "gene") <- gene_name
        }
        
        attr(etbl[[i]], "imputed") <- ensembl_name %in% sapply(ctwas_gene_res$id, function(x){unlist(strsplit(x, split="[.]"))[1]})
        
      }
      
      par(mar = c(0, 4.1, 0, 2.1))
      
      plot(NA, xlim = c(start, end), ylim = c(0, length(etbl)), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
      
      etbl <- etbl[order(-sapply(etbl, function(x){attr(x, "start")}))]
      
      for(tx_id in 1:length(etbl)){
        txtbl = etbl[[tx_id]]
        
        # if (attr(txtbl, "imputed")){
        #   exon_col = "#192a56"
        # } else {
        #   exon_col = "darkred"
        # }
        exon_col = "darkred"
        
        segments(x0 = attr(txtbl, "start"), y0 = tx_id-0.45, x1 = attr(txtbl, "end"), y1 = tx_id-0.45, col = exon_col, lwd = 1)
        
        if(is.na(attr(txtbl, "tx"))){
          text(x = start, y = tx_id-0.45, labels = paste0(attr(txtbl, "gene")), cex = 0.7, adj = 0, srt = 0, pos = 2, xpd = TRUE)
        } else {
          text(x = start, y = tx_id-0.45, labels = paste0(attr(txtbl, "tx"), " [", attr(txtbl, "gene"), "]"), cex = 0.7, adj = 0, srt = 0, pos = 2, xpd = TRUE)
        }
        
        rect(xleft = txtbl[[1]], ybottom = tx_id-0.75, xright = txtbl[[2]], ytop = tx_id-0.25, col = exon_col, border = NA)
        if(attr(txtbl, "strand") == "+"){
          dirat = pretty(x = c(min(txtbl[[1]]), max(txtbl[[2]])))
          dirat[1] = min(txtbl[[1]]) #Avoid drawing arrows outside gene length
          dirat[length(dirat)] = max(txtbl[[2]])
          points(x = dirat, y = rep(tx_id-0.45, length(dirat)), pch = ">", col = exon_col)
        }else{
          dirat = pretty(x = c(min(txtbl[[1]]), max(txtbl[[2]])))
          dirat[1] = min(txtbl[[1]]) #Avoid drawing arrows outside gene length
          dirat[length(dirat)] = max(txtbl[[2]])
          points(x = dirat, y = rep(tx_id-0.45, length(dirat)), pch = "<", col = exon_col)
        }
      }
    }
  }
  
}

