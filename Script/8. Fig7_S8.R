
# fig7, Characterization of the alteration of B cells and plasma cells undergoing azacytidine +nivolumab treatment in Relapsed/Refractory patients

# fig7A
  
  p1 <- DimPlot(object_b,cells.highlight = list(`Prior-treat`=WhichCells(object_b,expression = is_treatment == "Prior-Treat")),
                cols.highlight="#98df8a",cols="linen",sizes.highlight=0.05,raster=F)+
    theme(legend.text=element_text(size = 8),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm")) & NoAxes()
  p2 <- DimPlot(object_b,cells.highlight = list(`Post-treat`=WhichCells(object_b,expression = is_treatment == "Post-Treat")),
                cols.highlight="#2ca02c",cols="linen",sizes.highlight=0.05,raster=F)+
    theme(legend.text=element_text(size = 8),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm")) & NoAxes()
  p1+p2+plot_layout(nrow=1,guides="collect")
  

# fig7B
  
  df_pie <- object_b@meta.data
  df_pre <- df_pie[which(df_pie$is_treatment == "Prior-Treat"),]
  df_post <- df_pie[which(df_pie$is_treatment == "Post-Treat"),]
  pieChart(df_pre,group="celltype",title="Prior-Treat",color=celltype.colors,legend_name="Celltype")
  pieChart(df_post,group="celltype",title="Post-Treat",color=celltype.colors,legend_name="Celltype")
  
  
# fig7C
  
  load(file="Data/patient_obs.rda")
  patient_obs_rr <- patient_obs[patient_obs$sample_type == "RelRef",]

  clonalProportion(clone_data=patient_obs_rr,split = c(5,100,200,500),
                   group="is_treatment",group_order=names(is_treatment.colors))
  
  
# fig7D
  
  cloneDiversity(clone_data=patient_obs_rr,col_name="is_treatment",colors=is_treatment.colors,legend_name="Sample")
  
  
# fig7E
  
  patient_obs_rr$patient <- str_replace_all(patient_obs_rr$sample,c("A"="-Pre","B"="-Post","C"="-Post"))
  
  treat_clone <- countClones(patient_obs_rr,group=c("is_treatment","patient","Best_response"))
  pre_clone <- treat_clone[treat_clone$is_treatment =="Prior-Treat",]
  post_clone <- treat_clone[treat_clone$is_treatment =="Post-Treat",]
  treat_clone <- data.table(full_join(pre_clone,post_clone,by="clone_id"))
  treat_clone[is.na(treat_clone)] <- 0
  treat_clone[seq_freq.x == 0 & seq_freq.y > 0,clone_type:="Novel"]
  treat_clone[seq_freq.x > seq_freq.y & is.na(clone_type),clone_type:="Contracted"]
  treat_clone[seq_freq.x < seq_freq.y & is.na(clone_type),clone_type:="Expanded"]
  treat_clone[seq_freq.x == seq_freq.y & is.na(clone_type),clone_type:="Persitent"]
  treat_clone[Best_response.x == "0",response:=Best_response.y]
  treat_clone[is.na(response),response:=Best_response.x]
  treat_clone[patient.x == "0",patient:=patient.y]
  treat_clone[is.na(patient),patient:=patient.x]
  
  treat_clone$response <- factor(treat_clone$response,levels=names(Best_response.colors))
  treat_clone$clone_type <- factor(treat_clone$clone_type,levels=c("Novel","Expanded","Contracted"))
  phenoPropotion(data=as.data.frame(treat_clone),x_name="clone_type",y_name="response",
                 legend_name="Response",cols=Best_response.colors)  
  
  
# fig7F
  
  plasma_rr <- subset(object_b,subset=celltype == "Plasma" & sample_type == "RelRef")
  
  GeneCounts <- as.matrix(plasma_rr@assays$RNA@counts)
  iOrd <- rowSums(GeneCounts>0)
  GeneCounts <- GeneCounts[iOrd>10,]#only keep genes expressing in more than 10 cell
  plasma_rr_CytoTRACE <- CytoTRACE(GeneCounts, enableFast = TRUE,subsamplesize = 1000)
  plasma_rr$CytoTRACE <- plasma_rr_CytoTRACE$CytoTRACE[colnames(plasma_rr)]
  
  ggplot(plasma_rr@meta.data,aes(x=is_treatment,y=CytoTRACE))+
    geom_boxplot(aes(fill=is_treatment),width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    scale_fill_manual(values=is_treatment.colors)+
    stat_signif(comparisons = list(c("Post-Treat","Prior-Treat")),textsize=2.5,y_position=c(0.95),map_signif_level=T)+
    theme_classic()
  
  
# fig7G-I
  
  # gerenated from "4. Fig3_S5.R"
  # naive_mem <- readRDS(file.path("Data","object_naive_mem.rds"))
  naive_mem_rr <- subset(naive_mem,subset=sample_type == "RelRef" & celltype %in% c("Naive B","CD27+ memory B","Atypical memory B"))
  naive_mem_rr$is_treatment <- factor(naive_mem_rr$is_treatment,levels=c("Prior-Treat","Post-Treat"))
  
  k="AP.1.transcription.factor51"
  ggplot(naive_mem_rr@meta.data,aes_string(x="celltype",y=k,fill="is_treatment"))+
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    scale_fill_manual(values=is_treatment.colors,name="Sample")+        
    xlab("")+ylab("AP-1")+
    theme_classic()+
    stat_compare_means(aes(group = is_treatment),
                       label = "p.signif",size=2.5,
                       method = "wilcox.test",hide.ns = F)
  
  k="HALLMARK_TNFA_SIGNALING_VIA_NFKB45"
  ggplot(naive_mem_rr@meta.data,aes_string(x="celltype",y=k,fill="is_treatment"))+
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    scale_fill_manual(values=is_treatment.colors,name="Sample")+        
    xlab("")+ylab("TNFA signaling via NFKB")+
    theme_classic()+
    stat_compare_means(aes(group = is_treatment),
                       label = "p.signif",size=2.5,
                       method = "wilcox.test",hide.ns = F)
  
  k="Antigen.processing.and.presentation...Homo.sapiens..human.73"
  ggplot(naive_mem_rr@meta.data,aes_string(x="celltype",y=k,fill="is_treatment"))+
    geom_boxplot(width=0.5,color='gray40',outlier.shape=NA,size=0.3)+
    scale_fill_manual(values=is_treatment.colors,name="Sample")+        
    xlab("")+ylab("Antigen processing and presentation")+
    theme_classic()+
    stat_compare_means(aes(group = is_treatment),
                       label = "p.signif",size=2.5,
                       method = "wilcox.test",hide.ns = F)
  
  
# fig7J
  
  gseBP <- msigdbr(species = "Homo sapiens",
                   category = "C5", 
                   subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
  gseBP$gs_name <- gsub('GOBP_','',gseBP$gs_name)
  gseBP$gs_name <- gsub('_',' ',gseBP$gs_name)
  gseBP$gs_name <- str_to_title(gseBP$gs_name)
  
  mem_rr <- subset(object_b,subset=sample_type == "RelRef" & celltype %in% c("CD27+ memory B"))
  
  diff <- FindMarkers(mem_rr,ident.1="Post-Treat",ident.2="Prior-Treat",group.by="is_treatment",logfc.threshold =0)
  diff$gene <- row.names(diff)
  diff <- arrange(diff,desc(avg_log2FC))
  geneList <- diff$avg_log2FC
  names(geneList) <- diff$gene
  
  gseaRes <- fgsea(pathways = gseBP, 
                   stats = geneList,
                   minSize=5,
                   maxSize=500)
  Down <- gseaRes[NES < 0 & padj < 0.1][tail(order(desc(NES)), n=20), pathway]
  gseaSets_plot <- gseBP[Down]
  plotGseaTable(gseaSets_plot, geneList, gseaRes,gseaParam = 0.5)
  
  # write.xlsx(gseaRes[padj < 0.1][order(desc(NES))],file.path("Data/Tables","TableS8_CD27Mem_Post_vs_Prior_GSEA.xlsx"))
  
  
  
  
  
  
  
  
  
  
