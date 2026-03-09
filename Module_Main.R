# file name: PathFindeR_Modules
#
# Latest update: 27th September 2025
#
#-----------------------------------------------------------------------------
# List of disease identifiers for Open target queries
#-----------------------------------------------------------------------------
#
Disease_list<-yaml::read_yaml("mydiseases.yml")
# Trait_list<-yaml::read_yaml("mytraits.yml")
# Disease_list<-c(Disease_list,Trait_list)
#
#-----------------------------------------------------------------------------
# Save abbreviations separately
#-----------------------------------------------------------------------------
#
Disease_abbr<-sub(".*\\[(.*)\\].*", "\\1", names(Disease_list))
names(Disease_list)<-sub("\\s*\\[.*\\]$", "", names(Disease_list))
#
#-----------------------------------------------------------------------------
# Filters and parameters
#-----------------------------------------------------------------------------
#
L2G_co<-0.5       # Locus to Gene score lower limit
CG_co<-0.5        # Clingen score lower limit
GB_co<-0.5        # GeneBurden lower limit
sample_co<-0      # Sample size lower limit (for GWAS) 
STRING.co<-0.4    # Cut-off for PPI interactions on STRING
Nrand<-1000       # number of random disease modules
#
#-----------------------------------------------------------------------------
# Load functions and data
#-----------------------------------------------------------------------------
#
source("Module_Func.R")
#
#-------------------------------------------------------------------------------
# Create a directory for Modules if not present
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Modules")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
#-----------------------------------------------------------------------------
# Retrieve genes from GWAS and NGS for each disease
#-----------------------------------------------------------------------------
#
for (j in 1:length(Disease_list)) {
  #
  print(paste0("I am retrieving genes for ",names(Disease_list)[j]))
  #
  if (!file.exists(paste0("Modules/",names(Disease_list)[j],".csv"))) {
    #
    #-----------------------------------------------------------------------------
    # I have a custom module for Chronic fatigue syndrome
    #-----------------------------------------------------------------------------
    #
    if (Disease_list[[j]]=="customCFS") {
      mytargets<-fread("Data/all_genes_customCFS.tsv")
      write.table(mytargets,file=paste0("Modules/",names(Disease_list)[j],".csv"),
                  sep=",",row.names=FALSE,quote=FALSE)
    } else {
      #
      #-----------------------------------------------------------------------------
      # Retrieve targets and save
      #-----------------------------------------------------------------------------
      #
      success<-0
      while(success==0) {
        mytargets<-try(Targets4Disease(efo_id=Disease_list[[j]],L2G_cutoff=L2G_co,
                                       clingen_cutoff=CG_co,geneburden_cutoff=GB_co,
                                       sample_cutoff=sample_co),silent=T)
        if (!inherits(mytargets,"try-error")) {
          success<-1
        } else {
          print("Targets4Disease() failed and I'll try again")
        }
      }
      #
      if (nrow(mytargets)>0) {
        #
        #-----------------------------------------------------------------------------
        # Add STRING preferred name
        #-----------------------------------------------------------------------------
        #
        gene.wo.STRING<-0
        index<-c()
        mytargets$name<-rep(NA,nrow(mytargets))
        for (i in 1:nrow(mytargets)) {
          new.name<-STRING.name(mytargets$target.approvedSymbol[i])
          if (!is.na(new.name)) {
            mytargets$name[i]<-new.name
          } else if (is.na(new.name)) {
            print(paste(mytargets$target.approvedSymbol[i],"is not present in STRING's database"))
            gene.wo.STRING<-gene.wo.STRING+1
            index<-c(index,i)
          }
        }
        if (length(index>0)) mytargets<-mytargets[-index,] # remove genes not in STRING
        print(paste(gene.wo.STRING,"genes are not present in STRING's database"))
        #
        #-----------------------------------------------------------------------------
        # Add NCBI.id
        #-----------------------------------------------------------------------------
        #
        mytargets$NCBI.id<-rep(NA,nrow(mytargets))
        for (i in 1:nrow(mytargets)) {
          mytargets$NCBI.id[i]<-Symbol2NCBI.db(mytargets$target.approvedSymbol[i])
        }
        #
        #-----------------------------------------------------------------------------
        # Add list name (rare variant vs GWAS)
        #-----------------------------------------------------------------------------
        #
        mytargets$list.name<-rep(NA,nrow(mytargets))
        mytargets$list.count<-rep(NA,nrow(mytargets))
        #
        for (i in 1:nrow(mytargets)) {
          if (!is.na(mytargets$GWAS.id[i])&
              (!is.na(mytargets$clingen[i])|!is.na(mytargets$geneburden[i]))) {
            mytargets$list.name[i]<-"GWAS/Rare" 
            mytargets$list.count[i]<-2 
          } else if (!is.na(mytargets$GWAS.id[i])&
                     is.na(mytargets$clingen[i])&is.na(mytargets$geneburden[i])) {
            mytargets$list.name[i]<-"GWAS"
            mytargets$list.count[i]<-1
          } else if (is.na(mytargets$GWAS.id[i])&
                     (!is.na(mytargets$clingen[i])|!is.na(mytargets$geneburden[i]))) {
            mytargets$list.name[i]<-"Rare"
            mytargets$list.count[i]<-1
          }
        }
      }
    }
    #
    #-----------------------------------------------------------------------------
    # Save results
    #-----------------------------------------------------------------------------
    #
    write.table(mytargets,file=paste0("Modules/",names(Disease_list)[j],".csv"),
                sep=",",row.names=FALSE,quote=FALSE)
  }
}
#
#-----------------------------------------------------------------------------
# Build a catalog of disease proteins
#-----------------------------------------------------------------------------
#
myDiseaseGenes<-data.frame(NCBI.id=NA,name=NA)
#
for (j in 1:length(Disease_list)) {
  file.name<-paste0("Modules/",names(Disease_list)[j],".csv")
  if (file.exists(file.name)) {
    mytargets<-fread(file=file.name)  
    df<-subset.data.frame(mytargets,select=c("name","NCBI.id"))
    myDiseaseGenes<-rbind(myDiseaseGenes,df)
  }
}
myDiseaseGenes<-na.omit(myDiseaseGenes)
myDiseaseGenes<-unique(myDiseaseGenes)
#
# Save results
#
write.table(myDiseaseGenes,file="Modules/myDiseaseGenes.csv",sep=",",row.names=F,quote=F)
#
#-----------------------------------------------------------------------------
# Restrict STRING network to disease genes
#-----------------------------------------------------------------------------
#
myDiseaseGenes<-fread(file="Modules/myDiseaseGenes.csv",sep=",")
index<-which(STRING.names$preferred_name%in%myDiseaseGenes$name)
mynames<-STRING.names$string_protein_id[index]
index<-which(STRING.matrix$protein1%in%mynames)
STRING.matrix<-STRING.matrix[index,]
index<-which(STRING.matrix$protein2%in%mynames)
STRING.matrix<-STRING.matrix[index,]
#
#-----------------------------------------------------------------------------
# Build the gene network for the complete set of disease genes
#-----------------------------------------------------------------------------
#
if (!file.exists("Modules/myDiseaseGenes_net.rds")) {
  all.genes<-fread(file="Modules/myDiseaseGenes.csv",sep=",")
  gene.matrix<-GeneMatrix(all.genes)
  saveRDS(gene.matrix,file="Modules/myDiseaseGenes_net.rds")  
}
#
#-----------------------------------------------------------------------------
# Build the gene network for each disease
#-----------------------------------------------------------------------------
#
gene.matrix.full<-readRDS("Modules/myDiseaseGenes_net.rds") # read complete network
#
for (j in 1:length(Disease_list)) {
  #
  print(paste0("I am building the gene network for ",names(Disease_list)[j]))
  #
  if (!file.exists(file=paste0("Modules/",names(Disease_list)[j],".rds"))) {
    #
    mytargets<-fread(file=paste0("Modules/",names(Disease_list)[j],".csv"),sep=",")
    #
    #-----------------------------------------------------------------------------
    # Build the adjacency matrix associated with the merged gene list
    #-----------------------------------------------------------------------------
    #
    all.genes<-subset.data.frame(mytargets,select=c("name","list.name","list.count"))
    all.genes<-na.omit(all.genes)
    index<-which(rownames(gene.matrix.full)%in%all.genes$name)
    gene.matrix<-gene.matrix.full[index,index]
    saveRDS(gene.matrix,file=paste0("Modules/",names(Disease_list)[j],".rds"))
    #
    #-----------------------------------------------------------------------------
    # Build the graph associated with the Merged Gene List (MGL) and plot it
    #-----------------------------------------------------------------------------
    #
    graph<-graph_from_adjacency_matrix(gene.matrix,mode="undirected",weighted=TRUE)
    #
    # Color the genes according to the corresponding Expanded Gene List (EGL)
    #
    all.genes<-all.genes[match(rownames(gene.matrix),all.genes$name), ] # correct the order!
    colors<-hcl.colors(2,palette="Dark 3",alpha=1)
    vertex_colors<-c()
    NV<-length(vertex(graph)[[1]]) # number of verteces
    for (v in 1:NV) {
      index<-which(c("GWAS","Rare")==all.genes$list.name[v])
      if (length(index)>0) {
        vertex_colors[v]<-colors[index]
      }
    }
    index<-which(is.na(vertex_colors))
    vertex_colors[index]<-"white"
    #
    # Set the color of nodel labels
    #
    node.col<-rep("black",nrow(all.genes))
    #
    # Plot the image
    #
    tiff(paste0("Modules/",names(Disease_list)[j],".tiff"),width=10,height=10,units="in",res=600,compression="lzw")
    plot(graph,
         layout=layout_with_fr(graph,niter=30000,grid="nogrid",dim=2),
         vertex.size=300/NV,
         vertex.label.cex=20/NV,
         vertex.frame.color="black",
         vertex.color=vertex_colors,
         vertex.label.color=node.col,
         asp=1,
    )
    legend("topright",                                        # or "bottomleft", etc.  
           legend = c("GWAS","Rare","GWAS&Rare"),             # group names
           col = c(colors,"black"),                           # corresponding colors of circle
           pch = 21,                                          # filled circle
           pt.bg = c(colors,"white"),                         # fill color
           pt.cex = 1.5,                                      # size of points
           bty = "n")                                         # no box around legend
    title(main=Disease_list[[j]])
    dev.off()
    #
    # Save a jpeg version too
    #
    img<-image_read(paste0("Modules/",names(Disease_list)[j],".tiff"))
    image_write(img,path=paste0("Modules/",names(Disease_list)[j],".jpeg"),format="jpeg")
    #
    # Save graph for cytoscape
    #
    edges<-as.data.frame(as_edgelist(graph))
    weights<-E(graph)$weight
    edges_df<-data.frame(source=edges[,1],target=edges[,2],interaction="interacts_with",
                         weight=weights)
    file_name=paste0("Modules/",names(Disease_list)[j],"_cytoscape",".tsv")
    write.table(edges_df,file=file_name,sep="\t",row.names=FALSE,quote=FALSE)
  }
}
#
#-----------------------------------------------------------------------------
# Build gene networks for random diseases, Nrand random networks for each real disease
#-----------------------------------------------------------------------------
#
gene.matrix.full<-readRDS("Modules/myDiseaseGenes_net.rds") # read complete network
#
for (j in 1:length(Disease_list)) {
  file.name<-paste0("Modules/",names(Disease_list)[j],".csv")
  print(paste0("Random diseases for ",names(Disease_list)[j]))
  all.genes<-fread(file.name)  
  NG<-nrow(all.genes)
  list.matrix<-list()
  for (i in 1:Nrand) {
    #
    # Extract random genes
    #
    rows<-sample(seq(1:nrow(myDiseaseGenes)),NG)
    all.genes<-myDiseaseGenes[rows,]
    #
    # Build gene network
    #
    print(paste0("Gene network for random disease number ",i," of disease ",names(Disease_list)[j]))
    index<-which(rownames(gene.matrix.full)%in%all.genes$name)
    list.matrix[[i]]<-gene.matrix.full[index,index]
  }
  saveRDS(list.matrix,file=paste0("Random/",names(Disease_list)[j],"_",i,"_gene_matrix.rsd"))
}
#
#-----------------------------------------------------------------------------
# Study the properties of the disease modules and save them
#-----------------------------------------------------------------------------
#
myDisMod<-data.frame(Disease=c(names(Disease_list)))
for (j in 1:nrow(myDisMod)) {
  print(paste0("Working on ",names(Disease_list)[j]))
  #
  # Read disease matrix
  #
  gene.matrix<-readRDS(file=paste0("Modules/",names(Disease_list)[j],".rds"))
  #
  # Read corresponding random disease
  #
  list.matrix<-readRDS(file=paste0("Random/",names(Disease_list)[j],"_",1000,"_gene_matrix.rsd"))
  #
  # Number of genes
  #
  myDisMod$Vertices_Num[j]<-nrow(gene.matrix)
  #
  # Build the network for this disease
  #
  graph<-graph_from_adjacency_matrix(gene.matrix,mode="undirected",weighted=TRUE)
  #
  # Module Size (size of largest connected sub graph)
  #
  com<-igraph::components(graph)
  myDisMod$Module_Size[j]<-max(com$csize)
  myDisMod$Module_Size_Per[j]<-round(myDisMod$Module_Size[j]/myDisMod$Vertices_Num[j],2)
  #
  # Mean shortest distances within connected components
  #
  myDisMod$Mean_Short_Dist[j]<-mean_distance(graph,weights=E(graph)$weight,directed=F,
                                             unconnected=T)
  #
  # Mean degree
  #
  deg<-igraph::degree(graph)
  myDisMod$Mean_Degree[j]<-mean(deg)
  #
  # Mean strength (weighted degree)
  #
  wdeg<-strength(graph,mode="all",weights=E(graph)$weight)
  myDisMod$Mean_Strength[j]<-mean(wdeg)
  myDisMod$Mean_Strength_Rel[j]<-round(myDisMod$Mean_Strength[j]/myDisMod$Mean_Degree[j],2)
  #
  # Study of the corresponding random diseases
  #
  Module_Size<-c()
  Mean_Short_Dist<-c()
  Mean_Degree<-c()
  Mean_Strength_Rel<-c()
  #
  for (i in 1:Nrand) {
    #
    graph<-graph_from_adjacency_matrix(list.matrix[[i]],mode="undirected",weighted=TRUE)
    #
    # Module Size (size of largest connected sub graph)
    #
    com<-igraph::components(graph)
    Module_Size[i]<-max(com$csize)
    #
    # Mean shortest distances within connected components
    #
    Mean_Short_Dist[i]<-mean_distance(graph,weights=E(graph)$weight,directed=F)
    #
    # Mean degree
    #
    Mean_Degree[i]<-mean(igraph::degree(graph))
    #
    # Mean strength (weighted degree)
    #
    wdeg<-strength(graph,mode="all",weights=E(graph)$weight)
    Mean_Strength<-mean(wdeg)
    Mean_Strength_Rel[i]<-round(Mean_Strength/Mean_Degree[i],2)
    #
  }
  #
  # Run statistical tests
  #
  myDisMod$Module_Size_P[j]<-P_upper(myDisMod$Module_Size[j],Module_Size)
  myDisMod$Mean_Short_Dist_P[j]<-P_lower(myDisMod$Mean_Short_Dist[j],Mean_Short_Dist)
  myDisMod$Mean_Degree_P[j]<-P_upper(myDisMod$Mean_Degree[j],Mean_Degree)
  myDisMod$Mean_Strength_P[j]<-P_upper(myDisMod$Mean_Strength[j],Mean_Strength)
  myDisMod$Mean_Strength_Rel_P[j]<-P_upper(myDisMod$Mean_Strength_Rel[j],Mean_Strength_Rel)
}
#
# Save the data frame
#
write.table(myDisMod,file=paste0("Modules/Modules_analysis.csv"),
            sep=",",row.names=FALSE,quote=FALSE)
#
#-----------------------------------------------------------------------------
# Over-representation analysis (KEGG, Reactome, GO, Tissues)
#-----------------------------------------------------------------------------
#
for (j in 1:length(Disease_list)) {
  if (!file.exists(paste0("Modules/",names(Disease_list)[j],"_Tissue_ORA.tsv"))) {
    file.name<-paste0("Modules/",names(Disease_list)[j],".csv")
    if (file.exists(file.name)) {
      print(paste0("ORA for ",names(Disease_list)[j]))
      #
      # ORA on KEGG, Reactome, GOcc, and DO 
      #
      all.genes<-fread(file.name)
      ORA<-ORA.fun(all.genes)
      file.name<-paste0("Modules/",names(Disease_list)[j],"_ORA.tsv")
      write.table(ORA,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")  
      #
      # Tissue ORA
      #
      ORA.tissues<-Tissue.ORA(all.genes,pl=T)
      file.name<-paste0("Modules/",names(Disease_list)[j],"_Tissue_ORA.tsv")
      write.table(ORA.tissues,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
    }  
  }
}    
#
#-----------------------------------------------------------------------------
# Generate Random Diseases using the genes in myDiseaseGenes
# and compute ORA for each one of them
#-----------------------------------------------------------------------------
#
# Create output folders, if absent
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Random")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Calculate mean size of disease modules
#
myDisMod<-fread("Modules/Modules_analysis.csv")
N<-round(mean(myDisMod$Vertices_Num),0) # mean size of disease module
#
# Calculate ORA for random disease modules
#
for (i in 1:Nrand) {
  if (!file.exists(paste0("Random/",i,"_ORA.tsv"))) {
    #
    # Extract random disease module and save it 
    #
    rows<-sample(seq(1:nrow(myDiseaseGenes)),N)
    all.genes<-myDiseaseGenes[rows,]
    file.name<-paste0("Random/",i,"_genes.tsv")
    write.table(all.genes,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
    #
    # ORA on KEGG, Reactome, GO 
    #
    print(paste0("ORA for random disease number ",i))
    ORA<-ORA.fun(all.genes)
    file.name<-paste0("Random/",i,"_ORA.tsv")
    write.table(ORA,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")  
    #
    # Tissue ORA
    #
    ORA.tissues<-Tissue.ORA(all.genes,pl=F)
    file.name<-paste0("Random/",i,"_Tissue_ORA.tsv")
    write.table(ORA.tissues,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
  }
}
#
#-----------------------------------------------------------------------------
# Study similarity between diseases by gene overlap
# Jaccard index = intersection(geneA,geneB)/union(geneA,geneB)
#-----------------------------------------------------------------------------
#
# Create output folders, if absent
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Comparisons")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
current_dir<-getwd()
folder_path<-file.path(current_dir,"Comparisons/Jaccard")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Calculate similitude scores by gene overlap (Jaccard Index)
#
Score<-matrix(data=NA,nrow=length(Disease_list),ncol=length(Disease_list))
colnames(Score)<-names(Disease_list)
rownames(Score)<-names(Disease_list)
Overlap<-Score
LogP<-Score
for (i in 1:length(Disease_list)) {
  for (j in i:length(Disease_list)) {
    if (i!=j) {
      #
      # Select input files
      #
      module_1<-paste0("Modules/",names(Disease_list)[i],".csv")
      module_2<-paste0("Modules/",names(Disease_list)[j],".csv")
      if (file.exists(module_1)&file.exists(module_2)) {
        print(paste0("Comparing ",names(Disease_list)[i]," and ",names(Disease_list)[j]))
        #
        # Read files
        #
        mytargets_1<-fread(module_1)
        mytargets_1<-subset.data.frame(mytargets_1,select=c("name"))
        N1<-nrow(mytargets_1)
        mytargets_2<-fread(module_2)
        mytargets_2<-subset.data.frame(mytargets_2,select=c("name"))
        N2<-nrow(mytargets_2)
        #
        # Calculate overlap pval
        #
        NR<-length(which(mytargets_1$name%in%mytargets_2$name)) # overlap
        NB<-length(myuniverse)
        #
        LogP[i,j]<-as.numeric(phyper(q=NR-1,m=N1,n=NB-N1,k=N2,lower.tail=F)) # P(X>NR-1)=P(X>=NR)
        LogP[i,j]<--log10(LogP[i,j])
        if (is.infinite(LogP[i,j])) LogP[i,j]<-NA
        LogP[j,i]<-LogP[i,j] 
        #
        # Store gene overlap
        #
        Overlap[i,j]<-NR
        Overlap[j,i]<-NR
        #
        # Calculate Jaccard index
        #
        Score[i,j]<-NR/(N1+N2-NR)
        Score[j,i]<-NR/(N1+N2-NR)
      }
    } else {
      module_1<-paste0("Modules/",names(Disease_list)[i],".csv")
      mytargets_1<-fread(module_1)
      N1<-nrow(mytargets_1)
      Score[i,j]<-1
      Overlap[j,i]<-N1
    }
  } 
}
#
# Save results
#
Score.list<-list(Score,LogP,Overlap)
names(Score.list)<-c("Score","LogP","Overlap")
saveRDS(Score.list,file="Comparisons/Jaccard/Score.rds")
#
# Read results
#
Score.list<-readRDS("Comparisons/Jaccard/Score.rds")
Score<-Score.list[[1]]
LogP<-Score.list[[2]]
#
# Convert to distance
#
MaxScore<-max(Score,na.rm=T)
MinScore<-min(Score,na.rm=T)
Dist<-(MaxScore-Score)/(MaxScore-MinScore)
#
# Plot dendrogram
#
tiff("Comparisons/Jaccard/Tree_Jaccard.tiff",width=7,height=25,units="in",res=600,
     compression="lzw")
hc<-hclust(as.dist(Dist),method="average")
plot(hc)
dev.off()
#
# Save a jpeg version too
#
img<-image_read("Comparisons/Jaccard/Tree_Jaccard.tiff")
image_write(img,path="Comparisons/Jaccard/Tree_Jaccard.jpeg",format="jpeg")
#
# For each disease, build a plot -Log10P-Jaccard_index
#
distance<-"Jaccard"
alpha<-0.05
Distance_plot(Score,LogP,distance,alpha)
#
#-----------------------------------------------------------------------------
# Study similarity between diseases by regression of zScores from ORA
#-----------------------------------------------------------------------------
#
# Create output folders, if absent
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Comparisons/Correlation")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
current_dir<-getwd()
folder_path<-file.path(current_dir,"Comparisons/Correlation/Distributions")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Prepare a matrix for correlations
#
ORA_cor<-matrix(data=NA,nrow=length(Disease_list),ncol=Nrand)
row.names(ORA_cor)<-names(Disease_list)
colnames(ORA_cor)<-seq(1,Nrand,1)
#
# For each disease build a vector with regression coefficients with each random
# disease module
#
for (i in 1:length(Disease_list)) {
  #
  print(paste0("Comparing ",names(Disease_list)[i]," and random diseases by correlation of ORA terms"))
  #
  # Correlations of Disease i with random Diseases
  #
  for (k in 1:Nrand) {
    #
    file_1<-paste0("Modules/",names(Disease_list)[i],"_ORA.tsv")
    file_2<-paste0("Random/",k,"_ORA.tsv")
    if (file.exists(file_1)&file.exists(file_2)) {
      #
      # Read files
      #
      ORA_1<-fread(file_1)
      ORA_1<-subset.data.frame(ORA_1,select=c("ID","zScore"))
      ORA_2<-fread(file_2)
      ORA_2<-subset.data.frame(ORA_2,select=c("ID","zScore"))
      #
      # Merge
      #
      myORA<-merge(ORA_1,ORA_2,by=c("ID"),all=T)
      #
      # Scale
      #
      maxzS<-max(myORA$zScore.x,na.rm=T)
      minzS<-min(myORA$zScore.x,na.rm=T)
      myORA$zScore.x<-(myORA$zScore.x-minzS)/(maxzS-minzS)
      maxzS<-max(myORA$zScore.y,na.rm=T)
      minzS<-min(myORA$zScore.y,na.rm=T)
      myORA$zScore.y<-(myORA$zScore.y-minzS)/(maxzS-minzS)
    } 
    file_1<-paste0("Modules/",names(Disease_list)[i],"_Tissue_ORA.tsv")
    file_2<-paste0("Random/",k,"_Tissue_ORA.tsv")
    if (file.exists(file_1)&file.exists(file_2)) {
      #
      # Read files
      #
      ORA_1<-fread(file_1)
      ORA_1<-subset.data.frame(ORA_1,select=c("Tissue","zScore"))
      ORA_2<-fread(file_2)
      ORA_2<-subset.data.frame(ORA_2,select=c("Tissue","zScore"))
      #
      # Merge
      #
      myTissueORA<-merge(ORA_1,ORA_2,by=c("Tissue"),all=T)
      colnames(myTissueORA)<-c("ID","zScore.x","zScore.y")
      #
      # Scale
      #
      maxzS<-max(myTissueORA$zScore.x,na.rm=T)
      minzS<-min(myTissueORA$zScore.x,na.rm=T)
      myTissueORA$zScore.x<-(myTissueORA$zScore.x-minzS)/(maxzS-minzS)
      myTissueORAS<-max(myTissueORA$zScore.y,na.rm=T)
      minzS<-min(myTissueORA$zScore.y,na.rm=T)
      myTissueORA$zScore.y<-(myTissueORA$zScore.y-minzS)/(maxzS-minzS)
      #
      # Merge
      #
      myOra<-rbind(myORA,myTissueORA)
    }
    Random.cor[k]<-c(cor(myORA$zScore.x,myORA$zScore.y,method="spearman",
                         use="complete.obs"))
  }
  #
  # Save vector 
  #
  ORA_cor[i,]<-Random.cor
  saveRDS(ORA_cor,file="Comparisons/Correlation/ORA_cor.rds")
}
#
# Run correlation between two diseases and test for significance
#
Score<-matrix(data=NA,nrow=length(Disease_list),ncol=length(Disease_list))
colnames(Score)<-names(Disease_list)
rownames(Score)<-names(Disease_list)
LogP<-Score
ORA_cor<-readRDS("Comparisons/Correlation/ORA_cor.rds")
for (i in 1:length(Disease_list)) {
  for (j in i:length(Disease_list)) {
    if (i!=j) {
      #
      print(paste0("Comparing ",names(Disease_list)[i]," and ",names(Disease_list)[j]))
      #
      # Select input files for ORA (KEGG, GO, Reactome)
      #
      file_1<-paste0("Modules/",names(Disease_list)[i],"_ORA.tsv")
      file_2<-paste0("Modules/",names(Disease_list)[j],"_ORA.tsv")
      if (file.exists(file_1)&file.exists(file_2)) {
        #
        # Read files
        #
        ORA_1<-fread(file_1)
        ORA_1<-subset.data.frame(ORA_1,select=c("ID","zScore"))
        ORA_2<-fread(file_2)
        ORA_2<-subset.data.frame(ORA_2,select=c("ID","zScore"))
        #
        # Merge
        #
        myORA<-merge(ORA_1,ORA_2,by=c("ID"),all=T)
      } 
      #
      # Select input files for Tissue ORA
      #
      file_1<-paste0("Modules/",names(Disease_list)[i],"_Tissue_ORA.tsv")
      file_2<-paste0("Modules/",names(Disease_list)[j],"_Tissue_ORA.tsv")
      if (file.exists(file_1)&file.exists(file_2)) {
        #
        # Read files
        #
        ORA_1<-fread(file_1)
        ORA_1<-subset.data.frame(ORA_1,select=c("Tissue","zScore"))
        ORA_2<-fread(file_2)
        ORA_2<-subset.data.frame(ORA_2,select=c("Tissue","zScore"))
        #
        # Merge
        #
        myTissueORA<-merge(ORA_1,ORA_2,by=c("Tissue"),all=T)
        colnames(myTissueORA)<-c("ID","zScore.x","zScore.y")
        myOra<-rbind(myORA,myTissueORA)
      }
      Score[i,j]<-cor(myORA$zScore.x,myORA$zScore.y,method="spearman",use="complete.obs")
      Score[j,i]<-Score[i,j]
      #
      # Correlations of Disease i with random Diseases
      #
      Random.cor_i<-ORA_cor[i,]
      #
      # Correlations of Disease j with random Diseases
      #
      Random.cor_j<-ORA_cor[j,]
      #
      # Merge
      #
      Random.cor<-c(Random.cor_i,Random.cor_j)
      #
      # Calculate p value for correlation
      #
      Num<-length(which(Random.cor>=Score[i,j]))
      Den<-length(which(!is.na(Random.cor)))
      p_upper<-Num/Den
      #
      if (!is.infinite(log10(p_upper))) {
        LogP[i,j]<--log10(p_upper)
        LogP[j,i]<-LogP[i,j]
      }
      if (TRUE) {
        #
        # Plot null distribution
        #
        tiff(paste0("Comparisons/Correlation/Distributions/Hist_",names(Disease_list)[i],"_",
                    names(Disease_list)[j],".tiff"),width=7,height=7,units="in",
             res=600,compression="lzw")
        hist(Random.cor)
        abline(v=Score[i,j],col="red",lty=1)
        dev.off()
      }
    } 
  } 
}
#
# Save results
#
Score.list<-list(Score,LogP)
names(Score.list)<-c("Score","LogP")
saveRDS(Score.list,file="Comparisons/Correlation/Score.rds")
#
# Read results
#
Score.list<-readRDS("Comparisons/Correlation/Score.rds")
Score<-Score.list[[1]]
LogP<-Score.list[[2]]
#
# Convert to distance
#
MaxScore<-max(Score,na.rm=T)
MinScore<-min(Score,na.rm=T)
Dist<-(MaxScore-Score)/(MaxScore-MinScore)
v.ORA<-Dist[lower.tri(Dist)]
#
# Plot dendrogram
#
tiff("Comparisons/Correlation/Tree_ORA_cor.tiff",width=7,height=25,units="in",res=600,compression="lzw")
hc<-hclust(as.dist(Dist),method="average")
plot(hc)
dev.off()
#
# Save a jpeg version too
#
img<-image_read("Comparisons/Correlation/Tree_ORA_cor.tiff")
image_write(img,path="Comparisons/Correlation/Tree_ORA_cor.jpeg",format="jpeg")
#
# For each disease, build a plot -Log10P-Cor
#
distance<-"Correlation"
alpha<-0.05
Distance_plot(Score,LogP,distance,alpha)
#
#-----------------------------------------------------------------------------
# Study similarity between diseases by network-based separation between disease 
# A and B (SAB), according to Menche J 2015 
# (https://pmc.ncbi.nlm.nih.gov/articles/PMC4435741/)
# SAB = mean(dAB) - (mean(dAA)-mean(dBB))/2
#-----------------------------------------------------------------------------
#
# Create output folders, if absent
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Comparisons/Separation")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
current_dir<-getwd()
folder_path<-file.path(current_dir,"Comparisons/Separation/Distributions")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
file.name<-"Comparisons/Separation/SAB_rand.rds"
if (!file.exists(file.name)) {
  #
  # This list will store the SAB between each disease module and all the random
  # disease modules
  #
  SAB_rand.list<-list()
  #
  # We need the matrix of all the disease genes
  #
  gene.matrix.full<-readRDS("Modules/myDiseaseGenes_net.rds") # read complete network
  graph.full<-graph_from_adjacency_matrix(gene.matrix.full,mode="undirected",weighted=T)
  #
  # For each disease calculate all the shortest distances with random diseases. 
  # This task is performed in parallel
  #
  Ncores<-detectCores()-1
  cl<-makeCluster(Ncores) # set number of parallel workers
  clusterExport(cl,c("SAB_func","graph.full","Disease_list","Nrand"))
  clusterEvalQ(cl,library(igraph))
  #
  for (i in 1:length(Disease_list)) {
    #
    print(Sys.time())
    print(paste0("Comparing ",names(Disease_list)[i]," and random diseases"))
    #
    # Read disease matrix and build graph
    #
    gene.matrix.AA<-readRDS(file=paste0("Modules/",names(Disease_list)[i],".rds"))
    graph.AA<-graph_from_adjacency_matrix(gene.matrix.AA,mode="undirected",weighted=T)
    #
    # Prepare matrix for SABs 
    #
    SAB<-matrix(data=NA,nrow=length(Disease_list),ncol=Nrand)
    row.names(SAB)<-names(Disease_list)
    colnames(SAB)<-seq(1,Nrand,1)
    #
    clusterExport(cl,c("graph.AA"))
    #
    # start parallel computing
    #
    SAB_rows<-parLapply(cl, which(seq(1:length(Disease_list))!=i), function(j) {
      #
      # Read the gene matrices for the random disease generated for disease j
      #
      gene.matrix_rand<-readRDS(paste0("Random/",names(Disease_list)[j],"_1000_gene_matrix.rsd"))
      #
      # SAB with random diseases
      #
      row_j<-c() 
      for (k in 1:Nrand) {
        #
        # For random disease j-k
        #
        gene.matrix.BB<-gene.matrix_rand[[k]]
        graph.BB<-graph_from_adjacency_matrix(gene.matrix.BB,mode="undirected",weighted=T)
        #
        # Calculate SAB and store it
        #
        row_j[k]<-SAB_func(graph.AA,graph.BB,graph.full)
        #
      }
      return(list(row_j,j))
    })
    #
    # build SAB from the calculate rows
    #
    for (h in 1:length(SAB_rows)) {
      j<-SAB_rows[[h]][[2]]
      SAB[j,]<-SAB_rows[[h]][[1]]
    }
    #
    # save matrix SAB as element i of a list
    #
    SAB_rand.list[[i]]<-SAB
    #
    # Save partial results
    #
    saveRDS(SAB_rand.list,file="Comparisons/SAB_rand.rds")
  }      
}
#
# Calculate separation between two diseases and test for significance
#
SAB_rand.list<-readRDS(file="Comparisons/Separation/SAB_rand.rds")
gene.matrix.full<-readRDS("Modules/myDiseaseGenes_net.rds") # read complete network
graph.full<-graph_from_adjacency_matrix(gene.matrix.full,mode="undirected",weighted=T)
#
Score<-matrix(data=NA,nrow=length(Disease_list),ncol=length(Disease_list))
colnames(Score)<-names(Disease_list)
rownames(Score)<-names(Disease_list)
LogP<-Score
for (i in 1:length(Disease_list)) {
  for (j in i:length(Disease_list)) {
    if (i!=j) {
      #
      print(paste0("Comparing ",names(Disease_list)[i]," and ",names(Disease_list)[j]))
      #
      # Select disease modules 
      #
      gene.matrix.AA<-readRDS(file=paste0("Modules/",names(Disease_list)[i],".rds"))
      graph.AA<-graph_from_adjacency_matrix(gene.matrix.AA,mode="undirected",weighted=T)
      gene.matrix.BB<-readRDS(file=paste0("Modules/",names(Disease_list)[j],".rds"))
      graph.BB<-graph_from_adjacency_matrix(gene.matrix.BB,mode="undirected",weighted=T)
      #
      # Calculate SAB and store it
      #
      Score[i,j]<-SAB_func(graph.AA,graph.BB,graph.full)
      Score[j,i]<-Score[i,j]
      #
      # Separations of Disease i with random Diseases j
      #
      Random.S_i<-SAB_rand.list[[i]][j,]
      #
      # Correlations of Disease j with random Diseases
      #
      Random.S_j<-SAB_rand.list[[j]][i,]
      #
      # Merge
      #
      Random.S<-c(Random.S_i,Random.S_j)
      #
      # Calculate p value for separation
      #
      p<-P_lower(value=Score[i,j],distribution=Random.S)
      #
      if (!is.infinite(log10(p))) {
        LogP[i,j]<--log10(p)
        LogP[j,i]<-LogP[i,j]
      }
      if (TRUE) {
        #
        # Plot null distribution
        #
        tiff(paste0("Comparisons/Separation/Distributions/Hist_",names(Disease_list)[i],"_",
                    names(Disease_list)[j],".tiff"),width=7,height=7,units="in",
             res=600,compression="lzw")
        hist(Random.S)
        abline(v=Score[i,j],col="red",lty=1)
        dev.off()
      }
    } 
  } 
}
#
# Save results
#
Score.list<-list(Score,LogP)
names(Score.list)<-c("Score","LogP")
saveRDS(Score.list,file="Comparisons/Separation/Score.rds")
#
# Read results
#
Score.list<-readRDS("Comparisons/Separation/Score.rds")
Score<-Score.list[[1]]
LogP<-Score.list[[2]]
#
# Convert to distance
#
MaxDist<-max(Score,na.rm=T)
MinDist<-min(Score,na.rm=T)
Dist<-(Score-MinDist)/(MaxDist-MinDist) # normalization between zero and one
#
# Plot dendrogram
#
tiff("Comparisons/Separation/Tree_SAB.tiff",width=7,height=25,units="in",res=600,compression="lzw")
hc<-hclust(as.dist(Dist),method="average")
plot(hc)
dev.off()
#
# Save a jpeg version too
#
img<-image_read("Comparisons/Separation/Tree_SAB.tiff")
image_write(img,path="Comparisons/Separation/Tree_SAB.jpeg",format="jpeg")
#
# For each disease, build a plot -Log10P-Cor
#
distance<-"Separation"
alpha<-0.05
Distance_plot(Score,LogP,distance,alpha)
#
#-----------------------------------------------------------------------------
# Comparison between scores
#-----------------------------------------------------------------------------
#
# Read scores and put in a list
#
Scores<-list()
Score.list<-readRDS("Comparisons/Jaccard/Score.rds")
Scores[[1]]<-Score.list[[1]]
Score.list<-readRDS("Comparisons/Correlation/Score.rds")
Scores[[2]]<-Score.list[[1]]
Score.list<-readRDS("Comparisons/Separation/Score.rds")
Scores[[3]]<-Score.list[[1]]
names(Scores)<-c("Jaccard","Correlation","Separation")
#
# Convert into a vector 
#
for (i in 1:3) {
  Scores[[i]]<-Scores[[i]][lower.tri(Scores[[i]])]
}
#
# Pair-wise regressions and plot
#
for (p in 1:3) {
  Score_plot(Scores,i=2,j=1,ord=p)  
  Score_plot(Scores,i=3,j=1,ord=p)
  Score_plot(Scores,i=3,j=2,ord=p)  
}
#
#-----------------------------------------------------------------------------
# For each gene, retrieve degree in disease and degree in STRING
#-----------------------------------------------------------------------------
#
for (j in 1:length(Disease_list)) {
  #
  if (names(Disease_list)[j]!="Chronic Fatigue Syndrome") {
    print(paste0("I am studying a few network properties for ",names(Disease_list)[j]))
    #
    mytargets<-fread(file=paste0("Modules/",names(Disease_list)[j],".csv"),sep=",")
    gene.matrix<-readRDS(file=paste0("Modules/",names(Disease_list)[j],".rds"))
    graph<-graph_from_adjacency_matrix(gene.matrix,mode="undirected",weighted=TRUE)
    degrees<-igraph::degree(graph) # add degree
    #
    mytargets$N.PPI<-rep(NA,nrow(mytargets))
    mytargets$degree<-rep(NA,nrow(mytargets))
    for (i in 1:nrow(mytargets)) { # add number of interactions
      index<-which(names(degrees)==mytargets$target.approvedSymbol[i])
      if (length(index)==1) mytargets$degree[i]<-degrees[index]
      N<-nrow(STRING(mytargets$target.approvedSymbol[i])) # by STRING API
      if (length(N)==1) mytargets$N.PPI[i]<-N
    }
    #
    write.table(mytargets,file=paste0("Modules/",names(Disease_list)[j],".csv"),
                sep=",",row.names=FALSE,quote=FALSE)
  }
}
#
#-----------------------------------------------------------------------------
# Study the correlation between the degree of a gene and its degree in STRING
#-----------------------------------------------------------------------------
#
mydata<-data.frame(N.PPI=NA,degree=NA)
#
for (j in 1:length(Disease_list)) {
  file.name<-paste0("Modules/",names(Disease_list)[j],".csv")
  if (names(Disease_list)[j]!="Chronic Fatigue Syndrome") {
    if (file.exists(file.name)) {
      mytargets<-fread(file=file.name)  
      df<-subset.data.frame(mytargets,select=c("N.PPI","degree"))
      mydata<-rbind(mydata,df)
    }
  }
}
#
model<-lm(degree~N.PPI,data=mydata)
summary(model)
#
jpeg("Modules/interactors_degree.jpeg",width=800,height=600)
plot(mydata$N.PPI,mydata$degree,xlab="N.PPI",ylab="degree",pch=16)
abline(model,lwd=2,col="red")
dev.off()
