# file name: Module_Func
#
current_dir<-getwd() # current directory
set.seed(12345) # to make results reproducible 
#
#-------------------------------------------------------------------------------
# Packages and files
#-------------------------------------------------------------------------------
#
# Load libraries
#
library(dplyr) 
library(rentrez) 
library(httr)
library(jsonlite)
library(curl)
library(biomaRt)
library(stringr)
library(DOSE)
library(stats)
library(rstatix)
library(ReactomePA)
library(igraph)
library(MASS)
library(data.table)
library(clusterProfiler)
library(stats)
library(org.Hs.eg.db)
library(pathview)
library(enrichplot)
library(GOplot)
library(readxl) 
library(writexl)
library(TissueEnrich)
library(magick)
library(yaml)
library(calibrate)
library(Matrix)
library(parallel)
library(SNFtool)
library(mclust)
#
#-------------------------------------------------------------------------------
# Build (if necessary) and load STRING database
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/STRING")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Load data for all species, then keep only H. sapiens
#
url<-paste0("https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz")
destfile<-"9606.protein.links.v12.0.txt.gz"
file_path<-file.path(current_dir,"Data/STRING/",destfile)
if(!file.exists(file_path)) {
  print("Downloading STRING database...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}   
url<-paste0("https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz")
destfile<-"9606.protein.info.v12.0.txt.gz"
file_path<-file.path(current_dir,"Data/STRING/",destfile)
if(!file.exists(file_path)) {
  print("Downloading STRING database...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}   
#
print("Loading STRING database...")
#
file.name<-paste0(current_dir,"/Data/STRING/9606.protein.links.v12.0.txt.gz")
STRING.matrix<-fread(file.name,sep=" ") # PPI scores
STRING.matrix<-subset.data.frame(STRING.matrix,
                                 combined_score>=STRING.co*1000) # restrict to significant PPI
#
file.name<-paste0(current_dir,"/Data/STRING/9606.protein.info.v12.0.txt.gz")
STRING.names<-fread(file.name,sep="\t") # STRING preferred names
STRING.names<-STRING.names[,1:2] # we only need the first two columns
colnames(STRING.names)[1]<-"string_protein_id"
#
#-------------------------------------------------------------------------------
# Build (if necessary) and load NCBI database
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/NCBI")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
url<-paste0("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz")
destfile<-"gene_info.gz"
file_path<-file.path(current_dir,"Data/NCBI/",destfile)
if(!file.exists(file_path)&!file.exists("Data/NCBI/human_gene_info.gz")) {
  test<-0
  while(test==0) {
    print("Downloading NCBI database...")
    attempt<-try(curl_download(url,file_path,quiet=FALSE),silent=T)
    if (class(attempt)!="try-error") {
      test<-1
    }
  }
}   
#
file.name<-"Data/NCBI/human_gene_info.gz"
if (!file.exists(file.name)) {
  print("Editing NCBI database...")
  file.name<-"Data/NCBI/gene_info.gz"
  NCBI.names<-fread(file.name,sep="\t") 
  colnames(NCBI.names)[1]<-"tax_id"
  NCBI.names<-subset.data.frame(NCBI.names,tax_id==9606)
  NCBI.names<-subset.data.frame(NCBI.names,select=c("GeneID","Symbol","Synonyms","description"))
  colnames(NCBI.names)<-c("NCBI.id","name","Synonyms","description")
  file.name<-"Data/NCBI/human_gene_info.gz"
  fwrite(NCBI.names,file.name,sep="\t") 
  file.remove("Data/NCBI/gene_info.gz")
} else {
  print("Loading NCBI database...")
  NCBI.names<-fread(file.name,sep="\t")  
}
gc()
#
#-------------------------------------------------------------------------------
# Convert gene symbol to NCBI ID (when possible)
#-------------------------------------------------------------------------------
#
Symbol2NCBI<-function(gene.symbol) {
  #
  # print(paste("Seraching NCBI gene ID for",gene.symbol))
  retries<-100
  for (i in 1:retries) {
    # Use tryCatch to catch errors and warnings
    result<-tryCatch(
      {
        # Perform the search with entrez_search
        NCBI<-lapply(gene.symbol,function(gene) {
          search<-entrez_search(db="gene",term=paste(gene,"[Gene Name] AND human[Organism]"))
        })
      },
      error = function(e) {
        # Check if the error is HTTP 500
        if (grepl("HTTP failure: 500", e$message)) {
          message(paste("Attempt", i, "failed with HTTP 500. Retrying..."))
          return(NA)
        } else {
          stop(e)  # Stop the loop for other errors
        }
      }
    )
    # If result is not NULL, the search was successful
    if (!is.null(result)) {
      if (is.list(result$ids)) {
        return(NA)
      } else if (is.list(NCBI[[1]]$ids[1])) {
        return(NA)
      } else {
        NCBI<-NCBI[[1]]$ids[1]
        return(NCBI)
      }
    }
    # If we've reached this point, the search failed with HTTP 500.
    # Wait before retrying
    Sys.sleep(0.1)
  }
}
#
#-------------------------------------------------------------------------------
# Convert gene symbol to NCBI ID using local data base (not API)
#-------------------------------------------------------------------------------
#
Symbol2NCBI.db<-function(gene.symbol) {
  #
  # Search in "name" first
  #
  df<-subset.data.frame(NCBI.names,name==gene.symbol)
  if (nrow(df)>0) {
    return(as.character(df$NCBI.id[1]))
  } else { # Search among synonyms
    index<-which(grepl(paste0("\\b",gene.symbol,"\\b"),NCBI.names$Synonyms)) 
    if (length(index)==0) {
      return(NA)
    } else {
      return(as.character(NCBI.names$NCBI.id[index[1]]))
    }
  }
}
#
#-------------------------------------------------------------------------------
# Convert gene NCBI ID to gene symbol  (when possible)
#-------------------------------------------------------------------------------
#
NCBI2Symbol<-function(gene.NCBI.id) {
  #
  retries<-100
  for (i in 1:retries) {
    # Use tryCatch to catch errors and warnings
    result<-tryCatch(
      {
        # Perform the search with entrez_search
        entrez_summary(db ="gene",id=gene.NCBI.id)
      },
      error = function(e) {
        # Check if the error is HTTP 500
        if (grepl("HTTP failure: 500", e$message)) {
          message(paste("Attempt", i, "failed with HTTP 500. Retrying..."))
          return(NULL)
        } else {
          stop(e)  # Stop the loop for other errors
        }
      }
    )
    # If result is not NULL, the search was successful
    if (!is.null(result)) {
      if (length(result$name)==0) {
        return(NA)
      } else {
        return(result$name)
      }
    }
    # If we've reached this point, the search failed with HTTP 500.
    # Wait before retrying
    Sys.sleep(0)
  }
}  
#
#-------------------------------------------------------------------------------
# This function asks STRING API its preferred name for a gene
#-------------------------------------------------------------------------------
#
STRING.name<-function(gene) {
  #
  species_id<-9606 # Homo sapiens
  base_url<-"https://string-db.org/api/json/get_string_ids"
  identifiers<-paste(gene,collapse="%0d")
  response<-httr::GET(base_url,query=list(identifiers=identifiers,species=species_id)) 
  #
  if (status_code(response)==200) {
    response_content<-rawToChar(response$content)
    Encoding(response_content)<-"UTF-8"
    data<-jsonlite::fromJSON(response_content) # parse
    if (is.data.frame(data)) {
      preferred.name<-data$preferredName
    } else {
      preferred.name<-NA
    }
  } else {
    stop("Failed to retrieve data. Please check the NCBI IDs and your internet connection.")
  }
  return(preferred.name)
}
#
#-------------------------------------------------------------------------------
# This function finds genes that interact with the input gene, 
# according to STRING API
#-------------------------------------------------------------------------------
#
STRING<-function(gene) {
  #
  Sys.sleep(0.1)
  species_id<-9606 # Homo sapiens
  #
  # find STRING's identifier for the input gene
  #
  gene_name<-STRING.name(gene)
  if (is.na(gene_name)) gene_name<-"" # necessary, otherwise STRING search gene NA
  #
  # find interacting genes
  #
  base_url<-"https://string-db.org/api/json/network"
  response<-httr::GET(base_url,
                      query=list(identifiers=gene_name,species=species_id,
                                 required_score=STRING.co*1000,limit=500))
  #
  if (status_code(response)==200) {
    # 
    response_content<-rawToChar(response$content)
    Encoding(response_content)<-"UTF-8"
    data<-jsonlite::fromJSON(response_content) # parse
    #
    if (is.data.frame(data)>0) {
      data<-subset.data.frame(data,score>=STRING.co)
      data<-subset.data.frame(data,preferredName_A==gene|preferredName_B==gene)
      #
      # Force gene as protein1
      #
      index<-which(data$preferredName_A!=gene)
      data$preferredName_B[index]<-data$preferredName_A[index]
      data$preferredName_A[index]<-gene
      #
      # Kepp only gene B and score
      #
      data<-data.frame(name=data$preferredName_B,score=data$score)
      #
      # we keep for each gene only the entry with highest score
      #
      unique.genes<-unique(data$name) 
      interacting_genes<-data[1:length(unique.genes),]
      for (k in 1:length(unique.genes)) {
        temp.df<-subset.data.frame(data,name==unique.genes[k])
        temp.df<-temp.df[order(temp.df$score,decreasing=T),]
        interacting_genes[k,]<-temp.df[1,]
      }
      return(interacting_genes)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}
#
#-------------------------------------------------------------------------------
# This function finds genes that interact with the input gene, 
# using STRING database (not API)
#-------------------------------------------------------------------------------
#
STRING2<-function(gene) {
  #
  # Find all the interacting genes
  #
  df<-subset.data.frame(STRING.names,preferred_name==gene)
  gene<-df$string_protein_id
  df<-subset.data.frame(STRING.matrix,protein1==gene)
  df<-subset.data.frame(df,combined_score>=STRING.co*1000)
  #
  # edit the output
  #
  df<-df[,-1]
  colnames(df)<-c("name","score")
  #
  if (nrow(df)>0) {
    for (i in 1:nrow(df)) {
      temp.df<-subset.data.frame(STRING.names,string_protein_id==df$name[i])
      df$name[i]<-temp.df$preferred_name
    }
    df$score<-df$score/1000
    return(df)
  } else {
    return(NA)
  }
}
#
#-------------------------------------------------------------------------------
# This function find PPI score between two given genes 
# using STRING database (not API)
#-------------------------------------------------------------------------------
#
STRING.db<-function(gene1,gene2) {
  #
  # find PPI score between the two genes
  #
  if (!is.na(gene1)&!is.na(gene2)) {
    # find identifiers
    df<-subset.data.frame(STRING.names,preferred_name==gene1)
    gene1<-df$string_protein_id
    df<-subset.data.frame(STRING.names,preferred_name==gene2)
    gene2<-df$string_protein_id
    # find PPI score
    df<-subset.data.frame(STRING.matrix,protein1==gene1)
    df<-subset.data.frame(df,protein2==gene2)
    if (nrow(df)==1) {
      score<-df$combined_score/1000
    } else {
      score<-0
    }
    if (score<STRING.co) {
      score<-0
    } 
    return(score)
  } else {
    return(NA)
  }
}
#
#-------------------------------------------------------------------------------
# For genes that appears more than once, this function generate a single row with
# all the information and with the highest score associated to each gene
#-------------------------------------------------------------------------------
#  
ListCollapse<-function(all.genes.zero) {
  #
  # Indicate predicted genes and predictors
  #
  unique.genes<-unique(all.genes.zero$NCBI.id) # NCBI IDs of unique genes
  all.genes.unique<-all.genes.zero[1:length(unique.genes),]
  for (i in 1:length(unique.genes)) {
    temp<-subset.data.frame(all.genes.zero,NCBI.id==unique.genes[i])
    all.genes.unique$NCBI.id[i]<-temp$NCBI.id[1]
    all.genes.unique$name[i]<-temp$name[1]
    all.genes.unique$list.count[i]<-length(unique(temp$list.name))
    all.genes.unique$list.name[i]<-paste0(unique(temp$list.name),collapse="/")
  }
  return(all.genes.unique)
} 
#
#-------------------------------------------------------------------------------
# This function build adjacency matrix from a list of genes using STRING database 
# (not API). 
#-------------------------------------------------------------------------------
#
GeneMatrix<-function(all.genes) {
  #
  # build a STRING database with only genes of interest and with score above STRING.co
  #
  df1<-STRING.names
  colnames(df1)[2]<-"name"
  df.names<-merge(all.genes,df1,by="name",all.y=F)
  df.names<-subset.data.frame(df.names,select=c("name","string_protein_id"))
  dt.names<-as.data.table(df.names) # data tables should be faster
  remove(df1,df.names)
  #
  df.matrix<-subset.data.frame(STRING.matrix,combined_score>=STRING.co*1000)
  df.matrix<-subset.data.frame(df.matrix,protein1%in%dt.names$string_protein_id)
  df.matrix<-subset.data.frame(df.matrix,protein2%in%dt.names$string_protein_id)
  dt.matrix<-as.data.table(df.matrix) # data tables should be faster
  remove(df.matrix)
  #
  # build gene.matrix
  #
  nodes<-dt.names$string_protein_id
  i<-match(dt.matrix$protein1,nodes)
  j<-match(dt.matrix$protein2,nodes)
  #
  Msp<-sparseMatrix(
    i = i, j = j, x = dt.matrix$combined_score/1000,
    dims = c(length(nodes), length(nodes)),
    dimnames = list(nodes, nodes)
  )
  #
  gene.matrix<-as.matrix(Msp)
  #
  id_to_name <- dt.names$name
  names(id_to_name) <- dt.names$string_protein_id
  rownames(gene.matrix) <- id_to_name[rownames(gene.matrix)]
  colnames(gene.matrix) <- id_to_name[colnames(gene.matrix)]
  #
  # Test symmetry and end
  #
  test<-isSymmetric(gene.matrix)
  if (!test) {
    stop("The adjacency matrix is not symmetric; there is an error!")
  } else {
    return(gene.matrix)  
  }
}
#
#-------------------------------------------------------------------------------
# This function tests a sample of PPI score of a given gene.matrix using STRING 
# API
#-------------------------------------------------------------------------------
#
GeneTest<-function(gene.matrix) {
  N<-nrow(gene.matrix)
  genes<-rownames(gene.matrix)
  for (t in 1:300) {
    i<-sample(seq(1,N),1)
    j<-sample(seq(1,i-1),1)
    a<-STRING(genes[i])
    index<-which(a$name==genes[j])
    if (length(index)==0) {
      score<-0
    } else {
      score<-a$score[index]
    }
    # print(paste(t,score,gene.matrix[i,j]))
    if (score!=gene.matrix[i,j]) {
      stop("Your function GeneMatrix has a bug")
    }
  }
}
#
#-------------------------------------------------------------------------------
# This function performs Over-representation Analysis over several databases
#-------------------------------------------------------------------------------
#
ORA.fun<-function(all.genes) {
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis with KEGG (Kyoto Encyclopedia of Genes and Genomes - 
  # biological pathways)
  #-------------------------------------------------------------------------------
  #
  KEGG_results<-enrichKEGG(gene=all.genes$NCBI.id,organism='hsa',universe=myuniverse)
  KEGG.ORA<-KEGG_results@result
  KEGG.ORA<-KEGG.ORA[,3:ncol(KEGG.ORA)]
  KEGG.ORA$method<-rep("KEGG",nrow(KEGG.ORA))
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis with GO (Gene Ontology - functional role of genes within 
  # biological processes, molecular functions, and cellular components)
  #-------------------------------------------------------------------------------
  #
  GO_results<-enrichGO(gene=all.genes$NCBI.id,OrgDb=org.Hs.eg.db,universe=myuniverse,ont="ALL")
  GO.ORA<-GO_results@result
  GO.ORA$method<-rep("GO",nrow(GO.ORA))
  #
  #-------------------------------------------------------------------------------
  # Over-representation Analysis using Reactome 
  #-------------------------------------------------------------------------------
  #
  pathway_results<-enrichPathway(gene=all.genes$NCBI.id,organism="human",universe=myuniverse)
  Rea.ORA<-pathway_results@result
  Rea.ORA$method<-rep("Reactome",nrow(Rea.ORA))
  #
  #-------------------------------------------------------------------------------
  # Merge results from Over-representation Analysis and keep only significant ones
  #-------------------------------------------------------------------------------
  #
  if ("ONTOLOGY"%in%colnames(GO.ORA)) GO.ORA<-GO.ORA[,!names(GO.ORA) %in% "ONTOLOGY"]
  ORA<-rbind(KEGG.ORA,GO.ORA,Rea.ORA)
  #
  return(ORA)
}
#
#-------------------------------------------------------------------------------
# perform ORA with respect to tissue expression
#-------------------------------------------------------------------------------
# 
Tissue.ORA<-function(all.genes,pl) {
  gs<-GeneSet(geneIds=all.genes$name,organism='Homo Sapiens',geneIdType=SymbolIdentifier())
  bk<-GeneSet(geneIds=STRING.names$preferred_name,organism='Homo Sapiens',geneIdType=SymbolIdentifier())
  output<-teEnrichment(gs,rnaSeqDataset=1,backgroundGenes=bk,multiHypoCorrection=F) # 1 for HPA, 2 for GTEx
  seEnrichmentOutput<-output[[1]]
  enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names=rowData(seEnrichmentOutput)[,1]),colData(seEnrichmentOutput)[,1])
  enrichmentOutput$Tissue<-row.names(enrichmentOutput)
  enrichmentOutput$p.adjust<-10^-enrichmentOutput$Log10PValue # it is already adjusted, by default!
  enrichmentOutput$zScore<-qnorm(1 - enrichmentOutput$p.adjust/2) # unsigned zScore
  #
  # plot
  #
  if (pl) {
    tiff(paste0("Modules/",names(Disease_list)[j],"_Tissue_ORA.tiff"),width=10,
         height=6,units="in",res=600,compression="lzw")
    p<-ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,
                                   label = Tissue.Specific.Genes,fill = Tissue))+
      geom_bar(stat = 'identity')+
      geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
      labs(x='', y = '-LOG10(P-Value)')+
      theme_bw()+
      theme(legend.position='none')+
      theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
              element_text(size=15))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.grid.major= element_blank(),panel.grid.minor = element_blank())
    print(p)
    dev.off()
    #
    # Save a jpeg version too
    #
    img<-image_read(paste0("Modules/",names(Disease_list)[j],"_Tissue_ORA.tiff"))
    image_write(img,path=paste0("Modules/",names(Disease_list)[j],"_Tissue_ORA.jpeg"),format="jpeg")
  } 
  #
  return(enrichmentOutput)
}
#
#-------------------------------------------------------------------------------
# Prepare background genes
#-------------------------------------------------------------------------------
#
file.name<-"Data/My_Universe.tsv"
if (!file.exists(file.name)) {
  print("Building background for Over-representation analysis. This may take a while...")
  STRING.names$NCBI.id<-rep(NA,nrow(STRING.names))
  for (i in 1:nrow(STRING.names)) {
    STRING.names$NCBI.id[i]<-Symbol2NCBI.db(STRING.names$preferred_name[i])
  }
  myuniverse<-as.character(STRING.names$NCBI.id)  
  myuniverse<-na.omit(myuniverse)
  write.table(myuniverse,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
} else {
  myuniverse<-fread(file.name)
  myuniverse<-as.character(myuniverse$x)
}
#
#-------------------------------------------------------------------------------
# Retrieve GWAS study ID
#-------------------------------------------------------------------------------
#
GWAS4Gene<-function(efo_id,target_id,L2G_cutoff,sample_cutoff) {
  #
  # GraphQL endpoint of Open Target API
  #
  api_url<-"https://api.platform.opentargets.org/api/v4/graphql"
  #
  # Build the query
  #
  variables<-list(efoId=efo_id,ensemblId=target_id,size=100)
  graphql_query<- '
  query GwasCredibleSetsQuery($ensemblId: String!, $efoId: String!, $size: Int!) {
    target(ensemblId: $ensemblId) {
      approvedSymbol
    }
    disease(efoId: $efoId) {
      id
      name
      gwasCredibleSets: evidences(
        ensemblIds: [$ensemblId]
        enableIndirect: true
        datasourceIds: ["gwas_credible_sets"]
        size: $size
      ) {
        count
        rows {
          disease {
            id
            name
          }
          credibleSet {
            studyLocusId
            study {
              traitFromSource
              id
              projectId
              publicationFirstAuthor
              publicationDate
              pubmedId
              nSamples
            }
            variant {
              id
              chromosome
              position
              referenceAllele
              alternateAllele
            }
            pValueMantissa
            pValueExponent
            beta
            finemappingMethod
            confidence
          }
          score
        }
      }
    }
  }  
  '
  request_body<-list(query=graphql_query,variables=variables)
  #
  # Request data
  #
  response<-httr::POST(
    url = api_url,
    body = jsonlite::toJSON(request_body,auto_unbox=T),
    encode = "json",
    httr::content_type_json()
  )
  #
  # Verify the status of the request 
  #
  if (httr::status_code(response)!=200) {
    stop(paste("Errore nella richiesta API. Codice di stato:",httr::status_code(response)))
  }
  #
  # Parse JSON
  #
  json_content<-httr::content(response,"text",encoding="UTF-8")
  data<-jsonlite::fromJSON(json_content,flatten=T)
  gwas_list<-data$data$disease$gwasCredibleSets$rows
  #
  # Edit and filtering
  # 
  gwas_list<-subset.data.frame(gwas_list,credibleSet.pValueExponent<=-8) # keep only significant results
  gwas_list<-subset.data.frame(gwas_list,!(credibleSet.pValueExponent==-8&credibleSet.pValueMantissa>=0.5)) # keep only significant results
  gwas_list<-subset.data.frame(gwas_list,credibleSet.study.nSamples>=sample_cutoff) # keep only large studies
  gwas_list<-subset.data.frame(gwas_list,score>=L2G_cutoff) # keep only high confidence
  #
  return(gwas_list)
}
#
#-------------------------------------------------------------------------------
# Retrieve genetic targed for disease
#-------------------------------------------------------------------------------
#
Targets4Disease<-function(efo_id,L2G_cutoff,clingen_cutoff,geneburden_cutoff,sample_cutoff) {
  #
  # GraphQL endpoint of Open Target API
  #
  api_url<-"https://api.platform.opentargets.org/api/v4/graphql" 
  #
  # Build the query
  #
  variables<-list(efoId=efo_id)
  graphql_query <- '
  query TargetsForDisease($efoId: String!) {
    disease(efoId: $efoId) {
      name
      associatedTargets(page: {index: 0, size: 1000}) {
        count
        rows {
          target {
            id
            approvedSymbol
            biotype
          }
          score
          datasourceScores {
            id
            score
          }
        }
      }
    }
  }
  '
  request_body<-list(query=graphql_query,variables=variables)
  #
  # Request data
  #
  response<-httr::POST(
    url=api_url,
    body=jsonlite::toJSON(request_body,auto_unbox=T), # 'auto_unbox' è cruciale
    encode="json",
    httr::content_type_json()
  )
  #
  # Verify the status of the request 
  #
  if (httr::status_code(response)!=200) {
    stop(paste("Errore nella richiesta API. Codice di stato:",httr::status_code(response)))
  }
  #
  # Parse JSON
  #
  json_content<-httr::content(response,"text",encoding="UTF-8")
  data<-jsonlite::fromJSON(json_content,flatten=T)
  targets_list<-data$data$disease$associatedTargets$rows
  #
  # Retain data from GWAS, Clingen, and GeneBurden
  #  
  targets_results<-targets_list
  targets_results$gwas_credible_sets<-rep(NA,nrow(targets_results))
  targets_results$clingen<-rep(NA,nrow(targets_results))
  targets_results$geneburden<-rep(NA,nrow(targets_results))
  for (i in 1:nrow(targets_results)) {
    df<-as.data.frame(targets_list$datasourceScores[i])
    index<-which(df$id=="gwas_credible_sets")
    if (length(index)>0) {
      targets_results$gwas_credible_sets[i]<-df$score[index]
    }
    index<-which(df$id=="clingen")
    if (length(index)>0) {
      targets_results$clingen[i]<-df$score[index]
    }
    index<-which(df$id=="geneburden")
    if (length(index)>0) {
      targets_results$geneburden[i]<-df$score[index]
    }
  }
  #
  # Edit and filtering
  # 
  targets_results<-subset.data.frame(targets_results,select=c("score","target.id","target.approvedSymbol",
                                                              "target.biotype","gwas_credible_sets",
                                                              "clingen","geneburden"))
  targets_results<-subset.data.frame(targets_results,target.biotype=="protein_coding")
  targets_results<-subset.data.frame(targets_results,!is.na(gwas_credible_sets)|!is.na(clingen)|!is.na(geneburden))
  targets_results<-subset.data.frame(targets_results,(!is.na(gwas_credible_sets)&
                                                        gwas_credible_sets>=L2G_cutoff&
                                                        is.na(clingen))|
                                       !is.na(clingen)) 
  targets_results<-subset.data.frame(targets_results,is.na(clingen)|
                                       clingen>=clingen_cutoff) # only high confidence
  targets_results<-subset.data.frame(targets_results,is.na(geneburden)|
                                       geneburden>=geneburden_cutoff) # only high confidence
  #
  # Add GWAS Study id(s) 
  # 
  if (nrow(targets_results)>0) {
    targets_results$GWAS.id<-rep(NA,nrow(targets_results))
    targets_results$L2G<-rep(NA,nrow(targets_results))
    for (i in 1:nrow(targets_results)) {
      if (!is.na(targets_results$gwas_credible_sets[i])) {
        target_id<-targets_results$target.id[i]
        df<-GWAS4Gene(efo_id,target_id,L2G_cutoff,sample_cutoff)
        if (nrow(df)>0) {
          targets_results$GWAS.id[i]<-paste0(df$credibleSet.study.id,collapse="/")
          targets_results$L2G[i]<-paste0(df$score,collapse="/")  
        }
      }
    }
    #
    # Find unique GWAS studies
    # 
    unique_studies<-paste0(targets_results$GWAS.id,collapse="/")
    unique_studies<-unique(str_split(unique_studies,"/")[[1]]) 
    unique_studies<-na.omit(unique_studies)
    unique_studies<-unique_studies[!unique_studies%in%c("NA","")]
    #
    # Add a column for each unique GWAS study
    #
    for (j in 1:length(unique_studies)) {
      df<-data.frame(rep(NA,nrow(targets_results)))
      colnames(df)<-unique_studies[j]                  
      targets_results<-cbind(targets_results,df)
    }
    #
    # Fill column with scores
    #
    index<-c()
    for (i in 1:nrow(targets_results)) {
      if (!is.na(targets_results$GWAS.id[i])) {
        df<-data.frame(GWAS.id=str_split(targets_results$GWAS.id[i],"/")[[1]],
                       L2G=str_split(targets_results$L2G[i],"/")[[1]])
        for (j in 1:nrow(df)) {
          col<-which(colnames(targets_results)==df$GWAS.id[j])
          targets_results[i,col]<-round(as.numeric(df$L2G[j]),7)
        }
      } else if (is.na(targets_results$GWAS.id[i])&is.na(targets_results$clingen[i])) {
        index<-c(index,i)
      }
    }
    #
    # Remove rows for genes with no supporting evidence (using these filters for GWAS4Gene)
    #
    if (length(index)>0)targets_results<-targets_results[-index,]
    #
    # Count number of GWAS for each target
    # 
    for (i in 1:nrow(targets_results)) {
      count<-unique(str_split(targets_results$L2G[i],"/")[[1]])
      count(which(!is.na(count)))
      targets_results$GWAS.count[i]<-length(count)
    }
  }
  #
  return(targets_results)
}
#
#-----------------------------------------------------------------------------
# Download ME/CFS module form Zhang S. 2025 
# We use supplementary table 2.
#-----------------------------------------------------------------------------
#
url<-paste0("https://www.medrxiv.org/content/medrxiv/early/2025/05/11/2025.04.15.25325899/DC9/embed/media-9.xlsx")
destfile<-"media-9.xlsx"
file_path<-file.path(current_dir,"Data",destfile)
if(!file.exists(file_path)) {
  print("Downloading data from Zhang S. et al. 2025 (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/)...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}   
#
#-----------------------------------------------------------------------------
# Download ME/CFS module form Sardell JM 2025
# We use supplementary table 3.
#-----------------------------------------------------------------------------
#
url<-paste0("https://www.medrxiv.org/content/medrxiv/early/2025/12/03/2025.12.01.25341362/DC2/embed/media-2.xlsx")
destfile<-"media-2.xlsx"
file_path<-file.path(current_dir,"Data",destfile)
if(!file.exists(file_path)) {
  print("Downloading data from Sardell JM. et al. 2025 (https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v2)...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}    
#
#-----------------------------------------------------------------------------
# Read ME/CFS module form Zhang S. 2025 (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/)
# We use supplementary table 2.
#-----------------------------------------------------------------------------
#
Module<-read_xlsx("Data/media-9.xlsx") 
Module<-subset.data.frame(Module,q_value<0.02) # 115 genes associated with ME/CFS  
#
# Use STRING preferred name for genes
#
gene.wo.STRING<-0
index<-c()
Module$name<-rep(NA,nrow(Module))
for (i in 1:nrow(Module)) {
  new.name<-STRING.name(Module$Gene[i])
  if (!is.na(new.name)) {
    Module$name[i]<-new.name
  } else if (is.na(new.name)) {
    print(paste(Module$Gene[i],"is not present in STRING's database"))
    gene.wo.STRING<-gene.wo.STRING+1
    index<-c(index,i)
  }
}
if (length(index)>0) Module<-Module[-index,] # only genes in STRING
#
# Edit
#
Module$NCBI.id<-rep(NA,nrow(Module))
Module$list.name<-rep("Zhang",nrow(Module))
#
# Add NCBI id
#
for (i in 1:nrow(Module)) {
  Module$NCBI.id[i]<-Symbol2NCBI.db(Module$name[i])
}
Zhang_module<-as.data.frame(Module[,c("name","NCBI.id","list.name")])
#
#-----------------------------------------------------------------------------
# Read ME/CFS module form Sardell JM 2025 (https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v2)
# We use supplementary table 3.
#-----------------------------------------------------------------------------
#
Module<-read_xlsx("Data/media-2.xlsx",sheet=3,skip=2) 
Module<-data.frame(Gene=Module$`Gene name`)
#
# Use STRING preferred name for genes
#
gene.wo.STRING<-0
index<-c()
Module$name<-rep(NA,nrow(Module))
for (i in 1:nrow(Module)) {
  new.name<-STRING.name(Module$Gene[i])
  if (!is.na(new.name)) {
    Module$name[i]<-new.name
  } else if (is.na(new.name)) {
    print(paste(Module$Gene[i],"is not present in STRING's database"))
    gene.wo.STRING<-gene.wo.STRING+1
    index<-c(index,i)
  }
}
if (length(index)>0) Module<-Module[-index,] # only genes in STRING
#
# Edit
#
Module$NCBI.id<-rep(NA,nrow(Module))
Module$list.name<-rep("PL",nrow(Module))
#
# Add NCBI id
#
for (i in 1:nrow(Module)) {
  Module$NCBI.id[i]<-Symbol2NCBI.db(Module$name[i])
}
PL_module<-as.data.frame(Module[,c("name","NCBI.id","list.name")])
#
#-------------------------------------------------------------------------------
# Merge gene lists for ME/CFS
#-------------------------------------------------------------------------------
#
all.genes.zero<-rbind(Zhang_module,PL_module)
Exper_list<-all.genes.zero
#
#-------------------------------------------------------------------------------
# Build a list with a single row for each gene
#-------------------------------------------------------------------------------
#
all.genes<-ListCollapse(all.genes.zero)
all.genes<-all.genes[,c("name","NCBI.id","list.name","list.count")]
write.table(all.genes,file="Data/all_genes_customCFS.tsv",sep="\t",row.names=FALSE,quote=FALSE)
#
#-------------------------------------------------------------------------------
# P value for experimental distribution, right tail
#-------------------------------------------------------------------------------
#
P_upper<-function(value,distribution) {
  #
  Num<-length(which(distribution>=value))
  Den<-length(which(!is.na(distribution)))
  p_upper<-(Num+1)/(Den+1)
  #
  return(p_upper)
  #
}
#
#-------------------------------------------------------------------------------
# P value for experimental distribution, right tail
#-------------------------------------------------------------------------------
#
P_lower<-function(value,distribution) {
  #
  Num<-length(which(distribution<=value))
  Den<-length(which(!is.na(distribution)))
  p_upper<-(Num+1)/(Den+1)
  #
  return(p_upper)
  #
}
#
#-------------------------------------------------------------------------------
# P value for fitted parametric distribution, right tail
#-------------------------------------------------------------------------------
#
P_upper_param<-function(value, distribution) {
  #
  fit<-MASS::fitdistr(distribution[!is.na(distribution)],"normal")
  p <-pnorm(value, mean = fit$estimate["mean"],
               sd   = fit$estimate["sd"],
               lower.tail = FALSE)
  return(p)
}
#
#-------------------------------------------------------------------------------
# P value for fitted parametric distribution, left tail
#-------------------------------------------------------------------------------
#
P_lower_param<-function(value, distribution) {
  #
  fit<-MASS::fitdistr(distribution[!is.na(distribution)],"normal")
  p <-pnorm(value, mean = fit$estimate["mean"],
            sd   = fit$estimate["sd"],
            lower.tail = TRUE)
  return(p)
}
#
#-------------------------------------------------------------------------------
# Calculate SAB = mean(dAB) - (mean(dAA)-mean(dBB))/2
#-------------------------------------------------------------------------------
#
SAB_func<-function(graph.AA,graph.BB,graph.full) {
  #
  # Calculate mean(dAA) and mean(dBB)
  #
  dAA<-mean_distance(graph.AA,weights=E(graph.AA)$weight,directed=F)
  dBB<-mean_distance(graph.BB,weights=E(graph.BB)$weight,directed=F)
  #
  # Calculate mean(dAB)
  #
  index<-which(V(graph.full)$name%in%V(graph.AA)$name)
  V.AA<-V(graph.full)[index]
  index<-which(V(graph.full)$name%in%V(graph.BB)$name)
  V.BB<-V(graph.full)[index]
  dAB<-distances(graph.full,v=V.AA,to=V.BB,mode="all",
                 weights=NULL,
                 algorithm="dijkstra")
  #
  # Set zero for nodes that are present in both diseases
  #
  index_r<-which(rownames(dAB)%in%colnames(dAB))
  index_c<-which(colnames(dAB)%in%rownames(dAB))
  dAB[index_r,]<-0
  dAB[,index_c]<-0
  #
  # mean of dAB matrix (removing infinite values)
  #
  dAB<-mean(dAB[is.finite(dAB)])
  #
  # Calculate SAB 
  #
  SAB<-dAB-(dAA+dBB)/2
  return(SAB)
}
#
#-------------------------------------------------------------------------------
# Distance plot for each disease, with p value
#-------------------------------------------------------------------------------
#
Distance_plot<-function(Score,LogP,distance,alpha) {
  #
  # For each disease, build a plot -Log10P-SAB
  #
  for (j in 1:length(Disease_list)) {
    #
    # Retrieve Indexes, p values, and names
    #
    mydata<-data.frame(score=Score[,j],scoreLP=LogP[,j],
                       Disease=Disease_abbr)
    mydata<-mydata[-j,]
    mydata<-mydata[order(mydata$scoreLP,decreasing=T), ]
    mydata<-na.omit(mydata)
    #
    # Create folders
    #
    current_dir<-getwd()
    folder_path<-file.path(current_dir,paste0("Comparisons/",distance,"/TIFF"))  
    if(!dir.exists(folder_path)) {
      dir.create(folder_path) 
    }
    current_dir<-getwd()
    folder_path<-file.path(current_dir,paste0("Comparisons/",distance,"/JPEG"))  
    if(!dir.exists(folder_path)) {
      dir.create(folder_path) 
    }
    #
    # Plot
    #
    file.name<-paste0("Comparisons/",distance,"/TIFF/",names(Disease_list)[j],"_",distance,".tiff")
    tiff(file.name,width=7,height=7,units="in",res=600,compression="lzw")
    plot(mydata$score,mydata$scoreLP,xlab=distance,ylab="-Log10P",
         pch=21,main=names(Disease_list)[j])
    pco<--log10(alpha)
    abline(h=pco,col="red",lty=2)
    #
    sig<-mydata$scoreLP>=pco
    if (any(sig)) {
      text(mydata$score[sig],
           mydata$scoreLP[sig],
           labels=mydata$Disease[sig],
           pos=3,           # above the point
           cex=0.6,
           offset=0.4)
    } else {
      text(mydata$score[1:5],
           mydata$scoreLP[1:5],
           labels=mydata$Disease[1:5],
           pos=3,           # above the point
           cex=0.6,
           offset=0.4)
    }
    dev.off()
    #
    # Save a jpeg version too
    #
    img<-image_read(file.name)
    file.name<-paste0("Comparisons/",distance,"/JPEG/",names(Disease_list)[j],"_",distance,".jpeg")
    image_write(img,path=file.name,format="jpeg")
  }
}
#
#-------------------------------------------------------------------------------
# Plot regressions between scores
#-------------------------------------------------------------------------------
#
Score_plot<-function(Scores,i,j,ord) {
  #
  score_x<-Scores[[i]]
  score_y<-Scores[[j]]
  if (ord==1) {
    model<-lm(score_y ~ score_x)  
  } else if (ord==2) {
    model<-lm(score_y ~ score_x + I(score_x^2))
  } else if (ord==3) {
    model<-lm(score_y ~ score_x + I(score_x^2) + I(score_x^3))
  }
  R2adj<-round(summary(model)$adj.r.squared,2)
  TS<-summary(model)$fstatistic[1]
  df1<-summary(model)$fstatistic[2]
  df2<-summary(model)$fstatistic[3]
  pval<-format(pf(TS,df1=df1,df2=df2,lower.tail=F),digit=2)
  # 
  file.name<-paste0("Comparisons/",names(Scores)[i],"_",
                    names(Scores)[j],"_",ord,".tiff")
  tiff(file.name,width=5,height=5,units="in",res=600,compression="lzw")
  plot(
    score_x,score_y,
    xlab = names(Scores)[i],
    ylab = names(Scores)[j],
    pch = 16,
    main = paste0("R2adj=",R2adj,", pval=",pval)
  )
  xg  <- seq(min(score_x),max(score_x),length.out=200)
  #
  # Predict with confidence interval
  #
  pred <- predict(model, newdata=data.frame(score_x=xg), interval="prediction", level=0.95)
  yg   <- pred[,"fit"]
  ylo  <- pred[,"lwr"]
  yhi  <- pred[,"upr"]
  #
  # Shade confidence interval
  #
  polygon(c(xg, rev(xg)), c(yhi, rev(ylo)),
          col=adjustcolor("red", alpha.f=0.2),
          border=NA)
  #
  # Regression line
  #
  lines(xg,yg,col="red",lwd=3)
  dev.off()
  #
  # Save a jpeg version too
  #
  img<-image_read(file.name)
  file.name<-paste0("Comparisons/",names(Scores)[i],"_",
                    names(Scores)[j],"_",ord,".jpeg")
  image_write(img,path=file.name,format="jpeg")
  #
}
#-------------------------------------------------------------------------------
# Permutation test for Adjusted Rand Index
#-------------------------------------------------------------------------------
#
ARI_pvalue<-function(labels1,labels2,Nperm) {
  #
  # Observed ARI
  #
  observed_ARI<-round(adjustedRandIndex(labels1,labels2),2)
  #
  # Permutation distribution
  #
  perm_ARI<-rep(NA,Nperm)
  for (i in 1:Nperm) {
    perm_ARI[i]<-adjustedRandIndex(sample(labels1),labels2)
  }
  #
  # P value (right tail: is observed ARI higher than chance?)
  #
  p<-(sum(perm_ARI >= observed_ARI)+1)/Nperm
  p<-format(p,digit=2)
  #
  return(data.frame(ARI=observed_ARI,p.value=p))
}
#
#-------------------------------------------------------------------------------
# Hierarchical clustering with all linkages, adjusted Rand Index with p value,
# and Dendrograms
#-------------------------------------------------------------------------------
#
hierarchical_func<-function(dist,method) {
  #
  # Hierarchical clustering, dendrogram plots, Rand Index
  #
  hc_list<-list()
  linkage_list<-c("complete","single","average","centroid","median","mcquitty","ward.D","ward.D2")
  for (j in 1:length(linkage_list)) {
    hc_list[[j]]<-hclust(as.dist(Dist),method=linkage_list[j])  
  }
  #
  # Calculate adjusted Rand Index with p-value (permutation test)
  #
  ICD10<-fread("ICD10_classification.csv")
  N_cluster<-length(unique(ICD10$ICD10_category))
  myHC<-data.frame()
  Nperm<-20000
  for (j in 1:length(linkage_list)) {
    labels1<-cutree(hc_list[[j]],k=N_cluster)
    labels2<-as.integer(factor(ICD10$ICD10_category))
    adjRI<-ARI_pvalue(labels1,labels2,Nperm)
    adjRI$method<-linkage_list[j]
    myHC<-rbind(myHC,adjRI)
    ICD10[[linkage_list[j]]]<-labels1
  }
  #
  # Save
  #
  write.table(myHC,file=paste0("Comparisons/",method,"/Hierarchical_clustering_ARI.csv"),
              sep=",",row.names=F,quote=F)
  write.table(ICD10,file=paste0("Comparisons/",method,"/Hierarchical_clustering_Clusters.csv"),
              sep=",",row.names=F,quote=F)
  #
  # Plot dendrograms
  #
  for (j in 1:length(linkage_list)) {
    #
    # TIFF version
    #
    name<-paste0("Comparisons/",method,"/Tree_",method,"_",linkage_list[j])
    tiff(paste0(name,".tiff"),width=7,height=25,units="in",res=600,
         compression="lzw")
    hc<-hclust(as.dist(Dist),method=linkage_list[j])
    plot(hc)
    dev.off()
    #
    # JPEG version
    #
    img<-image_read(paste0(name,".tiff"))
    image_write(img,path=paste0(name,".jpeg"),format="jpeg")
  }
}
