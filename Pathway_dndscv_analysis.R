library(GSA)

library(dndscv)

# pathway level dndscv analysis
# the function is modified from the dndscv function of "dndscv" package
# there are only two modifications: (1)only output global results of globaldnds (no gene level results); (2) the original function does not output p-value, modified to be able to output p-value;p-value is extracted from the glm model summary (summary_model$coefficients[,4]). all of others including default parameters are as same as original dndscv function.
#' pathway level dndscv analysis
#'
#' @param mutations, a data frame includes 5 columns: Sample_ID, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, allele can be point mutations, deletions and insertions (sampleID, chr, pos, ref, alt)
#' @param gene_list, a vector containing included gene names
#' @param refdb
#' @param sm, Substitution model
#' @param kc
#' @param cv
#' @param max_muts_per_gene_per_sample
#' @param max_coding_muts_per_sample
#' @param use_indel_sites
#' @param min_indels
#' @param maxcovs
#' @param constrain_wnon_wspl
#' @param outp
#' @param numcode
#' @param outmats
#'
#' @return
#' @export
#'
#' @examples
dndscv_pathway = function(mutations, gene_list = NULL, refdb = "hg19", sm = "192r_3w", kc = "cgc81", cv = "hg19", max_muts_per_gene_per_sample = 3, max_coding_muts_per_sample = 3000, use_indel_sites = T, min_indels = 5, maxcovs = 20, constrain_wnon_wspl = T, outp = 3, numcode = 1, outmats = F) {

  ## 1. Environment
  message("[1] Loading the environment...")

  mutations = mutations[,1:5] # Restricting input matrix to first 5 columns
  mutations[,c(1,2,3,4,5)] = lapply(mutations[,c(1,2,3,4,5)], as.character) # Factors to character
  mutations[[3]] = as.numeric(mutations[[3]]) # Chromosome position as numeric
  mutations = mutations[mutations[,4]!=mutations[,5],] # Removing mutations with identical reference and mutant base
  colnames(mutations) = c("sampleID","chr","pos","ref","mut")

  # Removing NA entries from the input mutation table
  indna = which(is.na(mutations),arr.ind=T)
  if (nrow(indna)>0) {
    mutations = mutations[-unique(indna[,1]),] # Removing entries with an NA in any row
    warning(sprintf("%0.0f rows in the input table contained NA entries and have been removed. Please investigate.",length(unique(indna[,1]))))
  }

  # [Input] Reference database
  refdb_class = class(refdb)
  if ("character" %in% refdb_class) {
    if (refdb == "hg19") {
      data("refcds_hg19", package="dndscv")
      if (any(gene_list=="CDKN2A")) { # Replace CDKN2A in the input gene list with two isoforms
        gene_list = unique(c(setdiff(gene_list,"CDKN2A"),"CDKN2A.p14arf","CDKN2A.p16INK4a"))
      }
    } else {
      load(refdb)
    }
  } else if("array" %in% refdb_class) {
    # use the user-supplied RefCDS object
    RefCDS = refdb
  } else {
    stop("Expected refdb to be \"hg19\", a file path, or a RefCDS-formatted array object.")
  }

  # [Input] Gene list (The user can input a gene list as a character vector)
  if (is.null(gene_list)) {
    gene_list = sapply(RefCDS, function(x) x$gene_name) # All genes [default]
  } else { # Using only genes in the input gene list
    allg = sapply(RefCDS,function(x) x$gene_name)
    nonex = gene_list[!(gene_list %in% allg)]
    #if (length(nonex)>0) { stop(sprintf("The following input gene names are not in the RefCDS database: %s", paste(nonex,collapse=", "))) }
    RefCDS = RefCDS[allg %in% gene_list] # Only input genes
    gr_genes = gr_genes[gr_genes$names %in% gene_list] # Only input genes
  }

  # [Input] Covariates (The user can input a custom set of covariates as a matrix)
  if (is.character(cv)) {
    data(list=sprintf("covariates_%s",cv), package="dndscv")
  } else {
    covs = cv
  }

  # [Input] Known cancer genes (The user can input a gene list as a character vector)
  if (kc[1] %in% c("cgc81")) {
    data(list=sprintf("cancergenes_%s",kc), package="dndscv")
  } else {
    known_cancergenes = kc
  }

  # [Input] Substitution model (The user can also input a custom substitution model as a matrix)
  if (length(sm)==1) {
    data(list=sprintf("submod_%s",sm), package="dndscv")
  } else {
    substmodel = sm
  }

  # Expanding the reference sequences [for faster access]
  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
    RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
    RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
    if (!is.null(RefCDS[[j]]$seq_splice)) {
      RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
      RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
      RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
    }
  }


  ## 2. Mutation annotation
  message("[2] Annotating the mutations...")

  nt = c("A","C","G","T")
  trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
  trinucinds = setNames(1:64, trinucs)
  trinucsubs = NULL
  for (j in 1:length(trinucs)) {
    trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
  }
  trinucsubsind = setNames(1:192, trinucsubs)

  ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
  gr_genes_ind = ind[gr_genes$names]

  # Warning about possible unannotated dinucleotide substitutions
  if (any(diff(mutations$pos)==1)) {
    warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
  }

  # Warning about multiple instances of the same mutation in different sampleIDs
  if (nrow(unique(mutations[,2:5])) < nrow(mutations)) {
    warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
  }

  # Start and end position of each mutation
  mutations$end = mutations$start = mutations$pos
  l = nchar(mutations$ref)-1 # Deletions of multiple bases
  mutations$end = mutations$end + l
  ind = substr(mutations$ref,1,1)==substr(mutations$mut,1,1) & nchar(mutations$ref)>nchar(mutations$mut) # Position correction for deletions annotated in the previous base (e.g. CA>C)
  mutations$start = mutations$start + ind

  # Mapping mutations to genes
  gr_muts = GenomicRanges::GRanges(mutations$chr, IRanges::IRanges(mutations$start,mutations$end))
  ol = as.data.frame(GenomicRanges::findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  mutations = mutations[ol[,1],] # Duplicating subs if they hit more than one gene
  mutations$geneind = gr_genes_ind[ol[,2]]
  mutations$gene = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind]
  mutations = unique(mutations)

  # Optional: Excluding samples exceeding the limit of mutations/sample [see Default parameters]
  nsampl = sort(table(mutations$sampleID))
  exclsamples = NULL
  if (any(nsampl>max_coding_muts_per_sample)) {
    message(sprintf('    Note: %0.0f samples excluded for exceeding the limit of mutations per sample (see the max_coding_muts_per_sample argument in dndscv). %0.0f samples left after filtering.',sum(nsampl>max_coding_muts_per_sample),sum(nsampl<=max_coding_muts_per_sample)))
    exclsamples = names(nsampl[nsampl>max_coding_muts_per_sample])
    mutations = mutations[!(mutations$sampleID %in% names(nsampl[nsampl>max_coding_muts_per_sample])),]
  }

  # Optional: Limiting the number of mutations per gene per sample (to minimise the impact of unannotated kataegis and other mutation clusters) [see Default parameters]
  mutrank = ave(mutations$pos, paste(mutations$sampleID,mutations$gene), FUN = function(x) rank(x))
  exclmuts = NULL
  if (any(mutrank>max_muts_per_gene_per_sample)) {
    message(sprintf('    Note: %0.0f mutations removed for exceeding the limit of mutations per gene per sample (see the max_muts_per_gene_per_sample argument in dndscv)',sum(mutrank>max_muts_per_gene_per_sample)))
    exclmuts = mutations[mutrank>max_muts_per_gene_per_sample,]
    mutations = mutations[mutrank<=max_muts_per_gene_per_sample,]
  }

  # Additional annotation of substitutions

  mutations$strand = sapply(RefCDS,function(x) x$strand)[mutations$geneind]
  snv = (mutations$ref %in% nt & mutations$mut %in% nt)
  if (!any(snv)) { stop("Zero coding substitutions found in this dataset. Unable to run dndscv. Common causes for this error are inputting only indels or using chromosome names different to those in the reference database (e.g. chr1 vs 1)") }
  indels = mutations[!snv,]
  mutations = mutations[snv,]
  mutations$ref_cod = mutations$ref
  mutations$mut_cod = mutations$mut
  compnt = setNames(rev(nt), nt)
  isminus = (mutations$strand==-1)
  mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
  mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]

  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$N = array(0, dim=c(192,4)) # Initialising the N matrices
  }

  # Subfunction: obtaining the codon positions of a coding mutation given the exon intervals

  chr2cds = function(pos,cds_int,strand) {
    if (strand==1) {
      return(which(unlist(apply(cds_int, 1, function(x) x[1]:x[2])) %in% pos))
    } else if (strand==-1) {
      return(which(rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2]))) %in% pos))
    }
  }

  # Annotating the functional impact of each substitution and populating the N matrices

  ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = codonsub = array(NA, nrow(mutations))

  for (j in 1:nrow(mutations)) {

    geneind = mutations$geneind[j]
    pos = mutations$pos[j]

    if (any(pos == RefCDS[[geneind]]$intervals_splice)) { # Essential splice-site substitution

      impact[j] = "Essential_Splice"; impind = 4
      pos_ind = (pos==RefCDS[[geneind]]$intervals_splice)
      cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      aachange[j] = ntchange[j] = codonsub[j] = "."

    } else { # Coding substitution

      pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
      cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
      old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
      pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
      new_codon = old_codon; new_codon[pos_in_codon] = mutations$mut_cod[j]
      old_aa = seqinr::translate(old_codon, numcode = numcode)
      new_aa = seqinr::translate(new_codon, numcode = numcode)
      aachange[j] = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
      ntchange[j] = sprintf('%s%0.0f%s',mutations$ref_cod[j],pos_ind,mutations$mut_cod[j])
      codonsub[j] = sprintf('%s>%s',paste(old_codon,collapse=""),paste(new_codon,collapse=""))

      # Annotating the impact of the mutation
      if (new_aa == old_aa){
        impact[j] = "Synonymous"; impind = 1
      } else if (new_aa == "*"){
        impact[j] = "Nonsense"; impind = 3
      } else if (old_aa != "*"){
        impact[j] = "Missense"; impind = 2
      } else if (old_aa=="*") {
        impact[j] = "Stop_loss"; impind = NA
      }
    }

    if (mutations$ref_cod[j] != as.character(cdsnt)) { # Incorrect base annotation in the input mutation file (the mutation will be excluded with a warning)
      wrong_ref[j] = 1
    } else if (!is.na(impind)) { # Correct base annotation in the input mutation file
      trisub = trinucsubsind[ paste(ref3_cod[j], mut3_cod[j], sep=">") ]
      RefCDS[[geneind]]$N[trisub,impind] = RefCDS[[geneind]]$N[trisub,impind] + 1 # Adding the mutation to the N matrices
    }

    if (round(j/1e4)==(j/1e4)) { message(sprintf('    %0.3g%% ...', round(j/nrow(mutations),2)*100)) }
  }

  mutations$ref3_cod = ref3_cod
  mutations$mut3_cod = mut3_cod
  mutations$aachange = aachange
  mutations$ntchange = ntchange
  mutations$codonsub = codonsub
  mutations$impact = impact
  mutations$pid = sapply(RefCDS,function(x) x$protein_id)[mutations$geneind]

  if (any(!is.na(wrong_ref))) {
    if (mean(!is.na(wrong_ref)) < 0.1) { # If fewer than 10% of mutations have a wrong reference base, we warn the user
      warning(sprintf('%0.0f (%0.2g%%) mutations have a wrong reference base (see the affected mutations in dndsout$wrongmuts). Please identify the causes and rerun dNdScv.', sum(!is.na(wrong_ref)), 100*mean(!is.na(wrong_ref))))
    } else { # If more than 10% of mutations have a wrong reference base, we stop the execution (likely wrong assembly or a serious problem with the data)
      stop(sprintf('%0.0f (%0.2g%%) mutations have a wrong reference base. Please confirm that you are not running data from a different assembly or species.', sum(!is.na(wrong_ref)), 100*mean(!is.na(wrong_ref))))
    }
    wrong_refbase = mutations[!is.na(wrong_ref), 1:5]
    mutations = mutations[is.na(wrong_ref),]
  }

  if (any(nrow(indels))) { # If there are indels we concatenate the tables of subs and indels
    indels = cbind(indels, data.frame(ref_cod=".", mut_cod=".", ref3_cod=".", mut3_cod=".", aachange=".", ntchange=".", codonsub=".", impact="no-SNV", pid=sapply(RefCDS,function(x) x$protein_id)[indels$geneind]))

    # Annotation of indels
    ins = nchar(gsub("-","",indels$ref))<nchar(gsub("-","",indels$mut))
    del = nchar(gsub("-","",indels$ref))>nchar(gsub("-","",indels$mut))
    multisub = nchar(gsub("-","",indels$ref))==nchar(gsub("-","",indels$mut)) # Including dinucleotides
    l = nchar(gsub("-","",indels$ref))-nchar(gsub("-","",indels$mut))
    indelstr = rep(NA,nrow(indels))
    for (j in 1:nrow(indels)) {
      geneind = indels$geneind[j]
      pos = indels$start[j]:indels$end[j]
      if (ins[j]) { pos = c(pos-1,pos) } # Adding the upstream base for insertions
      pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
      if (length(pos_ind)>0) {
        inframe = (length(pos_ind) %% 3) == 0
        if (ins[j]) { # Insertion
          indelstr[j] = sprintf("%0.0f-%0.0f-ins%s",min(pos_ind),max(pos_ind),c("frshift","inframe")[inframe+1])
        } else if (del[j]) { # Deletion
          indelstr[j] = sprintf("%0.0f-%0.0f-del%s",min(pos_ind),max(pos_ind),c("frshift","inframe")[inframe+1])
        } else { # Dinucleotide and multinucleotide changes (MNVs)
          indelstr[j] = sprintf("%0.0f-%0.0f-mnv",min(pos_ind),max(pos_ind))
        }
      }
    }
    indels$ntchange = indelstr
    annot = rbind(mutations, indels)
  } else {
    annot = mutations
  }
  annot = annot[order(annot$sampleID, annot$chr, annot$pos),]


  ## 3. Estimation of the global rate and selection parameters
  message("[3] Estimating global rates...")

  Lall = array(sapply(RefCDS, function(x) x$L), dim=c(192,4,length(RefCDS)))
  Nall = array(sapply(RefCDS, function(x) x$N), dim=c(192,4,length(RefCDS)))
  L = apply(Lall, c(1,2), sum)
  N = apply(Nall, c(1,2), sum)

  # Subfunction: fitting substitution model

  fit_substmodel = function(N, L, substmodel) {

    l = c(L); n = c(N); r = c(substmodel)
    n = n[l!=0]; r = r[l!=0]; l = l[l!=0]

    params = unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
    indmat = as.data.frame(array(0, dim=c(length(r),length(params))))
    colnames(indmat) = params
    for (j in 1:length(r)) {
      indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1
    }

    model = glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
    summary_model<-summary(model)
    mle = exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
    ci = exp(confint.default(model)) # Wald confidence intervals
    par = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2],p_value=summary_model$coefficients[,4])
    return(list(par=par, model=model))
  }

  # Fitting all mutation rates and the 3 global selection parameters

  poissout = fit_substmodel(N, L, substmodel) # Original substitution model
  par = poissout$par
  poissmodel = poissout$model
  parmle =  setNames(par[,2], par[,1])
  mle_submodel = par
  rownames(mle_submodel) = NULL

  # Fitting models with 1 and 2 global selection parameters

  s1 = gsub("wmis","wall",gsub("wnon","wall",gsub("wspl","wall",substmodel)))
  par1 = fit_substmodel(N, L, s1)$par # Substitution model with 1 selection parameter
  s2 = gsub("wnon","wtru",gsub("wspl","wtru",substmodel))
  par2 = fit_substmodel(N, L, s2)$par # Substitution model with 2 selection parameter
  globaldnds = rbind(par, par1, par2)[c("wmis","wnon","wspl","wtru","wall"),]
  sel_loc = sel_cv = NULL

  return(globaldnds)
}

#' compare dndscv results of two related cohort. for example, primary mutations and metastatic mutations to see if the overall fitness of a specific pathway is changed from primary tumors to metastatic tumors.
#'
#' @param pathway_table, a data frame containing two columns: pathway, gene
#' @param mutation_table_1, mutation table of cohort 1
#' @param mutation_table_2, mutation table of cohort 2
#' @param name_mutation_table1,assign a name for cohort 1
#' @param name_mutation_table2, assign a name for cohort 2
#' @param mutation_gene1, genes contained in the muation_table_1, a vector
#' @param mutation_gene2, genes contained in the muation_table_2, a vector
#' @param cutoff, pathway must contain at least "cutoff" genes in both cohort 1 and cohort 2, otherwise, not performing analysis
#'
#' @return, melt table and dcast table
#' @export
#'
#' @examples
dndscv_compare<-function(pathway_table,mutation_table_1=ky744_mutation_primary,mutation_table_2=ky744_mutation_metastasis,name_mutation_table1="primary",name_mutation_table2="distance",mutation_gene1=primary_gene,mutation_gene2=metastasis_gene,cutoff=5){
  total_result_combine<-list()
  pair_result_combine<-list()
  for (pa in unique(pathway_table$pathway)){
    pa_gene<-as.character(pathway_table$gene[pathway_table$pathway==pa])
    intersect_gene_mutation1<-intersect(pa_gene,mutation_gene1)
    intersect_gene_mutation2<-intersect(pa_gene,mutation_gene2)
    if (length(intersect_gene_mutation1)>=cutoff & length(intersect_gene_mutation2)>=cutoff){
      pa_dndscv_1<-dndscv_pathway(mutations=mutation_table_1,sm="192r_3w",gene_list=pa_gene)
      result1<-as.data.frame(pa_dndscv_1)
      result1$overlap_num<-length(intersect_gene_mutation1)
      result1$group<-name_mutation_table1
      result1$pathway<-pa

      pa_dndscv_2<-dndscv_pathway(mutations=mutation_table_2,sm="192r_3w",gene_list=pa_gene)
      result2<-as.data.frame(pa_dndscv_2)
      result2$overlap_num<-length(intersect_gene_mutation2)
      result2$group<-name_mutation_table2
      result2$pathway<-pa

      total_result<-rbind(result1,result2)
      colnames(result1)[2:(ncol(result1)-1)]<-paste(colnames(result1)[2:(ncol(result1)-1)],name_mutation_table1,sep="_")
      colnames(result2)[2:(ncol(result2)-1)]<-paste(colnames(result2)[2:(ncol(result2)-1)],name_mutation_table2,sep="_")
      total_pair<-merge(result1,result2,by=c("name","pathway"))
      total_result_combine[[pa]]<-total_result
      pair_result_combine[[pa]]<-total_pair
    }

  }
  total_result_melt<-do.call(rbind,total_result_combine)
  total_result_pair<-do.call(rbind,pair_result_combine)
  final_result<-list()
  final_result[["melt_result"]]<-total_result_melt
  final_result[["pair_result"]]<-total_result_pair
  return(final_result)
}


# rectom<-GSA.read.gmt("c2.cp.reactome.v7.1.symbols.gmt")
# kegg<-GSA.read.gmt("c2.cp.kegg.v7.1.symbols.gmt")
# cytoband<-GSA.read.gmt("c1.all.v7.1.symbols.gmt")
# key_pathway<-GSA.read.gmt("h.all.v7.1.symbols.gmt")
# biocarta<-GSA.read.gmt("c2.cp.biocarta.v7.1.symbols.gmt")
# pathway_all<-GSA.read.gmt("c2.all.v7.4.symbols.gmt")
# oncogenic_pathway<-GSA.read.gmt("c6.all.v7.4.symbols.gmt")
# immune_signature<-GSA.read.gmt("c7.all.v7.4.symbols.gmt")
# Go_pathway<-GSA.read.gmt("c5.go.v7.4.symbols.gmt")

# Process DDR pathway table
kegg_pathway<-GSA.read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
kegg_pathway<-pathway_data_frame(kegg_pathway)
kegg_pathway$pathway<-gsub("KEGG_","",kegg_pathway$pathway)

# convert gmt_file to a data frame
pathway_data_frame<-function(gmt_file){
  total<-list()
  for (i in 1:length(gmt_file$geneset.names)){
    result<-data.frame(pathway=gmt_file$geneset.names[[i]],gene=gmt_file$genesets[[i]])
    total[[i]]<-result
  }
  final<-do.call(rbind,total)
  return(final)
}

#supplement Faconi anemia pathway
fanconi<-read.csv("fanconi.csv",stringsAsFactors = F)
fanconi<-fanconi[,c(2,1)]
kegg_pathway<-rbind(kegg_pathway,fanconi)
kegg_pathway_DNA_repair<-kegg_pathway[kegg_pathway$pathway %in% c("BASE_EXCISION_REPAIR","NUCLEOTIDE_EXCISION_REPAIR","MISMATCH_REPAIR","HOMOLOGOUS_RECOMBINATION","NON_HOMOLOGOUS_END_JOINING","P53_SIGNALING_PATHWAY","Fanconi_anemia_pathway"),]

# DDR contains 7 pathways
unique(kegg_pathway_DNA_repair$pathway)
DDR_pathway<-kegg_pathway_DNA_repair
DDR_pathway$pathway<-"DDR"
kegg_pathway_DNA_repair_total<-rbind(kegg_pathway_DNA_repair,DDR_pathway)
kegg_pathway_DNA_repair_cancer_gene<-kegg_pathway_DNA_repair_total[kegg_pathway_DNA_repair_total$gene %in% cancer_gene,]

# only including cancer genes from OncoKB and ref annotated genes
kegg_pathway_DNA_repair_cancer_gene_refgene<-kegg_pathway_DNA_repair_cancer_gene[kegg_pathway_DNA_repair_cancer_gene$gene %in% refcds_gene,]

# to maintain statistically comparable between primary and metastases, only paired samples are included
mutation_primary_paired<-mutation_primary[mutation_primary$Sample_ID %in%  pair_string_primary_distance_P49_rm$primary_sample_id,]
mutation_metastasis_paired<-mutation_metastasis[mutation_metastasis$Sample_ID %in%  pair_string_primary_distance_P49_rm$distance_sample_id,]

# study the fitness of DDR pathway in primary tumors and distant metastatic tumors

DDR_pathway_dndscv_pair<-dndscv_compare(pathway_table=kegg_pathway_DNA_repair_cancer_gene_refgene,mutation_table_1=mutation_primary_paired,mutation_table_2=mutation_metastasis_paired)
DDR_pathway_dndscv_pair$pair_result[DDR_pathway_dndscv_pair$pair_result$p_value_primary<=0.05 | DDR_pathway_dndscv_pair$pair_result$p_value_distance<=0.05,]
write.csv(DDR_pathway_dndscv_pair$pair_result,"DDR_pathway_dndscv_pair.csv")



