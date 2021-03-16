#!/usr/bin/R

############ READING FUNCTIONS #############

read_deseq2 <- function(res_fname, pthresh=0.05, downreg='AFR', upreg='EUR') {
  # downreg = "condition" in which DESeq2 log2FoldChange < 0
  res <- read.csv(res_fname, header=TRUE) %>% as_tibble()
  res <- res %>%
              dplyr::select(ensembl_gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,neglog10padj) %>%
              mutate(DE_significant=if_else(padj<pthresh,TRUE,FALSE),
                     DE_direction=if_else(log2FoldChange<0,downreg,upreg)) %>%
              dplyr::rename(DE_baseMean = baseMean,
                     DE_log2FoldChange = log2FoldChange,
                     DE_lfcSE = lfcSE,
                    DE_stat = stat,
                    DE_pvalue = pvalue,
                    DE_padj = padj,
                    DE_neglog10padj = neglog10padj)
  return(res)
}

read_deseq2_raw <- function(res_fname,
                            pthresh=0.05) {

    if(base::endsWith(res_fname,'.csv')) {
      sep=','
    } else {
      sep='\t'
    }

    # possible column names (doesn't strictly select)
    cnames <- c(
      'tissue',
      'dropped_cov',
      'ensembl_gene_id',
      'baseMean',
      'log2FoldChange',
      'lfcSE',
      'stat',
      'pvalue',
      'padj',
      'neglog10padj'
    )

    res <- read_delim(res_fname, sep) %>%
        dplyr::select(any_of(cnames)) %>%
        mutate(
            DE_significant=if_else(padj<pthresh,TRUE,FALSE)
        )


    return(res)

}

read_hg19_eqtl_granges_new <- function(eqtl_fname, ch_fname, gs_fname, fst=NULL, bf=NULL, de=NULL, lift_over=FALSE) {
  # add eQTLs from GTEx tissue tehaas format
  eqtl <- read_tsv(eqtl_fname, quote="") %>%
            separate(Gene_ID, c('ensembl_gene_id',NA), sep="\\.", remove=TRUE) %>%
            separate(SNP_ID, c('seqnames','start','ref','alt','build'), sep="_", remove=FALSE, convert=TRUE) %>%
            dplyr::rename(variant_id=SNP_ID) %>%
            dplyr::rename(eqtl_pval=`P-Value`) %>%
            as_granges(width=1)
  if (grepl("v8", eqtl_fname) | lift_over) {
    genome(eqtl) = "GRCh38"
    library(rtracklayer)
    ch = import.chain(ch_fname)
    seqlevelsStyle(eqtl) = "UCSC"  # necessary
    print("Lifting over eQTL coordinates from hg38 to hg19")
    eqtl_19 = liftOver(eqtl, ch)
    eqtl_19 = unlist(eqtl_19)
    genome(eqtl_19) = "hg19"
    print(paste(length(eqtl) - length(eqtl_19),"eQTLs dropped in liftOver"))
  } else {
    eqtl_19 <- eqtl
  }
  gs <- readLines(gs_fname)
  eqtl_19 <- eqtl_19 %>%
                  as_tibble() %>%
                  dplyr::select(c(seqnames,
                                 start,
                                 end,
                                 variant_id,
                                 ensembl_gene_id,
                                 eqtl_pval)) %>%
                  dplyr::rename(eqtl_id_new=variant_id,
                        eGene=ensembl_gene_id) %>%
                  mutate(eGene_in_set=if_else(eGene %in% gs,
                                       TRUE,
                                       FALSE),
                         eqtl_id=paste(seqnames,start,sep='_')) %>%
                  as_granges()

  if (!is.null(fst)) {
    if (class(fst)!="GRanges") {
      fst <- as_granges(fst)
    }
    print("Joining FST with eQTL...")
    eqtl_19 <- join_overlap_left(eqtl_19, fst) %>%
      dplyr::mutate(FST_eqtl=FST) %>%
      dplyr::select(-FST)
  }

  if (!is.null(bf)) {
    if (class(bf)!="GRanges") {
      bf <- as_granges(bf)
    }
    print("Joining BF with eQTL...")
    eqtl_19 <- join_overlap_left(eqtl_19, bf) %>%
      dplyr::mutate(BF_eqtl=median_BF) %>%
      dplyr::select(-median_BF)
  }

  if (!is.null(de)) {
      eqtl_19 <- eqtl_19 %>%
          as_tibble %>%
          merge(.,
                de %>%
                    dplyr::select(ensembl_gene_id, DE_log2FoldChange, DE_padj, DE_significant, DE_direction) %>%
                    dplyr::rename(eGene=ensembl_gene_id),
                by='eGene',
                all.x=TRUE)

  }

  return(eqtl_19)
}

read_caqtl_granges <- function(caqtl_fname, pop=1, fst=NULL, bf=NULL) {
  if (grepl('shared', caqtl_fname)) {
    qtls <- read_tsv(caqtl_fname,quote="")
    qtls <- qtls %>%
              mutate(seqnames=paste0("chr",chr)) %>%
              dplyr::select(-c(start,chr, peak_len)) %>%
              dplyr::rename(start=pos,
                            open_allele=open_best,
                            closed_allele=closed_best,
                            caqtl_pval=fishers_pval) %>%
              as_granges(width=1)
  } else {
    library(readxl)
    qtls <- read_excel(caqtl_fname, skip=0,sheet=pop,col_names=T) %>%
      dplyr::rename(seqnames=Chr,
                 start=position,
                 alt=ALTallele,
                 ref=REFallele,
                 open_allele=higherBindingAllele,
                 caqtl_pval=pvalue) %>%
      mutate(seqnames=paste0("chr",seqnames),
             closed_allele=if_else(open_allele==ref, alt, ref)) %>%
      as_granges(width=1)
  }

  qtls <- qtls %>%
    mutate(caqtl_id=paste(seqnames, start, sep="_"))

  if (!is.null(fst)) {
    if (class(fst)!="GRanges") {
      fst <- as_granges(fst)
    }
    print("Joining FST with caQTL...")
    qtls <- join_overlap_left(qtls, fst) %>%
      dplyr::mutate(FST_caqtl=FST) %>%
      dplyr::select(-FST)
  }

  if (!is.null(bf)) {
    if (class(bf)!="GRanges") {
      bf <- as_granges(bf)
    }
    print("Joining BF with caQTL...")
    qtls <- join_overlap_left(qtls, bf) %>%
      dplyr::mutate(BF_caqtl=median_BF) %>%
      dplyr::select(-median_BF)
  }

  return(qtls)

}

read_bed_noheader_granges <- function(bed) {
  cnames = c('seqnames','start','end')
  bedr <- read_tsv(bed, col_names=cnames) %>%
      as_granges()
  return(bedr)
}


################# REFORMATTING / SUBSETTING FUNCTIONS ################

liftover_grange <- function(df, ch_fname, from_b="GRCh38", to_b="hg19") {
    genome(df) = from_b
    library(rtracklayer)
    ch = import.chain(ch_fname)
    seqlevelsStyle(df) = "UCSC"  # necessary
    print(paste0("Lifting over eQTL coordinates from ", from_b, " to ", to_b, "..."))
    df_new = liftOver(df, ch)
    df_new = unlist(df_new)
    genome(df_new) = to_b
    nc <- length(df_new) - length(df)
    print(paste0("Ranges gained/lost in liftOver = ", if_else(nc>0,'+',''), nc))
    return(df_new)
}

hgnc2ensembl <- function (df) {

    hgncgenes <- df %>% pull(hgnc_symbol) %>% unique()
    hgncgenes <- hgncgenes[which(!is.na(hgncgenes))]
    merge_col <- "hgnc_symbol"

    library(EnsDb.Hsapiens.v79)
    gene_info <- ensembldb::select(EnsDb.Hsapiens.v79,
                                   keys = hgncgenes,
                                   keytype = "SYMBOL",
                                   columns = c("SYMBOL", "GENEID")) %>%
        as_tibble() %>%
        dplyr::rename(ensembl_gene_id = GENEID,
                      hgnc_symbol = SYMBOL) %>%
        dplyr::filter(grepl('ENSG', ensembl_gene_id)) # filter out alternate gene symbols

    print("Adding hgnc_symbols...")
    df <- merge(df, gene_info, by = merge_col, all.x = TRUE)
    return(df)

}

add_hgnc <- function (df, eGene = FALSE, biomart = FALSE, trans_eGene=FALSE) {
    if (eGene) {
        ensgenes <- df %>% pull(eGene) %>% unique()
        ensgenes <- ensgenes[which(!is.na(ensgenes))]
        merge_col <- "eGene"
    } else if (trans_eGene) {
        ensgenes <- df %>% pull(trans_eGene) %>% unique()
        ensgenes <- ensgenes[which(!is.na(ensgenes))]
        merge_col <- "trans_eGene"
    }
    else {
        ensgenes <- df %>% pull(ensembl_gene_id) %>% unique()
        ensgenes <- ensgenes[which(!is.na(ensgenes))]
        merge_col <- "ensembl_gene_id"
    }
    if (biomart) {
        library(biomaRt)
        ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
        mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
        gene_attributes = c("ensembl_gene_id", "hgnc_symbol")
        print("Getting hgnc_symbols from biomaRt...")
        gene_info = getBM(attributes = gene_attributes, filters = "ensembl_gene_id",
            values = ensgenes, mart = ensembl)
    }
    else {
        library(EnsDb.Hsapiens.v79)
        gene_info <- ensembldb::select(EnsDb.Hsapiens.v79, keys = ensgenes,
            keytype = "GENEID", columns = c("SYMBOL", "GENEID")) %>%
            as_tibble() %>% dplyr::rename(ensembl_gene_id = GENEID,
            hgnc_symbol = SYMBOL)
    }
    if (eGene) {
        gene_info <- gene_info %>% as_tibble() %>% dplyr::rename(eGene = ensembl_gene_id,
            eGene_hgnc = hgnc_symbol)
    }
    if (trans_eGene) {
        gene_info <- gene_info %>% as_tibble() %>% dplyr::rename(trans_eGene = ensembl_gene_id,
            trans_eGene_hgnc = hgnc_symbol)
    }

    print("Adding hgnc_symbols...")
    df <- merge(df, gene_info, by = merge_col, all.x = TRUE)
    return(df)
}

gprofiler_deseq2 <- function(df,
                               pthresh=0.05,
                               lfc_thresh=1,
                               ordered=FALSE,
                               forward=TRUE,
                               two_tailed=TRUE,
                               upreg=TRUE,
                               bg_all=TRUE,
                               lfc_bg=TRUE,
                               debug=FALSE,
                               outroot='DE') {

    # two_tailed = TRUE tests for GSE in up- and down-regulated genes relative to condition (bg_all=TRUE, upreg=FALSE)
    # upreg = TRUE, tests for GSE in up-regulated genes relative to condition
    # upreg = FALSE, tests for GSE in down-regulated genes relative to condition
    # bg_all = TRUE, all tested genes in background
    # bg_all = FALSE, tests against background of up- and down-regulated genes relative to condition
    # lfc_bg = TRUE, requires background genes to pass LFC threshold
    # lfc_bg = FALSE, when bg_all=FALSE background genes only have to pass pval threshold

    if ('tissue' %in% colnames(df)) {
      tissue_name <- df %>%
        pull(tissue) %>%
        unique()
    } else {
      tissue_name <- 'default'
    }

    library(gprofiler2)

    if (ordered) {

      if (forward) {

        outbase = paste(outroot,
                        'p',pthresh,
                        'lfc',lfc_thresh,
                        'orderedForward',
                        'bgALL',
                        sep="_")

        fg_genes <- df %>%
            dplyr::filter(!is.na(padj)) %>%
            arrange(desc(stat)) %>%
            pull(ensembl_gene_id) %>%
            unique()

      } else {

        outbase = paste(outroot,
                        'p',pthresh,
                        'lfc',lfc_thresh,
                        'orderedReverse',
                        'bgALL',
                        sep="_")

        fg_genes <- df %>%
            dplyr::filter(!is.na(padj)) %>%
            arrange(stat) %>%
            pull(ensembl_gene_id) %>%
            unique()

      }

      bg_genes <- df %>%
          dplyr::filter(!is.na(padj)) %>%
          pull(ensembl_gene_id) %>%
          unique()

    } else {

      dfs <- df %>%
          dplyr::filter(padj < pthresh,
                        abs(log2FoldChange) > lfc_thresh) %>%
          arrange(padj)

      if (two_tailed) {

          outbase = paste(outroot,
                          'p',pthresh,
                          'lfc',lfc_thresh,
                          'two-tailed',
                          'bgALL',
                          sep="_")

          fg_genes <- dfs %>%
              pull(ensembl_gene_id) %>%
              unique()
          bg_genes <- df %>%
              dplyr::filter(!is.na(padj)) %>%
              pull(ensembl_gene_id) %>%
              unique()

      } else {

          outbase = paste(outroot,
                          'p',pthresh,
                          'lfc',lfc_thresh,
                          if_else(upreg, 'upreg','downreg'),
                          if_else(bg_all, 'bgALL', 'bgDE'),
                          sep="_")

          if (upreg) {
              fg_genes <- dfs %>%
                  dplyr::filter(log2FoldChange > lfc_thresh) %>%
                  pull(ensembl_gene_id) %>%
                  unique()
          } else {
              fg_genes <- dfs %>%
                  dplyr::filter(log2FoldChange < -lfc_thresh) %>%
                  pull(ensembl_gene_id) %>%
                  unique()
          }
          if (bg_all) {
              bg_genes <- df %>%
                  dplyr::filter(!is.na(padj)) %>%
                  pull(ensembl_gene_id) %>%
                  unique()
          } else {
              if (lfc_bg) {
                  bg_genes <- dfs %>%
                      pull(ensembl_gene_id) %>%
                      unique()
              } else {
                  # only use p-value threshold for background, not LFC
                  bg_genes <- df %>%
                      dplyr::filter(padj < pthresh) %>%
                      pull(ensembl_gene_id) %>%
                      unique()
              }

          }
      }

    }



    print(paste0('Foreground genes at P<', pthresh,':'))
    print(length(fg_genes))
    print('Background genes:')
    print(length(bg_genes))

    custom_bg <- bg_genes

    # if (length(fg_genes) > 11000) {
    #
    #     print("Too many genes for g:Profiler to test quickly")
    #     write('To many genes for efficient testing',      paste0(outbase,'_too_many_genes_gprofiler.txt'))
    #     write('To many genes for efficient testing',      paste0(outbase,'_too_many_genes_gprofiler_summary.txt'))
    #     return(NULL)
    #
    # }

    print("Running g:Profiler query...")

    gostres <- tryCatch(
      {
        gostres <- gost(fg_genes, organism = "hsapiens", ordered_query = ordered,
                      multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                      measure_underrepresentation = FALSE, evcodes = TRUE,
                      user_threshold = 0.05, correction_method ="g_SCS",
                      domain_scope = "known", custom_bg = custom_bg,
                      numeric_ns = "", sources = NULL)
      },
      error=function(e) {
        message(e)
        return(NULL)
      }
    )

    if (is.null(gostres)) {
      print(paste(outbase, "bad request"))
      return(NULL)
    }


    if (debug) {
        return(gostres)
    }

    if (is.null(gostres$result)) {
      print("No results to show")
      write('No results to show',      paste0(outbase,'_no_sig_res_gprofiler.txt'))
      write('No results to show',      paste0(outbase,'_no_sig_res_gprofiler_summary.txt'))
    } else {
      gostdf <- gostres$result %>%
          as_tibble() %>%
          dplyr::select(-c(parents)) %>% # remove list column before writing
          mutate(tissue=tissue_name) %>%
          arrange(source,p_value)

      gostdf_summ <- gostdf %>%
          dplyr::select(tissue,
                        source,
                        term_name,
                        p_value,
                        query_size,
                      term_size,
                    intersection_size)


      print(gostdf_summ %>% dplyr::filter(source=='GO:BP'))

      write_tsv(gostdf,      paste0(outbase,'_gprofiler.txt'))
      write_tsv(gostdf_summ, paste0(outbase,'_gprofiler_summary.txt'))

      return(gostdf)
    }

}

################## PLOTTING FUNCTIONS #################

cor_label <- function(df, col1=2, col2=3, method='pearson') {
    require(grid)
    require(gridExtra)
    cres <- cor.test(df[[col1]], df[[col2]], method=method)
    groblab <- grobTree(
        textGrob(
            paste0(
                paste0(method," corr. = "),
                round(
                    cres$estimate,
                    4
                ),
                ",\nP = ",
                formatC(cres$p.value, format = "e", digits = 4)
            ),
            x = 0.3,
            y = 0.1,
            hjust = 0,
            gp = gpar(col = "gray",
                      fontsize = 20,
                      fontface = "italic")
        )
    )

    return(groblab)
}

################## trans_utils ####################

#!/usr/bin/R

bname <- function(f,pathstrip=TRUE) {
  b <- unlist(strsplit(f,'\\.'))
  b <- b[1] # grab element with all extensions removed
  if (length(unlist(strsplit(b,'\\.')))==1) {
    if (pathstrip) {
      b <- unlist(strsplit(b,'/'))
      return(b[length(b)])
    } else {
      return(b)
    }
  } else {
    return(bname(b,pathstrip=pathstrip))
  }
}

insert_suffix <- function(fname,suff='_subset') {
    fbase <- tools::file_path_sans_ext(fname)
    fext <- tools::file_ext(fname)
    fnew <- paste0(fbase,suff,'.',fext)
    return(fnew)
}

read_matrixEQTL_top_granges <- function(rfname) {
    if (str_detect(rfname,'fdr')) {
        df <- read_tsv(rfname) %>% as_granges
    } else {
       df <- read_tsv(rfname) %>%
        dplyr::rename(variant_id=SNP,
                     statistic=`t-stat`,
                     p_value=`p-value`) %>%
        separate(variant_id, c('seqnames','start','ref','alt','build'), sep="_", remove=FALSE, convert=TRUE) %>%
        as_granges(width=1)
    }

    return(df)
}

readAllGtexExpression <- function(samples = NULL, expressionFile = NULL) {

    sep <-  "\t"
    filecon <- gzfile(expressionFile, "r")
    commentLines <- readLines(filecon, n = 2)

    header <- unlist(strsplit(readLines(filecon, n = 1), sep))

    if(is.null(samples)) {
        sampleIndex <- rep(T, length(header))
        sampleIndex[2] <- F
    } else {
        sampleIndex <- header %in% samples
    }

    sampleIndex[1] <- T

    colClasses <- rep("NULL", length(header))
    colClasses[sampleIndex] <- "numeric"
    colClasses[1] <- "character"

    result <- read.table(filecon, header = F, sep = sep, colClasses = colClasses)
    colnames(result) <- header[sampleIndex]

    close(filecon)

    return(result)

}

get_mean_gs_expression_matrix <- function(e, key, samps=starts_with('GTEX')) {

    # mean expression for each named  set of genes provided
    # in the key (named list)

    i <- 0
    for (gsname in names(key)) {
        i <- i+1
        gs <- key[[gsname]]

        curr_gene <- e %>%
                          dplyr::filter(ensembl_gene_id %in% gs) %>%
                          dplyr::select(samps) %>%
                          colMeans() %>%
                          t() %>%
                          as_tibble() %>%
                          add_column(.,gene=gsname,.before=1)

        if (i==1) {
            gene <- curr_gene
        } else {
            gene <- rbind(gene,curr_gene)
        }

    }

    gene <- rename_all(gene, funs(gsub("^(GTEX-[^-]*)-.*-.*-.*$", "\\1", .)))

    rnames <- gene %>% pull(gene)
    gene <- as.matrix(gene[,-1])
    rownames(gene) <- rnames

    return(gene)

}

get_expression_matrix <- function(e, samps=dplyr::starts_with('GTEX'), genes=NULL) {

    # all expressed genes in matrix format for matrixEQTL
    esamps <- colnames(e)[-1]
    usamps <- dplyr::intersect(samps,esamps)
    e <- e %>%
        dplyr::select(gene,all_of(usamps)) %>%
        # remove any genes with NA values for the selected samples
        # e.g., if supplied expression df contained samples from multiple tissues
        drop_na(all_of(usamps))

    rnames <- e %>% pull(gene)
    gene <- e %>%
        dplyr::select(-gene) %>%
        dplyr::rename_all(., ~ gsub("^(GTEX-[^-]*)-.*-.*-.*$", "\\1", .))

    gene <- as.matrix(gene)
    rownames(gene) <- rnames

    if (!is.null(genes)) {
    	# subset to genes in expression set
    	genes <- intersect(genes, rnames)
      if (length(genes)==1) {
          gene <- t(
            as.matrix(
              gene[genes,]
            )
          )
          rownames(gene) <- genes
      } else {
        gene <- gene[genes,]
      }

    }

    return(gene)

}

transpose_tibble <- function(df, colnames_to="") {
    # transposes a tibble using the first column as the output column names
    # assumes numeric values
    # colnames_to = < column name for first column in the output df
    # storing the input column names >
    rnames <- colnames(df)[-1]
    cnames <- df %>% pull(1)
    df <- as_tibble(t(as.matrix(df[-1])))
    colnames(df) <- cnames
    df[colnames_to] <- rnames
    df <- df %>% dplyr::select(ncol(.), everything())
    coltypes <- paste0('c',paste(rep('d',ncol(df)-1), collapse=''))
    return(type_convert(df, col_types=coltypes))
}

readSelectGtex <- function(samples=NULL, expressionFile=NULL) {

  sep <-  "\t"
  filecon <- gzfile(expressionFile, "r")

  header <- unlist(strsplit(readLines(filecon, n = 1), sep))

  if(is.null(samples)) {
      sampleIndex <- rep(T, length(header))
      sampleIndex[2] <- F
  } else {
      sampleIndex <- header %in% samples
  }

  sampleIndex[1] <- T

  colClasses <- rep("NULL", length(header))
  colClasses[sampleIndex] <- "numeric"
  colClasses[1] <- "character"

  result <- read.table(filecon, header = F, sep = sep, colClasses = colClasses)
  colnames(result) <- header[sampleIndex]

  close(filecon)

  return(result)

}

read_qtl_pheno <- function(qp_fname, samp=NULL) {
  if (str_detect(qp_fname, '(DonovanEtAl2020|cibersort|Cibersort)')) {
    # parse qp_fname to feed into get_expression_matrix function
    # remove low cell type fractions and character pval column
    qp <- read_csv(qp_fname) %>%
      dplyr::select(-c(Fibroblast_3,Fibroblast_4,Lymphocyte_9,P.value))
    qp <- transpose_tibble(qp, colnames_to='gene')

  } else if (str_detect(qp_fname, '(RNASeQ|RNA-seq|RNAseq|tpm)')) {
    qp <- readAllGtexExpression(samples=samp, expressionFile=qp_fname) %>%
            as_tibble() %>%
            mutate(ensembl_gene_id=gsub("\\..*","",Name)) %>%
            dplyr::rename(gene=ensembl_gene_id) %>%
            dplyr::select(-c(Name)) %>%
            dplyr::select(gene, everything())
  } else if (str_detect(qp_fname, '(AFR|EUR|gene_set)')){
    qp <- read_tsv(qp_fname)
  } else if (str_detect(qp_fname, 'xCell')) {
    qp <- readSelectGtex(samples=samp, expressionFile=qp_fname) %>%
      as_tibble() %>%
      dplyr::rename(gene=cell_type)
  } else {
    stop(paste0("QTL phenotype filename not recognized: ", qp_fname))
  }

  # ensure subsetting to supplied samp ids
  if (!is.null(samp)) qp <- qp %>% dplyr::select(gene, any_of(samp))

  return(qp)

}

expression_modules <- function(gene, power=5, min_size=4, max_block=5000) {
  # uses WGCNA to get co-expressed gene modules from the the provided
  # gene-by-individual expression matrix represented as
  # eigengene-by-indivual loadings accessible as
  # net$MEs

  # also plots cluster dendrogram with module colors

  # code adapted from WGCNA tuturials:
  # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

  # library(WGCNA)
  options(stringsAsFactors = FALSE)

  # transpose to individual-by-gene matrix
  if (grepl('GTEX',colnames(gene)[1])) {
    gene <- t(gene)
  }

  print('Gene expression dimensions')
  print(dim(gene))

  if (ncol(gene) > 5000) {
    print("Enabling up to 32 threads")
    enableWGCNAThreads(nThreads=32) # hardcode 32 threads for now
  }

  gsg <- goodSamplesGenes(gene, verbose=3)

  # remove genes/samples with too many missing values or zero variance
  if (!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
    # printFlush(paste("Removing genes:", paste(colnames(gene)[!gsg$goodGenes], collapse = ", ")));
    printFlush(paste("Removing", length(colnames(gene)[!gsg$goodGenes]),"genes"));
    if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(gene)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    gene = gene[gsg$goodSamples, gsg$goodGenes]
  }
  print('power')
  print(power)
  print('min size')
  print(min_size)
  net = blockwiseModules(gene, power = power,
                        TOMType = "unsigned", minModuleSize = min_size,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        # maxBlockSize=max_block,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE,
                        saveTOMFileBase = "GTEx_GS_TOM",
                        verbose = 3)

  if (ncol(gene) < max_block) {
    # Convert labels to colors for plotting
    mergedColors = labels2colors(net$colors)
    # Plot the dendrogram and the module colors underneath
    plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  }

  return(net)

}

filter_expression <- function(gene) {
  gsg <- goodSamplesGenes(gene, verbose=3)

  # remove genes/samples with too many missing values or zero variance
  if (!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(colnames(gene)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(gene)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    gene = gene[gsg$goodSamples, gsg$goodGenes]
  }

  return(gene)
}

plot_sampletree <- function(gene) {
  # plot dendrogram of sample expression to detect outliers
  # code from:
  # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf

  sampleTree = hclust(dist(gene), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  # sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2)

}

plot_network_topology <- function(datExpr) {
  # makes plots for choosing the soft thresholding power beta
  # adjacency will be calculated by raising the co-expression
  # similarity to this power

  # code from:
  # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf

  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  # Plot the results:
  # sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
  xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
  main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

}

get_module_genes <- function(net, hgnc=TRUE, gene_set_fname=NULL) {
    # get dataframe of gene coexpression modules
    # with genes IDs corresponding to each module
    # for overlapping with corresponding tissue eGenes
    gene_mods <- as.data.frame(net$colors)
    colnames(gene_mods) <- 'module'
    gene_mods$ensembl_gene_id <- rownames(gene_mods)
    if (!is.null(gene_set_fname)) {
        # get modules that have at least one gene in gene set
        gs <- readLines(gs_fname)
        gs_mods <- gene_mods[gs,] %>%
            pull(module)
        # remove gene set genes removed prior to module creation
        # with no module (module number = 0)
        gs_mods <- gs_mods[!is.na(gs_mods) & gs_mods!=0] %>%
            unique()
        # filter gene mod df on modules containing gene set genes
        gene_mods <- gene_mods %>%
            filter(module %in% gs_mods)
    }

    if (hgnc) {
      gene_mods <- add_hgnc(gene_mods)
      summ <- gene_mods %>%
          as_tibble() %>%
          group_by(module) %>%
          summarize(n_genes=length(ensembl_gene_id),
                    ensembl_gene_ids=paste(ensembl_gene_id, collapse=';'),
                   hgnc_symbols=paste(hgnc_symbol, collapse=';')) %>%
          mutate(module=paste0('ME',module))
    } else {
      summ <- gene_mods %>%
          as_tibble() %>%
          group_by(module) %>%
          summarize(n_genes=length(ensembl_gene_id),
                    ensembl_gene_ids=paste(ensembl_gene_id, collapse=';')) %>%
          mutate(module=paste0('ME',module))
    }

    return(summ)
}

get_gs_mod_genes <- function(net, gene_set_fname=NULL) {
    # get dataframe of gene coexpression modules
    # with genes IDs corresponding to each module
    # for overlapping with corresponding tissue eGenes
    gene_mods <- as.data.frame(net$colors)
    colnames(gene_mods) <- 'module'
    gene_mods$ensembl_gene_id <- rownames(gene_mods)
    if (!is.null(gene_set_fname)) {
        # get modules that have at least one gene in gene set
        gs <- readLines(gs_fname)
        gs_mods <- gene_mods[gs,] %>%
            pull(module)
        # remove gene set genes removed prior to module creation
        # with no module (module number = 0)
        gs_mods <- gs_mods[!is.na(gs_mods) & gs_mods!=0] %>%
            unique()
        # filter gene mod df on modules containing gene set genes
        gene_mods <- gene_mods %>%
            filter(module %in% gs_mods)
    }

    return(gene_mods %>% as_tibble())

}

summarize_gsmod <- function(gene_set_mod_df) {
    # outputs df with eQTL and gene set numbers for
    # genes in each WGCNA module
    gene_set_mod_df %>%
        group_by(module, ensembl_gene_id) %>%
        arrange(eqtl_pval,.by_group=TRUE) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        group_by(module) %>%
        summarize(n_genes=length(ensembl_gene_id),
                  n_egenes=length(ensembl_gene_id[eqtl]),
                 n_genes_in_set=sum(in_gene_set),
                  n_egenes_in_set=sum(eqtl & in_gene_set),
                 prop_genes_in_set=n_genes_in_set/n_genes) %>%
        arrange(desc(prop_genes_in_set))

}

module_eqtl_egene_pairs <- function(gene_set_mod_df, mod_num=NULL, in_gene_set=FALSE, top_hit=FALSE) {
    rdf <- gene_set_mod_df %>%
                filter(eqtl)

    if (!is.null(mod_num)) {
      rdf <- rdf %>% filter(module==mod_num)
    }

    if(in_gene_set) {
        rdf <- rdf %>% filter(in_gene_set)
    }

    if (top_hit) {
        rdf <- rdf %>%
            group_by(ensembl_gene_id) %>%
            arrange(eqtl_pval, .by_group=TRUE) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup()
    }

    return(rdf)
}

get_ld_snps <- function(df) {

    if (!('rsID' %in% colnames(df))) {
        print("Adding rsIDs...")
        df <- add_rsIDs(df, format=TRUE)
    }

    rs_df <- df %>%
      dplyr::select(4,seqnames,start,rsID,variant_id) %>%
      distinct()
    rs_df$ld_snps <- apply(rs_df, MARGIN=1, FUN=ld_expand_df,
                            plink_bed_dir=plink_bed_dir,
                            plink_pre=plink_pre,
                            plink_dir=plink_dir,
                            rscol='variant_id',
                            sep=".")

    df <- merge(df, rs_df, by=c('module','seqnames','start','rsID','variant_id'))

    df <- unnest(df,ld_snps) %>%
      dplyr::select(-c(seqnames,start,end)) %>%
      separate(., ld_snps, c("seqnames","start"), sep = "_", convert=TRUE, remove=FALSE) %>%
      mutate(end=start) %>%
      dplyr::rename(tag_variant_id=variant_id)

    return(df)
}


subset_matrix <- function(m, rnames=NULL, cnames=NULL) {

  print("Number of rows to subset to...")
  print(length(rnames))
  print(rnames)
  print("Number of columns to subset to...")
  print(length(cnames))
  print(head(cnames))

  print("Dimensions of matrix...")
  print(dim(m))
  print(head(colnames(m)))
  print(rownames(m))

  if (is.null(rnames)) {
    if (is.null(rownames(m))) {
      rnames=c(1:nrow(m))
    } else {
      rnames=rownames(m)
    }
  }

  if (is.null(cnames)) {
    if (is.null(colnames(m))) {
      cnames=c(1:ncol(m))
    } else {
      cnames=colnames(m)
    }
  }

  # ensure rownames match
  rnames <- dplyr::intersect(rnames,rownames(m))

  if (length(rnames)==1 | length(cnames)==1) {

    msub <- t(
        as.matrix(
            m[rnames,cnames]
        )
    )
    rownames(msub) <- rnames

  } else {
      msub <- m[rnames,cnames]
  }

  return(msub)

}

lmMatFunction <- function(s,
                          g,
                          useModel = modelLINEAR,
                          pvalCutoff = 1,
                          outFile = tempfile(),
                          min.pv.by.genesnp = F,
                          noFDRsaveMemory = T,
                          cvrt = NULL,
                         verbose = T,
                         reformat = F){
  	library(MatrixEQTL)

  	sMat <- SlicedData$new()
  	gMat <- SlicedData$new()
  	cvrtMat <- SlicedData$new()

  	sMat$CreateFromMatrix(s)
  	gMat$CreateFromMatrix(g)

  	if(!is.null(cvrt)) {
        cvrtMat$CreateFromMatrix(cvrt)
        results <- Matrix_eQTL_engine(snps = sMat, gene = gMat, cvrt = cvrtMat, output_file_name = outFile,
  				      pvOutputThreshold = pvalCutoff, useModel = useModel, min.pv.by.genesnp = min.pv.by.genesnp,
                                  noFDRsaveMemory = noFDRsaveMemory, verbose = verbose)
    } else {
        results <- Matrix_eQTL_engine(snps = sMat, gene = gMat, output_file_name = outFile,
  				      pvOutputThreshold = pvalCutoff, useModel = useModel, min.pv.by.genesnp = min.pv.by.genesnp,
                                  noFDRsaveMemory = noFDRsaveMemory, verbose = verbose)
    }

    if (reformat) {

      rep <- results$all$eqtls %>%
        as_tibble() %>%
        dplyr::rename(variant_id=snps,
                      ensembl_gene_id=gene) %>%
     	  separate(., variant_id, c('seqnames','start',NA,NA,NA),
                 sep="_", convert=TRUE, remove=FALSE) %>%
        mutate(end=start,
               p_adj=pvalue*length(unique(variant_id))*length(unique(ensembl_gene_id))) %>%
        mutate(p_adj=if_else(p_adj>1, 1, p_adj)) %>%
        dplyr::select(variant_id,
                      seqnames,
  	                  start,
  	                  end,
                      ensembl_gene_id,
                      pvalue,
                      p_adj,
                      FDR,
                      beta,
                      statistic) %>%
        arrange(variant_id, FDR)

      return(rep)

    }


  	return(results)
}

# SMTSISCH = Total Ischemic time for a sample in sample attributes file
# AGE = age of donor in phenotypes file

get_SAMPID_DONORID <- function(pheno_fname,
                                tissues_fname,
                                attr_fname,
                                ancestry = c(2,3),
                                tissue_num = NULL,
                                admix_fname = NULL,
                                p_vars = c('SEX','ancestry','AGE'),
                                a_vars = c('SMTSISCH','SMRIN')) {

    p <- read_tsv(pheno_fname, quote = "", skip = 10)

    if (length(ancestry) == 1) {
        ancestry <- c(ancestry)
    }
    donorsub <- p %>%
        dplyr::rename(DONORID = SUBJID) %>%
        dplyr::filter(RACE %in% ancestry) %>%
        dplyr::rename(ancestry = RACE) %>%
        dplyr::select(DONORID, any_of(p_vars))
    a <- read_tsv(attr_fname, quote = "")
    if (!is.null(tissues_fname)) {
        tissues <- readLines(tissues_fname)
        tnum_df <- tibble(tissue_number=seq(1,length(tissues)),
                          SMTSD=tissues)
        a <- merge(tnum_df,a)
        if (!is.null(tissue_num)) {
            t_encode <- setNames(tissues, seq(1, length(tissues)))
            if (length(tissue_num) > 1) {
                ts <- t_encode[tissue_num] %>% as.vector
            }
            else {
                ts <- c(t_encode[[tissue_num]])
            }
            a <- a %>% dplyr::filter(SMTSD %in% ts)
        }

    }

    vars <- c('SAMPID','DONORID','SMTS','tissue_number','SMTSD', a_vars)

    a <- a %>%
        separate(.,
                 SAMPID,
                 c("PREFIX", "DONORID", NA, NA, NA),
                 sep = "-",
                 remove = FALSE) %>%
        mutate(DONORID = paste(PREFIX, DONORID, sep = "-")) %>%
        dplyr::select(any_of(vars)) %>%
        # add file for matching with filenames separated by tissue
        mutate(SMTSD_fileprefix = gsub(' ','_',gsub('\\(|\\)','',gsub(" - ","_",SMTSD))))

    print("Merging sex/ancestry info...")
    a <- merge(a, donorsub)
    if (!is.null(admix_fname)) {
        print("Loading 117 admixed DONORIDs from Gay et al.")
        ax <- read_excel(admix_fname, skip = 0, sheet = 2, col_names = c("DONORID")) %>%
            pull(DONORID)
        a <- a %>% filter(DONORID %in% ax)
    }

    return(
        a %>%
            dplyr::rename(ischemic_time=SMTSISCH,
                          RIN=SMRIN) %>%
            arrange(SMTSD,SAMPID)
    )

}

get_tissue_SAMPID <- function(attr_fname, tissues_fname, tissue_num) {
    tissues <- readLines(tissues_fname)
    t_encode <- setNames(tissues, seq(1, length(tissues)))
    if (length(tissue_num) > 1) {
        ts <- t_encode[tissue_num] %>% as.vector
    } else {
        ts <- c(t_encode[[tissue_num]])
    }
    a <- read_tsv(attr_fname, quote = "") %>%
        dplyr::filter(SMTSD %in% ts) %>%
        pull(SAMPID)
    return(a)
}

readGenoMatAll <- function(geno_fname, genotype="GT", samps=NULL) {

    gdf <- read_tsv(geno_fname, comment="##") %>%
        dplyr::select(ID, dplyr::starts_with('GTEX'))
    rnames <- gdf %>% pull(ID)
    genoMat <- as.matrix(gdf[,-1])
    rownames(genoMat) <- rnames

    if(genotype == "GL") {
       genoMat <- getGeno(genoMat)
    }else if(genotype == "GT") {
     genoMat <- getGenoGT(genoMat)
    }
    if (!is.null(samps)) {
        genoMat <- as.data.frame(genoMat) %>%
                        dplyr::select(samps) %>%
                        as.matrix()
    }

    if (nrow(genoMat) < 1){
        warning("No SNPs found for range")
        return(NULL)
    }

    return(genoMat)
}

## Gets all SNPs in a given region from a VCF file
regionSNPs <- function(vcfFile, chrom=NULL, startB=NULL, endB=NULL, grange = NULL, genotype = "GT", samps=NULL){ # Do not modify genotype or else function will not workd

	# Either explicit coordinates or granges
	if ( is.null(grange) &
      is.null(chrom) &
      is.null(startB) &
      is.null(endB) ) {

    param <- ScanVcfParam() # read everthing

  } else if (is.null(grange)) {

    param <-  GRanges(chrom, IRanges(startB, endB))

  } else {

    param <- grange

  }


  tab <- TabixFile(vcfFile, yieldSize=NA_integer_)
	open(tab)
  vcfInd <- samples(scanVcfHeader(vcfFile))
	genoMat <- readGeno(vcfFile, genotype, param=param)
  colnames(genoMat) <- vcfInd
  close(tab)

	if(genotype == "GL") {
	   genoMat <- getGeno(genoMat)
  }else if(genotype == "GT") {
     genoMat <- getGenoGT(genoMat)
  }
	if (!is.null(samps)) {
    	genoMat <- as.data.frame(genoMat) %>%
                		dplyr::select(samps) %>%
                		as.matrix()
	}

	if (nrow(genoMat) < 1){
		warning("No SNPs found for range")
		return(NULL)
	}

	return(genoMat)

}

# Given the genotype matrix "GL" of posterior probabilities
# this functions identifies the genotype with a posterior cutoff

getGeno <- function(genoMat, posterior = 0.9, parallel = F, ...){

	#Gets info from matrix
	sizeDim <- dim(genoMat)
	nPos <- sizeDim[1]
	nInd <- sizeDim[2]

	# passes the cutoff
	sharePos <- vector(mode = "numeric")
	genoMatLog <- genoMat >= posterior
	sharedGeno <- rowSums(genoMatLog) == nInd

	#Selects the matrix from those positions
	genoMat <- genoMat[sharedGeno,,]
	sizeDim <- dim(genoMat)
	nPos <- sizeDim[1]
	nInd <- sizeDim[2]

	# Gets the genotype for each position
	genotype <- matrix(0, nrow = nPos, ncol = nInd, dimnames= list(rownames(genoMat[,,1]), colnames(genoMat[,,1]) ) )
	for	(i in 1:nInd){
		#if (parallel){
		#	genotype[,i] <- do.call(c,bplapply(as.data.frame(t(genoMat[,i,])), function (x) which.max(x) , ...))
		#}else{
			genotype[,i] <- apply(genoMat[,i,], 1, function (x) which.max(x))
		#}
		#sharePos <- unique(c(sharePos, which(sharePosTemp)))
	}

	return(genotype)
}

# From genotype matrix GT gets genotypes in the form
# 0,1,2 from 0/0 0/1 1/1
getGenoGT <- function(x) {
    x <- as.matrix(x)
    allGeno <- c("0|0", "0|1", "1|0", "1|1")

    x[!x %in% allGeno] <- NA
    x[x == allGeno[1]] <- "0"
    x[x == allGeno[2]] <- "1"
    x[x == allGeno[3]] <- "1"
    x[x == allGeno[4]] <- "2"

    rowNames <- rownames(x)

    if(nrow(x) == 1) {
        x <- apply(rbind(x,x), 2, as.numeric)
        x <- x [-1, , drop=F]
    } else {
        x <- apply(x, 2, as.numeric)
    }

    rownames(x) <- rowNames
    return(x)

}

get_snp_matrix <- function(geno) {
	snps <- geno %>%
		dplyr::select(variant_id, starts_with('GTEX'))
	rnames <- snps %>% pull(variant_id)
	snps <- as.matrix(snps[,-1])
	rownames(snps) <- rnames
	return(snps)

}

get_covariate_matrix <- function(c) {

    # get covariate matrix with rows as covariates
    # columns as samples

    rnames <- c %>% pull(1)
    c <- c[,-1]
    c <- as.matrix(c)
    rownames(c) <- rnames
    return(c)

}

get_tissue_covariates_matrix <- function(cvrt_dir, tissue_num, use_samp, tissue=NULL) {

  if (!is.null(cvrt_dir)) {

    if (!is.null(tissue)) {

      cov_base <- gsub(' ','_',gsub('\\(|\\)','',gsub(" - ","_",tissue)))
      cov_fname <- paste0(
        cvrt_dir,
        '/',
        cov_base,
        '.v8.covariates.txt'
      )

    } else {

      if (tissue_num==46) {
        cov_fname <- paste0(
          cvrt_dir,
          '/',
          'Skin_Not_Sun_Exposed_Suprapubic.v8.covariates.txt'
        )
      } else if (tissue_num==47) {
        cov_fname <- paste0(
          cvrt_dir,
          '/',
          'Skin_Sun_Exposed_Lower_leg.v8.covariates.txt'
        )
      } else {
        stop(paste("Tissue number",tissue_num,"not currently supported."))
      }

    }

    cov <- tryCatch(
      {
        cov <- read_tsv(cov_fname)
      },
      error=function(e) {
        message(e)
        return(NULL)
      }
    )

    if (is.null(cov)) {
      cvrt <- NULL
    } else {
      print(paste(cov_fname, "successfully read."))
      cvrt <- get_covariate_matrix(cov)
      cvrt <- cvrt[,dplyr::intersect(use_samp,colnames(cvrt))]

      if (ncol(cvrt) < nrow(cvrt)) { # ensure number of covariates used is < number of samples
          cvrt <- cvrt[c('PC1','PC2','PC3','PC4','PC5','pcr','platform','sex'),]
      }
    }

  } else {
    cvrt <- NULL
  }

  return(cvrt)

}

get_tissue_coldata_deseq2 <- function(cvrt_dir,
                                         lookup_table,
                                         tissue_fileprefix=NULL,
                                         use_samp=NULL,
                                         debug=FALSE,
                                         cvars=c('SEX','AGE','ischemic_time','RIN','pcr','platform','ancestry')) {

  if (is.null(tissue_fileprefix)) {

      tissue_fileprefix <- lookup_table %>%
          pull(SMTSD_fileprefix) %>%
          unique

      if (length(tissue_fileprefix)>1) {
          stop("More than one tissue supplied in lookup_table.")
      }

  }

  if (is.null(use_samp)) {

      use_samp <- lookup_table %>%
          pull(DONORID) %>%
          unique

  }

  if (!is.null(cvrt_dir)) {

      cov_fname <- paste0(
        cvrt_dir,
        '/',
        tissue_fileprefix,
        '.v8.covariates.txt'
      )

    cov <- tryCatch(
      {
        cov <- read_tsv(cov_fname)
      },
      error=function(e) {
        message(e)
        return(NULL)
      }
    )

    if (debug) {
        return(cov)
    }

    if (is.null(cov)) {
      cvrt <- NULL
    } else {
      print(paste(cov_fname, "successfully read."))

      cvrt <- transpose_tibble(cov, colnames_to='DONORID') %>%
            dplyr::select(DONORID,pcr,platform,sex) %>%
            dplyr::filter(DONORID %in% use_samp)

      coldata <- merge(lookup_table,cvrt)

      lostsamps <- length(use_samp) - nrow(coldata)
      print(paste0("Samples missing GTEx covariate data = ",lostsamps))

      sexmatch <- coldata %>%
            dplyr::filter(SEX==sex) %>%
            nrow()
      if (sexmatch==nrow(coldata)) {
          print("Sexes from phenotypes and covariates have same encoding.")
      } else if (sexmatch==0) {
          print("Sexes from phenotypes and covariates have opposite encoding.")
      } else {
          stop("Sexes from phenotypes and covariates are mismatched!")
      }

      coldata <- coldata %>%
            dplyr::select(DONORID,SAMPID,any_of(cvars)) %>%
            mutate(SEXn = case_when(SEX==1 ~ 'male',
                                    SEX==2 ~ 'female'),
                   ancestryn = case_when(ancestry==2 ~ 'AFR',
                                         ancestry==3 ~ 'EUR'),
                   pcr = as.factor(pcr),
                   platform = as.factor(platform)) %>%
            dplyr::select(-c(SEX,ancestry)) %>%
            dplyr::rename(sex=SEXn, ancestry=ancestryn, age=AGE, rin=RIN) %>%
            dplyr::select(DONORID,SAMPID,str_to_lower(cvars))

    }

  } else {
    coldata <- NULL
  }

  return(
      type_convert(coldata, col_types='ccfdddfff')
  )

}

direction_agreement <- function(x) {
    if (length(x) > 1) {
        return( (all(x < 0) | all(x > 0)) )
    } else {
        return(NA)
    }
}



################ gene and variant formatting functions #################

get_gene_locs <- function(GRCh=38) {
    library(biomaRt)
    if (GRCh==37) {
        print("Getting GRCh37 coordinates...")
    } else {
        print("Getting GRCh38 coordinates...")
        GRCh=NULL
    }
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=GRCh)
    gene_attributes = c("ensembl_gene_id",
                        "hgnc_symbol",
                       "chromosome_name",
                       "start_position",
                       "end_position",
                       "strand")
    gene_info = getBM(attributes=gene_attributes,
                      mart=ensembl)

    gene_info <- gene_info %>%
        as_tibble() %>%
        dplyr::filter(chromosome_name %in% c(1:22,'X','Y'))
    return(gene_info)

}

add_gene_locs <- function(df,GRCh=38,grange_format=FALSE) {
  gene_locs_bm <- get_gene_locs(GRCh=GRCh)

  if (grange_format) {

    gene_locs_tomerge <- gene_locs_bm %>%
        mutate(seqnames=paste0('chr',chromosome_name),
               start=start_position,
               end=end_position) %>%
        dplyr::select(ensembl_gene_id, seqnames, start, end, strand)

  } else {
    gene_locs_tomerge <- gene_locs_bm %>%
        mutate(gene_loc=paste0('chr',
                               chromosome_name,
                               '_',
                               start_position,
                               '_',
                              end_position)) %>%
        dplyr::select(ensembl_gene_id, gene_loc)

  }

  df <- merge(df, gene_locs_tomerge,
              by='ensembl_gene_id', all.x=TRUE, sort=FALSE)

  return(df)

}

add_hgnc <- function(df, eGene=FALSE) {
  library(biomaRt)

  if (eGene) {
    ensgenes <- df %>%
        pull(eGene) %>%
        unique()
    merge_col <- 'eGene'
  } else {
    ensgenes <- df %>%
        pull(ensembl_gene_id) %>%
        unique()
    merge_col <- 'ensembl_gene_id'
  }

  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  gene_attributes = c("ensembl_gene_id","hgnc_symbol")
  print("Getting hgnc_symbols from biomaRt...")
  gene_info = getBM(attributes=gene_attributes,filters="ensembl_gene_id",values=ensgenes,mart=ensembl)
  if (eGene) {
    gene_info <- gene_info %>%
      as_tibble() %>%
      dplyr::rename(eGene=ensembl_gene_id,
                   eGene_hgnc=hgnc_symbol)
  }
  print("Adding hgnc_symbols...")
  df <- merge(df, gene_info, by=merge_col, all.x=TRUE)

  return(df)
}

get_hgnc <- function(ensgenes) {
    # returns vector of hgnc symbols corresponding
    # to input ensembl gene id vector

    library(biomaRt)
    ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    gene_attributes = c("ensembl_gene_id", "hgnc_symbol")
    print("Getting hgnc_symbols from biomaRt...")
    gene_info = getBM(attributes = gene_attributes, filters = "ensembl_gene_id",
        values = ensgenes, mart = ensembl)
    return(gene_info %>% pull(hgnc_symbol))
}

add_rsIDs <- function(df, format=FALSE) {
    # adds rsID column to granges-formatted df
    if (format) {
      df <- df %>%
        separate(variant_id,
                 c('seqnames','start','ref','alt','build'),
                 sep="_",
                 remove=FALSE,
                 convert=TRUE) %>%
        mutate(end=start)
    }

    library(SNPlocs.Hsapiens.dbSNP151.GRCh38)

    h2rs <- df %>%
        mutate(seqnames=gsub('chr','',seqnames)) %>%
        as_granges()

    # print('Getting rsIDs with SNPloc...')
    rsids <- as_tibble(snpsByOverlaps(SNPlocs.Hsapiens.dbSNP151.GRCh38, h2rs)) %>%
                dplyr::rename(start=pos) %>%
                mutate(seqnames=paste0('chr',seqnames),
                      end=start) %>%
                dplyr::rename(rsID=RefSNP_id) %>%
                dplyr::select(-c(strand,alleles_as_ambig))

    # print('Removing SNPs with no rsID...')
    dfrs <- as_tibble(merge(df, rsids))

    print(paste0(nrow(df)-nrow(dfrs)," variants dropped due to no rsID."))

    return(dfrs)
}

ld_expand_df <- function(r, plink_bed_dir, plink_pre,
                         chrcol='seqnames',
                         poscol='start',
                         rscol='rsID',
                         plink_dir="/home/trevor/Programs",
                         sep="") {

  cursnp = r[rscol]

  curchr = as.character(r[chrcol])
  curpos = as.integer(r[poscol])

  # print(
  #   paste(curchr,curpos,sep='_')
  # )

  curchrname = paste(plink_pre,curchr,sep=sep)
  # print(curchrname)
  syscomm <- paste0(
      plink_dir,
      "/plink --bfile ",
      plink_bed_dir,
      "/",
      curchrname,
      " --r2 --ld-snp ",
      cursnp,
      " --ld-window 9999 --ld-window-kb 500 --ld-window-r2 .8 --out /tmp/snpld",
      cursnp
  )
  # print(syscomm)

  cursnpcall = system(syscomm,ignore.stdout=TRUE,ignore.stderr=TRUE)

  curreadfile = paste0("/tmp/snpld",cursnp,".ld")

  if(file.exists(curreadfile)) {

    cursnpdatac <- tryCatch(
      {
        cursnpdata <- read.table(curreadfile,header=TRUE) %>%
            as_tibble() %>%
            mutate(coord=paste0('chr',CHR_B,'_',BP_B))

        cursnpdatac = as.character(cursnpdata$coord)
      },
      error=function(e) {
        message(e)
        return("novalidvar")
      }
    )

  } else {
      print('No LD file found')
  }

  if(!file.exists(curreadfile)) {

    cursnpdatac = "novalidvar"

  }

  #deal with rsIDs separated by semicolon that point to same coordinates
  snpIDs = unlist(strsplit(cursnpdatac,";"))

  if(is.null(snpIDs)) {
    snpIDs = cursnp
  }

  if(snpIDs == "novalidvar") {
    snpIDs = cursnp
  }

  return(as.list(snpIDs))

}
