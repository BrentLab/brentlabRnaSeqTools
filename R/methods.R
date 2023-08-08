# brentlabRnaSeqSet methods ----------------------------------------------------

## coercetoDds

#' coerce brentlabRnaSeqSet to DESeqDataSet
#'
#' @rdname coerceToDds
#'
#' @description coerce brentlabRnaSeqSet to DESeqDataSet
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix sizeFactors
#'
#' @param x a brentlabRnaSeqSet
#'
#' @export
setMethod("coerceToDds", "brentlabRnaSeqSet", function(x) {

  size_factors = sizeFactors(x)
  design = design(x)
  colData = extractColData(x)
  counts = counts(x)

  dds = DESeqDataSetFromMatrix(
    colData = colData,
    countData = counts,
    design = design
  )

  if(!is.null(size_factors)){
    sizeFactors(dds) = size_factors
  }

  dds
})

## createExperimentSet ----

#'
#' Create a pre-defined experiment set
#'
#' @description for pre-defined experiment definitions, see the github repo
#'   R/ExperimentSetFunctions.R
#'
#' @rdname createExperimentSet
#'
#' @param x a brentlabRnaSeqSet object
#' @inheritParams brentlabRnaSeqExperiment
#'
#' @export
setMethod("createExperimentSet", "brentlabRnaSeqSet", function(x, set_name) {

  subset_blrs = -1

  if(set_name == "ninetyMin_2016Grant"){
    subset_blrs = createNinetyMinuteInduction_2016grant(x)
  }
  else if(set_name == "ninetyMin_2016GrantWithDoubles"){
    subset_blrs = createNinetyMinuteInduction_2016grantWithDoubles(x)
  }
  else if(set_name == "ninetyMin_non2016grant"){
    subset_blrs = createNinetyMinuteInduction_non2016grant(x)
  }
  else if(set_name == "ninetyMin_all"){
    subset_blrs = createNinetyMinuteInduction_all(x)
  }
  else if(set_name == "envPert_epWT"){
    subset_blrs = createEnvPert_epWT(x)
  }
  else if(set_name == "envPert_perturbed"){
    subset_blrs = createEnvPertSet_perturbed(x)
  }
  else if(set_name == "envPert_titrationWT"){
    subset_blrs = createEnvPertSet_titrationWT(x)
  }

  if(is.numeric(subset_blrs)){
    stop(paste0("set_name '", set_name,
                "' is not recognized. see ?createExperimentSet for
                acceptable set_names"))
  }

  brentlabRnaSeqExperiment(subset_blrs, set_name)
})

## qaFilter(brentlabRnaSeqSet) ----

#'
#' Quality filter a brentlabRnaSeqSet object
#'
#' Filter for manual passes (overrides auto fail),
#'   automatic passes (unless auto failed), and optionally RLE IQR
#'
#' @note on the brentlabRnaSeqSet (as opposed to the brentlabRnaSeqExperiment),
#'   this only filters on manualAudit and autoAudit. no rle filtering
#'
#' @rdname qaFilter
#'
#' @return a quality assessment filtered brentlabRnaSeqSet object
#'
#' @export
setMethod("qaFilter",
          "brentlabRnaSeqSet",
          function(x) {


  passing_fastqFileNumbers = as_tibble(colData(x)) %>%
    filter(manualAudit == FALSE |
             (is.na(manualAudit) & autoAudit == FALSE)) %>%
    pull(fastqFileNumber)

  # create a sample quality filtered brentlabRnaSeqSet
  blrs_qc_fltr = x[,colData(x)$fastqFileNumber %in% passing_fastqFileNumbers]

  blrs_qc_fltr
})

## splitProtocolGroups ----

#' Split a brentlabRnaSeqSet into two by protocol
#'
#' @description split an object into 'SolexaPrep' and 'E7420L'
#'
#' @rdname splitProtocolGroups
#'
#' @param x a brentlabRnaSeqSet
#'
#' @return a list 'SolexaPrep' and 'E7420L' group. We call 'SolexaPrep' "old"
#'
#' @export
setMethod("splitProtocolGroups", "brentlabRnaSeqSet", function(x) {

  solexaprep = x[,colData(x)$libraryProtocol == 'SolexaPrep']
  colData(solexaprep)$libraryProtocol =
    droplevels(as.factor(colData(solexaprep)$libraryProtocol))

  # helper function to perform split and clean up relevant factored columns
  # TODO: think about which columns to factor in R/Database.R, and which others
  # to clean here
  cleanSplit = function(obj, protocol){
    split_obj = obj[,colData(obj)$libraryProtocol == as.character(protocol)]

    colData(split_obj)$libraryProtocol =
      droplevels(as.factor(colData(split_obj)$libraryProtocol))

    colData(split_obj)$libraryDate =
      droplevels(colData(split_obj)$libraryDate)

    split_obj
  }

  list(
    'SolexaPrep' = cleanSplit(x, 'SolexaPrep'),
    'E7420L' = cleanSplit(x, 'E7420L')
  )
})

## extractColData ----

#' extract colData as tibble
#'
#' @importFrom tibble as_tibble
#'
#' @rdname extractColData
#'
#' @param x a brentlabRnaSeqSet
#'
#' @return colData of the object as a tibble
#'
#' @export
setMethod("extractColData", "brentlabRnaSeqSet", function(x) {
  as_tibble(colData(x))
})

#' extract the model.matrix using the design and colData
#'
#' @rdname extractDesignMatrix
#'
#' @importFrom rlang is_formula
#'
#' @param x a brentlabRnaSeqSet
#'
#' @return the model matrix, using the design and coldata
#'
#' @export
setMethod("extractDesignMatrix", "brentlabRnaSeqSet", function(x) {

  if(is_formula(design(x))){
    message(paste0("Model Matrix for design: ", as.character(design(x))))
    model.matrix(design(x), extractColData(x))
  } else if(is.matrix(design(x))){
    design(x)
  } else{
    stop("data type of the object is not recognized")
  }

})

# brentlabRnaSeqExperiment methods ---------------------------------------------

## qaFilter(brentlabRnaSeqExperiment) ----

#'
#' Quality filter a brentlabRnaSeqExperiment object
#'
#' Filter for manual passes (overrides auto fail),
#'   automatic passes (unless auto failed), and optionally RLE IQR
#'
#' @note If there are less than 3 replicates, the RLE stats will be NA. These
#'   are returned, also. Filter them out of the returned object if you wish to
#'   have only those samples with more than 3 replicates.
#'
#' @rdname qaFilter
#'
#' @importFrom dplyr filter as_tibble
#'
#' @param x a brentlabRnaSeqSet object
#' @param rle_iqr_threshold default is NA, which means no rle_iqr_threshold filtering.
#'   Set a value to filter based on RLE IQR. Note: RLE stats are all calculated
#'   after removing the libraryDate effect.
#' @param iqr_colname default is NA, which creates the iqr col name from the
#'   set name in the object. Pass the col name itself, eg "my_set_iqr", to
#'   specify a iqr column name manually
#'
#' @return a quality assessment filtered brentlabRnaSeqExperiment object
#'
#' @export
setMethod("qaFilter",
          "brentlabRnaSeqExperiment",
          function(x, rle_iqr_threshold = NA, iqr_colname = NA) {

    # note: I tried the filter commented out below first, but had trouble with
    #       handling NAs. dplyr filter works. Used the same method for rle
    #       filtering below. This should be fixed for clarity/speed/reducing
    #       dependencies eventually

    # # create a filter on manualAudit and autoAudit
    # qc1_fltr = (colData(x)$manualAudit == FALSE |
    #              is.na(colData(x)$manualAudit)) &
    #             colData(x)$autoAudit == FALSE

    passing_fastqFileNumbers = as_tibble(colData(x)) %>%
      filter(manualAudit == FALSE |
               (is.na(manualAudit) & autoAudit == FALSE)) %>%
      pull(fastqFileNumber)

     # create a sample quality filtered brentlabRnaSeqSet
     blrs_qc_fltr = x[,colData(x)$fastqFileNumber %in% passing_fastqFileNumbers]

     # if rle_iqr_threshold is set, filter based on IQR
     if(!is.na(rle_iqr_threshold) & is.numeric(rle_iqr_threshold)){


       iqr_colname =
        if(is.na(iqr_colname)){
          iqr_colname = paste0(x@set_name, "_iqr")
         } else{
           iqr_colname
         }

       if(!iqr_colname %in% colnames(colData(blrs_qc_fltr))){
         stop("iqr column, '", iqr_colname, "', is not in the metadata. This
              is a problem in the definition of the ExperimentSet object --
              please issue a report on github. Paste this error, and the name
              of the experiment set that you are using.")
       }

       iqr_passing_fastqFileNumbers = as_tibble(colData(blrs_qc_fltr)) %>%
         filter(!!rlang::sym(iqr_colname) < rle_iqr_threshold |
                  is.na(!!rlang::sym(iqr_colname))) %>%
         pull(fastqFileNumber)

       blrs_qc_fltr =  blrs_qc_fltr[,colData(blrs_qc_fltr)$fastqFileNumber %in%
                                      iqr_passing_fastqFileNumbers]
     }

     blrs_qc_fltr
})

## estimateSizeFactorsByProtocol ----

#'
#' Calculate Size Factors within protocol groups
#'
#' Size Factors will be calculated for the set libraryProtocol 'SolexaPrep'
#'   and the complement separately
#'
#' @rdname estimateSizeFactorsByProtocol
#'
#' @importFrom DESeq2 estimateSizeFactors
#' @importFrom SummarizedExperiment cbind
#'
#' @param x a brentlabRnaSeqSet object
#'
#' @return a brentlabRnaSeqExperiment with sizeFactors estimated from within
#'   protocol groups (SolexaPrep and not SolexaPrep)
#'
#' @export
setMethod("estimateSizeFactorsByProtocol",
          "brentlabRnaSeqExperiment",
          function(x) {

  old_dds = estimateSizeFactors(x[, colData(x)$libraryProtocol == 'SolexaPrep'])
  new_dds = estimateSizeFactors(x[, colData(x)$libraryProtocol != 'SolexaPrep'])

  SummarizedExperiment::cbind(old_dds, new_dds)

})

## replicateByProtocolTally ----

#'
#' experimental replicates by libraryProtocol tally
#'
#' @description tally experimental replicate columns against the library prep
#'   protocol
#'
#' @rdname replicateByProtocolTally
#'
#' @importFrom tibble as_tibble
#'
#' @param x a brentlabRnaSeqExperiment
#'
#' @export
setMethod("replicateByProtocolTally",
          "brentlabRnaSeqExperiment",
          function(x) {
  if(x@set_name == "ninetyMin_2016Grant"){
    # see R/ExperimentSetFunctions.R
    replicateByProtocolTally_90min(as_tibble(colData(x)))
  } else{
    message(paste0('There is no replicateByProtocolTally for set: ', x@set_name))
  }

})

## test_train_partition ----

#' create test train set
#'
#' @description For brentlabRnaSeqExperiment objects, this offers methods to
#'   create train test sets
#'
#' @rdname test_train_partition
#'
#' @param x a brentlabRnaSeqExperiment object
#' @param min_set_size the minimum replicate set size to consider for hold outs
#'
#' @return a list with slots train and test, each with brentlabRnaSeqSetExperiments
#'
#' @export
setMethod("test_train_partition",
          "brentlabRnaSeqExperiment",
          function(x, min_set_size) {

    if(x@set_name == 'ninetyMin_2016Grant'){

      genotype_list = as_tibble(colData(x)) %>%
        group_by(genotype1) %>%
        tally() %>%
        filter(genotype1 != "CKF44_00000" & n > min_set_size) %>%
        pull(genotype1)

      test_ffn = as_tibble(colData(x)) %>%
        filter(genotype1 %in% c(genotype_list)) %>%
        group_by(genotype1) %>%
        sample_n(1) %>%
        pull(fastqFileNumber)

      test = x[,colData(x)$fastqFileNumber %in% test_ffn]
      colData(test) = colData(test,droplevels(as_tibble(colData(test))))

      # note -- this will include all wildtypes, also
      train = x[,!colData(x)$fastqFileNumber %in% test_ffn]
      colData(train) = colData(train,droplevels(as_tibble(colData(train))))

    } else{

      stop(paste0("There is no defined method for set: ", x@set_name, ". If you need one,
           post an issues report along with a completed explanation of how you
           would like it done."))

    }

    list(
      train = train,
      test = test
    )
})


#' @rdname VariantExplorer
#' @param x a VariantExplorer object
#' @inheritParams brentlabRnaSeqExperiment
setMethod("igv_script", "VariantExplorer", function(x,...){

  # check igv_genome attribute
  if (identical(x@igv_genome, character(0))){
    stop("The igv_genome slot must exist")
  } else if(!file.exists(x@igv_genome)){
    stop("The path to the igv genome is not valid")
  }

  # create the locus granges
  granges = x@gff[x@gff$gene == locus & x@gff$type == "gene"]
  if (length(granges) == 0){
    stop(sprintf("Gene: %s is not recognized -- no ranges present in the gff",
                 locus))
  }

  # create output directories
  script_output_dir = file.path(output_dir, "igv", "igv_scripts")
  dir.create(script_output_dir, recursive = TRUE)
  browser_output_dir = file.path(output_dir, "igv", "browser_shots")
  dir.create(browser_output_dir, recursive=TRUE)


  # create output paths
  if (image_basename == ""){
    script_output_path = file.path(script_output_dir,paste0(locus,".txt"))
    browser_output_filepath = file.path(browser_output_dir,paste0(locus,".png"))
  } else{
    script_output_path = file.path(script_output_dir,paste0(image_basename,".txt"))
    browser_output_filepath = file.path(browser_output_dir,paste0(image_basename,".png"))
  }

  # get bam list
  sample_id_list = x@variants %>%
    filter(gene_id == locus) %>%
    pull(sample) %>%
    unique()

  bam_list = x@metadata %>%
    filter(sample_id %in% sample_id_list) %>%
    pull(bam) %>%
    unique()

  # write igv batchscript lines
  load_samples = paste0(map(bam_list, ~sprintf("\tload %s\n", .)), collapse=" ")
  parsed_range = paste(as.character(granges@seqnames),
                       as.character(granges@ranges), sep=":")
  batch_script = paste0("new\nsnapshotDirectory %s\nmaxPanelHeight %s\ngenome %s\n",
                        load_samples, "goto %s", collapse = " ")
  batch_script = sprintf(batch_script,
                         path.expand(output_dir),
                         maxPanelHeight,
                         path.expand(x@igv_genome),
                         parsed_range)
  if(exit_browser){
    batch_script = paste0(batch_script, "\nsnapshot %s\nexit", collapse = " ")
    batch_script = sprintf(batch_script, browser_output_filepath)
  }

  cat(batch_script, file = script_output_path)
})

#' @export
#' @rdname VariantExplorer
#' @param object a VariantExplorer instance
setMethod('summary', signature(object = 'VariantExplorer'), function(object){
  split_variants = object@variants %>%
    ungroup() %>%
    group_by(effect,impact) %>%
    tally() %>%
    dplyr::rename(num_genes = n) %>%
    ungroup() %>%
    group_by(impact) %>%
    arrange(impact, desc(num_genes), .by_group = TRUE) %>%
    droplevels() %>%
    group_split()


   names(split_variants) = unlist(map(split_variants, ~unique(pull(.,impact))))

   pretty_print = function(impact_table, impact_level){

       message(sprintf("Impact level: %s", impact_level))
     dplyr::select(impact_table,-impact) %>% print()
   }

   for(i in names(split_variants)){
     pretty_print(split_variants[[i]], i)
   }


})

#' Visualize coverage at a locus
#' @rdname VariantExplorer
#' @param x a VariantExplorer object
#' @param gene the name of the gene you wish to visualize
#' @param sample the sample_id you wish to see
#' @param plot_type either 'pileup' or 'coverage'. Defaults to 'pileup'
#'
#' @importFrom withr local_options
#' @importFrom Gviz AlignmentsTrack AnnotationTrack GenomeAxisTrack SequenceTrack plotTracks
#' @importFrom dplyr filter pull
#'
#' @export
setMethod("visualize", "VariantExplorer",
          function(x, gene, sample, plot_type = 'pileup'){
            # set options locally -- this will revert to outer when function completes
            withr::local_options(list(ucscChromosomeNames = FALSE))

            gene_granges = x@gff[x@gff$gene == gene]
            bamfile = filter(ve@metadata, sample_id == sample) %>%
              pull(bam)


            track_list = list(
              axis_track = GenomeAxisTrack(),
              annotation = AnnotationTrack(gene_granges, group = gene_granges$ID),
              alignment = AlignmentsTrack(bamfile,
                                          isPaired = TRUE,
                                          showMismatches = TRUE,
                                          showIndels = TRUE,
                                          type = plot_type),
              sequence  = SequenceTrack(ve@bsgenome,chromosome = unique(gene_granges$seqnames))
            )


            plotTracks(track_list,
                       chromosome = unique(gene_granges$seqnames),
                       start = gene_granges[gene_granges$type == 'gene']$start,
                       end = gene_granges[gene_granges$type == 'gene']$end,
                       add53 = TRUE,
                       showId = TRUE,
                       labelPos = 'below')
})
