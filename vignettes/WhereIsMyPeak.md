---
title: "WhereIsMyPeak"
author: "Thomas Schwarzl"
date: "24 March 2015"
output: html_document
---

# Where is my peak? 

WhereIsMyPeak does genomic annotation with priorisation what gene types you are interested. This package was inspired
by ChIPseeker, however overcomes limited functionlity in the genomic annotation.


## Install this package via Github

The package is Open Source and can be downloaded here
- [Github](https://github.com/Distue/WhereIsMyPeak "Where is my peak")

In R you can directly install it with this command

<div class="chunk" id="unnamed-chunk-1"><div class="rcode"><div class="source"><pre class="knitr r">library(devtools)
devtools::install_github("distue/WhereIsMyPeak")
</pre></div>
<div class="message"><pre class="knitr r">## Downloading github repo distue/WhereIsMyPeak@master
## Installing WhereIsMyPeak
## '/Library/Frameworks/R.framework/Resources/bin/R' --vanilla CMD INSTALL  \
##   '/private/var/folders/x6/5mjddpkd3y327ppd7njqzvkc0000gn/T/RtmphDgRT3/devtoolsf47b7cd4de87/Distue-WhereIsMyPeak-6c57b40'  \
##   --library='/Library/Frameworks/R.framework/Versions/3.1/Resources/library'  \
##   --install-tests 
## 
## Reloading installed WhereIsMyPeak
</pre></div>
</div></div>


Load the package
<div class="chunk" id="unnamed-chunk-2"><div class="rcode"><div class="source"><pre class="knitr r">suppressPackageStartupMessages(library(WhereIsMyPeak))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(biomaRt))
</pre></div>
</div></div>

## Load peaks

Peaks or Intervals are stored as GRanges object. 

Here we load a small example data set
<div class="chunk" id="unnamed-chunk-3"><div class="rcode"><div class="source"><pre class="knitr r">data("peaks", package="WhereIsMyPeak")
</pre></div>
</div></div>

Let's preview the peaks
<div class="chunk" id="unnamed-chunk-4"><div class="rcode"><div class="source"><pre class="knitr r">head(peaks)
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 6 ranges and 2 metadata columns:
##       seqnames                 ranges strand |      name       score
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt; | &lt;integer&gt;   &lt;numeric&gt;
##   [1]        1 [153643717, 153643723]      + |      5039 9.657284435
##   [2]       12 [ 62997806,  62998017]      + |     15618 8.840176975
##   [3]        1 [204380047, 204380178]      - |      6622 8.663731336
##   [4]       12 [125348232, 125348415]      - |     17293 8.381262242
##   [5]        1 [205686650, 205686793]      - |      6738 8.362087961
##   [6]       11 [118888683, 118888771]      - |     13136 8.354382931
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>

## Load annotation

We can load the annotation using a TxDB Ensembl Build from biomaRt. 
<div class="chunk" id="unnamed-chunk-5"><div class="rcode"><div class="source"><pre class="knitr r">txdb <- makeTranscriptDbFromBiomart()
</pre></div>
</div></div>

We also can load an specific annotation via the ensembl Archive page:

<div class="chunk" id="unnamed-chunk-6"><div class="rcode"><div class="source"><pre class="knitr r"># ENSEMBL GENES 75
txdb <- makeTranscriptDbFromBiomart( biomart         = "ENSEMBL_MART_ENSEMBL",
                                     dataset         = "hsapiens_gene_ensembl",
                                     transcript_ids  = NULL,
                                     filters         = "",
                                     id_prefix       = "ensembl_",
                                     host            = "feb2014.archive.ensembl.org",
                                     port            = 80,
                                     miRBaseBuild    = NA )
</pre></div>
</div></div>


We load additional information about gene biotypes and gene description which is not available in the TxDB
(Note: This package is currently under development, this could be included in a custom TxDB build in future)

For this we can load the latest ensembl annotation
<div class="chunk" id="unnamed-chunk-7"><div class="rcode"><div class="source"><pre class="knitr r">ensembl.mart <- useMart(biomart="ensembl")
</pre></div>
</div></div>

or we can load a specific from the archive
<div class="chunk" id="unnamed-chunk-8"><div class="rcode"><div class="source"><pre class="knitr r">ensembl.mart <- useMart(host="feb2014.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL")        
</pre></div>
</div></div>

Then we load the dataset 

<div class="chunk" id="unnamed-chunk-9"><div class="rcode"><div class="source"><pre class="knitr r">ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl.mart) 
</pre></div>
</div></div>

and load the necessary tables for biotype and gene description
<div class="chunk" id="unnamed-chunk-10"><div class="rcode"><div class="source"><pre class="knitr r">mappingTable <- getBM(c("ensembl_gene_id", "gene_biotype") , filters="", values=T, uniqueRows = TRUE, mart=ensembl)  
rownames(mappingTable) <- mappingTable[,1]  

namingTable <- getBM(c("ensembl_gene_id", "external_gene_id", "description") , filters="", values=T, uniqueRows = TRUE, mart=ensembl)  
rownames(namingTable) <- namingTable[,1]
namingTable <- namingTable[,-1]
</pre></div>
</div></div>

First we proprocess the txdb since it is not stable, loses connection and serialization is often broken when saved to disk. Also, when annotating multiple batches of peaks (GRanges objects), we do not have to proprocess every single time

<div class="chunk" id="unnamed-chunk-11"><div class="rcode"><div class="source"><pre class="knitr r">txdbParsed <- initRegions(txdb, mappingTable, namingTable)
</pre></div>
</div></div>

## Annoation Priority

Now we choose an annotation priority.

The first position is always prioritized. 

When a peak hits multiple features, there are decisions made:

1) We prefer Exon annotation over Intron annotation
<div class="chunk" id="unnamed-chunk-12"><div class="rcode"><div class="source"><pre class="knitr r">annotationOrder = c("Exon", "Intron")
</pre></div>
</div></div>

2) We decide on Exon annation priority
In Exon annotation, we prefer, CDS over 5', and 5' over 3' annotation. The rest is annotated as Exon. Note that non-coding RNAs for example do not have a CDS, therfore they still remain as Exon. 

(Note: "rest" Exon over everything else will be implemented)  

<div class="chunk" id="unnamed-chunk-13"><div class="rcode"><div class="source"><pre class="knitr r">exonAnnotationOrder = c("CDS", "5' UTR", "3' UTR")
</pre></div>
</div></div>

3) We decide on types
At last, we set types priority.  

<div class="chunk" id="unnamed-chunk-14"><div class="rcode"><div class="source"><pre class="knitr r">typeOrder = c( "non_coding",
                 "known_ncrna",
                 "rRNA",
                 "snoRNA",
                 "snRNA",
                 "miRNA", 
                 "lincRNA", 
                 "Mt_rRNA", 
                 "Mt_tRNA",
                 "sense_intronic",
                 "transcribed_processed_pseudogene",
                 "transcribed_unitary_pseudogene",   
                 "transcribed_unprocessed_pseudogene",
                 "translated_processed_pseudogene",
                 "translated_unprocessed_pseudogene",
                 "unitary_pseudogene",
                 "unprocessed_pseudogene",
                 "sense_overlapping",
                 "3prime_overlapping_ncrna",   
                 "antisense",
                 "IG_C_gene",
                 "IG_D_gene",
                 "IG_J_gene",
                 "IG_V_gene",
                 "IG_C_pseudogene", 
                 "IG_J_pseudogene",
                 "IG_V_pseudogene", 
                 "polymorphic_pseudogene",
                 "misc_RNA",
                 "TEC",
                 "TR_C_gene",
                 "TR_D_gene",
                 "TR_J_gene",
                 "TR_J_pseudogene",
                 "TR_V_gene",
                 "TR_V_pseudogene",
                 "pseudogene",
                 "processed_pseudogene", 
                 "processed_transcript",
                 "LRG_gene", 
                 "protein_coding")
</pre></div>
</div></div>

## Annotate

Then simply with one line we can annotate all peaks 

<div class="chunk" id="unnamed-chunk-15"><div class="rcode"><div class="source"><pre class="knitr r">annotatedPeaks <- WhereIsMyPeak(peaks, txdbParsed, annotationOrder, exonAnnotationOrder, typeOrder)
</pre></div>
</div></div>




You will get back a GRanges object with the annotations stored in the elementMetadata columns
<div class="chunk" id="unnamed-chunk-17"><div class="rcode"><div class="source"><pre class="knitr r">head(annotatedPeaks) 
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 6 ranges and 8 metadata columns:
##       seqnames                 ranges strand |      name     score
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt; | &lt;integer&gt; &lt;numeric&gt;
##   [1]        1 [153643717, 153643723]      + |      5039  9.657284
##   [2]       12 [ 62997806,  62998017]      + |     15618  8.840177
##   [3]        1 [204380047, 204380178]      - |      6622  8.663731
##   [4]       12 [125348232, 125348415]      - |     17293  8.381262
##   [5]        1 [205686650, 205686793]      - |      6738  8.362088
##   [6]       11 [118888683, 118888771]      - |     13136  8.354383
##                    ID            V2
##              &lt;factor&gt;      &lt;factor&gt;
##   [1]              NA          &lt;NA&gt;
##   [2] ENSG00000257354 RP11-631N16.2
##   [3] ENSG00000158615      PPP1R15B
##   [4] ENSG00000073060        SCARB1
##   [5] ENSG00000069275        NUCKS1
##   [6] ENSG00000118181         RPS25
##                                                                                                 V3
##                                                                                           &lt;factor&gt;
##   [1]                                                                                         &lt;NA&gt;
##   [2]                                                                                             
##   [3]                 protein phosphatase 1, regulatory subunit 15B [Source:HGNC Symbol;Acc:14951]
##   [4]                           scavenger receptor class B, member 1 [Source:HGNC Symbol;Acc:1664]
##   [5] nuclear casein kinase and cyclin-dependent kinase substrate 1 [Source:HGNC Symbol;Acc:29923]
##   [6]                                         ribosomal protein S25 [Source:HGNC Symbol;Acc:10413]
##         region           type combinedRegionType
##       &lt;factor&gt;       &lt;factor&gt;           &lt;factor&gt;
##   [1]       NA             NA                 NA
##   [2]   Intron      antisense             Intron
##   [3]      CDS protein_coding                CDS
##   [4]      CDS protein_coding                CDS
##   [5]   3' UTR protein_coding             3' UTR
##   [6]      CDS protein_coding                CDS
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>

To export this to a table, simple convert it to a data.frame
<div class="chunk" id="unnamed-chunk-18"><div class="rcode"><div class="source"><pre class="knitr r">head(as.data.frame(annotatedPeaks))
</pre></div>
<div class="output"><pre class="knitr r">##   seqnames     start       end width strand  name    score              ID
## 1        1 153643717 153643723     7      +  5039 9.657284              NA
## 2       12  62997806  62998017   212      + 15618 8.840177 ENSG00000257354
## 3        1 204380047 204380178   132      -  6622 8.663731 ENSG00000158615
## 4       12 125348232 125348415   184      - 17293 8.381262 ENSG00000073060
## 5        1 205686650 205686793   144      -  6738 8.362088 ENSG00000069275
## 6       11 118888683 118888771    89      - 13136 8.354383 ENSG00000118181
##              V2
## 1          &lt;NA&gt;
## 2 RP11-631N16.2
## 3      PPP1R15B
## 4        SCARB1
## 5        NUCKS1
## 6         RPS25
##                                                                                             V3
## 1                                                                                         &lt;NA&gt;
## 2                                                                                             
## 3                 protein phosphatase 1, regulatory subunit 15B [Source:HGNC Symbol;Acc:14951]
## 4                           scavenger receptor class B, member 1 [Source:HGNC Symbol;Acc:1664]
## 5 nuclear casein kinase and cyclin-dependent kinase substrate 1 [Source:HGNC Symbol;Acc:29923]
## 6                                         ribosomal protein S25 [Source:HGNC Symbol;Acc:10413]
##   region           type combinedRegionType
## 1     NA             NA                 NA
## 2 Intron      antisense             Intron
## 3    CDS protein_coding                CDS
## 4    CDS protein_coding                CDS
## 5 3' UTR protein_coding             3' UTR
## 6    CDS protein_coding                CDS
</pre></div>
</div></div>

If you want to continue with annotated peaks you can also load them here

<div class="chunk" id="unnamed-chunk-19"><div class="rcode"><div class="source"><pre class="knitr r">data("annotatedPeaks", package="WhereIsMyPeak")
</pre></div>
</div></div>

## Shuffle

This function shuffles annotated peaks, so that each annotated peak will return an unbound region on the same genomic information. For example, if the peak is on an exon, the shuffled peak will be on an exon of the same gene. If the peak was on an intron it will return a region on the same gene in the intron. It tries to preserve the width of the peak, if there is not enough 'rest - annotation' left, it will return the rest. If the peak was not annotated, it does not shuffle this peak because we do not know on what kind of feature it was.

First lets shuffle only one feature and do not care if the feature is overlapping with other peaks

<div class="chunk" id="unnamed-chunk-20"><div class="rcode"><div class="source"><pre class="knitr r">shuffleRegion(annotatedPeaks[3], txdbParsed)
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 1 range and 0 metadata columns:
##       seqnames                 ranges strand
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt;
##   [1]        1 [204379830, 204379962]      -
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>

If we don't want it to overlap with the same peak, we just put it as background

<div class="chunk" id="unnamed-chunk-21"><div class="rcode"><div class="source"><pre class="knitr r">shuffleRegion(annotatedPeaks[3], txdbParsed, background = annotatedPeaks[3])
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 1 range and 0 metadata columns:
##       seqnames                 ranges strand
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt;
##   [1]        1 [204374729, 204374861]      -
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>

If we do not want the shuffled peak to overlap with any of our peak list then we put in the whole peak list in the background


<div class="chunk" id="unnamed-chunk-22"><div class="rcode"><div class="source"><pre class="knitr r">shuffleRegion(annotatedPeaks[3], txdbParsed, background = annotatedPeaks)
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 1 range and 0 metadata columns:
##       seqnames                 ranges strand
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt;
##   [1]        1 [204374615, 204374747]      -
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>


If we want to apply shuffling as described above on the whole peak list, the command shuffle will do this for a whole GRanges object.

<div class="chunk" id="unnamed-chunk-23"><div class="rcode"><div class="source"><pre class="knitr r">shuffledPeaks <- shuffle(annotatedPeaks, txdbParsed, background = annotatedPeaks)
</pre></div>
</div></div>

You can also parallelize the call using parallel = T

<div class="chunk" id="unnamed-chunk-24"><div class="rcode"><div class="source"><pre class="knitr r">shuffledPeaks <- shuffle(annotatedPeaks, txdbParsed, background = annotatedPeaks, parallel = T)
</pre></div>
</div></div>

Display the shuffled peaks
<div class="chunk" id="unnamed-chunk-25"><div class="rcode"><div class="source"><pre class="knitr r">head(shuffledPeaks)
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 6 ranges and 0 metadata columns:
##       seqnames                 ranges strand
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt;
##   [1]       12 [ 62999022,  62999234]      +
##   [2]        1 [204379291, 204379423]      -
##   [3]       12 [125270006, 125270190]      -
##   [4]        1 [205685156, 205685300]      -
##   [5]       11 [118886555, 118886644]      -
##   [6]       10 [121336124, 121336328]      -
##   -------
##   seqinfo: 690 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>

Let's test if they overlap
<div class="chunk" id="unnamed-chunk-26"><div class="rcode"><div class="source"><pre class="knitr r">sum(countOverlaps(shuffledPeaks, annotatedPeaks))
</pre></div>
<div class="output"><pre class="knitr r">## [1] 0
</pre></div>
</div></div>




You can load a set of shuffled peaks here
<div class="chunk" id="unnamed-chunk-28"><div class="rcode"><div class="source"><pre class="knitr r">data("annotatedPeaks", package="WhereIsMyPeak")
</pre></div>
</div></div>



##  Benchmark
To benchmark this package, do this:

<div class="chunk" id="unnamed-chunk-29"><div class="rcode"><div class="source"><pre class="knitr r">devtools::install_github("hadley/lineprof")
library(lineprof)
x <- lineprof(WhereIsMyPeak(peaks, txdbParsed, annotationOrder, exonAnnotationOrder, typeOrder))
shine(x)
</pre></div>
</div></div>

## Remove Chr

If your peaks are not in the ensembl format, you want to remove the "chr" prefix
and you want to rewrite the mitochondria information.

(Note: this is not perfect and can be certainly improved)

<div class="chunk" id="unnamed-chunk-30"><div class="rcode"><div class="source"><pre class="knitr r">ensemblizePeaks <- function(x) {
      x <- as.data.frame(x)
      x[,1] <- gsub("chr", "", x[,1])
      x[,1][ x[,1] == "M" ] = "MT"
      x <- x[,-4]
      return(makeGRangesFromDataFrame(x, keep.extra.columns = T))
}
</pre></div>
</div></div>


## Session Info
<div class="chunk" id="unnamed-chunk-31"><div class="rcode"><div class="source"><pre class="knitr r">sessionInfo()
</pre></div>
<div class="output"><pre class="knitr r">## R version 3.1.2 (2014-10-31)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] WhereIsMyPeak_0.1      devtools_1.7.0         knitrBootstrap_0.9.0  
##  [4] biomaRt_2.22.0         BiocParallel_1.0.2     GenomicFeatures_1.18.3
##  [7] AnnotationDbi_1.28.1   Biobase_2.26.0         GenomicRanges_1.18.4  
## [10] GenomeInfoDb_1.2.4     IRanges_2.0.1          S4Vectors_0.4.0       
## [13] BiocGenerics_0.12.1   
## 
## loaded via a namespace (and not attached):
##  [1] base64enc_0.1-2         BatchJobs_1.5          
##  [3] BBmisc_1.8              Biostrings_2.34.1      
##  [5] bitops_1.0-6            brew_1.0-6             
##  [7] checkmate_1.5.1         codetools_0.2-10       
##  [9] DBI_0.3.1               digest_0.6.8           
## [11] evaluate_0.5.5          fail_1.2               
## [13] foreach_1.4.2           formatR_1.0            
## [15] GenomicAlignments_1.2.1 highr_0.4              
## [17] httr_0.6.1              iterators_1.0.7        
## [19] knitr_1.9               markdown_0.7.4         
## [21] RCurl_1.95-4.5          Rsamtools_1.18.2       
## [23] RSQLite_1.0.0           rtracklayer_1.26.2     
## [25] sendmailR_1.2-1         stringr_0.6.2          
## [27] tools_3.1.2             XML_3.98-1.1           
## [29] XVector_0.6.0           zlibbioc_1.12.0
</pre></div>
</div></div>

## Compile this page
<div class="chunk" id="unnamed-chunk-32"><div class="rcode"><div class="source"><pre class="knitr r">library(knitrBootstrap)
knit_bootstrap("WhereIsMyPeak.Rmd", show_code=TRUE, boot_style="Readable")
</pre></div>
</div></div>



