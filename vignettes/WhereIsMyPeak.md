---
title: "WhereIsMyPeak"
author: "Thomas Schwarzl"
date: "24 March 2015"
output: html_document
---

# Where is my peak? 

WhereIsMyPeak does genomic annotation with priorisation what gene types you are interested. This package was inspired
by ChIPseeker, however overcomes limited functionlity in the genomic annotation.

Load the package
<div class="chunk" id="unnamed-chunk-1"><div class="rcode"><div class="source"><pre class="knitr r">suppressPackageStartupMessages(library(WhereIsMyPeak))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(biomaRt))
</pre></div>
</div></div>

## Load peaks

Peaks or Intervals are stored as GRanges object. 

Here we load a small example data set
<div class="chunk" id="unnamed-chunk-2"><div class="rcode"><div class="source"><pre class="knitr r">data("peaks", package="WhereIsMyPeak")
peaks
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 99 ranges and 2 metadata columns:
##        seqnames                 ranges strand   |      name       score
##           &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt;   | &lt;integer&gt;   &lt;numeric&gt;
##    [1]        1 [153643717, 153643723]      +   |      5039 9.657284435
##    [2]       12 [ 62997806,  62998017]      +   |     15618 8.840176975
##    [3]        1 [204380047, 204380178]      -   |      6622 8.663731336
##    [4]       12 [125348232, 125348415]      -   |     17293 8.381262242
##    [5]        1 [205686650, 205686793]      -   |      6738 8.362087961
##    ...      ...                    ...    ... ...       ...         ...
##   [95]        1 [ 87572174,  87572278]      +   |      3051  7.41561665
##   [96]        1 [101443010, 101443125]      +   |      3653  7.41561665
##   [97]        1 [197412116, 197412221]      +   |      6450 7.415546471
##   [98]       11 [102268991, 102269101]      -   |     12890  7.41175604
##   [99]        1 [ 77556444,  77556546]      -   |      2695 7.408157717
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>

## Load annotation

Here we load load the annotation using a TxDB Ensembl Build from biomaRt. 

<div class="chunk" id="unnamed-chunk-3"><div class="rcode"><div class="source"><pre class="knitr r"># ENSEMBL GENES 75
txdb <- makeTranscriptDbFromBiomart( biomart         = "ENSEMBL_MART_ENSEMBL",
                                     dataset         = "hsapiens_gene_ensembl",
                                     transcript_ids  = NULL,
                                     filters         = "",
                                     id_prefix       = "ensembl_",
                                     host            = "feb2014.archive.ensembl.org",
                                     port            = 80,
                                     miRBaseBuild    = NA )
</pre></div>
<div class="message"><pre class="knitr r">## Download and preprocess the 'transcripts' data frame ... OK
## Download and preprocess the 'splicings' data frame ... OK
## Download and preprocess the 'genes' data frame ... OK
## Prepare the 'metadata' data frame ... OK
## Make the TxDb object ...
</pre></div>
<div class="warning"><pre class="knitr r">## Warning in .normargChrominfo(chrominfo, transcripts$tx_chrom,
## splicings$exon_chrom): chromosome lengths and circularity flags are not
## available for this TxDb object
</pre></div>
<div class="warning"><pre class="knitr r">## Warning: 'dbBeginTransaction' is deprecated.
## Use 'dbBegin' instead.
## See help(&quot;Deprecated&quot;)
</pre></div>
<div class="warning"><pre class="knitr r">## Warning: 'dbBeginTransaction' is deprecated.
## Use 'dbBegin' instead.
## See help(&quot;Deprecated&quot;)
</pre></div>
<div class="warning"><pre class="knitr r">## Warning: 'dbBeginTransaction' is deprecated.
## Use 'dbBegin' instead.
## See help(&quot;Deprecated&quot;)
</pre></div>
<div class="warning"><pre class="knitr r">## Warning: 'dbBeginTransaction' is deprecated.
## Use 'dbBegin' instead.
## See help(&quot;Deprecated&quot;)
</pre></div>
<div class="warning"><pre class="knitr r">## Warning: 'dbBeginTransaction' is deprecated.
## Use 'dbBegin' instead.
## See help(&quot;Deprecated&quot;)
</pre></div>
<div class="warning"><pre class="knitr r">## Warning: 'dbBeginTransaction' is deprecated.
## Use 'dbBegin' instead.
## See help(&quot;Deprecated&quot;)
</pre></div>
<div class="message"><pre class="knitr r">## OK
</pre></div>
<div class="source"><pre class="knitr r">ensembl75 <- useMart(host="feb2014.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")        
</pre></div>
</div></div>


We load additional information about gene biotypes and gene description which is not available in the TxDB
(Note: This package is currently under development, this could be included in a custom TxDB build in future)
<div class="chunk" id="unnamed-chunk-4"><div class="rcode"><div class="source"><pre class="knitr r">ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl75) 

mappingTable <- getBM(c("ensembl_gene_id", "gene_biotype") , filters="", values=T, uniqueRows = TRUE, mart=ensembl)  
rownames(mappingTable) <- mappingTable[,1]  

namingTable <- getBM(c("ensembl_gene_id", "external_gene_id", "description") , filters="", values=T, uniqueRows = TRUE, mart=ensembl)  
rownames(namingTable) <- namingTable[,1]
namingTable <- namingTable[,-1]
</pre></div>
</div></div>

First we proprocess the txdb since it is not stable, loses connection and serialization is often broken when saved to disk. Also, when annotating multiple batches of peaks (GRanges objects), we do not have to proprocess every single time

<div class="chunk" id="unnamed-chunk-5"><div class="rcode"><div class="source"><pre class="knitr r">txdbParsed <- initRegions(txdb, mappingTable, namingTable)
</pre></div>
</div></div>

## Annoation Priority

Now we choose an annotation priority.

The first position is always prioritized. 

When a peak hits multiple features, there are decisions made:

1) We prefer Exon annotation over Intron annotation
<div class="chunk" id="unnamed-chunk-6"><div class="rcode"><div class="source"><pre class="knitr r">annotationOrder = c("Exon", "Intron")
</pre></div>
</div></div>

2) We decide on Exon annation priority
In Exon annotation, we prefer, CDS over 5', and 5' over 3' annotation. The rest is annotated as Exon. Note that non-coding RNAs for example do not have a CDS, therfore they still remain as Exon. 

(Note: "rest" Exon over everything else will be implemented)  

<div class="chunk" id="unnamed-chunk-7"><div class="rcode"><div class="source"><pre class="knitr r">exonAnnotationOrder = c("CDS", "5' UTR", "3' UTR")
</pre></div>
</div></div>

3) We decide on types
At last, we set types priority.  

<div class="chunk" id="unnamed-chunk-8"><div class="rcode"><div class="source"><pre class="knitr r">typeOrder = c( "non_coding",
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

<div class="chunk" id="unnamed-chunk-9"><div class="rcode"><div class="source"><pre class="knitr r">annotatedPeaks <- WhereIsMyPeak(peaks, txdbParsed, annotationOrder, exonAnnotationOrder, typeOrder)
</pre></div>
</div></div>




You will get back a GRanges object with the annotations stored in the elementMetadata columns
<div class="chunk" id="unnamed-chunk-11"><div class="rcode"><div class="source"><pre class="knitr r">head(annotatedPeaks) 
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
<div class="chunk" id="unnamed-chunk-12"><div class="rcode"><div class="source"><pre class="knitr r">as.data.frame(annotatedPeaks)
</pre></div>
<div class="output"><pre class="knitr r">##    seqnames     start       end width strand  name    score
## 1         1 153643717 153643723     7      +  5039 9.657284
## 2        12  62997806  62998017   212      + 15618 8.840177
## 3         1 204380047 204380178   132      -  6622 8.663731
## 4        12 125348232 125348415   184      - 17293 8.381262
## 5         1 205686650 205686793   144      -  6738 8.362088
## 6        11 118888683 118888771    89      - 13136 8.354383
## 7        10 121356141 121356344   204      -  9878 8.342655
## 8        10 129903023 129903148   126      - 10026 8.203284
## 9        11  93468208  93468245    38      - 12717 8.190256
## 10        1  54293720  54293861   142      -  1826 8.155124
## 11        1 207968203 207968429   227      +  6851 8.132756
## 12        1  45987494  45987558    65      -  1563 8.132713
## 13        1 155196923 155197065   143      -  5343 8.100723
## 14        1 207498972 207499055    84      +  6822 8.096089
## 15        1 110091260 110091430   171      +  3936 8.086775
## 16        1 150208339 150208415    77      -  4794 8.072690
## 17       12  88887255  88887368   114      - 16066 8.048950
## 18        1 100503275 100503575   301      -  3556 8.048607
## 19       10  98287809  98287930   122      -  9152 8.034489
## 20        1  87810606  87810722   117      +  3085 8.029682
## 21       11  85374762  85374838    77      - 12605 8.024767
## 22        1  70705122  70705174    53      +  2491 7.980237
## 23       10 102106835 102106868    34      +  9272 7.970189
## 24       12 121746284 121746409   126      - 16981 7.965067
## 25       12 104326055 104326147    93      + 16344 7.959764
## 26       12  51325164  51325204    41      + 14894 7.954913
## 27        1 149858465 149858524    60      +  4755 7.949784
## 28       11  47488400  47488480    81      - 11163 7.913415
## 29       11  34112065  34112152    88      + 10951 7.908137
## 30        1  93620042  93620134    93      -  3355 7.892185
## 31       10 120817615 120817713    99      -  9828 7.892185
## 32       11   8008881   8008954    74      + 10336 7.886879
## 33        1 109506005 109506103    99      -  3837 7.881503
## 34       11  47491174  47491250    77      - 11175 7.876106
## 35       12 125397135 125397171    37      - 17296 7.876106
## 36        1 115277099 115277147    49      -  4266 7.859794
## 37        1 109235279 109235397   119      +  3776 7.854315
## 38        1  71872398  71872507   110      -  2581 7.848816
## 39        1 117530233 117530345   113      +  4367 7.843347
## 40       12 107366147 107366234    88      + 16468 7.826659
## 41        1 207940395 207940509   115      +  6837 7.821053
## 42       12 132413806 132413917   112      + 17405 7.815425
## 43        1  70695612  70695684    73      +  2475 7.798408
## 44        1  70687462  70687679   218      +  2456 7.793130
## 45        1 169080613 169080711    99      +  5930 7.775566
## 46        1  92147722  92147793    72      -  3275 7.752134
## 47       11  61568469  61568545    77      - 11586 7.740302
## 48        1 154556485 154556639   155      -  5202 7.734349
## 49        1 205682964 205683066   103      -  6714 7.722370
## 50        1 100488820 100489039   220      +  3530 7.716684
## 51       13  41634888  41634914    27      - 17720 7.716456
## 52       11  10819341  10819419    79      - 10482 7.716343
## 53        1 155277895 155277956    62      -  5377 7.710290
## 54        1 167758119 167758131    13      +  5862 7.710290
## 55        1 150238824 150238942   119      -  4801 7.692094
## 56       11  18429003  18429100    98      + 10735 7.667255
## 57       12 125411964 125411972     9      - 17308 7.651687
## 58        1  77556887  77557021   135      -  2696 7.648445
## 59        1  76253985  76254023    39      +  2660 7.648385
## 60       12  56705069  56705166    98      - 15440 7.642040
## 61        1  76253316  76253383    68      +  2658 7.622896
## 62       11  93464334  93464379    46      - 12707 7.609951
## 63        1  87573461  87573610   150      +  3059 7.603495
## 64       10  17277286  17277327    42      +  7914 7.603434
## 65       12  12715128  12715298   171      - 13910 7.592776
## 66        1 115250181 115250220    40      -  4248 7.590311
## 67        1  76256114  76256150    37      +  2671 7.577131
## 68       10 129904025 129904107    83      - 10030 7.577131
## 69       12  56122774  56122799    26      - 15307 7.563766
## 70        1 154570337 154570425    89      -  5218 7.557036
## 71       11 111896223 111896316    94      + 12989 7.557036
## 72       11  61183644  61183820   177      - 11551 7.550339
## 73       12  69652676  69652796   121      + 15829 7.543547
## 74       11  93465578  93465625    48      - 12709 7.529327
## 75       10 119042842 119043044   203      -  9752 7.528499
## 76       10  69524254  69524327    74      +  8402 7.527522
## 77        1 179813893 179813983    91      -  6194 7.516053
## 78        1  89299080  89299212   133      +  3119 7.509097
## 79       13  80126691  80126801   111      + 17898 7.509097
## 80       13 113539234 113539339   106      + 18192 7.495083
## 81       12  20836675  20836746    72      + 14077 7.488092
## 82       12  20522431  20522540   110      + 14034 7.480999
## 83        1 214803883 214803963    81      +  6989 7.473872
## 84       10 129907388 129907488   101      - 10050 7.459510
## 85       11 122931887 122931924    38      - 13271 7.452275
## 86        1 169100530 169100646   117      +  5932 7.449044
## 87        1  76252915  76252949    35      +  2656 7.445073
## 88        1   1326899   1326967    69      -    70 7.430420
## 89        1  70716796  70716870    75      +  2497 7.430350
## 90        1  87572712  87572785    74      +  3054 7.430350
## 91        1 155279493 155279566    74      +  5388 7.430350
## 92       10  79797156  79797214    59      +  8780 7.423037
## 93        1  28835696  28835728    33      +   777 7.420126
## 94        1  62916787  62916905   119      +  2109 7.419012
## 95        1  87572174  87572278   105      +  3051 7.415617
## 96        1 101443010 101443125   116      +  3653 7.415617
## 97        1 197412116 197412221   106      +  6450 7.415546
## 98       11 102268991 102269101   111      - 12890 7.411756
## 99        1  77556444  77556546   103      -  2695 7.408158
##                 ID             V2
## 1               NA           &lt;NA&gt;
## 2  ENSG00000257354  RP11-631N16.2
## 3  ENSG00000158615       PPP1R15B
## 4  ENSG00000073060         SCARB1
## 5  ENSG00000069275         NUCKS1
## 6  ENSG00000118181          RPS25
## 7  ENSG00000151923          TIAL1
## 8  ENSG00000148773          MKI67
## 9  ENSG00000166012          TAF1D
## 10 ENSG00000058804           NDC1
## 11 ENSG00000117335           CD46
## 12 ENSG00000117450          PRDX1
## 13 ENSG00000160766          GBAP1
## 14 ENSG00000196352           CD55
## 15 ENSG00000065135          GNAI3
## 16 ENSG00000143401         ANP32E
## 17 ENSG00000049130          KITLG
## 18              NA           &lt;NA&gt;
## 19 ENSG00000077147         TM9SF3
## 20 ENSG00000143013           LMO4
## 21 ENSG00000137504         CREBZF
## 22 ENSG00000116754         SRSF11
## 23              NA           &lt;NA&gt;
## 24 ENSG00000089053         ANAPC5
## 25 ENSG00000166598        HSP90B1
## 26 ENSG00000185432        METTL7A
## 27              NA           &lt;NA&gt;
## 28 ENSG00000149187          CELF1
## 29 ENSG00000135387        CAPRIN1
## 30 ENSG00000117500          TMED5
## 31 ENSG00000107581          EIF3A
## 32 ENSG00000175390          EIF3F
## 33 ENSG00000162641         AKNAD1
## 34 ENSG00000149187          CELF1
## 35 ENSG00000150991            UBC
## 36 ENSG00000009307          CSDE1
## 37 ENSG00000134186        PRPF38B
## 38 ENSG00000172260          NEGR1
## 39 ENSG00000134247         PTGFRN
## 40 ENSG00000151135       C12orf23
## 41 ENSG00000117335           CD46
## 42 ENSG00000177192           PUS1
## 43 ENSG00000116754         SRSF11
## 44 ENSG00000116754         SRSF11
## 45 ENSG00000143153         ATP1B1
## 46 ENSG00000069702         TGFBR3
## 47 ENSG00000149485          FADS1
## 48 ENSG00000160710           ADAR
## 49 ENSG00000069275         NUCKS1
## 50 ENSG00000117620        SLC35A3
## 51 ENSG00000120690           ELF1
## 52 ENSG00000110321         EIF4G2
## 53              NA           &lt;NA&gt;
## 54 ENSG00000197965          MPZL1
## 55 ENSG00000117362          APH1A
## 56 ENSG00000134333           LDHA
## 57              NA           &lt;NA&gt;
## 58 ENSG00000142892           PIGK
## 59 ENSG00000137955        RABGGTB
## 60 ENSG00000144785 RP11-977G19.10
## 61 ENSG00000137955        RABGGTB
## 62 ENSG00000166012          TAF1D
## 63 ENSG00000153936         HS2ST1
## 64 ENSG00000026025            VIM
## 65 ENSG00000111266         DUSP16
## 66 ENSG00000213281           NRAS
## 67 ENSG00000137955        RABGGTB
## 68 ENSG00000148773          MKI67
## 69 ENSG00000135404           CD63
## 70 ENSG00000160710           ADAR
## 71 ENSG00000150768           DLAT
## 72 ENSG00000149532          CPSF7
## 73 ENSG00000111605          CPSF6
## 74 ENSG00000207304         SNORA8
## 75 ENSG00000165650          PDZD8
## 76              NA           &lt;NA&gt;
## 77 ENSG00000169905       TOR1AIP2
## 78 ENSG00000065243           PKN2
## 79 ENSG00000102471         NDFIP2
## 80 ENSG00000068650         ATP11A
## 81 ENSG00000172572          PDE3A
## 82 ENSG00000172572          PDE3A
## 83 ENSG00000117724          CENPF
## 84 ENSG00000148773          MKI67
## 85 ENSG00000109971          HSPA8
## 86 ENSG00000143153         ATP1B1
## 87 ENSG00000137955        RABGGTB
## 88 ENSG00000221978          CCNL2
## 89 ENSG00000116754         SRSF11
## 90 ENSG00000153936         HS2ST1
## 91 ENSG00000160752           FDPS
## 92 ENSG00000138326          RPS24
## 93 ENSG00000242125          SNHG3
## 94 ENSG00000162607           USP1
## 95 ENSG00000153936         HS2ST1
## 96 ENSG00000162695        SLC30A7
## 97 ENSG00000134376           CRB1
## 98 ENSG00000152558        TMEM123
## 99 ENSG00000142892           PIGK
##                                                                                                                       V3
## 1                                                                                                                   &lt;NA&gt;
## 2                                                                                                                       
## 3                                           protein phosphatase 1, regulatory subunit 15B [Source:HGNC Symbol;Acc:14951]
## 4                                                     scavenger receptor class B, member 1 [Source:HGNC Symbol;Acc:1664]
## 5                           nuclear casein kinase and cyclin-dependent kinase substrate 1 [Source:HGNC Symbol;Acc:29923]
## 6                                                                   ribosomal protein S25 [Source:HGNC Symbol;Acc:10413]
## 7                            TIA1 cytotoxic granule-associated RNA binding protein-like 1 [Source:HGNC Symbol;Acc:11804]
## 8                                                            marker of proliferation Ki-67 [Source:HGNC Symbol;Acc:7107]
## 9            TATA box binding protein (TBP)-associated factor, RNA polymerase I, D, 41kDa [Source:HGNC Symbol;Acc:28759]
## 10                                                         NDC1 transmembrane nucleoporin [Source:HGNC Symbol;Acc:25525]
## 11                                            CD46 molecule, complement regulatory protein [Source:HGNC Symbol;Acc:6953]
## 12                                                                         peroxiredoxin 1 [Source:HGNC Symbol;Acc:9352]
## 13                                                    glucosidase, beta, acid pseudogene 1 [Source:HGNC Symbol;Acc:4178]
## 14            CD55 molecule, decay accelerating factor for complement (Cromer blood group) [Source:HGNC Symbol;Acc:2665]
## 15 guanine nucleotide binding protein (G protein), alpha inhibiting activity polypeptide 3 [Source:HGNC Symbol;Acc:4387]
## 16                       acidic (leucine-rich) nuclear phosphoprotein 32 family, member E [Source:HGNC Symbol;Acc:16673]
## 17                                                                              KIT ligand [Source:HGNC Symbol;Acc:6343]
## 18                                                                                                                  &lt;NA&gt;
## 19                                                   transmembrane 9 superfamily member 3 [Source:HGNC Symbol;Acc:21529]
## 20                                                                       LIM domain only 4 [Source:HGNC Symbol;Acc:6644]
## 21                                                     CREB/ATF bZIP transcription factor [Source:HGNC Symbol;Acc:24905]
## 22                                                serine/arginine-rich splicing factor 11 [Source:HGNC Symbol;Acc:10782]
## 23                                                                                                                  &lt;NA&gt;
## 24                                                   anaphase promoting complex subunit 5 [Source:HGNC Symbol;Acc:15713]
## 25                                        heat shock protein 90kDa beta (Grp94), member 1 [Source:HGNC Symbol;Acc:12028]
## 26                                                              methyltransferase like 7A [Source:HGNC Symbol;Acc:24550]
## 27                                                                                                                  &lt;NA&gt;
## 28                                                        CUGBP, Elav-like family member 1 [Source:HGNC Symbol;Acc:2549]
## 29                                                         cell cycle associated protein 1 [Source:HGNC Symbol;Acc:6743]
## 30                              transmembrane emp24 protein transport domain containing 5 [Source:HGNC Symbol;Acc:24251]
## 31                                   eukaryotic translation initiation factor 3, subunit A [Source:HGNC Symbol;Acc:3271]
## 32                                   eukaryotic translation initiation factor 3, subunit F [Source:HGNC Symbol;Acc:3275]
## 33                                                               AKNA domain containing 1 [Source:HGNC Symbol;Acc:28398]
## 34                                                        CUGBP, Elav-like family member 1 [Source:HGNC Symbol;Acc:2549]
## 35                                                                            ubiquitin C [Source:HGNC Symbol;Acc:12468]
## 36                                           cold shock domain containing E1, RNA-binding [Source:HGNC Symbol;Acc:29905]
## 37                                                         pre-mRNA processing factor 38B [Source:HGNC Symbol;Acc:25512]
## 38                                                            neuronal growth regulator 1 [Source:HGNC Symbol;Acc:17302]
## 39                                                     prostaglandin F2 receptor inhibitor [Source:HGNC Symbol;Acc:9601]
## 40                                                    chromosome 12 open reading frame 23 [Source:HGNC Symbol;Acc:28281]
## 41                                            CD46 molecule, complement regulatory protein [Source:HGNC Symbol;Acc:6953]
## 42                                                             pseudouridylate synthase 1 [Source:HGNC Symbol;Acc:15508]
## 43                                                serine/arginine-rich splicing factor 11 [Source:HGNC Symbol;Acc:10782]
## 44                                                serine/arginine-rich splicing factor 11 [Source:HGNC Symbol;Acc:10782]
## 45                                          ATPase, Na+/K+ transporting, beta 1 polypeptide [Source:HGNC Symbol;Acc:804]
## 46                                          transforming growth factor, beta receptor III [Source:HGNC Symbol;Acc:11774]
## 47                                                                 fatty acid desaturase 1 [Source:HGNC Symbol;Acc:3574]
## 48                                                        adenosine deaminase, RNA-specific [Source:HGNC Symbol;Acc:225]
## 49                          nuclear casein kinase and cyclin-dependent kinase substrate 1 [Source:HGNC Symbol;Acc:29923]
## 50 solute carrier family 35 (UDP-N-acetylglucosamine (UDP-GlcNAc) transporter), member A3 [Source:HGNC Symbol;Acc:11023]
## 51                                     E74-like factor 1 (ets domain transcription factor) [Source:HGNC Symbol;Acc:3316]
## 52                                     eukaryotic translation initiation factor 4 gamma, 2 [Source:HGNC Symbol;Acc:3297]
## 53                                                                                                                  &lt;NA&gt;
## 54                                                              myelin protein zero-like 1 [Source:HGNC Symbol;Acc:7226]
## 55                                                          APH1A gamma secretase subunit [Source:HGNC Symbol;Acc:29509]
## 56                                                                 lactate dehydrogenase A [Source:HGNC Symbol;Acc:6535]
## 57                                                                                                                  &lt;NA&gt;
## 58                                phosphatidylinositol glycan anchor biosynthesis, class K [Source:HGNC Symbol;Acc:8965]
## 59                                             Rab geranylgeranyltransferase, beta subunit [Source:HGNC Symbol;Acc:9796]
## 60                                                         Uncharacterized protein  [Source:UniProtKB/TrEMBL;Acc:F8VP03]
## 61                                             Rab geranylgeranyltransferase, beta subunit [Source:HGNC Symbol;Acc:9796]
## 62           TATA box binding protein (TBP)-associated factor, RNA polymerase I, D, 41kDa [Source:HGNC Symbol;Acc:28759]
## 63                                                  heparan sulfate 2-O-sulfotransferase 1 [Source:HGNC Symbol;Acc:5193]
## 64                                                                               vimentin [Source:HGNC Symbol;Acc:12692]
## 65                                                        dual specificity phosphatase 16 [Source:HGNC Symbol;Acc:17909]
## 66                                        neuroblastoma RAS viral (v-ras) oncogene homolog [Source:HGNC Symbol;Acc:7989]
## 67                                             Rab geranylgeranyltransferase, beta subunit [Source:HGNC Symbol;Acc:9796]
## 68                                                           marker of proliferation Ki-67 [Source:HGNC Symbol;Acc:7107]
## 69                                                                           CD63 molecule [Source:HGNC Symbol;Acc:1692]
## 70                                                        adenosine deaminase, RNA-specific [Source:HGNC Symbol;Acc:225]
## 71                                                    dihydrolipoamide S-acetyltransferase [Source:HGNC Symbol;Acc:2896]
## 72                                  cleavage and polyadenylation specific factor 7, 59kDa [Source:HGNC Symbol;Acc:30098]
## 73                                  cleavage and polyadenylation specific factor 6, 68kDa [Source:HGNC Symbol;Acc:13871]
## 74                                                       small nucleolar RNA, H/ACA box 8 [Source:HGNC Symbol;Acc:32596]
## 75                                                                PDZ domain containing 8 [Source:HGNC Symbol;Acc:26974]
## 76                                                                                                                  &lt;NA&gt;
## 77                                                         torsin A interacting protein 2 [Source:HGNC Symbol;Acc:24055]
## 78                                                                       protein kinase N2 [Source:HGNC Symbol;Acc:9406]
## 79                                                     Nedd4 family interacting protein 2 [Source:HGNC Symbol;Acc:18537]
## 80                                                             ATPase, class VI, type 11A [Source:HGNC Symbol;Acc:13552]
## 81                                                    phosphodiesterase 3A, cGMP-inhibited [Source:HGNC Symbol;Acc:8778]
## 82                                                    phosphodiesterase 3A, cGMP-inhibited [Source:HGNC Symbol;Acc:8778]
## 83                                                        centromere protein F, 350/400kDa [Source:HGNC Symbol;Acc:1857]
## 84                                                           marker of proliferation Ki-67 [Source:HGNC Symbol;Acc:7107]
## 85                                                              heat shock 70kDa protein 8 [Source:HGNC Symbol;Acc:5241]
## 86                                          ATPase, Na+/K+ transporting, beta 1 polypeptide [Source:HGNC Symbol;Acc:804]
## 87                                             Rab geranylgeranyltransferase, beta subunit [Source:HGNC Symbol;Acc:9796]
## 88                                                                              cyclin L2 [Source:HGNC Symbol;Acc:20570]
## 89                                                serine/arginine-rich splicing factor 11 [Source:HGNC Symbol;Acc:10782]
## 90                                                  heparan sulfate 2-O-sulfotransferase 1 [Source:HGNC Symbol;Acc:5193]
## 91                                                           farnesyl diphosphate synthase [Source:HGNC Symbol;Acc:3631]
## 92                                                                  ribosomal protein S24 [Source:HGNC Symbol;Acc:10411]
## 93                                   small nucleolar RNA host gene 3 (non-protein coding) [Source:HGNC Symbol;Acc:10118]
## 94                                                         ubiquitin specific peptidase 1 [Source:HGNC Symbol;Acc:12607]
## 95                                                  heparan sulfate 2-O-sulfotransferase 1 [Source:HGNC Symbol;Acc:5193]
## 96                                  solute carrier family 30 (zinc transporter), member 7 [Source:HGNC Symbol;Acc:19306]
## 97                                                           crumbs homolog 1 (Drosophila) [Source:HGNC Symbol;Acc:2343]
## 98                                                              transmembrane protein 123 [Source:HGNC Symbol;Acc:30138]
## 99                                phosphatidylinositol glycan anchor biosynthesis, class K [Source:HGNC Symbol;Acc:8965]
##    region           type combinedRegionType
## 1      NA             NA                 NA
## 2  Intron      antisense             Intron
## 3     CDS protein_coding                CDS
## 4     CDS protein_coding                CDS
## 5  3' UTR protein_coding             3' UTR
## 6     CDS protein_coding                CDS
## 7  5' UTR protein_coding             5' UTR
## 8     CDS protein_coding                CDS
## 9    Exon protein_coding               Exon
## 10    CDS protein_coding                CDS
## 11 3' UTR protein_coding             3' UTR
## 12 5' UTR protein_coding             5' UTR
## 13 Intron     pseudogene             Intron
## 14    CDS protein_coding                CDS
## 15    CDS protein_coding                CDS
## 16 5' UTR protein_coding             5' UTR
## 17   Exon protein_coding               Exon
## 18     NA             NA                 NA
## 19    CDS protein_coding                CDS
## 20 3' UTR protein_coding             3' UTR
## 21 3' UTR protein_coding             3' UTR
## 22    CDS protein_coding                CDS
## 23     NA             NA                 NA
## 24    CDS protein_coding                CDS
## 25    CDS protein_coding                CDS
## 26 3' UTR protein_coding             3' UTR
## 27     NA             NA                 NA
## 28 3' UTR protein_coding             3' UTR
## 29    CDS protein_coding                CDS
## 30 3' UTR protein_coding             3' UTR
## 31    CDS protein_coding                CDS
## 32    CDS protein_coding                CDS
## 33    CDS protein_coding                CDS
## 34 3' UTR protein_coding             3' UTR
## 35    CDS protein_coding                CDS
## 36    CDS protein_coding                CDS
## 37    CDS protein_coding                CDS
## 38 3' UTR protein_coding             3' UTR
## 39 3' UTR protein_coding             3' UTR
## 40 3' UTR protein_coding             3' UTR
## 41    CDS protein_coding                CDS
## 42 5' UTR protein_coding             5' UTR
## 43 Intron protein_coding             Intron
## 44    CDS protein_coding                CDS
## 45    CDS protein_coding                CDS
## 46 3' UTR protein_coding             3' UTR
## 47 3' UTR protein_coding             3' UTR
## 48 3' UTR protein_coding             3' UTR
## 49 3' UTR protein_coding             3' UTR
## 50 3' UTR protein_coding             3' UTR
## 51 Intron protein_coding             Intron
## 52    CDS protein_coding                CDS
## 53     NA             NA                 NA
## 54 3' UTR protein_coding             3' UTR
## 55    CDS protein_coding                CDS
## 56 3' UTR protein_coding             3' UTR
## 57     NA             NA                 NA
## 58 3' UTR protein_coding             3' UTR
## 59   Exon protein_coding               Exon
## 60    CDS protein_coding                CDS
## 61   Exon protein_coding               Exon
## 62 3' UTR protein_coding             3' UTR
## 63 3' UTR protein_coding             3' UTR
## 64    CDS protein_coding                CDS
## 65 5' UTR protein_coding             5' UTR
## 66 3' UTR protein_coding             3' UTR
## 67 5' UTR protein_coding             5' UTR
## 68    CDS protein_coding                CDS
## 69 5' UTR protein_coding             5' UTR
## 70    CDS protein_coding                CDS
## 71    CDS protein_coding                CDS
## 72    CDS protein_coding                CDS
## 73    CDS protein_coding                CDS
## 74   Exon         snoRNA               Exon
## 75    CDS protein_coding                CDS
## 76     NA             NA                 NA
## 77 3' UTR protein_coding             3' UTR
## 78    CDS protein_coding                CDS
## 79 3' UTR protein_coding             3' UTR
## 80 3' UTR protein_coding             3' UTR
## 81 3' UTR protein_coding             3' UTR
## 82    CDS protein_coding                CDS
## 83    CDS protein_coding                CDS
## 84    CDS protein_coding                CDS
## 85    CDS protein_coding                CDS
## 86    CDS protein_coding                CDS
## 87   Exon protein_coding               Exon
## 88    CDS protein_coding                CDS
## 89 3' UTR protein_coding             3' UTR
## 90 3' UTR protein_coding             3' UTR
## 91   Exon protein_coding               Exon
## 92   Exon protein_coding               Exon
## 93   Exon sense_intronic               Exon
## 94 3' UTR protein_coding             3' UTR
## 95 3' UTR protein_coding             3' UTR
## 96 3' UTR protein_coding             3' UTR
## 97 3' UTR protein_coding             3' UTR
## 98 3' UTR protein_coding             3' UTR
## 99 3' UTR protein_coding             3' UTR
</pre></div>
</div></div>

If you want to continue with annotated peaks you can also load them here

<div class="chunk" id="unnamed-chunk-13"><div class="rcode"><div class="source"><pre class="knitr r">data("annotatedPeaks", package="WhereIsMyPeak")
</pre></div>
</div></div>

## Shuffle

This function shuffles annotated peaks, so that each annotated peak will return an unbound region on the same genomic information. For example, if the peak is on an exon, the shuffled peak will be on an exon of the same gene. If the peak was on an intron it will return a region on the same gene in the intron. It tries to preserve the width of the peak, if there is not enough 'rest - annotation' left, it will return the rest. If the peak was not annotated, it does not shuffle this peak because we do not know on what kind of feature it was.

First lets shuffle only one feature and do not care if the feature is overlapping with other peaks

<div class="chunk" id="unnamed-chunk-14"><div class="rcode"><div class="source"><pre class="knitr r">shuffleRegion(annotatedPeaks[3], txdbParsed)
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 1 range and 0 metadata columns:
##       seqnames                 ranges strand
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt;
##   [1]        1 [204372745, 204372877]      -
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>

If we don't want it to overlap with the same peak, we just put it as background

<div class="chunk" id="unnamed-chunk-15"><div class="rcode"><div class="source"><pre class="knitr r">shuffleRegion(annotatedPeaks[3], txdbParsed, background = annotatedPeaks[3])
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 1 range and 0 metadata columns:
##       seqnames                 ranges strand
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt;
##   [1]        1 [204375016, 204375148]      -
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>

If we do not want the shuffled peak to overlap with any of our peak list then we put in the whole peak list in the background


<div class="chunk" id="unnamed-chunk-16"><div class="rcode"><div class="source"><pre class="knitr r">shuffleRegion(annotatedPeaks[3], txdbParsed, background = annotatedPeaks)
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 1 range and 0 metadata columns:
##       seqnames                 ranges strand
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt;
##   [1]        1 [204379292, 204379424]      -
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>


If we want to apply shuffling as described above on the whole peak list, the command shuffle will do this for a whole GRanges object.

<div class="chunk" id="unnamed-chunk-17"><div class="rcode"><div class="source"><pre class="knitr r">shuffledPeaks <- shuffle(annotatedPeaks, txdbParsed, background = annotatedPeaks)
</pre></div>
</div></div>

You can also parallelize the call using parallel = T

<div class="chunk" id="unnamed-chunk-18"><div class="rcode"><div class="source"><pre class="knitr r">shuffledPeaks <- shuffle(annotatedPeaks, txdbParsed, background = annotatedPeaks, parallel = T)
</pre></div>
</div></div>

Display the shuffled peaks
<div class="chunk" id="unnamed-chunk-19"><div class="rcode"><div class="source"><pre class="knitr r">head(shuffledPeaks)
</pre></div>
<div class="output"><pre class="knitr r">## GRanges object with 6 ranges and 0 metadata columns:
##       seqnames                 ranges strand
##          &lt;Rle&gt;              &lt;IRanges&gt;  &lt;Rle&gt;
##   [1]       12 [ 63001914,  63002126]      +
##   [2]        1 [204372594, 204372726]      -
##   [3]       12 [125270536, 125270720]      -
##   [4]        1 [205682592, 205682736]      -
##   [5]       11 [118886422, 118886468]      -
##   [6]       10 [121336284, 121336488]      -
##   -------
##   seqinfo: 690 sequences from an unspecified genome; no seqlengths
</pre></div>
</div></div>

Let's test if they overlap
<div class="chunk" id="unnamed-chunk-20"><div class="rcode"><div class="source"><pre class="knitr r">sum(countOverlaps(shuffledPeaks, annotatedPeaks))
</pre></div>
<div class="output"><pre class="knitr r">## [1] 0
</pre></div>
</div></div>




You can load a set of shuffled peaks here
<div class="chunk" id="unnamed-chunk-22"><div class="rcode"><div class="source"><pre class="knitr r">data("annotatedPeaks", package="WhereIsMyPeak")
</pre></div>
</div></div>



##  Benchmark
To benchmark this package, do this:

<div class="chunk" id="unnamed-chunk-23"><div class="rcode"><div class="source"><pre class="knitr r">devtools::install_github("hadley/lineprof")
library(lineprof)
x <- lineprof(WhereIsMyPeak(peaks, txdbParsed, annotationOrder, exonAnnotationOrder, typeOrder))
shine(x)
</pre></div>
</div></div>

## Remove Chr

If your peaks are not in the ensembl format, you want to remove the "chr" prefix
and you want to rewrite the mitochondria information.

<div class="chunk" id="unnamed-chunk-24"><div class="rcode"><div class="source"><pre class="knitr r">ensemblizePeaks <- function(x) {
      x <- as.data.frame(x)
      x[,1] <- gsub("chr", "", x[,1])
      x[,1][ x[,1] == "M" ] = "MT"
      x <- x[,-4]
      return(makeGRangesFromDataFrame(x, keep.extra.columns = T))
}
</pre></div>
</div></div>


## Session Info
<div class="chunk" id="unnamed-chunk-25"><div class="rcode"><div class="source"><pre class="knitr r">sessionInfo()
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
##  [1] knitrBootstrap_0.9.0   WhereIsMyPeak_0.1      biomaRt_2.22.0        
##  [4] BiocParallel_1.0.2     GenomicFeatures_1.18.3 AnnotationDbi_1.28.1  
##  [7] Biobase_2.26.0         GenomicRanges_1.18.4   GenomeInfoDb_1.2.4    
## [10] IRanges_2.0.1          S4Vectors_0.4.0        BiocGenerics_0.12.1   
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
## [17] iterators_1.0.7         knitr_1.9              
## [19] markdown_0.7.4          RCurl_1.95-4.5         
## [21] Rsamtools_1.18.2        RSQLite_1.0.0          
## [23] rtracklayer_1.26.2      sendmailR_1.2-1        
## [25] stringr_0.6.2           tools_3.1.2            
## [27] XML_3.98-1.1            XVector_0.6.0          
## [29] zlibbioc_1.12.0
</pre></div>
</div></div>



