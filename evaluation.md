% Reproduction of "Human Gut Microbiome Viewed Across Age and Geography"
% Hermann Pauly
% \today

Abstract
======

Soon to come...

Introduction
========

In science we try to produce reliable and unambigous results, working hard to avoid
common pitfalls and sloppy work, and we expect the same from our fellow scientists.
Nevertheless, an unacceptable number of publications *(numbers)* contain grave errors
*(quote needed)*, which has lead to waning trust in science, among scientists as well
as in the general public. One feature of trustworthy scientific work is
reproducability. If a second person with proper knowledge analyses a researchers' data
and finds the same results, this is an indicator that the first scientist did work
correctly. *(or both are wrong ;-) )*  
Here I work on the data Yasunenko et al. published in 2012 along with their publication
"Human Gut Microbiome Viewed Across Age and Geography" in Nature magazine.
They sequenced genomic data of fecal samples from participants hailing from three
different countries. Using ANOVA post hoc testing they found significant diversities
among ages, geographic locations and families (see [@Yats12]).
I recreated calculations, data processing, and data visualisation following the
descriptions given in the original paper to check if their findings and conclusions
can be reproduced knowing nothing more than the published sequencer data.

Methods and results
==============

<!--
## System information


```r
require(cluster)
```

```
## Loading required package: cluster
```

```r
sessionInfo()
```

```
## R version 3.1.0 (2014-04-10)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## locale:
## [1] de_DE.UTF-8/de_DE.UTF-8/de_DE.UTF-8/C/de_DE.UTF-8/C
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] cluster_1.15.2 knitr_1.6     
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.5 formatR_0.10   stringr_0.6.2  tools_3.1.0
```
-->

## Data aquisition ##

I downloaded the raw data and metadata files from MG-RAST[@MG-RAST] using the
UNIX command line tool *wget*. They can be found under project ID
numbers 98 (whole genome shotgun sequences) and 401 (16S rRNA V4 sequences).

## Preprocessing ##

### Metadata mapping files ###

To categorize and compare genome data, a link between sequenced sample
files and their respective metadata is required. The MG-RAST database [@MG-RAST]
provides a
mapping file between sample data files and unique IDs with incomplete
metadata, while Yatsunenko et al. [@Yats12] provide a complete metadata file without
mapping to the sample data files. I applied the "calc" module of the
*LibreOffice* suite to combine the contents of both files to a
complete mapping file. For this I sorted both provided files' contents
by their sample ID strings and copy-pasted missing rows from the
MG-RAST mapping file into the metadata table. Two entries from the
mapping file were not found in the metadata table, and Yatsunenko et al. also
stated that two samples could not be used in the analysis, so I
removed them from the metadata-mapping file. The result was saved as a
tab-separated .csv file. On my system with German environment settings
I had to convert decimal values from comma-separated values to
international dot-separated format using the command line tool *sed*:

```
sed 's/,/\./g' 16s_mapping.csv > 16s_mapping_decimaldot.csv
```

The unique sequencer IDs differed from the sample IDs only by a
numerical appendix, so I used a Python script to check if all samples
were matched correctly:

TODO: move out to .py file

~~~python
def checkSampleIDmapping(fileName):
    print "checking", fileName, "-->",
    csv = open(fileName, "r")
    count = 0
    errors = 0
    csv.readline() # skip header

    for line in csv:
        line_split = line.split('\t')
        try:
            longID = line_split[0]
            shortID = line_split[1]
            if not longID.startswith(shortID):
                print "error in line", count, shortID, longID
                errors += 1
        except:
            pass
        count += 1

    print "read", count, "samples,", errors, "errors"

for fileName in ["16s_mapping_decimaldot.csv", "wgs_mapping_decimaldot.csv"]:
    checkSampleIDmapping(fileName)
~~~

Calling the script produced the following output:

    checking 16s_mapping_decimaldot.csv --> read 528 samples, 0 errors
	checking wgs_mapping_decimaldot.csv --> read 110 samples, 0 errors  

### $\alpha$- and $\beta$ diversity ###

Using the supported data and the manually created metadata mapping
table, I followed the QIIME workflow to create tables of OTU
abundances and analyse them for $\alpha$- and $\beta$ diversity.

The first step required is to categorize the sequenced microorganisms
into groups by phylogenetic distance. This step is called "picking". To
do this, I used the QIIME tool *pick\_closed\_reference.py*, which
compares the samples to an existing similarity tree. As picking
reference I used the GreenGenes database from 2011-02-04 and a
similarity threshold of 97%, as did Yatsunenko et al. [@Yats12].
The picking was done on
each sample individually, ten samples in parallel at any time, by
using the custom Python script *picking.py*.

The QIIME tools for $\alpha$-rarefaction and
$\beta$-diversity-analysis require a single biom table as input file,
so I combined the picking results. A try to combine all 528 sample
tables exceeded the computational power available, so I used the
custom Python script *combine.py* which creates a shell script that
calls the QIIME tool *merge\_otu\_tables.py* to iteratively add
sample tables to a common table. A copy of this table was converted to
a tab-separated file (\texttt{combined.csv}) for easy use with *R*.

The following steps were all done with QIIME tools. To assess species
richness from the samples, the "rarefaction" technique is used. A
given population (here: microbiome inside a fecal sample) is
subsampled to calculate the overall species richness in the population
while keeping the sample size as small as possible and as big as
necessary. I used *alpha\_rarefaction.py* to create an overview of
suitable rarefaction depths.  
(TODO: nice rarefaction image here)
The read numbers per sample ranged from 305,631 to 5,826,936, with a mean
of 1,932,291 and a median of 1,884,081. To guarantee that all samples are
represented and each one is subsampled I chose a rarefaction depth of
290,000 (as compared to 290,603 in the original publication), which is
less than the smallest sample's read count. I repeated the
rarefactioning ten times, using the tool
*multiple\_rarefaction\_even\_depth.py*.

From those rarefied OTU tables I calculated the $\alpha$- and $\beta$-diversity.
I calculated $\alpha$-diversity using the tool *alpha\_diversity.py*
with the number of observed species as a disance measure and merged
the results into a single table (\texttt{observed\_species.csv}) with the tool *collate\_alpha.py*.
To measure $\beta$-diversity I used UniFrac distance, the percentage
of shared branches on the phylogenetic tree, compared against the
GreenGenes reference tree. With the tool *beta\_diversity.py* I
calculated $\beta$-diversity distance matrices for all ten rarefied
tables, each with weighted and unweighted UniFrac distances as a
measure. The resulting files are called
\texttt{(unweighted\_)unifrac\_rarefaction\_290000\_N.txt}.

### Load data from preprocessed files ###

I loaded the data files created by the previous steps 3.2.1 and 3.2.2
and created convenience access methods.


```r
theCountries <- c("Malawi", "USA", "Venezuela")
theColors <- c("red", "blue", "green")
names(theColors) <- theCountries
alphaTable <- read.delim("local_copy/observed_species.csv")
betaTable <- read.delim("local_copy/unweighted_unifrac_rarefaction_290000_4.txt")
rownames(betaTable) <- betaTable[,1]; betaTable <- betaTable[,-1]
theMetadata <- read.delim("16s_mapping_decimaldot.csv")
rownames(theMetadata) <- theMetadata$X.SampleID
colnames(theMetadata)
```

```
##  [1] "X.SampleID"      "IDstring"        "Country"        
##  [4] "Region"          "FamilyID"        "FamilyMember"   
##  [7] "Age"             "Zygosity"        "Gender"         
## [10] "BMI"             "BreastFed"       "SequenceCount"  
## [13] "BarcodeSequence" "LinkerSequence"  "RunPrefix"      
## [16] "FlowcellLane"
```

```r
 # get sample ids with label of specific value
getGroupIDs <- function(label, value) {
  theMetadata[theMetadata[,label] %in% value,]$X.SampleID
}
 # get sample ids with label in numeric range
getRangeIDs <- function(label, lower, upper) {
  values <- theMetadata[,label]
  theMetadata[values >= lower & values <= upper,]$X.SampleID
}
 # get beta variance data for samples with given ids
getBetaGroup <- function(label=NULL, value=NULL, range=FALSE, ids=NULL) {
  if (length(ids) == 0) {
    if (any(is.na(c(label,value)))) { # no input at all
	  print("either label and value or ids required")
	  return(c())
	}
	if (!range)	ids <- getGroupIDs(label, value)
	else ids <- getRangeIDs(label, value[1], value[2])
  }
  betaTable[which(rownames(betaTable) %in% ids), which(colnames(betaTable) %in% ids)]
}
```

In this document, only one of 20 possible beta diversity matrices is
loaded. The others were also processed, and the produced clusterings
were very similar to the one shown here (see 3.2.5).

<!-- Usable sample sizes: 188517 or lower -->

### UniFrac distance variation with age ###

Yatsunenko et al. observed the change of microbiome composition from
infant-specific to adult configuration by comparing the composition of
each child's microbiome against the microbiome composition of all adults from the
same country [@Yats12]. As it is not completely clear how the distance to all
adults was calculated I chose to use the mean of distances of the
child to each adult.


```r
plot(NULL, NULL, xlim=c(0, 18), ylim=c(0.35, 0.85), xlab="age", ylab="UniFrac distance")
country <- list(
    "USA" = getGroupIDs("Country", "USA"),
    "Malawi" = getGroupIDs("Country", "Malawi"),
    "Venezuela" = getGroupIDs("Country", "Venezuela")
    )
adults <- list(
    "USA" = intersect(country[["USA"]], getRangeIDs("Age", 19, 99)),
    "Malawi" = intersect(country[["Malawi"]], getRangeIDs("Age", 19, 99)),
    "Venezuela" = intersect(country[["Venezuela"]], getRangeIDs("Age", 19, 99))
    )
children <- list(
    "USA" = intersect(country[["USA"]], getRangeIDs("Age", 0, 18)),
    "Malawi" = intersect(country[["Malawi"]], getRangeIDs("Age", 0, 18)),
    "Venezuela" = intersect(country[["Venezuela"]], getRangeIDs("Age", 0, 18))
    )
for (country in theCountries) {
    for (child in children[[country]]) {
        betadiv <- getBetaGroup(ids=c(child, adults[[country]]))
        x <- theMetadata[theMetadata$X.SampleID==child,]$Age
        y <- sum(betadiv[child,]) / (nrow(betadiv) - 1)
        points(x, y, col=theColors[[country]], pch=20)
    }
}
legend(x=14, y=1.0, legend=theCountries, text.col=theColors)
```

![plot of chunk DiversityChildToAdult](figure/DiversityChildToAdult.png) 

The UniFrac distances start with high values at early ages, show a
strong decline until approximately three years of age and stay
steadily low until adulthood (see fig. 1).

### PCoA analysis of betadiversity

Partitioning Around Medoids (PAM) is a method of clustering data
points into a given number k of groups. I initially assigns a data
point as center for each of the groups and then minimizes a global distance
function by iteratively swapping points and centers and assigning all
data points to the nearest new center.
Yatsunenko et al. used the $\beta$-distance matrix as a one-dimensional
dissimilarity measure for the PAM algorithm and chose k=3 clusters to
refind the samples' countries of origin in the microbiome composition [@Yats12].
I used the *pam* function from the *R* package *cluster* to repeat the
process.


```r
clu <- pam(betaTable, diss=TRUE, k=3, keep.diss=TRUE)
clusplot(clu, col.p=c("red","blue","green")[theMetadata[names(clu$clustering),]$Country], main="")
legend(x=0.1, y=-0.2, legend=theCountries, text.col=theColors)
```

![plot of chunk Clustering](figure/Clustering.png) 

```r
 # note: color order is different, because the pamobject$clustering vector has a different order
contTable <- table(clu$clustering, theMetadata[names(clu$clustering),]$Country)
contTable <- cbind(contTable[,2], contTable[,1], contTable[,3]) # correct order
pam.result <- sum(diag(contTable) / sum(contTable))
print(pam.result)
```

```
## [1] 0.8407
```

Figure 2: Colors follow the real countries of origin, while symbols denote the
country-clusters predicted by PAM.
The plot shows a separated, banana-shaped distribution of US Americans
and intermingled groups of Malawians and Venezuelans with
concentrations for each subgroups. Using different
rarefaction tables showed the same results of circa 84 per cent correct
assignment.

### SVM classificator analysis of $\beta$-diversity ###

The PCoA plot suggested that discrimination of samples by
microbiome diversity is possible. To improve the assignment I used the
implementation of support vector machines (SVM) in the *R* package
*e1071* [@e1071]. Support vector machines classify multi-dimensional datasets
by finding a hyperplane that separates the classes among (a subset
of) their features. It solves the problem of seemingly inseperable
datasets by applying a "kernel function" that adds additional
dimensions to the data's used features.

I trained a SVM to discriminate the dataset by different subsets of the
OTUs that showed the greatest variance and compared the predicted
samples with their countries of origin, using 20-fold cross validation.


```r
require(e1071)
species <- read.delim("local_copy/combined.csv", sep="\t", header=TRUE)
rownames(species) <- species[,1]
species <- species[,-1] # species level OTU table
species <- t(species)   # now rows = samples, cols = OTU counts
reps <- 20              # number of repetitions for cross validation
L <- nrow(species)
N <- L / reps           # number of samples per validation run
featureSizes <- c(1:5, seq(10, 500, len=10), ncol(species)) # 
strength <- matrix(0, nrow=reps, ncol=length(featureSizes))
for (j in 1:length(featureSizes)) {
	allOfThem <- 1:L    # available samples
	for (i in 1:reps) {
		testSet <- sample(allOfThem, N)
		allOfThem <- allOfThem[-testSet]
		training <- c(1:L)[-testSet]
		variances <- apply(species[training,], 2, var)
		topVariables <- order(variances, decreasing=TRUE)
		topVariables <- topVariables[1:featureSizes[j]]
		trnNames <- rownames(species)[training]
		tstNames <- rownames(species)[testSet]
		model <- svm(x=species[training,topVariables], y=theMetadata[trnNames,]$Country);
		prediction <- predict(model, species[testSet,topVariables])
		contingency <- table(prediction, theMetadata[tstNames,]$Country)
		strength[i,j] <- sum(diag(contingency)) / sum(contingency)
	}
}
```

```
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X332204' and 'X457371' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X283864' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X584933' and 'X312833' and 'X548603' and 'X296798' and 'X77904' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X297964' and 'X594119' and 'X563315' and 'X2879' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X33939' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X230798' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X209031' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X326170' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X197600' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X551992' and 'X326985' and 'X16058' and 'X300443' and 'X55666' and 'X319751' and 'X293754' and 'X37302' and 'X442949' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X226596' and 'X131845' and 'X142695' and 'X227521' and 'X306635' and 'X539496' and 'X159436' and 'X225899' and 'X331980' and 'X152545' and 'X513591' and 'X32254' and 'X240018' and 'X541215' and 'X149878' and 'X227060' and 'X244699' and 'X534921' and 'X255515' and 'X366191' and 'X156986' and 'X257237' and 'X566337' and 'X360826' and 'X200652' and 'X237783' and 'X531775' and 'X185766' and 'X518998' and 'X311656' and 'X178988' and 'X574890' and 'X351877' and 'X494253' and 'X562807' and 'X21491' and 'X249926' and 'X90192' and 'X548020' and 'X559239' and 'X241718' and 'X366814' and 'X226516' and 'X161090' and 'X214306' and 'X574911' and 'X53533' and 'X574265' and 'X109654' and 'X11396' and 'X593659' and 'X222810' and 'X269541' and 'X110434' and 'X396217' and 'X313816' and 'X90311' and 'X353143' and 'X288882' and 'X210920' and 'X326517' and 'X227404' and 'X262760' and 'X263805' and 'X294909' and 'X214392' and 'X296723' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X237034' and 'X561204' and 'X341939' and 'X356218' and 'X288346' and 'X237073' and 'X417573' and 'X93749' and 'X251691' and 'X240933' and 'X540317' and 'X494140' and 'X282746' and 'X236797' and 'X219729' and 'X185444' and 'X273400' and 'X103895' and 'X350772' and 'X248695' and 'X9061' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X178106' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X516411' and 'X566717' and 'X255719' and 'X150338' and 'X312275' and 'X296079' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X51161' and 'X141739' and 'X14282' and 'X576671' and 'X12591' and 'X120648' and 'X156429' and 'X518471' and 'X165020' and 'X197240' and 'X145402' and 'X77708' and 'X113319' and 'X550070' and 'X297865' and 'X292202' and 'X527641' and 'X110229' and 'X302511' and 'X341908' and 'X158223' and 'X558344' and 'X565814' and 'X163340' and 'X91416' and 'X219876' and 'X139040' and 'X587111' and 'X240189' and 'X558389' and 'X581638' and 'X361428' and 'X305547' and 'X279442' and 'X209377' and 'X169379' and 'X375076' and 'X179541' and 'X549650' and 'X553041' and 'X97947' and 'X25255' and 'X255423' and 'X347891' and 'X333448' and 'X182248' and 'X511876' and 'X470688' and 'X261434' and 'X76251' and 'X277303' and 'X593420' and 'X72406' and 'X591347' and 'X327001' and 'X467388' and 'X364600' and 'X14122' and 'X64874' and 'X593891' and 'X242915' and 'X125235' and 'X211106' and 'X15851' and 'X63942' and 'X8067' and 'X78716' and 'X550048' and 'X14059' and 'X263005' and 'X362188' and 'X139922' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X297847' and 'X270447' and 'X289202' and 'X259498' and 'X208044' and 'X70671' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X589380' and 'X44758' and 'X301405' and 'X252146' and 'X336333' and 'X171239' and 'X532384' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X14507' and 'X108282' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X457371' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X150728' and 'X318733' and 'X227017' and 'X231593' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X328422' and 'X579422' and 'X141548' and 'X288438' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X301852' and 'X34261' and 'X221369' and 'X436548' and 'X226964' and 'X151886' and 'X264974' and 'X527879' and 'X230970' and 'X111201' and 'X158223' and 'X28461' and 'X293780' and 'X515774' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X15828' and 'X206803' and 'X138206' and 'X132661' and 'X307909' and 'X144291' and 'X27088' and 'X538707' and 'X428705' and 'X46393' and 'X201491' and 'X525091' and 'X12131' and 'X261010' and 'X276258' and 'X130529' and 'X135577' and 'X250036' and 'X101262' and 'X143369' and 'X162964' and 'X545997' and 'X315609' and 'X177473' and 'X101249' and 'X272478' and 'X307781' and 'X291838' and 'X127012' and 'X244063' and 'X589570' and 'X255658' and 'X142794' and 'X104872' and 'X358745' and 'X588682' and 'X208701' and 'X565780' and 'X221621' and 'X563358' and 'X553759' and 'X273309' and 'X168547' and 'X113624' and 'X357330' and 'X523807' and 'X93483' and 'X257150' and 'X543154' and 'X513639' and 'X531410' and 'X136995' and 'X305064' and 'X249770' and 'X141539' and 'X534686' and 'X235693' and 'X469179' and 'X558503' and 'X322897' and 'X113347' and 'X563703' and 'X230108' and 'X81666' and 'X162705' and 'X105915' and 'X213030' and 'X32571' and 'X90697' and 'X541511' and 'X113599' and 'X152024' and 'X200086' and 'X565027' and 'X152039' and 'X222566' and 'X168491' and 'X398839' and 'X19216' and 'X177392' and 'X261435' and 'X27938' and 'X102227' and 'X106027' and 'X514663' and 'X190531' and 'X342774' and 'X272809' and 'X312486' and 'X332547' and 'X282458' and 'X175586' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X208044' and 'X319810' and 'X302376' and 'X293997' and 'X44758' and 'X259498' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X317218' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X100680' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X155703' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X14507' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X223318' and 'X43339' and 'X324275' and 'X108861' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X100068' and 'X265716' and 'X594352' and 'X326693' and 'X342719' and 'X181683' and 'X518565' and 'X270263' and 'X28165' and 'X148303' and 'X106422' and 'X243412' and 'X273437' and 'X537859' and 'X141548' and 'X134788' and 'X522744' and 'X108748' and 'X294604' and 'X319438' and 'X514611' and 'X303508' and 'X144685' and 'X210958' and 'X159398' and 'X231028' and 'X96001' and 'X532771' and 'X140034' and 'X330247' and 'X284456' and 'X141366' and 'X534964' and 'X383716' and 'X244390' and 'X14284' and 'X223706' and 'X326600' and 'X369119' and 'X319241' and 'X338749' and 'X108266' and 'X109334' and 'X60771' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289202' and 'X289208' and 'X97369' and 'X42642' and 'X355836' and 'X206709' and 'X32016' and 'X540672' and 'X143227' and 'X583489' and 'X314357' and 'X265971' and 'X524595' and 'X32247' and 'X244056' and 'X166076' and 'X512939' and 'X114159' and 'X143097' and 'X548602' and 'X562772' and 'X45520' and 'X130481' and 'X235519' and 'X243468' and 'X257403' and 'X255454' and 'X218287' and 'X568948' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X323206' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X70671' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X589380' and 'X44758' and 'X171239' and 'X532384' and 'X268454' and 'X351744' and 'X301405' and 'X252146' and 'X336333' and 'X190976' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X50040' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X301606' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X15059' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X347159' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X108861' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X347809' and 'X276889' and 'X11113' and 'X93761' and 'X100068' and 'X265716' and 'X328695' and 'X10302' and 'X230933' and 'X21434' and 'X539398' and 'X237892' and 'X591519' and 'X223183' and 'X134301' and 'X463007' and 'X306735' and 'X5598' and 'X577152' and 'X274257' and 'X274560' and 'X586378' and 'X161584' and 'X261254' and 'X552887' and 'X236325' and 'X557255' and 'X9458' and 'X339847' and 'X16301' and 'X243412' and 'X344500' and 'X276258' and 'X130529' and 'X155884' and 'X269994' and 'X71248' and 'X209866' and 'X168506' and 'X135577' and 'X250036' and 'X101262' and 'X143369' and 'X101323' and 'X127012' and 'X244063' and 'X589570' and 'X255658' and 'X142794' and 'X104872' and 'X358745' and 'X588682' and 'X153137' and 'X208701' and 'X565780' and 'X221621' and 'X563358' and 'X553759' and 'X273309' and 'X168547' and 'X113624' and 'X357330' and 'X523807' and 'X93483' and 'X257150' and 'X543154' and 'X144150' and 'X513639' and 'X531410' and 'X136995' and 'X305064' and 'X249770' and 'X141539' and 'X534686' and 'X235693' and 'X469179' and 'X509899' and 'X558503' and 'X322897' and 'X113347' and 'X563703' and 'X230108' and 'X81666' and 'X162705' and 'X105915' and 'X213030' and 'X32571' and 'X90697' and 'X541511' and 'X113599' and 'X152024' and 'X200086' and 'X565027' and 'X152039' and 'X332976' and 'X14265' and 'X206709' and 'X32016' and 'X540672' and 'X143227' and 'X583489' and 'X314357' and 'X265971' and 'X524595' and 'X32247' and 'X244056' and 'X166076' and 'X512939' and 'X114159' and 'X143097' and 'X548602' and 'X562772' and 'X45520' and 'X130481' and 'X235519' and 'X243468' and 'X257403' and 'X255454' and 'X218287' and 'X568948' and 'X262446' and 'X188416' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X258376' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X342642' and 'X349634' and 'X578119' and 'X494937' and 'X429838' and 'X279277' and 'X560064' and 'X520283' and 'X140377' and 'X193229' and 'X463007' and 'X306735' and 'X5598' and 'X577152' and 'X529472' and 'X589923' and 'X228065' and 'X150663' and 'X251317' and 'X264560' and 'X167112' and 'X142457' and 'X292384' and 'X91691' and 'X46239' and 'X543710' and 'X541375' and 'X26896' and 'X311254' and 'X12296' and 'X102242' and 'X210760' and 'X217695' and 'X4530' and 'X571641' and 'X202014' and 'X525417' and 'X71543' and 'X344620' and 'X109263' and 'X511542' and 'X311032' and 'X140772' and 'X99274' and 'X206237' and 'X572372' and 'X185650' and 'X34798' and 'X321220' and 'X105823' and 'X554372' and 'X593617' and 'X511844' and 'X288689' and 'X512790' and 'X291747' and 'X338371' and 'X166153' and 'X231606' and 'X292053' and 'X6849' and 'X414219' and 'X262611' and 'X147459' and 'X338821' and 'X150' and 'X518474' and 'X142117' and 'X551463' and 'X319318' and 'X217779' and 'X132318' and 'X114803' and 'X565936' and 'X5497' and 'X147281' and 'X291934' and 'X131852' and 'X6144' and 'X534993' and 'X271528' and 'X125257' and 'X74125' and 'X77568' and 'X242049' and 'X185640' and 'X267323' and 'X583979' and 'X246762' and 'X14559' and 'X64356' and 'X246819' and 'X584366' and 'X77057' and 'X329045' and 'X11428' and 'X52053' and 'X583887' and 'X533172' and 'X238291' and 'X266522' and 'X311345' and 'X77712' and 'X586961' and 'X508084' and 'X247938' and 'X263760' and 'X99041' and 'X564501' and 'X517399' and 'X136867' and 'X322801' and 'X153617' and 'X247875' and 'X343296' and 'X554990' and 'X138031' and 'X329326' and 'X354851' and 'X426848' and 'X72805' and 'X102610' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X100680' and 'X302376' and 'X44758' and 'X171239' and 'X293997' and 'X336333' and 'X532384' and 'X319810' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X535573' and 'X592633' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X288438' and 'X14507' and 'X301405' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X289341' and 'X296079' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X294146' and 'X584326' and 'X352502' and 'X138036' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X190646' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X11113' and 'X93761' and 'X35138' and 'X255169' and 'X141313' and 'X554916' and 'X217993' and 'X383023' and 'X270739' and 'X340680' and 'X133472' and 'X554797' and 'X184695' and 'X66635' and 'X15719' and 'X548876' and 'X173100' and 'X158930' and 'X220686' and 'X92717' and 'X252937' and 'X296038' and 'X528107' and 'X12486' and 'X296516' and 'X169815' and 'X128261' and 'X355777' and 'X577176' and 'X113319' and 'X550070' and 'X297865' and 'X292202' and 'X527641' and 'X302511' and 'X560064' and 'X140377' and 'X558344' and 'X180141' and 'X529472' and 'X589923' and 'X228065' and 'X150663' and 'X251317' and 'X264560' and 'X167112' and 'X142457' and 'X292384' and 'X215972' and 'X274257' and 'X274560' and 'X59529' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X42877' and 'X573052' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X155703' and 'X210958' and 'X159398' and 'X231028' and 'X96001' and 'X532771' and 'X140034' and 'X330247' and 'X284456' and 'X141366' and 'X534964' and 'X383716' and 'X244390' and 'X14284' and 'X223706' and 'X326600' and 'X369119' and 'X333033' and 'X225511' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X324283' and 'X366048' and 'X147031' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X353409' and 'X532384' and 'X171239' and 'X301405' and 'X544894' and 'X293092' and 'X457990' and 'X299900' and 'X262741' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X3861' and 'X584933' and 'X332204' and 'X457371' and 'X305211' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X242868' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X556604' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X584296' and 'X110034' and 'X533537' and 'X594119' and 'X563315' and 'X537729' and 'X2879' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X237063' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X163461' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X550996' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X12984' and 'X356039' and 'X179040' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X103842' and 'X255234' and 'X177331' and 'X296146' and 'X141548' and 'X537651' and 'X152568' and 'X7553' and 'X579047' and 'X563303' and 'X160679' and 'X288049' and 'X14540' and 'X148392' and 'X299370' and 'X519353' and 'X183478' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X409404' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X527507' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X298447' and 'X528107' and 'X12486' and 'X296516' and 'X169815' and 'X128261' and 'X355777' and 'X577176' and 'X9565' and 'X204226' and 'X558344' and 'X274290' and 'X135960' and 'X163406' and 'X12530' and 'X137961' and 'X528709' and 'X315979' and 'X200178' and 'X512361' and 'X65413' and 'X6236' and 'X9568' and 'X112277' and 'X132196' and 'X546951' and 'X6997' and 'X6733' and 'X539360' and 'X324234' and 'X83477' and 'X530869' and 'X435512' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X243412' and 'X277460' and 'X272752' and 'X101144' and 'X355736' and 'X561950' and 'X327286' and 'X71049' and 'X108408' and 'X114384' and 'X539301' and 'X206221' and 'X103476' and 'X353702' and 'X540454' and 'X33643' and 'X252471' and 'X446356' and 'X242798' and 'X241214' and 'X240850' and 'X220305' and 'X138741' and 'X77610' and 'X142028' and 'X230396' and 'X236550' and 'X338628' and 'X167066' and 'X138602' and 'X3538' and 'X14641' and 'X138316' and 'X269218' and 'X163551' and 'X239951' and 'X200217' and 'X578443' and 'X248784' and 'X320812' and 'X247439' and 'X114526' and 'X594040' and 'X81554' and 'X199105' and 'X4416' and 'X379604' and 'X358044' and 'X319537' and 'X513982' and 'X340207' and 'X486153' and 'X151804' and 'X534213' and 'X187236' and 'X383501' and 'X128201' and 'X20974' and 'X264010' and 'X248199' and 'X509171' and 'X369064' and 'X541286' and 'X100639' and 'X304693' and 'X12364' and 'X256308' and 'X345664' and 'X137631' and 'X365101' and 'X139166' and 'X25041' and 'X1422' and 'X277641' and 'X309862' and 'X264742' and 'X354510' and 'X19356' and 'X285281' and 'X233682' and 'X295539' and 'X593970' and 'X185752' and 'X309361' and 'X17434' and 'X349548' and 'X414219' and 'X262611' and 'X247268' and 'X242096' and 'X168183' and 'X249191' and 'X571956' and 'X164506' and 'X526198' and 'X308932' and 'X484437' and 'X358429' and 'X593325' and 'X252984' and 'X415943' and 'X566333' and 'X557691' and 'X101880' and 'X321738' and 'X71891' and 'X287452' and 'X327712' and 'X222507' and 'X494503' and 'X532885' and 'X254162' and 'X509986' and 'X587498' and 'X346479' and 'X193662' and 'X593400' and 'X557465' and 'X311867' and 'X218620' and 'X265444' and 'X157893' and 'X166436' and 'X101225' and 'X256598' and 'X513265' and 'X115015' and 'X8717' and 'X402651' and 'X136243' and 'X67654' and 'X6953' and 'X189488' and 'X589725' and 'X561739' and 'X62675' and 'X343091' and 'X11301' and 'X259544' and 'X593016' and 'X31050' and 'X278562' and 'X167215' and 'X222902' and 'X551031' and 'X366378' and 'X237395' and 'X554177' and 'X248158' and 'X42141' and 'X14642' and 'X257520' and 'X15022' and 'X193304' and 'X571004' and 'X15021' and 'X588865' and 'X14990' and 'X387035' and 'X269212' and 'X298787' and 'X300849' and 'X359303' and 'X101657' and 'X17956' and 'X540549' and 'X277506' and 'X153253' and 'X104789' and 'X553112' and 'X100099' and 'X287776' and 'X17599' and 'X578586' and 'X587323' and 'X247009' and 'X180400' and 'X187366' and 'X48785' and 'X544351' and 'X88595' and 'X289958' and 'X130998' and 'X328205' and 'X547395' and 'X183736' and 'X471352' and 'X575425' and 'X559213' and 'X188247' and 'X100562' and 'X242915' and 'X125235' and 'X198983' and 'X266941' and 'X133105' and 'X32167' and 'X583557' and 'X91993' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X140393' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X179131' and 'X335755' and 'X6721' and 'X143735' and 'X252834' and 'X34222' and 'X528107' and 'X12486' and 'X296516' and 'X169815' and 'X128261' and 'X355777' and 'X577176' and 'X279277' and 'X515662' and 'X206143' and 'X261270' and 'X9838' and 'X91270' and 'X44695' and 'X312781' and 'X47454' and 'X288379' and 'X570916' and 'X164160' and 'X28461' and 'X293780' and 'X515774' and 'X217695' and 'X4530' and 'X571641' and 'X229171' and 'X320378' and 'X305758' and 'X189891' and 'X271865' and 'X99071' and 'X181996' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X361428' and 'X87601' and 'X131680' and 'X312294' and 'X350384' and 'X42483' and 'X549656' and 'X99155' and 'X343424' and 'X562486' and 'X99304' and 'X366545' and 'X414219' and 'X262611' and 'X168491' and 'X320835' and 'X143178' and 'X247045' and 'X202045' and 'X159766' and 'X195654' and 'X367606' and 'X533060' and 'X397321' and 'X225760' and 'X217413' and 'X290559' and 'X87249' and 'X251112' and 'X213161' and 'X223672' and 'X335630' and 'X107195' and 'X41720' and 'X21795' and 'X239601' and 'X532312' and 'X292243' and 'X266173' and 'X369734' and 'X69993' and 'X327592' and 'X291491' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X208044' and 'X319810' and 'X302376' and 'X293997' and 'X44758' and 'X259498' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X341286' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X100680' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X516411' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X155703' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X14507' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X327694' and 'X271891' and 'X183753' and 'X570082' and 'X12357' and 'X73352' and 'X293243' and 'X546078' and 'X582179' and 'X254195' and 'X559587' and 'X11113' and 'X93761' and 'X582466' and 'X494937' and 'X429838' and 'X10302' and 'X594352' and 'X149284' and 'X326693' and 'X560064' and 'X140377' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X312781' and 'X47454' and 'X288379' and 'X570916' and 'X164160' and 'X45848' and 'X578429' and 'X250514' and 'X255802' and 'X210792' and 'X384373' and 'X575450' and 'X332574' and 'X556036' and 'X532912' and 'X533602' and 'X548891' and 'X99509' and 'X123738' and 'X543676' and 'X552729' and 'X487233' and 'X237879' and 'X238675' and 'X9994' and 'X453756' and 'X346926' and 'X325677' and 'X336984' and 'X112916' and 'X217403' and 'X89541' and 'X141394' and 'X219106' and 'X148662' and 'X154113' and 'X217798' and 'X324704' and 'X241908' and 'X554341' and 'X548878' and 'X218845' and 'X73544' and 'X549839' and 'X178178' and 'X334226' and 'X7598' and 'X239195' and 'X410855' and 'X251499' and 'X29284' and 'X227361' and 'X39461' and 'X207957' and 'X130709' and 'X240983' and 'X154220' and 'X316541' and 'X537483' and 'X3094' and 'X160115' and 'X242284' and 'X511637' and 'X111394' and 'X91588' and 'X244970' and 'X55946' and 'X5163' and 'X237583' and 'X591751' and 'X58634' and 'X193867' and 'X517065' and 'X261978' and 'X525666' and 'X585097' and 'X534483' and 'X113942' and 'X239741' and 'X142304' and 'X76549' and 'X255914' and 'X515421' and 'X97778' and 'X576689' and 'X140788' and 'X557604' and 'X154684' and 'X584480' and 'X98983' and 'X529793' and 'X164098' and 'X583293' and 'X239095' and 'X244548' and 'X362116' and 'X237095' and 'X329486' and 'X137809' and 'X319977' and 'X277776' and 'X250987' and 'X147665' and 'X273437' and 'X537859' and 'X141110' and 'X558297' and 'X225887' and 'X142173' and 'X143400' and 'X144142' and 'X410399' and 'X148779' and 'X410272' and 'X147943' and 'X149631' and 'X105233' and 'X314203' and 'X309588' and 'X221399' and 'X225933' and 'X18245' and 'X565593' and 'X110332' and 'X164828' and 'X207298' and 'X86557' and 'X166153' and 'X366545' and 'X210958' and 'X159398' and 'X231028' and 'X96001' and 'X532771' and 'X140034' and 'X330247' and 'X284456' and 'X141366' and 'X338075' and 'X239441' and 'X246316' and 'X352244' and 'X241392' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289202' and 'X289208' and 'X97369' and 'X42642' and 'X224877' and 'X78716' and 'X550048' and 'X14059' and 'X263005' and 'X362188' and 'X144186' and 'X318862' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X302308' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X457990' and 'X190976' and 'X299900' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X3861' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X237188' and 'X248972' and 'X255719' and 'X150338' and 'X312275' and 'X289341' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X533537' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X163461' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X177331' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X14540' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X527507' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X108843' and 'X323595' and 'X574338' and 'X147308' and 'X296402' and 'X104585' and 'X171639' and 'X19867' and 'X230933' and 'X21434' and 'X539398' and 'X237892' and 'X591519' and 'X223183' and 'X181573' and 'X576232' and 'X198202' and 'X7303' and 'X525159' and 'X211927' and 'X2588' and 'X291173' and 'X243412' and 'X362169' and 'X525091' and 'X12131' and 'X146138' and 'X230801' and 'X244892' and 'X550215' and 'X152715' and 'X6365' and 'X351023' and 'X240686' and 'X63343' and 'X191576' and 'X312051' and 'X471122' and 'X6849' and 'X534964' and 'X383716' and 'X244390' and 'X14284' and 'X223706' and 'X326600' and 'X369119' and 'X52593' and 'X168491' and 'X247268' and 'X242096' and 'X168183' and 'X249191' and 'X571956' and 'X164506' and 'X526198' and 'X308932' and 'X484437' and 'X358429' and 'X593325' and 'X252984' and 'X415943' and 'X566333' and 'X557691' and 'X101880' and 'X321738' and 'X71891' and 'X287452' and 'X327712' and 'X222507' and 'X494503' and 'X532885' and 'X254162' and 'X509986' and 'X587498' and 'X346479' and 'X193662' and 'X593400' and 'X557465' and 'X311867' and 'X218620' and 'X157893' and 'X166436' and 'X256598' and 'X513265' and 'X115015' and 'X8717' and 'X8895' and 'X402651' and 'X136243' and 'X67654' and 'X6953' and 'X189488' and 'X589725' and 'X561739' and 'X62675' and 'X343091' and 'X177392' and 'X213152' and 'X547346' and 'X207832' and 'X347680' and 'X320835' and 'X225511' and 'X264539' and 'X143178' and 'X247045' and 'X202045' and 'X159766' and 'X195654' and 'X367606' and 'X533060' and 'X397321' and 'X225760' and 'X217413' and 'X290559' and 'X87249' and 'X251112' and 'X213161' and 'X223672' and 'X335630' and 'X107195' and 'X41720' and 'X21795' and 'X239601' and 'X532312' and 'X292243' and 'X266173' and 'X369734' and 'X519627' and 'X213877' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X524900' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X38061' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X544546' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X114748' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X103120' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X89402' and 'X105304' and 'X560064' and 'X140377' and 'X251213' and 'X193229' and 'X14648' and 'X245864' and 'X291317' and 'X263804' and 'X257022' and 'X78216' and 'X273134' and 'X356889' and 'X223711' and 'X169060' and 'X200998' and 'X287758' and 'X107742' and 'X131214' and 'X319790' and 'X543389' and 'X109915' and 'X350331' and 'X114455' and 'X435512' and 'X181573' and 'X576232' and 'X198202' and 'X7303' and 'X524433' and 'X128098' and 'X256695' and 'X242252' and 'X77638' and 'X251102' and 'X16001' and 'X259903' and 'X70502' and 'X123281' and 'X315994' and 'X102595' and 'X87601' and 'X131680' and 'X312294' and 'X350384' and 'X9245' and 'X322850' and 'X303387' and 'X111945' and 'X547346' and 'X207832' and 'X347680' and 'X259544' and 'X593016' and 'X278562' and 'X560114' and 'X110474' and 'X560962' and 'X581039' and 'X247473' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X297847' and 'X270447' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X544894' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X584933' and 'X332204' and 'X2879' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X9349' and 'X216695' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X457371' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X550996' and 'X179040' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X155703' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X110737' and 'X255719' and 'X150338' and 'X312275' and 'X184872' and 'X312808' and 'X321950' and 'X206189' and 'X310913' and 'X579184' and 'X141618' and 'X572519' and 'X169164' and 'X140787' and 'X161410' and 'X564867' and 'X42640' and 'X470178' and 'X126508' and 'X594352' and 'X326693' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X530224' and 'X134038' and 'X572970' and 'X134081' and 'X74034' and 'X150235' and 'X262102' and 'X293284' and 'X539124' and 'X342821' and 'X535088' and 'X346430' and 'X149934' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X109798' and 'X299309' and 'X15053' and 'X9136' and 'X583608' and 'X83477' and 'X530869' and 'X70247' and 'X25925' and 'X525091' and 'X12131' and 'X514611' and 'X244892' and 'X550215' and 'X152715' and 'X210958' and 'X159398' and 'X231028' and 'X96001' and 'X532771' and 'X140034' and 'X330247' and 'X284456' and 'X141366' and 'X11301' and 'X223975' and 'X245248' and 'X328536' and 'X574241' and 'X319219' and 'X267557' and 'X150955' and 'X293015' and 'X91993' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X100680' and 'X302376' and 'X235381' and 'X44758' and 'X171239' and 'X293997' and 'X336333' and 'X532384' and 'X319810' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X288438' and 'X14507' and 'X301405' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X323964' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X578230' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X42877' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X108843' and 'X323595' and 'X574338' and 'X450645' and 'X539260' and 'X153307' and 'X585106' and 'X279277' and 'X463007' and 'X306735' and 'X5598' and 'X577152' and 'X529472' and 'X589923' and 'X228065' and 'X150663' and 'X251317' and 'X264560' and 'X167112' and 'X142457' and 'X292384' and 'X270263' and 'X28165' and 'X109798' and 'X299309' and 'X274560' and 'X28461' and 'X293780' and 'X515774' and 'X525417' and 'X71543' and 'X344620' and 'X109263' and 'X511542' and 'X311032' and 'X140772' and 'X99274' and 'X206237' and 'X572372' and 'X185650' and 'X340706' and 'X329497' and 'X139424' and 'X15828' and 'X206803' and 'X138206' and 'X132661' and 'X307909' and 'X144291' and 'X27088' and 'X250563' and 'X511130' and 'X136529' and 'X546063' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X70346' and 'X203489' and 'X560552' and 'X231606' and 'X292053' and 'X88727' and 'X518444' and 'X176117' and 'X178223' and 'X33773' and 'X24308' and 'X172303' and 'X220551' and 'X129425' and 'X270244' and 'X224877' and 'X519627' and 'X213877' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X524900' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X532384' and 'X202911' and 'X301405' and 'X336333' and 'X171239' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X535573' and 'X592633' and 'X9349' and 'X584933' and 'X332204' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X457371' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X190646' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X103842' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X155703' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X436548' and 'X226964' and 'X151886' and 'X264974' and 'X527879' and 'X230970' and 'X111201' and 'X383023' and 'X270739' and 'X340680' and 'X554797' and 'X184695' and 'X66635' and 'X494937' and 'X429838' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X560064' and 'X140377' and 'X592917' and 'X249280' and 'X241853' and 'X255109' and 'X557292' and 'X195370' and 'X238944' and 'X572591' and 'X349012' and 'X87984' and 'X514449' and 'X69411' and 'X100108' and 'X588461' and 'X277460' and 'X272752' and 'X101144' and 'X355736' and 'X561950' and 'X71049' and 'X108408' and 'X114384' and 'X539301' and 'X206221' and 'X103476' and 'X353702' and 'X540454' and 'X33643' and 'X252471' and 'X446356' and 'X242798' and 'X273437' and 'X537859' and 'X162964' and 'X545997' and 'X315609' and 'X240612' and 'X136553' and 'X545988' and 'X33152' and 'X143890' and 'X114895' and 'X34670' and 'X230376' and 'X27536' and 'X552651' and 'X210958' and 'X159398' and 'X231028' and 'X96001' and 'X532771' and 'X140034' and 'X330247' and 'X284456' and 'X141366' and 'X319241' and 'X338749' and 'X333033' and 'X223975' and 'X245248' and 'X328536' and 'X574241' and 'X560114' and 'X110474' and 'X560962' and 'X581039' and 'X247473' and 'X198983' and 'X266941' and 'X133105' and 'X32167' and 'X583557' and 'X519627' and 'X213877' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X265676' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X190646' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X179131' and 'X270840' and 'X368350' and 'X144163' and 'X328695' and 'X143735' and 'X252834' and 'X34222' and 'X10302' and 'X312978' and 'X328981' and 'X568580' and 'X518407' and 'X312767' and 'X276703' and 'X109798' and 'X299309' and 'X421948' and 'X273437' and 'X537859' and 'X6849' and 'X33773' and 'X24308' and 'X172303' and 'X333033' and 'X255423' and 'X347891' and 'X333448' and 'X182248' and 'X264539' and 'X198983' and 'X266941' and 'X133105' and 'X32167' and 'X583557' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X453159' and 'X270447' and 'X297847' and 'X208044' and 'X319810' and 'X302376' and 'X293997' and 'X44758' and 'X259498' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X341286' and 'X584933' and 'X332204' and 'X457371' and 'X178106' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X100680' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X323964' and 'X255719' and 'X150338' and 'X312275' and 'X324596' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X297379' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X546128' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X14507' and 'X136707' and 'X315609' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X450645' and 'X539260' and 'X585106' and 'X528107' and 'X12486' and 'X296516' and 'X169815' and 'X128261' and 'X355777' and 'X577176' and 'X328981' and 'X568580' and 'X518407' and 'X346596' and 'X533016' and 'X25337' and 'X255109' and 'X557292' and 'X195370' and 'X238944' and 'X572591' and 'X349012' and 'X87984' and 'X514449' and 'X69411' and 'X100108' and 'X588461' and 'X250514' and 'X255802' and 'X210792' and 'X384373' and 'X575450' and 'X332574' and 'X556036' and 'X532912' and 'X533602' and 'X548891' and 'X99509' and 'X123738' and 'X543676' and 'X552729' and 'X487233' and 'X237879' and 'X238675' and 'X453756' and 'X346926' and 'X325677' and 'X336984' and 'X112916' and 'X217403' and 'X89541' and 'X141394' and 'X219106' and 'X148662' and 'X154113' and 'X217798' and 'X324704' and 'X241908' and 'X554341' and 'X548878' and 'X218845' and 'X73544' and 'X549839' and 'X178178' and 'X334226' and 'X7598' and 'X239195' and 'X410855' and 'X251499' and 'X29284' and 'X227361' and 'X39461' and 'X207957' and 'X130709' and 'X240983' and 'X154220' and 'X316541' and 'X537483' and 'X3094' and 'X160115' and 'X242284' and 'X511637' and 'X111394' and 'X91588' and 'X244970' and 'X55946' and 'X5163' and 'X237583' and 'X591751' and 'X58634' and 'X193867' and 'X517065' and 'X261978' and 'X525666' and 'X585097' and 'X534483' and 'X113942' and 'X239741' and 'X142304' and 'X76549' and 'X255914' and 'X515421' and 'X97778' and 'X576689' and 'X140788' and 'X154684' and 'X584480' and 'X98983' and 'X529793' and 'X164098' and 'X583293' and 'X239095' and 'X244548' and 'X362116' and 'X237095' and 'X329486' and 'X137809' and 'X319977' and 'X277776' and 'X250987' and 'X147665' and 'X16001' and 'X259903' and 'X70502' and 'X421948' and 'X7553' and 'X563303' and 'X511844' and 'X231606' and 'X292053' and 'X579129' and 'X563533' and 'X101106' and 'X14281' and 'X590820' and 'X111945' and 'X104744' and 'X368332' and 'X529429' and 'X547346' and 'X207832' and 'X347680' and 'X320835' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289202' and 'X289208' and 'X97369' and 'X42642' and 'X64874' and 'X593891' and 'X519627' and 'X213877' and 'X101707' and 'X255964' and 'X108062' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X453159' and 'X323206' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X302376' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X584933' and 'X332204' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X546128' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X15059' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X103120' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X42877' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X183478' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X11113' and 'X93761' and 'X328695' and 'X342642' and 'X349634' and 'X578119' and 'X450645' and 'X539260' and 'X585106' and 'X323964' and 'X327084' and 'X362128' and 'X556604' and 'X312833' and 'X548603' and 'X296798' and 'X301924' and 'X262102' and 'X293284' and 'X539124' and 'X180141' and 'X524433' and 'X128098' and 'X256695' and 'X339847' and 'X16301' and 'X514449' and 'X69411' and 'X100108' and 'X588461' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X421948' and 'X123281' and 'X315994' and 'X102595' and 'X535998' and 'X544412' and 'X327858' and 'X10993' and 'X244892' and 'X550215' and 'X152715' and 'X593970' and 'X185752' and 'X309361' and 'X17434' and 'X349548' and 'X135577' and 'X250036' and 'X101262' and 'X143369' and 'X147459' and 'X338821' and 'X150' and 'X518474' and 'X142117' and 'X551463' and 'X319318' and 'X217779' and 'X132318' and 'X114803' and 'X565936' and 'X5497' and 'X104744' and 'X368332' and 'X529429' and 'X108266' and 'X109334' and 'X60771' and 'X355836' and 'X198983' and 'X266941' and 'X133105' and 'X32167' and 'X583557' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X453159' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X293997' and 'X302376' and 'X100680' and 'X44758' and 'X171239' and 'X532384' and 'X202911' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X299008' and 'X558276' and 'X50040' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X271680' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X516411' and 'X255719' and 'X9349' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X550996' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X546128' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X579047' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X582466' and 'X368350' and 'X144163' and 'X342642' and 'X349634' and 'X578119' and 'X134301' and 'X592917' and 'X249280' and 'X150172' and 'X73191' and 'X435512' and 'X68348' and 'X15238' and 'X292979' and 'X79242' and 'X268529' and 'X174572' and 'X162913' and 'X270519' and 'X249360' and 'X533148' and 'X277967' and 'X356164' and 'X587574' and 'X290178' and 'X288689' and 'X512790' and 'X291747' and 'X338371' and 'X104744' and 'X368332' and 'X562084' and 'X529429' and 'X167215' and 'X222902' and 'X551031' and 'X366378' and 'X76251' and 'X277303' and 'X264539' and 'X332976' and 'X14265' and 'X310510' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X224877' and 'X78716' and 'X550048' and 'X14059' and 'X263005' and 'X173071' and 'X362188' and 'X198983' and 'X266941' and 'X133105' and 'X32167' and 'X583557' and 'X91993' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X293997' and 'X302376' and 'X302308' and 'X100680' and 'X532384' and 'X171239' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X299008' and 'X558276' and 'X50040' and 'X584933' and 'X332204' and 'X457371' and 'X550996' and 'X232356' and 'X42877' and 'X579422' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X353109' and 'X323964' and 'X255719' and 'X9349' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X297379' and 'X594119' and 'X563315' and 'X2879' and 'X19141' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X328422' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X91993' and 'X108449' and 'X110737' and 'X322112' and 'X450645' and 'X539260' and 'X585106' and 'X15719' and 'X548876' and 'X158930' and 'X220686' and 'X92717' and 'X252937' and 'X296038' and 'X193229' and 'X222685' and 'X167766' and 'X549261' and 'X181573' and 'X576232' and 'X198202' and 'X7303' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X75333' and 'X44758' and 'X338567' and 'X14577' and 'X340706' and 'X329497' and 'X139424' and 'X538707' and 'X428705' and 'X46393' and 'X201491' and 'X273437' and 'X537859' and 'X68348' and 'X15238' and 'X79242' and 'X244892' and 'X550215' and 'X152715' and 'X579129' and 'X563533' and 'X101106' and 'X14281' and 'X436019' and 'X590820' and 'X52593' and 'X320835' and 'X310510' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' constant. Cannot scale data.
```

```r
boxplot(strength, xlab="number of features", ylab="prediction success", axes=FALSE)
abline(h=pam.result, col="blue")
axis(side=1, at=1:length(featureSizes), labels=as.integer(featureSizes), las=2)
axis(side=2, at=seq(from=0, to=1.2, by=0.2))
```

![plot of chunk Classifiers](figure/Classifiers.png) 
<!-- all hail to the cache=TRUE parameter! -->

Figure 4 shows barplots of SVM-classifier prediction of the 16S dataset. X-axis
describes the number of features used. The features were selected by greatest
variance. A blue horizontal line shows the prediction success of PAM clustering for
comparison.

TODO: Draw evenly from all populations according to their
size. Improve feature sizes, include 90, as stated in the paper. Maybe
use p-values instead of variances (will be even slower).

### Bacterial diversity with age ###

Yatsunenko et al. found that the number of OTUs inside the fecal samples
increased with age [@Yats12]. To verify this I calculated the means of OTU
counts found in the ten rarefaction repetitions for each sample and plotted
them against each samples' age.


```r
rarefaction <- 188517
alpha <-alphaTable[alphaTable$sequences.per.sample==rarefaction,][,4:ncol(alphaTable)]
counts <- colMeans(data.matrix(alpha))
names(counts) <- colnames(alpha)
x <- theMetadata[names(counts),]$Age
cols <- c("green", "blue", "red")[theMetadata[names(counts),]$Country]
plot(x, counts, col=cols, pch=20, ylab="Number of OTUs", xlab="Age", lab=c(20, 6, 7))
legend(x=55, y=500, legend=theCountries, text.col=theColors)
```

![plot of chunk BacterialDiversityWithAge](figure/BacterialDiversityWithAge.png) 

TODO: Recompile on linux box, where the data doesn't "drop to the
floor". Find out, why these conversion errors appear on the Mac only.

The plot shows an increasing number of OTUs with increasing age (fig. 4).

### Processing of whole genome shotgun sequence data ###

*   Paper: "preprocessing was done using custom Perl scripts and publicly available
    software tools"
*   Filtering for degenerate sequences, duplicates: possible with custom script, but
    huge amounts of cpu time needed
*   Further preprocessing unknown

Decisions

*   created a *filter\_wgs\_reads.py* and then *faster\_filter.py* for
    proof of concept, but didn't run them on the whole dataset 
*   downloaded data processed by automatized qiime workflow from
    MG-RAST (.stats files)
*   try to recreate findings with this (different?) data
*   found that even improved mapping files could not identify samples
    ->  decided to ignore that part

\clearpage

Discussion
==========

Pro:

*   16S data analysis can be reproduced very well, plots virtually identical
*   conclusions very plausible
*   vast supplementary data available, so recreation of sample mapping possible

Con:

*   some minor vaguenesses in the description (eg "distance to all adults")
*   missing mapping files: this analysis relies very heavily on the mapping
    between samples and their metadata
*   no meaningful reproduction of WGS data possible due to lacking mapping file
    (Gigabytes of data useless because of one metafile)
*   switch from good description and documentation of 16S data to general,
	useless statements for WGS analysis ("custom scripts and publicly available
	tools")
*   PAM clustering of beta-div not very meaningful, SVM proved a lot better

Software used
=============

TODO: Integrate properly into methods

The metadata mapping files were created using *LibreOffice* software
suite version 4.2.4.2.  
The QIIME workflow was followed using the *QIIME* suite of software
tools, version 1.8.0 as described by Caporaso et al. [@QIIME].  
Data analysis and visualization was done using the *R* statistical sofware
version 3.1.0 [@R], the package *cluster* by Maechler et al. [@cluster] for
partitioning around medoids
and the package *e1071* by Meyer et al. [@e1071] for support vector machine
classifiers.  
The custom Python scripts mentioned in the scope of this document can
be found on
[my github page](https://github.com/hermann-p/yatsunenko-2012-microbiome).

References
==========

<!-- auto-filled by pandoc-citeproc -->
