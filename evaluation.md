% Reproduction of "Human Gut Microbiome Viewed Across Age and Geography"
% Hermann Pauly
% \today

Abstract
======

Soon to come...

Introduction
========

In science we try to produce reliable and unambigous results, working hard to avoid
common pitfalls and sloppy word, and we expect the same from our fellow scientists.
Nevertheless, an unacceptable number of publications *(numbers)* contain grave errors
*(quote needed)*, which has lead to waning trust in science, among scientists as well
as in the general public. One feature of trustworthy scientific work is
reproducability. If a second person with proper knowledge analyses a researchers' data
and finds the same results, this is an indicator that the first scientist did work
correctly. *(or both are wrong ;-) )*  
Here I work on the data Yasunenko et al published in 2012 along with their publication
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
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=C              
##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
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
UNIX command line tool *wget*.

## Preprocessing ##

### Metadata mapping files ###

To categorize and compare genome data, a link between sequenced sample
files and their respective metadata is required. The MG-RAST database [@MG-RAST]
provides a
mapping file between sample data files and unique IDs with incomplete
metadata, while Yatsunenko et al [@Yats12] provide a complete metadata file without
mapping to the sample data files. I applied the "calc" module of the
*LibreOffice* suite to combine the contents of both files to a
complete mapping file. For this I sorted both provided files' contents
by their sample ID strings and copy-pasted missing rows from the
MG-RAST mapping file into the metadata table. Two entries from the
mapping file were not found in the metadata table, and Yatsunenko et al also
stated that two samples could not be used in the analysis, so I
removed them from the metadata-mapping file. The result was saved as a
tab-seperated .csv file. On my system with German environment settings
I had to convert decimal values from comma-seperated values to
international dot-seperated format using the command line tool *sed*:

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
similarity threshold of 97%, as did Yatsunenko et al [@Yats12].
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
a tab-seperated file (\texttt{combined.csv}) for easy use with *R*.

The following steps were all done with QIIME tools. To assess species
richness from the samples, the "rarefaction" technique is used. A
given population (here: microbiome inside a fecal sample) is
subsampled to calculate the overall species richness in the population
while keeping the sample size as small as possible and as big as
necessary. I used *alpha\_rarefaction.py* to create an overview of
suitable rarefaction depths.  
(TODO: nice rarefaction image here)
The read numbers per sample ranged from 305631 to 5826936, with a mean
of 1932291 and a median of 1884081. To guarantee that all samples are
represented and each one is subsampled I chose a rarefaction depth of
290000 (as compared to 290603 in the original publication), which is
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
<!-- usable sample sizes: 188517 377024 -->

### UniFrac distance variation with age ###

[@Yats12] observed the change of microbiome composition from
infant-specific to adult configuration by comparing the composition of
each child's microbiome against the microbiome composition of all adults from the
same country. As it is not completely clear how the distance to all
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
points into a given number k of groups. It minimizes a global distance
function by randomly assigning one data point as starting center for
each group, iteratively swapping points and centers and assigning all
data points to the nearest center.
[@Yats12] used the $\beta$-distance matrix as a one-dimensional
dissimilarity measure for the PAM algorithm and chose k=3 clusters to
refind the samples' countries of origin in the microbiome composition.
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

The plot shows a seperated, banana-shaped distribution of US Americans
and intermingled groups of Malawians and Venezuelans with
concentrations for each subgroups (see fig. 2). Using different
rarefaction tables showed similar results of circa 84 per cent correct
assignment.

### SVM classificator analysis of $\beta$-diversity ###

The PCoA plot (fig. 2) suggested that discrimination of samples by
microbiome diversity is possible. To improve the assignment I used the
implementation of support vector machines (SVM) in the *R* package
*e1071*[@e1071]. Support vector machines classify multi-dimensional datasets
by finding a hyperplane that seperates the classes among (a subset
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
## Warning: Variable(s) 'X288706' and 'X323206' and 'X270447' and 'X297847' and 'X260290' and 'X289202' and 'X259498' and 'X105272' and 'X208044' and 'X302376' and 'X100680' and 'X293997' and 'X319810' and 'X235381' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X544894' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X578230' and 'X110034' and 'X594119' and 'X563315' and 'X2879' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X533807' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X15059' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X11113' and 'X93761' and 'X218435' and 'X49369' and 'X251180' and 'X327084' and 'X200337' and 'X362128' and 'X549284' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X228694' and 'X220924' and 'X152276' and 'X307692' and 'X156895' and 'X139346' and 'X360415' and 'X330873' and 'X253054' and 'X95477' and 'X239678' and 'X529472' and 'X589923' and 'X228065' and 'X150663' and 'X251317' and 'X264560' and 'X167112' and 'X142457' and 'X292384' and 'X29618' and 'X34959' and 'X180930' and 'X16185' and 'X571847' and 'X294146' and 'X584326' and 'X352502' and 'X83477' and 'X530869' and 'X339847' and 'X252964' and 'X16301' and 'X556185' and 'X269172' and 'X114706' and 'X256977' and 'X188183' and 'X191528' and 'X562530' and 'X70346' and 'X203489' and 'X560552' and 'X281004' and 'X279351' and 'X326879' and 'X150918' and 'X270215' and 'X295364' and 'X138760' and 'X223975' and 'X245248' and 'X144084' and 'X328536' and 'X574241' and 'X264539' and 'X355836' and 'X9502' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X297847' and 'X524900' and 'X270447' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X457990' and 'X299900' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X265676' and 'X558276' and 'X50040' and 'X9349' and 'X3861' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X298531' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X533537' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X163461' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X550996' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X138036' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X177331' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X14540' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X527507' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X270840' and 'X559587' and 'X150021' and 'X2605' and 'X22928' and 'X583673' and 'X270263' and 'X28165' and 'X64871' and 'X291843' and 'X265983' and 'X109798' and 'X299309' and 'X59529' and 'X277106' and 'X243412' and 'X137367' and 'X250563' and 'X511130' and 'X136529' and 'X546063' and 'X261931' and 'X358964' and 'X533078' and 'X35642' and 'X201190' and 'X221705' and 'X233334' and 'X19846' and 'X141191' and 'X142181' and 'X352269' and 'X225145' and 'X436019' and 'X111945' and 'X534964' and 'X383716' and 'X244390' and 'X14284' and 'X223706' and 'X326600' and 'X369119' and 'X247268' and 'X242096' and 'X168183' and 'X249191' and 'X571956' and 'X164506' and 'X526198' and 'X308932' and 'X484437' and 'X358429' and 'X593325' and 'X252984' and 'X415943' and 'X566333' and 'X557691' and 'X101880' and 'X321738' and 'X71891' and 'X287452' and 'X327712' and 'X222507' and 'X494503' and 'X532885' and 'X254162' and 'X509986' and 'X587498' and 'X346479' and 'X193662' and 'X593400' and 'X557465' and 'X311867' and 'X218620' and 'X157893' and 'X166436' and 'X256598' and 'X513265' and 'X115015' and 'X8717' and 'X402651' and 'X136243' and 'X67654' and 'X6953' and 'X189488' and 'X589725' and 'X561739' and 'X62675' and 'X343091' and 'X533011' and 'X560114' and 'X110474' and 'X560962' and 'X581039' and 'X247473' and 'X426848' and 'X72805' and 'X91993' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X208044' and 'X302376' and 'X100680' and 'X44758' and 'X259498' and 'X268454' and 'X293997' and 'X336333' and 'X171239' and 'X532384' and 'X319810' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X2879' and 'X550996' and 'X232356' and 'X579422' and 'X14507' and 'X301405' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X323964' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X556604' and 'X312833' and 'X548603' and 'X296798' and 'X457371' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X179040' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X42877' and 'X573052' and 'X1879' and 'X255234' and 'X296146' and 'X328422' and 'X141548' and 'X288438' and 'X7553' and 'X563303' and 'X519353' and 'X183478' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X347809' and 'X276889' and 'X450645' and 'X539260' and 'X585106' and 'X105304' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X560064' and 'X140377' and 'X229171' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X201726' and 'X70804' and 'X112746' and 'X421948' and 'X538707' and 'X428705' and 'X46393' and 'X201491' and 'X537651' and 'X269541' and 'X514611' and 'X160679' and 'X347159' and 'X288049' and 'X148392' and 'X299370' and 'X155884' and 'X269994' and 'X71248' and 'X209866' and 'X593970' and 'X185752' and 'X309361' and 'X17434' and 'X349548' and 'X162964' and 'X545997' and 'X315609' and 'X222566' and 'X11301' and 'X529406' and 'X561309' and 'X303244' and 'X114813' and 'X147281' and 'X291934' and 'X131852' and 'X198983' and 'X266941' and 'X133105' and 'X32167' and 'X583557' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X336333' and 'X301405' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X558276' and 'X50040' and 'X535573' and 'X592633' and 'X9349' and 'X584933' and 'X298531' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X579422' and 'X288438' and 'X14507' and 'X299008' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X13850' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X244792' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X584202' and 'X296146' and 'X328422' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X108861' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X22981' and 'X515116' and 'X303242' and 'X383023' and 'X270739' and 'X340680' and 'X554797' and 'X184695' and 'X405046' and 'X538856' and 'X128236' and 'X339015' and 'X14021' and 'X105304' and 'X2599' and 'X193229' and 'X262102' and 'X293284' and 'X539124' and 'X71267' and 'X359367' and 'X108677' and 'X240912' and 'X120989' and 'X223578' and 'X261745' and 'X594373' and 'X243570' and 'X255713' and 'X219963' and 'X212652' and 'X558390' and 'X336584' and 'X201417' and 'X161224' and 'X536865' and 'X65457' and 'X316721' and 'X579256' and 'X571700' and 'X179995' and 'X564025' and 'X215375' and 'X33613' and 'X249078' and 'X593232' and 'X250043' and 'X509259' and 'X217748' and 'X569476' and 'X304694' and 'X76700' and 'X251955' and 'X238490' and 'X351128' and 'X64871' and 'X291843' and 'X265983' and 'X320378' and 'X305758' and 'X189891' and 'X271865' and 'X99071' and 'X181996' and 'X538707' and 'X428705' and 'X46393' and 'X201491' and 'X325633' and 'X42372' and 'X360054' and 'X311974' and 'X554372' and 'X593617' and 'X250563' and 'X511130' and 'X136529' and 'X546063' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X242915' and 'X125235' and 'X206709' and 'X32016' and 'X540672' and 'X143227' and 'X583489' and 'X314357' and 'X265971' and 'X524595' and 'X32247' and 'X244056' and 'X166076' and 'X512939' and 'X114159' and 'X143097' and 'X548602' and 'X562772' and 'X45520' and 'X130481' and 'X235519' and 'X243468' and 'X257403' and 'X255454' and 'X218287' and 'X568948' and 'X191316' and 'X179854' and 'X44311' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X208044' and 'X319810' and 'X293997' and 'X44758' and 'X259498' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X302376' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X100680' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X342685' and 'X584933' and 'X298531' and 'X332204' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X114748' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X179040' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X190646' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X584202' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X14507' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X143735' and 'X252834' and 'X34222' and 'X105304' and 'X528107' and 'X12486' and 'X296516' and 'X169815' and 'X128261' and 'X355777' and 'X577176' and 'X310781' and 'X171914' and 'X544546' and 'X312833' and 'X548603' and 'X296798' and 'X2599' and 'X14648' and 'X245864' and 'X291317' and 'X263804' and 'X257022' and 'X78216' and 'X273134' and 'X356889' and 'X223711' and 'X169060' and 'X200998' and 'X287758' and 'X107742' and 'X131214' and 'X91691' and 'X46239' and 'X543710' and 'X541375' and 'X26896' and 'X311254' and 'X12296' and 'X102242' and 'X210760' and 'X64871' and 'X291843' and 'X265983' and 'X319790' and 'X543389' and 'X109915' and 'X350331' and 'X114455' and 'X567966' and 'X217815' and 'X325633' and 'X42372' and 'X360054' and 'X311974' and 'X276258' and 'X130529' and 'X70671' and 'X572313' and 'X531776' and 'X398839' and 'X19216' and 'X333033' and 'X11301' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289202' and 'X289208' and 'X97369' and 'X42642' and 'X242915' and 'X125235' and 'X179854' and 'X44311' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X77904' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X297964' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X138036' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X1879' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X110737' and 'X51161' and 'X141739' and 'X14282' and 'X576671' and 'X301924' and 'X149199' and 'X539762' and 'X262102' and 'X293284' and 'X539124' and 'X226596' and 'X131845' and 'X230342' and 'X101394' and 'X278947' and 'X2955' and 'X269416' and 'X236335' and 'X246297' and 'X59529' and 'X435512' and 'X201726' and 'X70804' and 'X112746' and 'X255234' and 'X123281' and 'X315994' and 'X102595' and 'X146324' and 'X570614' and 'X329501' and 'X514611' and 'X110434' and 'X396217' and 'X313816' and 'X322850' and 'X303387' and 'X534964' and 'X383716' and 'X244390' and 'X14284' and 'X223706' and 'X326600' and 'X369119' and 'X319241' and 'X338749' and 'X26459' and 'X555365' and 'X2256' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X100680' and 'X263881' and 'X44758' and 'X171239' and 'X293997' and 'X336333' and 'X532384' and 'X319810' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X9349' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X288438' and 'X14507' and 'X301405' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X50040' and 'X255719' and 'X150338' and 'X312275' and 'X216695' and 'X541232' and 'X538085' and 'X584933' and 'X340514' and 'X302376' and 'X332204' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X297964' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X42877' and 'X573052' and 'X255234' and 'X537651' and 'X7553' and 'X563303' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X159833' and 'X29693' and 'X582883' and 'X184872' and 'X312808' and 'X321950' and 'X206189' and 'X310913' and 'X579184' and 'X141618' and 'X572519' and 'X169164' and 'X140787' and 'X161410' and 'X564867' and 'X42640' and 'X470178' and 'X126508' and 'X312833' and 'X548603' and 'X296798' and 'X77904' and 'X262102' and 'X293284' and 'X539124' and 'X180141' and 'X592917' and 'X249280' and 'X167872' and 'X109798' and 'X299309' and 'X514449' and 'X69411' and 'X100108' and 'X588461' and 'X145887' and 'X296146' and 'X141548' and 'X554372' and 'X593617' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X110434' and 'X396217' and 'X313816' and 'X558226' and 'X303552' and 'X73836' and 'X561616' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X298531' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X138036' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X4772' and 'X16057' and 'X35138' and 'X255169' and 'X141313' and 'X342642' and 'X349634' and 'X578119' and 'X143735' and 'X252834' and 'X34222' and 'X470273' and 'X144258' and 'X226596' and 'X131845' and 'X64871' and 'X291843' and 'X265983' and 'X59529' and 'X362169' and 'X255234' and 'X141442' and 'X570614' and 'X329501' and 'X244445' and 'X294007' and 'X101323' and 'X33773' and 'X24308' and 'X172303' and 'X108266' and 'X109334' and 'X60771' and 'X511876' and 'X470688' and 'X261434' and 'X26459' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X208044' and 'X319810' and 'X302376' and 'X293997' and 'X44758' and 'X259498' and 'X532384' and 'X301405' and 'X3838' and 'X336333' and 'X171239' and 'X575799' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X100680' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X50040' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X556604' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X158775' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X183478' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X14507' and 'X136707' and 'X231606' and 'X292053' and 'X266447' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X257151' and 'X423769' and 'X380753' and 'X274299' and 'X292672' and 'X532925' and 'X288119' and 'X68999' and 'X144522' and 'X207007' and 'X251213' and 'X328981' and 'X568580' and 'X518407' and 'X342821' and 'X535088' and 'X103303' and 'X346430' and 'X331156' and 'X11164' and 'X119236' and 'X282134' and 'X295440' and 'X592917' and 'X249280' and 'X521165' and 'X323234' and 'X296496' and 'X134498' and 'X34798' and 'X321220' and 'X15828' and 'X206803' and 'X138206' and 'X132661' and 'X307909' and 'X144291' and 'X27088' and 'X16001' and 'X259903' and 'X70502' and 'X325633' and 'X42372' and 'X360054' and 'X311974' and 'X276037' and 'X517044' and 'X7553' and 'X563303' and 'X146138' and 'X230801' and 'X510484' and 'X332585' and 'X240168' and 'X278889' and 'X276563' and 'X114375' and 'X302115' and 'X590326' and 'X369685' and 'X593970' and 'X185752' and 'X309361' and 'X17434' and 'X349548' and 'X162964' and 'X545997' and 'X315609' and 'X33773' and 'X24308' and 'X172303' and 'X569672' and 'X13837' and 'X173574' and 'X47316' and 'X77837' and 'X220693' and 'X143469' and 'X148281' and 'X172538' and 'X111923' and 'X137761' and 'X13420' and 'X343089' and 'X227919' and 'X221818' and 'X201637' and 'X25868' and 'X152486' and 'X142828' and 'X242857' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289202' and 'X289208' and 'X97369' and 'X42642' and 'X114813' and 'X211106' and 'X15851' and 'X63942' and 'X8067' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X575799' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X3838' and 'X336333' and 'X301405' and 'X279387' and 'X541325' and 'X160031' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X535573' and 'X592633' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X178106' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X50040' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X293006' and 'X312833' and 'X548603' and 'X296798' and 'X77904' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X297964' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X158775' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X115051' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X278501' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X571667' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X266447' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X22981' and 'X515116' and 'X303242' and 'X335755' and 'X6721' and 'X383023' and 'X270739' and 'X340680' and 'X554797' and 'X184695' and 'X128236' and 'X339015' and 'X14021' and 'X295188' and 'X592858' and 'X249748' and 'X15888' and 'X463007' and 'X306735' and 'X5598' and 'X577152' and 'X360102' and 'X594313' and 'X131000' and 'X112161' and 'X509038' and 'X181573' and 'X576232' and 'X198202' and 'X7303' and 'X15828' and 'X206803' and 'X138206' and 'X132661' and 'X307909' and 'X144291' and 'X27088' and 'X362169' and 'X241214' and 'X240850' and 'X220305' and 'X138741' and 'X142028' and 'X230396' and 'X236550' and 'X338628' and 'X167066' and 'X138602' and 'X3538' and 'X14641' and 'X138316' and 'X269218' and 'X163551' and 'X239951' and 'X200217' and 'X578443' and 'X248784' and 'X320812' and 'X247439' and 'X114526' and 'X594040' and 'X81554' and 'X199105' and 'X4416' and 'X379604' and 'X358044' and 'X319537' and 'X513982' and 'X340207' and 'X486153' and 'X151804' and 'X534213' and 'X33768' and 'X187236' and 'X383501' and 'X128201' and 'X574943' and 'X20974' and 'X264010' and 'X509171' and 'X369064' and 'X541286' and 'X304693' and 'X12364' and 'X256308' and 'X345664' and 'X137631' and 'X365101' and 'X139166' and 'X25041' and 'X1422' and 'X277641' and 'X309862' and 'X264742' and 'X99512' and 'X86707' and 'X297307' and 'X327035' and 'X273126' and 'X110434' and 'X396217' and 'X313816' and 'X107665' and 'X8250' and 'X277421' and 'X590820' and 'X533011' and 'X558789' and 'X569672' and 'X13837' and 'X173574' and 'X47316' and 'X77837' and 'X220693' and 'X143469' and 'X148281' and 'X172538' and 'X111923' and 'X137761' and 'X13420' and 'X343089' and 'X227919' and 'X221818' and 'X201637' and 'X25868' and 'X152486' and 'X142828' and 'X242857' and 'X417573' and 'X93749' and 'X251691' and 'X240933' and 'X540317' and 'X64874' and 'X593891' and 'X105192' and 'X393026' and 'X342774' and 'X272809' and 'X312486' and 'X332547' and 'X282458' and 'X175586' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X110080' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X278501' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X327694' and 'X271891' and 'X183753' and 'X570082' and 'X12357' and 'X73352' and 'X293243' and 'X546078' and 'X582179' and 'X254195' and 'X110737' and 'X554916' and 'X217993' and 'X298447' and 'X328695' and 'X18932' and 'X248337' and 'X168054' and 'X131539' and 'X560064' and 'X140377' and 'X293189' and 'X322090' and 'X100110' and 'X518565' and 'X342245' and 'X229452' and 'X287880' and 'X291314' and 'X238970' and 'X335120' and 'X301613' and 'X512820' and 'X6539' and 'X547303' and 'X25532' and 'X142695' and 'X150728' and 'X325633' and 'X42372' and 'X360054' and 'X311974' and 'X570614' and 'X329501' and 'X88727' and 'X518444' and 'X176117' and 'X178223' and 'X220551' and 'X129425' and 'X270244' and 'X114813' and 'X541327' and 'X342790' and 'X137870' and 'X105192' and 'X393026' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X532384' and 'X301405' and 'X336333' and 'X171239' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X114591' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X573426' and 'X537767' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X540458' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X74554' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X174315' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X110737' and 'X15719' and 'X548876' and 'X158930' and 'X220686' and 'X249818' and 'X92717' and 'X252937' and 'X296038' and 'X2489' and 'X292684' and 'X243150' and 'X293189' and 'X322090' and 'X100110' and 'X534168' and 'X559986' and 'X553559' and 'X253774' and 'X341708' and 'X255859' and 'X316413' and 'X559200' and 'X160173' and 'X206963' and 'X22673' and 'X557467' and 'X136472' and 'X217401' and 'X513189' and 'X534321' and 'X525670' and 'X262768' and 'X358635' and 'X368658' and 'X245246' and 'X244867' and 'X551005' and 'X581900' and 'X251806' and 'X250277' and 'X257891' and 'X247445' and 'X9062' and 'X356027' and 'X261663' and 'X244127' and 'X131989' and 'X8146' and 'X221896' and 'X316335' and 'X584477' and 'X305400' and 'X146541' and 'X100737' and 'X111264' and 'X204919' and 'X307581' and 'X352830' and 'X100174' and 'X558023' and 'X254262' and 'X312767' and 'X276703' and 'X521165' and 'X323234' and 'X296496' and 'X134498' and 'X114455' and 'X278947' and 'X2955' and 'X269416' and 'X236335' and 'X246297' and 'X243412' and 'X150728' and 'X16001' and 'X259903' and 'X70502' and 'X255234' and 'X162964' and 'X545997' and 'X315609' and 'X223219' and 'X149938' and 'X261034' and 'X273018' and 'X511346' and 'X265137' and 'X243161' and 'X242145' and 'X334587' and 'X225301' and 'X573041' and 'X19369' and 'X88144' and 'X158376' and 'X223155' and 'X167052' and 'X235601' and 'X512264' and 'X178912' and 'X259537' and 'X74727' and 'X547810' and 'X527168' and 'X213225' and 'X510745' and 'X97664' and 'X205713' and 'X235345' and 'X237616' and 'X236904' and 'X321764' and 'X79280' and 'X359824' and 'X179479' and 'X417385' and 'X226508' and 'X93202' and 'X221405' and 'X236131' and 'X203959' and 'X163173' and 'X562639' and 'X204906' and 'X542078' and 'X270711' and 'X162144' and 'X146126' and 'X216594' and 'X584085' and 'X284643' and 'X217221' and 'X563370' and 'X222512' and 'X235795' and 'X247201' and 'X267825' and 'X8084' and 'X154633' and 'X544615' and 'X122054' and 'X60733' and 'X236657' and 'X548516' and 'X29077' and 'X218840' and 'X192248' and 'X331818' and 'X334500' and 'X340065' and 'X429841' and 'X214806' and 'X254353' and 'X512755' and 'X44239' and 'X21532' and 'X249252' and 'X238020' and 'X549166' and 'X513801' and 'X171640' and 'X318753' and 'X560846' and 'X206488' and 'X128216' and 'X226410' and 'X333608' and 'X91919' and 'X553534' and 'X244549' and 'X292764' and 'X215846' and 'X562144' and 'X140918' and 'X247121' and 'X5207' and 'X243247' and 'X244840' and 'X317483' and 'X86480' and 'X194964' and 'X256945' and 'X349315' and 'X218242' and 'X537673' and 'X510547' and 'X251135' and 'X564413' and 'X137607' and 'X250637' and 'X309465' and 'X141407' and 'X325796' and 'X146531' and 'X172332' and 'X584129' and 'X579734' and 'X100187' and 'X593259' and 'X140763' and 'X554086' and 'X82843' and 'X12970' and 'X547624' and 'X137132' and 'X465713' and 'X511558' and 'X557941' and 'X320892' and 'X493116' and 'X100282' and 'X262115' and 'X358280' and 'X515258' and 'X152567' and 'X302479' and 'X254609' and 'X248044' and 'X510823' and 'X99452' and 'X67358' and 'X77547' and 'X258809' and 'X179461' and 'X147459' and 'X338821' and 'X150' and 'X518474' and 'X142117' and 'X551463' and 'X319318' and 'X217779' and 'X132318' and 'X114803' and 'X565936' and 'X5497' and 'X143178' and 'X247045' and 'X202045' and 'X159766' and 'X195654' and 'X367606' and 'X533060' and 'X397321' and 'X225760' and 'X217413' and 'X290559' and 'X87249' and 'X251112' and 'X213161' and 'X223672' and 'X335630' and 'X107195' and 'X41720' and 'X21795' and 'X239601' and 'X532312' and 'X292243' and 'X266173' and 'X369734' and 'X198983' and 'X266941' and 'X133105' and 'X32167' and 'X583557' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X262741' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X112062' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X103120' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X35138' and 'X255169' and 'X141313' and 'X10302' and 'X101041' and 'X329012' and 'X9565' and 'X312767' and 'X276703' and 'X360102' and 'X594313' and 'X131000' and 'X112161' and 'X28461' and 'X293780' and 'X515774' and 'X142695' and 'X524433' and 'X128098' and 'X256695' and 'X586378' and 'X190023' and 'X161584' and 'X261254' and 'X552887' and 'X236325' and 'X557255' and 'X9458' and 'X15828' and 'X206803' and 'X141707' and 'X138206' and 'X132661' and 'X307909' and 'X144291' and 'X27088' and 'X12618' and 'X149265' and 'X115211' and 'X224427' and 'X49771' and 'X222549' and 'X588053' and 'X539763' and 'X305625' and 'X562121' and 'X249700' and 'X33773' and 'X24308' and 'X172303' and 'X177392' and 'X558226' and 'X144186' and 'X318862' and 'X342774' and 'X272809' and 'X312486' and 'X332547' and 'X282458' and 'X175586' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X297847' and 'X270447' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X310510' and 'X84248' and 'X299008' and 'X558276' and 'X50040' and 'X535573' and 'X592633' and 'X9349' and 'X245157' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X205987' and 'X274257' and 'X509038' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X594359' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X1879' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X571667' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X383023' and 'X270739' and 'X340680' and 'X554797' and 'X184695' and 'X15719' and 'X548876' and 'X158930' and 'X220686' and 'X92717' and 'X252937' and 'X296038' and 'X158223' and 'X301924' and 'X262102' and 'X293284' and 'X539124' and 'X565814' and 'X163340' and 'X91416' and 'X219876' and 'X139040' and 'X587111' and 'X240189' and 'X558389' and 'X581638' and 'X114455' and 'X106394' and 'X550936' and 'X27055' and 'X453710' and 'X100828' and 'X10867' and 'X163140' and 'X589528' and 'X132728' and 'X461524' and 'X256583' and 'X70804' and 'X112746' and 'X471399' and 'X144388' and 'X111455' and 'X421948' and 'X554372' and 'X593617' and 'X269541' and 'X540416' and 'X270553' and 'X554336' and 'X510458' and 'X165982' and 'X335081' and 'X530805' and 'X330208' and 'X572313' and 'X33773' and 'X24308' and 'X172303' and 'X558789' and 'X170079' and 'X560811' and 'X224745' and 'X557840' and 'X139833' and 'X520731' and 'X300862' and 'X555365' and 'X2256' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X351744' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X584933' and 'X298531' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X9349' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X108843' and 'X323595' and 'X574338' and 'X335755' and 'X6721' and 'X16086' and 'X368350' and 'X144163' and 'X255719' and 'X150338' and 'X312275' and 'X101041' and 'X329012' and 'X342719' and 'X181683' and 'X293189' and 'X322090' and 'X100110' and 'X558344' and 'X180141' and 'X310662' and 'X344129' and 'X16915' and 'X217821' and 'X212713' and 'X348047' and 'X255053' and 'X210945' and 'X592917' and 'X249280' and 'X64871' and 'X291843' and 'X265983' and 'X109798' and 'X299309' and 'X229281' and 'X344007' and 'X276037' and 'X517044' and 'X356657' and 'X549030' and 'X305547' and 'X279442' and 'X209377' and 'X169379' and 'X375076' and 'X170555' and 'X179541' and 'X224427' and 'X49771' and 'X222549' and 'X588053' and 'X539763' and 'X305625' and 'X562121' and 'X249700' and 'X147459' and 'X338821' and 'X150' and 'X518474' and 'X142117' and 'X551463' and 'X319318' and 'X217779' and 'X132318' and 'X114803' and 'X565936' and 'X5497' and 'X239946' and 'X352042' and 'X417573' and 'X93749' and 'X251691' and 'X240933' and 'X540317' and 'X342774' and 'X272809' and 'X312486' and 'X332547' and 'X282458' and 'X175586' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X202911' and 'X301405' and 'X336333' and 'X211171' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X341286' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X566717' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X51161' and 'X510182' and 'X141739' and 'X14282' and 'X576671' and 'X12591' and 'X120648' and 'X156429' and 'X518471' and 'X310519' and 'X89580' and 'X298447' and 'X128236' and 'X339015' and 'X14021' and 'X15719' and 'X548876' and 'X158930' and 'X220686' and 'X92717' and 'X252937' and 'X296038' and 'X528107' and 'X12486' and 'X296516' and 'X169815' and 'X128261' and 'X355777' and 'X577176' and 'X565814' and 'X163340' and 'X91416' and 'X219876' and 'X139040' and 'X587111' and 'X240189' and 'X558389' and 'X581638' and 'X75475' and 'X78150' and 'X331156' and 'X11164' and 'X119236' and 'X282134' and 'X295440' and 'X142695' and 'X250514' and 'X255802' and 'X210792' and 'X384373' and 'X575450' and 'X332574' and 'X556036' and 'X532912' and 'X533602' and 'X548891' and 'X99509' and 'X123738' and 'X543676' and 'X552729' and 'X487233' and 'X237879' and 'X238675' and 'X9994' and 'X453756' and 'X346926' and 'X325677' and 'X336984' and 'X112916' and 'X217403' and 'X89541' and 'X141394' and 'X219106' and 'X148662' and 'X154113' and 'X217798' and 'X324704' and 'X241908' and 'X554341' and 'X548878' and 'X218845' and 'X73544' and 'X549839' and 'X178178' and 'X334226' and 'X7598' and 'X239195' and 'X410855' and 'X251499' and 'X29284' and 'X227361' and 'X39461' and 'X207957' and 'X130709' and 'X240983' and 'X154220' and 'X316541' and 'X537483' and 'X3094' and 'X160115' and 'X242284' and 'X511637' and 'X111394' and 'X91588' and 'X244970' and 'X55946' and 'X5163' and 'X237583' and 'X591751' and 'X58634' and 'X193867' and 'X517065' and 'X261978' and 'X525666' and 'X585097' and 'X534483' and 'X113942' and 'X239741' and 'X142304' and 'X76549' and 'X255914' and 'X515421' and 'X97778' and 'X576689' and 'X140788' and 'X154684' and 'X584480' and 'X98983' and 'X529793' and 'X164098' and 'X583293' and 'X239095' and 'X244548' and 'X362116' and 'X237095' and 'X329486' and 'X137809' and 'X319977' and 'X277776' and 'X250987' and 'X147665' and 'X141110' and 'X558297' and 'X225887' and 'X142173' and 'X143400' and 'X144142' and 'X410399' and 'X148779' and 'X410272' and 'X147943' and 'X149631' and 'X105233' and 'X314203' and 'X309588' and 'X221399' and 'X225933' and 'X18245' and 'X565593' and 'X110332' and 'X164828' and 'X207298' and 'X86557' and 'X322850' and 'X303387' and 'X107241' and 'X270993' and 'X246533' and 'X547346' and 'X207832' and 'X347680' and 'X222171' and 'X138500' and 'X541327' and 'X342790' and 'X137870' and 'X6144' and 'X534993' and 'X271528' and 'X125257' and 'X74125' and 'X77568' and 'X242049' and 'X185640' and 'X267323' and 'X583979' and 'X246762' and 'X14559' and 'X64356' and 'X246819' and 'X584366' and 'X77057' and 'X329045' and 'X11428' and 'X52053' and 'X583887' and 'X533172' and 'X238291' and 'X266522' and 'X311345' and 'X77712' and 'X586961' and 'X508084' and 'X247938' and 'X263760' and 'X99041' and 'X564501' and 'X517399' and 'X136867' and 'X322801' and 'X153617' and 'X247875' and 'X343296' and 'X554990' and 'X138031' and 'X329326' and 'X354851' and 'X103895' and 'X350772' and 'X248695' and 'X9061' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X293997' and 'X302376' and 'X100680' and 'X44758' and 'X171239' and 'X532384' and 'X202911' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X299008' and 'X558276' and 'X50040' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X19141' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X1879' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X91993' and 'X255719' and 'X9349' and 'X150338' and 'X312275' and 'X193229' and 'X328981' and 'X568580' and 'X518407' and 'X310662' and 'X344129' and 'X16915' and 'X217821' and 'X212713' and 'X348047' and 'X255053' and 'X217815' and 'X586378' and 'X161584' and 'X261254' and 'X552887' and 'X236325' and 'X557255' and 'X9458' and 'X243412' and 'X70804' and 'X112746' and 'X87601' and 'X131680' and 'X312294' and 'X350384' and 'X68348' and 'X15238' and 'X79242' and 'X269541' and 'X261931' and 'X358964' and 'X533078' and 'X35642' and 'X201190' and 'X221705' and 'X233334' and 'X19846' and 'X141191' and 'X142181' and 'X352269' and 'X225145' and 'X289071' and 'X164832' and 'X66924' and 'X239441' and 'X246316' and 'X352244' and 'X241392' and 'X107241' and 'X270993' and 'X246533' and 'X555365' and 'X2256' and 'X310510' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X202911' and 'X301405' and 'X336333' and 'X9249' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X9349' and 'X584933' and 'X332204' and 'X2879' and 'X550996' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X558276' and 'X50040' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X457371' and 'X594119' and 'X563315' and 'X573426' and 'X537767' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X231593' and 'X277624' and 'X510498' and 'X573052' and 'X1879' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X335755' and 'X6721' and 'X154459' and 'X89195' and 'X316587' and 'X574741' and 'X558408' and 'X298447' and 'X128236' and 'X339015' and 'X14021' and 'X143735' and 'X252834' and 'X34222' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X587804' and 'X109305' and 'X4551' and 'X557327' and 'X247484' and 'X564264' and 'X103425' and 'X534168' and 'X559986' and 'X553559' and 'X253774' and 'X341708' and 'X255859' and 'X316413' and 'X559200' and 'X160173' and 'X206963' and 'X22673' and 'X557467' and 'X136472' and 'X217401' and 'X513189' and 'X534321' and 'X525670' and 'X262768' and 'X358635' and 'X368658' and 'X245246' and 'X244867' and 'X551005' and 'X251806' and 'X250277' and 'X9062' and 'X356027' and 'X261663' and 'X244127' and 'X131989' and 'X8146' and 'X316335' and 'X584477' and 'X305400' and 'X146541' and 'X100737' and 'X111264' and 'X204919' and 'X307581' and 'X352830' and 'X100174' and 'X558023' and 'X254262' and 'X10390' and 'X207629' and 'X109798' and 'X299309' and 'X278947' and 'X2955' and 'X269416' and 'X236335' and 'X246297' and 'X16001' and 'X259903' and 'X70502' and 'X361428' and 'X70804' and 'X112746' and 'X229281' and 'X344007' and 'X87601' and 'X131680' and 'X312294' and 'X350384' and 'X269541' and 'X239441' and 'X246316' and 'X352244' and 'X241392' and 'X239946' and 'X352042' and 'X190531' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X270447' and 'X297847' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X532384' and 'X336333' and 'X171239' and 'X301405' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X341286' and 'X584933' and 'X332204' and 'X457371' and 'X2879' and 'X550996' and 'X140904' and 'X232356' and 'X42877' and 'X288438' and 'X14507' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X111824' and 'X104752' and 'X594359' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X231593' and 'X277624' and 'X510498' and 'X573052' and 'X255234' and 'X296146' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X22981' and 'X515116' and 'X303242' and 'X89580' and 'X298447' and 'X15719' and 'X548876' and 'X158930' and 'X220686' and 'X92717' and 'X252937' and 'X296038' and 'X143735' and 'X252834' and 'X34222' and 'X560064' and 'X140377' and 'X301816' and 'X360102' and 'X594313' and 'X131000' and 'X112161' and 'X586378' and 'X161584' and 'X261254' and 'X552887' and 'X236325' and 'X557255' and 'X9458' and 'X250514' and 'X255802' and 'X210792' and 'X384373' and 'X575450' and 'X332574' and 'X556036' and 'X532912' and 'X533602' and 'X548891' and 'X99509' and 'X123738' and 'X543676' and 'X552729' and 'X487233' and 'X237879' and 'X238675' and 'X453756' and 'X346926' and 'X325677' and 'X336984' and 'X112916' and 'X217403' and 'X89541' and 'X141394' and 'X219106' and 'X148662' and 'X154113' and 'X217798' and 'X324704' and 'X241908' and 'X554341' and 'X548878' and 'X218845' and 'X73544' and 'X549839' and 'X178178' and 'X334226' and 'X7598' and 'X239195' and 'X410855' and 'X251499' and 'X29284' and 'X227361' and 'X39461' and 'X207957' and 'X130709' and 'X240983' and 'X154220' and 'X316541' and 'X537483' and 'X3094' and 'X160115' and 'X242284' and 'X511637' and 'X111394' and 'X91588' and 'X244970' and 'X55946' and 'X5163' and 'X237583' and 'X591751' and 'X58634' and 'X193867' and 'X517065' and 'X261978' and 'X525666' and 'X585097' and 'X534483' and 'X113942' and 'X239741' and 'X142304' and 'X76549' and 'X255914' and 'X515421' and 'X97778' and 'X576689' and 'X140788' and 'X154684' and 'X584480' and 'X98983' and 'X529793' and 'X164098' and 'X583293' and 'X239095' and 'X244548' and 'X362116' and 'X237095' and 'X329486' and 'X137809' and 'X319977' and 'X277776' and 'X250987' and 'X147665' and 'X15828' and 'X206803' and 'X138206' and 'X132661' and 'X307909' and 'X144291' and 'X27088' and 'X540416' and 'X270553' and 'X554336' and 'X510458' and 'X165982' and 'X335081' and 'X530805' and 'X330208' and 'X572313' and 'X162964' and 'X545997' and 'X315609' and 'X239946' and 'X352042' and 'X533011' and 'X417573' and 'X93749' and 'X251691' and 'X240933' and 'X540317' and 'X555365' and 'X2256' and 'X541327' and 'X342790' and 'X137870' and 'X190531' constant. Cannot scale data.
## Warning: Variable(s) 'X288706' and 'X297847' and 'X270447' and 'X524900' and 'X289202' and 'X259498' and 'X208044' and 'X319810' and 'X302376' and 'X100680' and 'X293997' and 'X44758' and 'X171239' and 'X532384' and 'X301405' and 'X336333' and 'X279387' and 'X541325' and 'X75333' and 'X14577' and 'X162964' and 'X545997' and 'X307781' and 'X14974' and 'X300862' and 'X310510' and 'X299008' and 'X558276' and 'X50040' and 'X9349' and 'X341286' and 'X584933' and 'X332204' and 'X2879' and 'X550996' and 'X42877' and 'X288438' and 'X14507' and 'X22981' and 'X515116' and 'X303242' and 'X89402' and 'X316587' and 'X574741' and 'X558408' and 'X255719' and 'X150338' and 'X312275' and 'X541232' and 'X538085' and 'X312833' and 'X548603' and 'X296798' and 'X457371' and 'X594119' and 'X563315' and 'X511480' and 'X12561' and 'X86772' and 'X539421' and 'X222685' and 'X167766' and 'X549261' and 'X302793' and 'X355912' and 'X349669' and 'X246206' and 'X6697' and 'X307844' and 'X222843' and 'X234716' and 'X336615' and 'X296409' and 'X293147' and 'X137539' and 'X198966' and 'X344208' and 'X73120' and 'X303173' and 'X296671' and 'X32126' and 'X356039' and 'X274257' and 'X294146' and 'X584326' and 'X352502' and 'X313991' and 'X291304' and 'X585970' and 'X289851' and 'X188772' and 'X381694' and 'X591496' and 'X587586' and 'X555754' and 'X149723' and 'X151714' and 'X338567' and 'X232356' and 'X207813' and 'X103103' and 'X333505' and 'X104511' and 'X150728' and 'X318733' and 'X227017' and 'X277624' and 'X510498' and 'X573052' and 'X103842' and 'X255234' and 'X141548' and 'X537651' and 'X7553' and 'X563303' and 'X160679' and 'X288049' and 'X148392' and 'X299370' and 'X519353' and 'X244196' and 'X536014' and 'X320540' and 'X295095' and 'X556428' and 'X266368' and 'X12538' and 'X542417' and 'X31334' and 'X351663' and 'X306508' and 'X548855' and 'X325318' and 'X58790' and 'X136707' and 'X315609' and 'X231606' and 'X292053' and 'X101323' and 'X34513' and 'X565390' and 'X584810' and 'X557184' and 'X131721' and 'X266381' and 'X531776' and 'X177473' and 'X101249' and 'X272478' and 'X291838' and 'X46015' and 'X78125' and 'X146014' and 'X551873' and 'X245596' and 'X232231' and 'X379056' and 'X355680' and 'X353125' and 'X512892' and 'X31669' and 'X548932' and 'X294610' and 'X289208' and 'X97369' and 'X42642' and 'X223318' and 'X43339' and 'X324275' and 'X33780' and 'X318963' and 'X289424' and 'X91993' and 'X144753' and 'X278545' and 'X562002' and 'X275085' and 'X16026' and 'X295250' and 'X110034' and 'X560064' and 'X140377' and 'X558344' and 'X310662' and 'X344129' and 'X16915' and 'X217821' and 'X212713' and 'X348047' and 'X255053' and 'X463007' and 'X306735' and 'X5598' and 'X577152' and 'X229171' and 'X111824' and 'X104752' and 'X243412' and 'X250514' and 'X255802' and 'X210792' and 'X384373' and 'X575450' and 'X332574' and 'X556036' and 'X532912' and 'X533602' and 'X548891' and 'X99509' and 'X123738' and 'X543676' and 'X348570' and 'X552729' and 'X487233' and 'X237879' and 'X238675' and 'X151439' and 'X453756' and 'X346926' and 'X325677' and 'X336984' and 'X112916' and 'X217403' and 'X89541' and 'X141394' and 'X219106' and 'X148662' and 'X154113' and 'X217798' and 'X324704' and 'X241908' and 'X554341' and 'X548878' and 'X218845' and 'X73544' and 'X549839' and 'X178178' and 'X334226' and 'X7598' and 'X239195' and 'X410855' and 'X251499' and 'X218626' and 'X29284' and 'X227361' and 'X39461' and 'X207957' and 'X130709' and 'X240983' and 'X154220' and 'X316541' and 'X537483' and 'X3094' and 'X160115' and 'X242284' and 'X511637' and 'X111394' and 'X91588' and 'X244970' and 'X55946' and 'X5163' and 'X237583' and 'X591751' and 'X58634' and 'X193867' and 'X517065' and 'X261978' and 'X525666' and 'X585097' and 'X534483' and 'X113942' and 'X239741' and 'X142304' and 'X76549' and 'X255914' and 'X515421' and 'X97778' and 'X256584' and 'X576689' and 'X140788' and 'X154684' and 'X584480' and 'X98983' and 'X529793' and 'X164098' and 'X583293' and 'X239095' and 'X244548' and 'X362116' and 'X237095' and 'X329486' and 'X137809' and 'X319977' and 'X277776' and 'X230389' and 'X250987' and 'X147665' and 'X277460' and 'X272752' and 'X101144' and 'X355736' and 'X561950' and 'X71049' and 'X108408' and 'X114384' and 'X539301' and 'X206221' and 'X103476' and 'X353702' and 'X540454' and 'X33643' and 'X252471' and 'X446356' and 'X242798' and 'X296146' and 'X241214' and 'X240850' and 'X220305' and 'X138741' and 'X142028' and 'X230396' and 'X236550' and 'X338628' and 'X167066' and 'X138602' and 'X3538' and 'X14641' and 'X138316' and 'X269218' and 'X163551' and 'X239951' and 'X200217' and 'X578443' and 'X248784' and 'X320812' and 'X247439' and 'X114526' and 'X594040' and 'X81554' and 'X199105' and 'X4416' and 'X379604' and 'X358044' and 'X319537' and 'X513982' and 'X340207' and 'X486153' and 'X151804' and 'X534213' and 'X187236' and 'X383501' and 'X128201' and 'X20974' and 'X264010' and 'X509171' and 'X369064' and 'X541286' and 'X304693' and 'X12364' and 'X256308' and 'X345664' and 'X137631' and 'X365101' and 'X139166' and 'X25041' and 'X1422' and 'X277641' and 'X309862' and 'X264742' and 'X146324' and 'X322850' and 'X303387' and 'X263805' and 'X239946' and 'X352042' and 'X107241' and 'X270993' and 'X246533' and 'X560114' and 'X110474' and 'X560962' and 'X581039' and 'X247473' and 'X100104' and 'X134126' and 'X250119' and 'X169856' and 'X221734' and 'X587796' and 'X270852' constant. Cannot scale data.
```

```r
boxplot(strength, xlab="number of features", ylab="prediction success", axes=FALSE)
abline(h=pam.result, col="blue")
axis(side=1, at=1:length(featureSizes), labels=as.integer(featureSizes), las=2)
axis(side=2, at=seq(from=0, to=1.2, by=0.2))
```

![plot of chunk Classifiers](figure/Classifiers.png) 
<!-- all hail the cache=TRUE parameter! -->

Figure 4 shows barplots of SVM-classifier prediction of the 16S dataset. X-axis
describes the number of features used. The features were selected by greatest
variance. A blue horizontal line shows the prediction success of PAM clustering for
comparison.

### Bacterial diversity with age ###

Yatsunenko et al found that the number of OTUs inside the fecal samples
increased with age[@Yats12]. To verify this I calculated the means of OTU
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

Discussion
==========

Feel free to speak your mind.

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
tools, version 1.8.0 as described by Caporaso et al[@QIIME].  
Data analysis and visualization was done using the *R* statistical sofware
version 3.1.0 [-@R], the package *cluster* by Maechler et al[@cluster] for
partitioning around medoids
and the package *e1071* by Meyer et al[@e1071] for support vector machine
classifiers.  
The custom Python scripts mentioned in the scope of this document can
be found on
[my github page](https://github.com/hermann-p/yatsunenko-2012-microbiome).

References
==========

<!-- auto-filled by pandoc-citeproc -->
