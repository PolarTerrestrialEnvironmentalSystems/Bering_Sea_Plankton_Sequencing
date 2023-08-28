setwd("")

########################################################################################################################################################
####################### 1 Filtering ObiTools Output  ###################################################################################################
########################################################################################################################################################

rm(list = ls())

dir.create("Output")
dir.create("Output/01_Filtering")

library("tidyverse")
library("readxl")
library("openxlsx")

# after running the ObiTools pipeline, the metabarcoding data are filtered.
# several consecutive filter steps are applied:
# (1) ASVs are assigned a taxonomic name based on 96-100% similarity to an entry in the reference database,
# (2) ASVs are represented with at least 100 read counts in total and
# (3) ASVs have at least 10 read counts per PCR-product,
# (4) ASVs are present at least 3 times among all PCR-products,
# (5) ASVs show taxonomic resolution at least to phylum level "Bacillariophyta",
# (6) ASVs are tagged as "internal" by obiclean in less than 50% of the different replicates per sample.

# load data
diatoms <- read.xlsx("APMG-41_embl143_rbcL.xlsx")

# overview dataset
names(diatoms)
names(diatoms)[1:10]

diatoms_onlycounts <- as.data.frame(sapply(diatoms[, 4:205], as.numeric)) # counts numeric

# create a data frame with counts in numeric, rest in characters
diatoms_classes <- cbind(diatoms[, 1:3], diatoms_onlycounts, diatoms[, 206:ncol(diatoms)])

names(diatoms)
ASV_sum <- rowSums(diatoms_classes[, c(4:205)], na.rm = TRUE)
raw_sum <- sum(ASV_sum)
raw_sum # total counts prior to filtering

### reads in samples and blanks
samples_blanks <- diatoms_classes[, c(4:205)]
samples_blanks$ASV <- diatoms_classes$NUC_SEQ
ncol(samples_blanks)
samples_blanks <- samples_blanks[, c(203, 1:202)]
samples_blanks <- as.data.frame(t(samples_blanks))
colnames(samples_blanks) <- samples_blanks[1, ]
samples_blanks <- samples_blanks[-1,]
samples_blanks$sample_name <- rownames(samples_blanks)
# load sample_names
sample_names <- read.xlsx("Sample_names.xlsx")
samples_blanks <- merge(samples_blanks, sample_names, by = "sample_name")
# raw counts samples
samples_counts <- subset(samples_blanks, subset = !age %in% c("EB", "NTC"))
ncol(samples_counts)
samples_counts <- samples_counts[, 2:12124]
samples_counts <- as.data.frame(sapply(samples_counts, as.numeric)) # counts numeric
samples_counts[is.na(samples_counts)] <- 0 # replace NAs by 0
samples_counts <- colSums(samples_counts)
tot_sample_counts <- sum(samples_counts) # total raw reads in samples
(tot_sample_counts/raw_sum)*100

# raw counts extraction blanks
EB_counts <- subset(samples_blanks, subset = age %in% "EB")
ncol(EB_counts)
EB_counts <- EB_counts[, 2:12124]
EB_counts <- as.data.frame(sapply(EB_counts, as.numeric)) # counts numeric
EB_counts[is.na(EB_counts)] <- 0 # replace NAs by 0
EB_counts <- colSums(EB_counts)
tot_EB_counts <- sum(EB_counts) # total raw reads in samples
(tot_EB_counts/raw_sum)*100

# raw counts PCR blanks
NTC_counts <- subset(samples_blanks, subset = age %in% "NTC")
ncol(NTC_counts)
NTC_counts <- NTC_counts[, 2:12124]
NTC_counts <- as.data.frame(sapply(NTC_counts, as.numeric)) # counts numeric
NTC_counts[is.na(NTC_counts)] <- 0 # replace NAs by 0
NTC_counts <- colSums(NTC_counts)
tot_NTC_counts <- sum(NTC_counts) # total raw reads in samples
(tot_NTC_counts/raw_sum)*100

#############################################
###       Filter 1: Best Indentity        ###
#############################################

# filter for best identity >= 0.96
id96_filter1 <- subset(diatoms_classes, subset = BEST_IDENTITY >=0.96)
ASV_sum_filtered1 <- rowSums(id96_filter1[, c(4:205)], na.rm = TRUE)
filtered1_sum <- sum(ASV_sum_filtered1)

# calculate percent of filtered counts to total counts of raw data
(filtered1_sum/raw_sum)*100


#############################################
###        Filter 2: Total counts         ###
#############################################

# counts of ASV over all samples >= 100
id96_filter2 <- subset(id96_filter1, subset = rowSums(id96_filter1[, c(4:205)], na.rm = TRUE) >= 100)

# percent to raw data
ASV_sum_filtered2 <- rowSums(id96_filter2[, c(4:205)], na.rm = TRUE)
filtered2_sum <- sum(ASV_sum_filtered2)
(filtered2_sum/raw_sum)*100
# percent to filter 1
(filtered2_sum/filtered1_sum)*100


###################################################
### Filter 3: Minimum Read Counts per replicate ###
###################################################

# minimum read count per replicate = 10
# otherwise replace by 0

# replace NA by 0
id96_filter3 <- id96_filter2
id96_filter3[is.na(id96_filter3)] <- 0
# replace counts < 10 by 0 
id96_filter3[,c(4:205)][id96_filter3[,c(4:205)] < 10] <- 0

# replace 0 by NA again (important for next filtering step)
id96_filter3[id96_filter3 == 0] <- NA

# percent to raw data
ASV_sum_filtered3 <- rowSums(id96_filter3[, c(4:205)], na.rm = TRUE)
filtered3_sum <- sum(ASV_sum_filtered3)
(filtered3_sum/raw_sum)*100
# percent to filter 1
(filtered3_sum/filtered1_sum)*100
# percent to filter 2
(filtered3_sum/filtered2_sum)*100


##################################################
### Filter 4: Presence over all PCR replicates ###
##################################################

# number of samples
ncol(id96_filter3[, c(4:205)]) # = 202
# --> in at least 198 samples no NA

# occurrence in at least 3 of all replicates 
id96_filter4 <- subset(id96_filter3, subset = rowSums(is.na(id96_filter3[, c(4:205)])) < 199)

# percent to raw data
ASV_sum_filtered4 <- rowSums(id96_filter4[, c(4:205)], na.rm = TRUE)
filtered4_sum <- sum(ASV_sum_filtered4)
(filtered4_sum/raw_sum)*100
# percent to filter 1
(filtered4_sum/filtered1_sum)*100
# percent to filter 2
(filtered4_sum/filtered2_sum)*100
# percent to filter 3
(filtered4_sum/filtered3_sum)*100


##################################################
### Filter 5: Classification to Bacillariophta ###
##################################################

# taxonomic assignment to at least level to division Bacillariophyta

#  a list with all taxa that occur in the data set after filter 3
# and do not belong to the Bacillariophyta is loaded

# load in list with non-Bacillariophyta taxa
nonbac <- read_delim("Not_Bacillariophyceae.txt", delim = "\t", col_names = F)
nonbac_v <- unlist(nonbac) # make vector
is.vector(nonbac_v) # check if it is a vector

id96_filter5 <- subset(id96_filter4, subset = !SCIENTIFIC_NAME %in% nonbac_v)

# percent to raw data
ASV_sum_filtered5 <- rowSums(id96_filter5[, c(4:205)], na.rm = TRUE)
filtered5_sum <- sum(ASV_sum_filtered5)
(filtered5_sum/raw_sum)*100
# percent to filter 1
(filtered5_sum/filtered1_sum)*100
# percent to filter 2
(filtered5_sum/filtered2_sum)*100
# percent to filter 3
(filtered5_sum/filtered3_sum)*100
# percent to filter 4
(filtered5_sum/filtered4_sum)*100


#############################################
###   Filter 6: Not only internal reads   ###
#############################################

# filter for all ASVs that don't occur as headcounts at all
#and where the internal counts make up more than 50% of the total counts
# that means: no headcounts, more internal counts than singleton counts
id96_filter5$obiclean_headcount <- as.numeric(id96_filter5$obiclean_headcount)
id96_filter5$obiclean_internalcount <- as.numeric(id96_filter5$obiclean_internalcount)
id96_filter5$obiclean_singletoncount <- as.numeric(id96_filter5$obiclean_singletoncount)

temp <- subset(id96_filter5, subset = is.na(obiclean_headcount) & obiclean_internalcount >= obiclean_singletoncount)
temp_seq <- as.vector(temp$NUC_SEQ) # make ASVs that should be filtered out if the data set a vetor

id96_filter6 <- subset(id96_filter5, subset = !NUC_SEQ %in% temp_seq)

# percent to raw data
ASV_sum_filtered6 <- rowSums(id96_filter6[, c(4:205)], na.rm = TRUE)
filtered6_sum <- sum(ASV_sum_filtered6)
(filtered6_sum/raw_sum)*100
# percent to filter 1
(filtered6_sum/filtered1_sum)*100
# percent to filter 2
(filtered6_sum/filtered2_sum)*100
# percent to filter 3
(filtered6_sum/filtered3_sum)*100
# percent to filter 4
(filtered6_sum/filtered4_sum)*100
# percent to filter 5
(filtered6_sum/filtered5_sum)*100


#############################################
###           OVERVIEW DATASET            ###
#############################################

# family level
family <- subset(id96_filter6, subset = !is.na(family_name)) # only sequence types assigned at least up to family level
nrow(family) # unique ASVs
ASV_sum_family <- rowSums(family[, c(4:205)], na.rm = TRUE)
family_sum <- sum(ASV_sum_family)
family_sum # counts assigned up to family level
(family_sum/raw_sum)*100 # % of raw data assigned at least to family level 

# genus level
genus <- subset(id96_filter6, subset = !is.na(genus_name)) # only sequence types assigned at least up to genus level
nrow(genus) # unique ASVs
ASV_sum_genus <- rowSums(genus[, c(4:205)], na.rm = TRUE)
genus_sum <- sum(ASV_sum_genus)
genus_sum # counts assigned up to genus level
(genus_sum/raw_sum)*100 # % of raw data assigned at least to genus level 

# species level
species <- subset(id96_filter6, subset = !is.na(species_name)) # only sequence types assigned up to species level
nrow(species) # unique ASVs
ASV_sum_species <- rowSums(species[, c(4:205)], na.rm = TRUE)
species_sum <- sum(ASV_sum_species)
species_sum # counts assigned to species level
(species_sum/raw_sum)*100 # % of raw data assigned at least to species level 

# export
write.xlsx(genus, file = "Output/01_Filtering/rbcL_filtered_6.xlsx")

# unique taxonomic levels
# unique families
unique(id96_filter6$family_name)
length(unique(id96_filter6$family_name)) -1 # -1 because of NA
# unique genera
unique(id96_filter6$genus_name)
length(unique(id96_filter6$genus_name)) -1
# unique species
unique(id96_filter6$species_name)
length(unique(id96_filter6$species_name)) -1

# counts on family level
ASV_sum_family
data_per_family <- data.frame(family$family_name, ASV_sum_family)
names(data_per_family) <- c("family_name", "ASV_sum_family")

sumcounts_family <- data_per_family %>%
  filter(!is.na(family_name)) %>%
  group_by(family_name) %>%
  summarise(sumcount_sample = sum(ASV_sum_family))

# counts on genus level
ASV_sum_genus
data_per_genus <- data.frame(genus$genus_name, ASV_sum_genus)
names(data_per_genus) <- c("genus_name", "ASV_sum_genus")

sumcounts_genus <- data_per_genus %>%
  filter(!is.na(genus_name)) %>%
  group_by(genus_name) %>%
  summarise(sumcount_sample = sum(ASV_sum_genus))


########################################################################################################################################################
####################### 2 Counts per replicate  ########################################################################################################
########################################################################################################################################################

dir.create("Output/02_Replicates")

rm(list = ls())


library("tidyverse")
library("readxl")
library("openxlsx")

# check the similarity of metabarcoding replicates with an NMDS
# based on presence/absence of a taxon in the filtered dataset
#and its assignment to a "scientific name"
# (= lowest taxonomic unit possible to assign, not restricted to a taxonomic level)

#load data
diatoms <- read.xlsx("Output/01_Filtering/rbcL_filtered_6.xlsx")

# overview dataset
names(diatoms)
names(diatoms)[1:10]
nuc_seq <- diatoms$NUC_SEQ # get vector with ASVs

diatoms_onlycounts <- as.data.frame(sapply(diatoms[, 4:205], as.numeric)) # counts numeric

names(diatoms_onlycounts)

diatoms_onlycounts <- as.data.frame(t(diatoms_onlycounts))
colnames(diatoms_onlycounts) <- nuc_seq

rownames(diatoms_onlycounts)

diatoms_onlycounts$sample_name <- rownames(diatoms_onlycounts)

# get sample ages
sample_info <- read.xlsx("Sample_names.xlsx")
replicates <- diatoms_onlycounts
replicates$sample_name <- rownames(replicates)
replicates <- merge(replicates, sample_info, by = "sample_name")
ncol(replicates)
replicates <- replicates[, c(477, 1:476)]

replicates <- replicates[!grepl(c("EB", "NTC"), replicates$age),] # remove blanks
ncol(replicates)
replicates <- replicates[, -2] # remove sample_name

replicates$age <- as.numeric(replicates$age)/1000 # age in ka BP
replicates <- arrange(replicates, age) # sort by age

# preparation for export
replicates <- as.data.frame(t(replicates))
age <- replicates[1, ]
colnames(replicates) <- age # age as column names
replicates <- replicates[-1, ]
replicates[is.na(replicates)] <- 0 # replace NAs by 0
scientific_name <- rownames(replicates)
replicates <- as.data.frame(sapply(replicates, as.numeric)) # counts numeric
replicates$scientific_name <- scientific_name
ncol(replicates)
replicates <- replicates[, c(163, 1:162)]

# export
write.xlsx(replicates, file = "Output/02_Replicates/ASVs_per_replicate_pre_rarefaction.xlsx")


########################################################################################################################################################
####################### 3 NMDS  ########################################################################################################################
########################################################################################################################################################

dir.create("Output/03_NMDS")

rm(list = ls())

library("readxl")
library("openxlsx")
library("vegan")
library("ggplot2")
library("tidypaleo")

theme_set(theme_paleo(8))

#############################################
###         (A)  Data preparation         ###
#############################################

# load count data
counts <- read.xlsx("Output/02_Replicates/ASVs_per_replicate_pre_rarefaction.xlsx")

scientific_name <- counts[, 1]

# Helliger transformation
counts_hell <-decostand(counts[, 2:ncol(counts)], method = "hellinger")
rownames(counts_hell) <- scientific_name
counts_hell <- as.data.frame(t(counts_hell))

# Check the Axis lengths of DCA1, use CCA if >4
decorana(counts_hell) # error, there is a row without any counts

rowSums(counts_hell)
which(rowSums(counts_hell) == 0) # see which is row without counts

#remove that row
rownames(counts_hell)
counts_hell_without_zero_row <- counts_hell[-138,]

counts_hell <- counts_hell_without_zero_row # overwrite original seq_data_hell

# for presence/absence
counts_hell[counts_hell > 0] <- 1
decorana(counts_hell)


#############################################
###       (B)  NMDS on name level         ###
#############################################

# # what nmds does:
# 1.Define the original positions of communities in multidimensional space.
# 2.Specify the number of reduced dimensions (typically 2).
# 3.Construct an initial configuration of the samples in 2-dimensions.
# 4.Regress distances in this initial configuration against the observed (measured) distances.
# 5.Determine the stress, or the disagreement between 2-D configuration and predicted values from the regression.
# stress < 0.05 provides an excellent representation in reduced dimensions, 
# stress < 0.1 is great, 
# strss < 0.2 is good/ok, and 
# stress > 0.3 provides a poor representation.
# high stress is bad, low stress is good!

set.seed(123)
nmds_name <- metaMDS(counts_hell, k=2, trymax = 100)
nmds_name
plot(nmds_name)

# extract position of replicates
age_position <- as.data.frame(nmds_name$points)


# import age periods
age_period <- read.xlsx("Age_period.xlsx")
age_period <- age_period[-138, ]
age_position$age <- age_period$age # add age column
names(age_position)
age_position$period <- age_period$period # add period to ages
age_position$period <- as.factor(age_position$period) # make period a factor
levels(age_position$period) # check factor levels

# get list with colors for plot later;
# note that colors are given to groups in alphabetical order:
# 1 Bolling-Allerod
# 2 Eemian
#3 Glacial
# Holocene
# LGM
# Younger Dryaas

col_list <- c("gold",  "brown2", "steelblue1",   "orange2",  "blue1", "turquoise1")


NMDS_replicates <- ggplot(age_position, aes(x = MDS1, y = MDS2)) +
  theme(axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 15, colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  geom_polygon(data = age_position, aes(group = age, color = period, fill = period), alpha = 0.3) +
  scale_color_manual(values = col_list) +
  scale_fill_manual(values = col_list) +
  geom_point(data = age_position, aes(x = MDS1, y = MDS2), pch = 20, size = 2, stroke = 0.8, colour = adjustcolor("gray26", alpha.f = 0.3)) +
  geom_text(data = age_position, aes(label=age),hjust=0, vjust=0, size = 3) + # add all family names
  annotate("text", x = -0.8, y = 0.95, label = paste0("stress: ", format(nmds_name$stress, digits = 4)), size = 6)

NMDS_replicates

ggsave(NMDS_replicates, file = "Output/03_NMDS/NMDS_rbcL_replicates_scientific_name_level_presence_absence.pdf", width = 23, height = 18, units = "cm")


########################################################################################################################################################
####################### 4 Sum of replicates  ###########################################################################################################
########################################################################################################################################################

dir.create("Output/04_Sum_of_Replicates")

rm(list = ls())

library("readxl")
library("openxlsx")
library("tidyverse")

replicates <- read.xlsx("Output/02_Replicates/ASVs_per_replicate_pre_rarefaction.xlsx")

# import ages and add age column
age_period <- read.xlsx("Age_period.xlsx")
names(age_period)
replicates <- as.data.frame(t(replicates))
colnames(replicates) <- replicates[1, ]
replicates <- replicates[-1, ]
replicates$age <- age_period$age
ncol(replicates)
replicates <- replicates[, c(476, 1:475)]

# sum of replicates
replicates$age <- as.factor(replicates$age)
levels(replicates$age)

replicates_long <- pivot_longer(replicates, # turn into long format
                            cols = !c(age),
                            names_to = c("NUC_SEQ"),
                            values_to = "counts")

replicates_long$NUC_SEQ <- as.factor(replicates_long$NUC_SEQ)
replicates_long$counts <- as.numeric(replicates_long$counts)

#summarize read counts over family even if lower taxa is known
sample <- replicates_long %>%
  group_by(NUC_SEQ, age) %>%
  summarise(sumcount_ASV = sum(counts))


###turn to wide format
sample_wide <- sample %>%
  pivot_wider(names_from = age,
              values_from = sumcount_ASV)

# export
write.xlsx(sample_wide, file = "Output/04_Sum_of_Replicates/ASV_per_sample_pre_rarefaction.xlsx")


########################################################################################################################################################
####################### 5 Rarefaction ##################################################################################################################
########################################################################################################################################################

rm(list = ls())

dir.create("Output/05_Rarefaction")

options(stringsAsFactors=FALSE)

library("openxlsx")

#############################################
###         PREPARE DATA FRAME            ###
#############################################

# read data and check visually if the df is as requested
# requirements: (1) blanks deleted
#               (2) "x" infront of the age, which are the column; otherwise they are not recognized as the column names

# import data
data <- read.xlsx("Output/04_Sum_of_Replicates/ASV_per_sample_pre_rarefaction.xlsx")

# add the "x" infront of the column names (age)
names(data)
ncol(data)
data_temp <- data[, 2:55] # remove columns where no "x" will be added (ASV and blank columns)
names(data_temp)
ages <- names(data_temp) # extract ages
ages_mod <- paste("x", ages, sep = "") # add the "x"
ages_mod
colnames(data_temp) <- ages_mod # add modified ages to data frame
data_temp$family <- data$NUC_SEQ # add ASV again to data frame (name the column "family" because Stefan's script uses the column like this)
ncol(data_temp)
data_temp <- data_temp[, c(55, 1:54)] # rearrange so ASV is first column


#####################################################################################
###  Rarefaction_process_data_family_level_Stefan's codes ###########################
#####################################################################################

# number of repeats
repeatnumber=100 # 100 standard

specseq <- data_temp
str(specseq)

# cols with data
# in case of other structure of your df adapt numbers
colstart=2
colend=as.numeric(dim(specseq)[2]) # sets dimensions of df
ncol(specseq)

# determine minimum number to rarefy to
# automatically from df, you can manipulate it here (even to different values e.g. c(100,200,mincounts) to check also the others)
mincounts=sort(apply(specseq[2:55], 2, sum))[1]
mincounts# print and check

# family and taxa is merged for faster computation --> not possible here as only family information exist (no species given)
specseq_final_name=paste(specseq$family)
famorig=specseq$family

# resampling
specseq=as.data.frame(specseq)
nsampleff=mincounts
genrare=list()		# empty container for data
yr=NULL			# sample name in case of numbers the leading X is deleted
missing=0			# needed for if columns without data exist, they are removed
for(yrcoli in colstart:colend)
{
  print(yrcoli)
  
  allspec=specseq_final_name[which(specseq[,yrcoli]>0)]
  allspec_counts=specseq[which(specseq[,yrcoli]>0),yrcoli]
  
  if(length(allspec)>0)
  {
    #yr=c(yr, round(as.numeric(gsub("X","",names(specseq)[yrcoli]))))
    yr=c(yr, names(specseq)[yrcoli])
    
    sampleeffort=list()
    for(nsampleeffi in nsampleff)
    {
      repeatsample=list()
      for(repi in 1:repeatnumber)
      {
        repeatsample[[repi]]=sample(allspec,nsampleeffi,replace=TRUE, prob=allspec_counts/sum(allspec_counts))
      }
      sampleeffort[[which(nsampleff==nsampleeffi)]]=repeatsample
    }
    genrare[[yrcoli-missing-(colstart-1)]]=list(allspec,sampleeffort)
  } else
  {
    missing=missing+1
  }
}
names(genrare)=yr
# save for later use
save(genrare, file="Output/05_Rarefaction/rarefaction_metabarcoding.RDATA")


# analyses
# ... count total species number 
# ... ... and total species per family

par(mar=c(10,2,0,0));barplot(table(famorig),las=2)
familylevels=names(rev(sort(table(famorig))))
totspec=NULL
totfam=NULL
famreads=NULL
specreads=NULL
for(li in 1:length(genrare))
{
  for(li2 in 1:length(genrare[[li]][[2]]))
  {
    print(paste0(li," - ",li2))
    
    spectot=NULL
    spectot4fam=NULL
    famread=NULL
    speccounts=NULL
    for(repi in 1:100)
    {
      pei=unique(genrare[[li]][[2]][[li2]][[repi]])
      spectot=c(spectot,length(pei))
      spectot4fam=rbind(spectot4fam,table(factor(unlist(lapply(strsplit(split=" ",pei),function(x)return(x[1]))),levels=familylevels)))
      
      famread=rbind(famread,table(factor(unlist(lapply(strsplit(split=" ",genrare[[li]][[2]][[li2]][[repi]]),function(x)return(x[1]))),levels=familylevels)))
      speccounts=rbind(speccounts,table(factor(genrare[[li]][[2]][[li2]][[repi]],levels=specseq_final_name)))
    }
    totspec=rbind(totspec, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),Nspecies=spectot))
    totfam=rbind(totfam, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),spectot4fam))
    famreads=rbind(famreads, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),famread))
    specreads=rbind(specreads, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),speccounts))
    
  }
}
str(totspec)	# count family number per sample/year
#write.csv2(totspec, "Rarefaction_pelagic_FamiliesPerSample.csv", row.names=FALSE)
write.xlsx(totspec, file = "Output/05_Rarefaction/Rarefaction_Diatoms_ASVsPerSample.xlsx")

str(famreads)	# count counts per families and sample/year

write.xlsx(famreads, file = "Output/05_Rarefaction/Rarefaction_Diatoms_ASVCounts.xlsx")

################# aggregate mean and confidence interval 95% for families ###################
famdf=NULL
for(sampleidi in unique(totfam$T))
{
  sampleidi
  totspecsub=totspec[totspec$T==sampleidi,"Nspecies"]
  famdf=rbind(famdf,data.frame(Sample=sampleidi, Family="ALL", NSpecMean=mean(totspecsub), NSpecCI95=sd(totspecsub)*1.96))
  
  for(familyi in names(totfam)[3:dim(totfam)[2]])
  {
    totfamsub=totfam[totfam$T==sampleidi,familyi]	
    famdf=rbind(famdf,data.frame(Sample=sampleidi, Family=familyi, NSpecMean=mean(totfamsub), NSpecCI95=sd(totfamsub)*1.96))
  }
}
str(famdf)
write.xlsx(famdf, "Output/05_Rarefaction/Rarefaction_Diatoms_ASVsCI95.xlsx")


#############################################
###       MEAN FAMILIES PER SAMPLE        ###
#############################################

rm(list = ls())

library("readxl")
library("openxlsx")

# load data coming from Rarefaction Script
diatoms_rar <- read_excel("Output/05_Rarefaction/Rarefaction_Diatoms_ASVsPerSample.xlsx")
colnames(diatoms_rar) <- c("age", "SampleEff", "NASVs")

diatoms_rar$age <- sub("x", "", diatoms_rar$age) # removes the x infront of the age
diatoms_rar$age <- factor(diatoms_rar$age, levels = unique(diatoms_rar$age)) # age as factor
is.factor(diatoms_rar$age)

# calculate mean for each year (100 repeats for rarefaction)
nasv_mean <- data.frame(c(1:length(levels(diatoms_rar$age))))
for(i in c(1:nrow(diatoms_rar))) {
  nasv_mean[i, ] <- mean(diatoms_rar$NASVs[diatoms_rar$age == levels(diatoms_rar$age)[i]])
}
nasv_mean <- na.omit(nasv_mean) # delete rows with NA

nasv_mean$age <- levels(diatoms_rar$age)
nasv_mean <- nasv_mean[, c(2, 1)]
colnames(nasv_mean) <- c("age", "NASVs_mean")

write.xlsx(nasv_mean, file = "Output/05_Rarefaction/Rarefaction_Diatoms_ASVsPerSample_mean.xlsx")


#############################################
###   MEAN COUNTS PER SAMPLE AND FAMILY   ###
#############################################

rm(list = ls())

library("readxl")
library("openxlsx")

# load data coming from Rarefaction Script
diatoms <- read_excel("Output/05_Rarefaction/Rarefaction_Diatoms_ASVCounts.xlsx")

names(diatoms)
colnames(diatoms)[1] <- "age" # rename age column
diatoms <- diatoms[, -2] # removes sampleEff column

diatoms$age <- factor(diatoms$age, levels = unique(diatoms$age)) # make age a factor
is.factor(diatoms$age)
length(levels(diatoms$age)) # how many age levels (= samples)

# get subset of all resampling repeats per sample age (100 repeats for rarefaction)
a_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[1])
b_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[2])
c_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[3])
d_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[4])
e_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[5])
f_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[6])
g_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[7])
h_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[8])
i_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[9])
j_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[10])
k_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[11])
l_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[12])
m_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[13])
n_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[14])
o_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[15])
p_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[16])
q_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[17])
r_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[18])
s_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[19])
t_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[20])
u_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[21])
v_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[22])
w_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[23])
x_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[24])
y_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[25])
z_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[26])
aa_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[27])
ab_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[28])
ac_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[29])
ad_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[30])
ae_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[31])
af_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[32])
ag_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[33])
ah_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[34])
ai_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[35])
aj_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[36])
ak_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[37])
al_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[38])
am_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[39])
an_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[40])
ao_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[41])
ap_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[42])
aq_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[43])
ar_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[44])
as_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[45])
at_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[46])
au_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[47])
av_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[48])
aw_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[49])
ax_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[50])
ay_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[51])
az_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[52])
ba_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[53])
bb_diatoms <- subset(diatoms, subset = age %in% levels(diatoms$age)[54])

# calculate mean for each sampling age
a <- colMeans(a_diatoms[, 2:ncol(a_diatoms)])
b <- colMeans(b_diatoms[, 2:ncol(b_diatoms)])
c <- colMeans(c_diatoms[, 2:ncol(c_diatoms)])
d <- colMeans(d_diatoms[, 2:ncol(d_diatoms)])
e <- colMeans(e_diatoms[, 2:ncol(e_diatoms)])
f <- colMeans(f_diatoms[, 2:ncol(f_diatoms)])
g <- colMeans(g_diatoms[, 2:ncol(g_diatoms)])
h <- colMeans(h_diatoms[, 2:ncol(h_diatoms)])
i <- colMeans(i_diatoms[, 2:ncol(i_diatoms)])
j <- colMeans(j_diatoms[, 2:ncol(j_diatoms)])
k <- colMeans(k_diatoms[, 2:ncol(k_diatoms)])
l <- colMeans(l_diatoms[, 2:ncol(l_diatoms)])
m <- colMeans(m_diatoms[, 2:ncol(m_diatoms)])
n <- colMeans(n_diatoms[, 2:ncol(n_diatoms)])
o <- colMeans(o_diatoms[, 2:ncol(o_diatoms)])
p <- colMeans(p_diatoms[, 2:ncol(p_diatoms)])
q <- colMeans(q_diatoms[, 2:ncol(q_diatoms)])
r <- colMeans(r_diatoms[, 2:ncol(r_diatoms)])
s <- colMeans(s_diatoms[, 2:ncol(s_diatoms)])
t <- colMeans(t_diatoms[, 2:ncol(t_diatoms)])
u <- colMeans(u_diatoms[, 2:ncol(u_diatoms)])
v <- colMeans(v_diatoms[, 2:ncol(v_diatoms)])
w <- colMeans(w_diatoms[, 2:ncol(w_diatoms)])
x <- colMeans(x_diatoms[, 2:ncol(x_diatoms)])
y <- colMeans(y_diatoms[, 2:ncol(y_diatoms)])
z <- colMeans(z_diatoms[, 2:ncol(z_diatoms)])
aa <- colMeans(aa_diatoms[, 2:ncol(aa_diatoms)])
ab <- colMeans(ab_diatoms[, 2:ncol(ab_diatoms)])
ac <- colMeans(ac_diatoms[, 2:ncol(ac_diatoms)])
ad <- colMeans(ad_diatoms[, 2:ncol(ad_diatoms)])
ae <- colMeans(ae_diatoms[, 2:ncol(ae_diatoms)])
af <- colMeans(af_diatoms[, 2:ncol(af_diatoms)])
ag <- colMeans(ag_diatoms[, 2:ncol(ag_diatoms)])
ah <- colMeans(ah_diatoms[, 2:ncol(ah_diatoms)])
ai <- colMeans(ai_diatoms[, 2:ncol(ai_diatoms)])
aj <- colMeans(aj_diatoms[, 2:ncol(aj_diatoms)])
ak <- colMeans(ak_diatoms[, 2:ncol(ak_diatoms)])
al <- colMeans(al_diatoms[, 2:ncol(al_diatoms)])
am <- colMeans(am_diatoms[, 2:ncol(am_diatoms)])
an <- colMeans(an_diatoms[, 2:ncol(an_diatoms)])
ao <- colMeans(ao_diatoms[, 2:ncol(ao_diatoms)])
ap <- colMeans(ap_diatoms[, 2:ncol(ap_diatoms)])
aq <- colMeans(aq_diatoms[, 2:ncol(aq_diatoms)])
ar <- colMeans(ar_diatoms[, 2:ncol(ar_diatoms)])
as <- colMeans(as_diatoms[, 2:ncol(as_diatoms)])
at <- colMeans(at_diatoms[, 2:ncol(at_diatoms)])
au <- colMeans(au_diatoms[, 2:ncol(au_diatoms)])
av <- colMeans(av_diatoms[, 2:ncol(av_diatoms)])
aw <- colMeans(aw_diatoms[, 2:ncol(aw_diatoms)])
ax <- colMeans(ax_diatoms[, 2:ncol(ax_diatoms)])
ay <- colMeans(ay_diatoms[, 2:ncol(ay_diatoms)])
az <- colMeans(az_diatoms[, 2:ncol(az_diatoms)])
ba <- colMeans(ba_diatoms[, 2:ncol(ba_diatoms)])
bb <- colMeans(bb_diatoms[, 2:ncol(bb_diatoms)])

# combine means per age in one data frame
countASV_mean <- as.data.frame(rbind(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab, ac, ad, ae, af, ag, ah, ai, aj, ak, al, am, an, ao, ap, aq, ar, as, at, au, av, aw, ax, ay, az, ba, bb))

ASVs <- names(countASV_mean) # save vector with families for later

age <- levels(diatoms$age)
countASV_mean$age <- age
ncol(countASV_mean)
countASV_mean <- countASV_mean[, c(476, 1:475)]
rownames(countASV_mean) <- NULL

countASV_mean$age <- sub("x", "", countASV_mean$age) # removes the x infront of the age
countASV_mean <- as.data.frame(t(countASV_mean)) # transposing so age is in columns
age <- countASV_mean[1, ] # vector with ages
colnames(countASV_mean) <- age # ages to column names
countASV_mean <- countASV_mean[-1, ] # remove row with ages

countASV_mean <- as.data.frame(sapply(countASV_mean, as.numeric)) # counts numeric

countASV_mean$NUC_SEQ <- ASVs  # add column with ASV
ncol(countASV_mean)
countASV_mean <- countASV_mean[, c(55, 1:54)] # make family the first column
rownames(countASV_mean) <- NULL # remove row names

write.xlsx(countASV_mean, file = "Output/05_Rarefaction/Rarefaction_Diatoms_ASVCounts_means.xlsx")


#############################################
###         TAXONOMIC INFORMATION         ###
#############################################

rm(list = ls())

# load data
diatoms <- read.xlsx("Output/05_Rarefaction/Rarefaction_Diatoms_ASVCounts_means.xlsx")
# load taxonomic indormation
tax_info <- read.xlsx("AVS_to_Genus_info.xlsx")

diatoms <- merge(diatoms, tax_info, by = "NUC_SEQ") # merge
ncol(diatoms)
diatoms <- diatoms[, c(56, 1:55)]
diatoms <- diatoms[, -2] # remove NUC_SEQ column

# export
write.xlsx(diatoms, file = "Output/05_Rarefaction/Rarefaction_Diatoms_ASVCounts_means_genera.xlsx")


########################################################################################################################################################
####################### 6 Stratigrams ##################################################################################################################
########################################################################################################################################################

rm(list = ls())

dir.create("Output/06_Stratigram")

library("tidypaleo")
library("readxl")
library("openxlsx")
library("tidyverse")
library("ggplot2")
library("ggh4x")

theme_set(theme_paleo(8))

#############################################
###         (A) DATA PREPARARION          ###
#############################################

# read in data
diatoms <- read_excel("Output/05_Rarefaction/Rarefaction_Diatoms_ASVCounts_means_genera.xlsx")

# create data frame with only columns with counts
ncol(diatoms)
names(diatoms)
diatoms_onlycounts <- diatoms[, -1] # remove columns with taxonomic information, so only counts are left

# create data frame with relative counts per column
diatoms_rel <- data.frame(c(1:nrow(diatoms_onlycounts)))
for( i in c(1:ncol(diatoms_onlycounts))) {
  diatoms_rel[, i] <- data.frame(diatoms_onlycounts[, i]) %>%
    transmute(rel_ab = (diatoms_onlycounts[, i]/sum((diatoms_onlycounts[, i])))) # calculation for relative abundance for each cell of the data frame
}
diatoms_rel <- diatoms_rel*100 # for making relative abundance in % (somehow it makes the column names weird if yoh do it right in the for-loop)

# add taxonomic information again
diatoms_rel <- cbind(diatoms_rel, diatoms$ASV) # add ASV column and class column
names(diatoms_rel)
ncol(diatoms_rel)
colnames(diatoms_rel)[55] <- "ASV" # rename ASV column correctly
diatoms_rel <- diatoms_rel[, c(55, 1:54)]

# export
write.xlsx(diatoms_rel, file = "Output/06_Stratigram/ASVs_all_rel_abund.xlsx")


#############################################
###             (B) FILTERING             ###
#############################################

#change to long format (needed for plotting)
diatoms_long <- pivot_longer(diatoms_rel,
                             cols = !c(ASV),
                             names_to = c("age"),
                             values_to = "rel_abund")

# make age numeric
diatoms_long$age <- as.numeric(diatoms_long$age)

# For plotting I want to keep only ASVs that occur in at least 20 samples
#and in at least 1 sample with a relative abundance > 3%

# 1. remove all rows (in long format) of ASVs occurring less than 20 times
long_biggerzero <-  filter(diatoms_long, rel_abund > 0) # subset of long format data with only relative abundance > 0
fam_biggerzero <- long_biggerzero$ASV # only ASV names
fam_counts <- as.data.frame(table(fam_biggerzero)) # counts how often ASVs occur with relative abundance > 0
threshold_occurrence <- filter(fam_counts, Freq >= 30)$fam_biggerzero # subset of ASVs that occur more than 10 times
threshold_occurrence <- as.character(threshold_occurrence) # make it a vector of characters


diatoms_filter1 <- subset(diatoms_long, subset = ASV %in% threshold_occurrence)

# 2. remove all rows (in long format) of ASVs that don't occur in at least 1 sample with min. 3%
long_biggerpro <-  filter(diatoms_filter1, rel_abund >= 3) # subset of already filtered data set to select rows that have a relative >= 1% 
fam_biggerpro <- long_biggerpro$ASV # only ASV names 
fam_countspro <- as.data.frame(table(fam_biggerpro)) # counts how often ASVs occur with relative abundance > 5
threshold_relabund <- filter(fam_countspro, Freq >= 1)$fam_biggerpro # subset of ASVs that occur in at least 1 sample with relative abundance of 1%
threshold_relabund <- as.character(threshold_relabund) # make it a vector of characters

diatoms_filter2 <- subset(diatoms_filter1, subset = ASV %in% threshold_relabund)

# for rel abund of main filtered ASVs
diatoms_main <- pivot_wider(diatoms_filter2,
                            names_from = "age",
                            values_from = "rel_abund")

names(diatoms_main)
nrow(diatoms_main)

write.xlsx(diatoms_main, file = "Output/06_Stratigram/ASVs_main_rel_abund.xlsx")


#add shape infos
shape <- read.xlsx("Shape_info.xlsx")
diatoms_main <- merge(diatoms_main, shape, by = "ASV")
names(diatoms_main)
ncol(diatoms_main)
diatoms_main <- diatoms_main[, c(1, 56, 2:55)]

# get total pennate and total centric counts
tot_pennate <- subset(diatoms_main, subset = shape %in% "pennate")
tot_centric <- subset(diatoms_main, subset = shape %in% "centric")
tot_pennate_v <- colSums(tot_pennate[, 3:ncol(tot_pennate)]) # vector with sums of pennate
tot_pennate_v <- as.vector(c("total pennate diatoms", "pennate", tot_pennate_v))
tot_centric_v <- colSums(tot_centric[, 3:ncol(tot_centric)])
tot_centric_v <- as.vector(c("total centric diatoms", "centric", tot_centric_v))

diatoms_main <- rbind(tot_centric,
                      tot_centric_v,
                      tot_pennate,
                      tot_pennate_v)

# transfer shape information to filter2
diatoms_filter2 <- pivot_longer(diatoms_main,
                                cols = !c(ASV, shape),
                                names_to = c("age"),
                                values_to = "rel_abund")


#########################################
###      (C) PLOTTING MAIN ASVs       ###
#########################################

# preparation of data
diatoms_filter2$shape <- factor(diatoms_filter2$shape, levels = unique(diatoms_filter2$shape)) # make shape a factor

# order
names(diatoms_filter2)

diatoms_filter2$rel_abund <- as.numeric(diatoms_filter2$rel_abund) # make re_abund numeric
diatoms_filter2$age <- as.numeric(diatoms_filter2$age) # make age numeric
diatoms_filter2$ASV <- factor(diatoms_filter2$ASV, levels = unique(diatoms_filter2$ASV)) # make ASV a factor with order it has right now
diatoms_filter2$shape <- factor(diatoms_filter2$shape, levels = unique(diatoms_filter2$shape)) # make shape a factor with order it has right now


# diatoms plot:
diatoms_plot <- ggplot(diatoms_filter2, aes(x = age, y = rel_abund)) +
  geom_area(aes(fill = shape)) +
  geom_segment(aes(xend = age, yend = 0), lwd = 0.5, color = "black") +
  expand_limits(x = 0, y = 0) + # axis start with zero
  scale_y_continuous(n.breaks = 3, # numbers mark y axis
                     expand = c(0,0),
                     name = "Relative abundance [%] \nbased on total reads of ASVs") + 
  facet_abundanceh(vars(ASV), rotate_facet_labels = 70, space = "fixed", scales = "free_x") + # shapes plots by ASV, makes every plot the same size (different scale)
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, color = c(NA, "black")), # rotates x axis by 90 degree, labels in middle of tick (vjust)
        axis.title.x = element_text(vjust = 0.5, size = 14),
        axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14, colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(colour = "black"),
        panel.grid = element_blank(), # removes background panel
        panel.border = element_blank(),
        strip.text = element_text(size = 14, vjust = 1, face = "bold"),
        strip.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  geom_hline(yintercept = 0, color = "black", size = 0.4) +
  scale_x_reverse(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120), # ticks on x axis 
                  expand = c(0,1),
                  name = "Age [ka BP]") +
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") + # reverse axes
  force_panelsizes(rows = unit(9, "cm"),
                   cols = unit(0.8, "cm"))

diatoms_plot
ggsave(diatoms_plot, file = "Output/06_Stratigram/ASVs_area.pdf", width = 33, height = 20, units = "cm")


########################################################################################################################################################
####################### 7 PCA/RDA ######################################################################################################################
########################################################################################################################################################

rm(list = ls())

dir.create("Output/07_PCA_RDA")

library("openxlsx")
library("vegan")
library("ggplot2")
library("tidypaleo")
library("openxlsx")
library("olsrr")
library("ggrepel")

theme_set(theme_paleo(8))


#########################################
###      (A) ENVIRONMENTAL DATA       ###
#########################################

# NGRIP and zooplankton data (Calanidae and heterotrophic zooplankton)
# is the same as for the Shotgun analysis.
# Therefore, the same data are imported.
# For the creation of the dataframe, see the Shotgun R script.

envi_data <- read.xlsx("Envi_for_RDA.xlsx")


#########################################
###       (B) SEQUENCING DATA         ###
#########################################

# # import relative abundance of the main ASVs (same as in stratigram)
# and create vector with ASVs for PCA
seq_data_ASVs <-read.xlsx("Output/06_Stratigram/ASVs_main_rel_abund.xlsx")
names(seq_data_ASVs)
seq_ASVs <- as.vector(seq_data_ASVs$ASV)

# import full relative abundance dataset for Hellinger transformation
seq_data <-read.xlsx("Output/06_Stratigram/ASVs_all_rel_abund.xlsx")
names(seq_data)

age_long <- names(seq_data[-1]) # get ages (= all metabarcoding samples)
length(age_long)

seq_data <- as.data.frame(t(seq_data))
names(seq_data) <- seq_data[1,]
seq_data <- seq_data[-1,]
seq_data <- as.data.frame(lapply(seq_data ,as.numeric))
seq_data$age <- age_long

# because Zooplankton is an enviromental facgor in the RDA, but zooplankton data
# derive from the shotgun approach, not all metabarcoding samples can be used for
# the PCA and RDA (explaining variables must have the same number of samples)
samples <- read.xlsx("Age_Period_Group.xlsx") # import ages (= samples) that were used in Shotgun as well
sample_ages <- as.vector(samples$age)
class(sample_ages)
class(seq_data$age)
seq_data$age <- as.numeric(seq_data$age)
class(seq_data$age)

# filter for samples analyszed in both metabarcoding and shotgun
seq_data_filtered <- subset(seq_data, subset = age %in% sample_ages)
length(sample_ages)
nrow(seq_data_filtered)

names(seq_data_filtered)
ncol(seq_data_filtered)
seq_data_filtered <- seq_data_filtered[, c(476, 1:475)] # change column order (age = col 1)
names(seq_data_filtered)


#########################################
###   (C) HELLINGER TRANSFORMATION    ###
#########################################

# Hellinger transformation = square root of relative abundance data
# here: count data already transformed to relative abundances,
# therefore only squre root transformation on seq_data needed for Hellinger transformed data

seq_data_hell <- sqrt(seq_data_filtered[, 2:ncol(seq_data_filtered)])
rownames(seq_data_hell) <- seq_data_filtered$age

#filter for main ASVs
seq_data_hell <- as.data.frame(t(seq_data_hell))
seq_data_hell$ASV <- rownames(seq_data_hell)
seq_data_filtered <- subset(seq_data_hell, subset = ASV %in% seq_ASVs) # filter
names(seq_data_filtered)
ncol(seq_data_filtered)
seq_data_filtered <- seq_data_filtered[, -43] # remove ASV column again
seq_data_filtered <- as.data.frame(t(seq_data_filtered))

seq_data_filtered <- as.data.frame(sapply(seq_data_filtered, as.numeric)) # counts numeric
rownames(seq_data_filtered) <- names(seq_data_hell[1:42])


#########################################
###            (D) PCA/RDA            ###
#########################################

# run a  pca
pca <- rda(seq_data_filtered, scale=TRUE)
plot(pca)

summary <- summary(pca)

# check collinearity (on untransformed data):
model <- lm(age ~ NGRIP + Calanidae + het_pro, data = envi_data)
ols_vif_tol(model)
# VIF should be <10, but better is <4

#https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
# with all variables
rownames(envi_data) <- envi_data$age # make age rownames
names(envi_data)
selected_envi_data <- envi_data[, -1] # remove age column
names(selected_envi_data)

data_envi_rda_all <- rda(seq_data_filtered ~., data = selected_envi_data)
plot(data_envi_rda_all)
adjR2.all <- RsquareAdj (data_envi_rda_all)$adj.r.squared
adjR2.all


data_envi_rda0 <- rda(seq_data_filtered ~ 1, data = selected_envi_data) # "empty" model with intercept only
data_envi_rda_all <- rda(seq_data_filtered ~., data = selected_envi_data) # full model
sel.osR2 <- ordiR2step (data_envi_rda0, scope = formula (data_envi_rda_all), R2scope = adjR2.all, direction = 'forward', permutations = 999)
sel.osR2
sel.osR2$anova


sel.osR2_adj <- sel.osR2
sel.osR2_adj$anova$`Pr(>F)` <- p.adjust (sel.osR2$anova$`Pr(>F)`, method = 'holm', n = ncol (selected_envi_data))
sel.osR2_adj$anova


# check backward selection (gives the same result as forward selection)
sel.osR2_back <- ordiR2step (data_envi_rda0, scope = formula (data_envi_rda_all), R2scope = adjR2.all, direction = 'backward', permutations = 999)
sel.osR2_back
sel.osR2_back$anova


sel.osR2_adj_back <- sel.osR2_back
sel.osR2_adj_back$anova$`Pr(>F)` <- p.adjust (sel.osR2_back$anova$`Pr(>F)`, method = 'holm', n = ncol (selected_envi_data))
sel.osR2_adj_back$anova


#########################################
###         (E) SUMMARY RDA           ###
#########################################

#https://mb3is.megx.net/gustame/constrained-analyses/rda for interpretation

# final RDA
rda_final <- rda(seq_data_filtered ~., data = selected_envi_data)
plot(rda_final)

adjR2.final <- RsquareAdj (rda_final1)$adj.r.squared
adjR2.final

summary(rda_final)

screeplot(rda_final)

#########################################
###           (F) PLOTTING            ###
#########################################

plot(pca)

# incorporate environmental factors as restricting variables
fit_all <- envfit(pca ~ . , selected_envi_data, perm = 999)
# perm = to test significance by permutations
fit_all

# as multiple environmental variables are used, the p-values are corrected by Bonferroni correction
# ("The function envfit executed three permutation tests, all three with highly significant results.
# When the number of tested variables increases, correction for multiple testing is desirable.
# In the case above, we may extract the significance values from ef object and apply Bonferroni
# correction using the function p.adjust";
# from https://www.davidzeleny.net/anadat-r/doku.php/en:suppl_vars_examples)
fit.adj <- fit_all
pvals.adj <- p.adjust (fit_all$vectors$pvals, method = 'bonferroni')
fit.adj$vectors$pvals <- pvals.adj
fit.adj

# extract positions of age, families and environmental parameter vectors
summary <- summary(pca)

# extract age information
age_position_pca <- as.data.frame(summary$sites[, 1:2])
age_position_pca$age <- as.numeric(rownames(age_position_pca)) # make column with age

age_info <- read.xlsx("Age_Period_Group.xlsx") # import age_period table
age_position_pca <- merge(age_position_pca, age_info, by = "age")
age_position_pca$pch <- as.factor(age_position_pca$pch)

# extract ASV information
ASV_position_pca <- as.data.frame(summary$species[, 1:2])
# type I scaling for families:
ASV_scaling <- plot(pca, scaling = 1)
ASV_position_scaleII <- as.data.frame(ASV_scaling$species)

# add size information in colour of the data point
shape_info <- read.xlsx("Shape_info.xlsx")
ASV_temp <- ASV_position_scaleII
ASV_temp$ASV <- rownames(ASV_temp)
ASV_temp <- merge(ASV_temp, shape_info, by = "ASV")
names(ASV_temp)

rownames(ASV_temp) <- ASV_temp$ASV

ASV_position_final <- ASV_temp
names(ASV_position_final)
ASV_position_final <- ASV_position_final[, -1] # remove ASV column

# extract environmental parameter information
envi_position_pca <- as.data.frame(fit_all$vectors$arrows)
envi_position_pca
NGRIP_PC1 <- envi_position_pca["NGRIP", "PC1"]
NGRIP_PC2 <- envi_position_pca["NGRIP", "PC2"]
Cala_PC1 <- envi_position_pca["Calanidae", "PC1"]
Cala_PC2 <- envi_position_pca["Calanidae", "PC2"]
het_PC1 <- envi_position_pca["het_pro", "PC1"]
het_PC2 <- envi_position_pca["het_pro", "PC2"]

# check summary$cont and extract explained proportion for PC1 and PC2
axis_explained <- as.data.frame(summary$cont$importance[, 1:2])

PC1_explain <- paste(format(round(axis_explained$PC1[2]*100, 2), nsmall = 2), "%", sep = "")
PC1_explain

PC2_explain <- paste(format(round(axis_explained$PC2[2]*100, 2), nsmall = 2), "%", sep = "")
PC2_explain


PCA_envi <- ggplot(ASV_position_final, aes(x = PC1, y = PC2)) +
  geom_point(data = ASV_position_final, aes(x = PC1, y = PC2), size = 1, stroke = 0.8)+
  
  geom_segment(aes(x = 0, y = 0, xend = NGRIP_PC1 , yend = NGRIP_PC2), colour = "gray14", # NGRIP
               arrow = arrow(length = unit(0.35, "cm")), size = 0.5) +
  annotate("text", x = NGRIP_PC1+0.03, y = NGRIP_PC2-0.03, label = expression(paste(delta^"18", "O NGRIP")), colour = "gray14", angle = 0, size = 5, fontface = "bold") +
  
  geom_segment(aes(x = 0, y = 0, xend = Cala_PC1, yend = Cala_PC2), colour = "gray14", # Calanidae
               arrow = arrow(length = unit(0.35, "cm")),  size = 0.5) +
  annotate("text", x = Cala_PC1-0.03, y = Cala_PC2+0.03, label = "Calanidae",colour = "gray14", angle = 0, size = 5, fontface = "bold") +
  
  geom_segment(aes(x = 0, y = 0, xend = het_PC1 , yend = het_PC2), colour = "gray14", # Salpingoecidae
               arrow = arrow(length = unit(0.35, "cm")),  size = 0.5) +
  
  annotate("text", x = het_PC1+0.03, y = het_PC2-0.03, label = "heterotrophic protists",colour = "gray14", angle = 0, size = 5, fontface = "bold") +
  geom_point(data = age_position_pca, aes(x = PC1, y = PC2, shape = group),size = 4) + # age
  geom_text_repel(data = ASV_position_final,
                  aes(label=rownames(ASV_position_final),
                      point.size = 10,
                      colour = shape,
                      segment.color = "black"),
                  min.segment.length = 0.1, seed = 45, box.padding = 0.1)+
  scale_color_manual(values = c("#FF996FFF", "#a0ceFFFF")) + # colors
  scale_shape_manual(values=c(18, 1, 13)) + 
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20, colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_x_continuous(name = PC1_explain) +
  expand_limits(x = c(0.8, 0.6)) +
  scale_y_continuous(name = PC2_explain) +
  geom_vline(xintercept = 0, size = 0.4, linetype = "dashed", colour = adjustcolor("gray26", alpha.f = 0.7)) +
  geom_hline(yintercept = 0, size = 0.4, linetype = "dashed", colour = adjustcolor("gray26", alpha.f = 0.7))
PCA_envi

ggsave(PCA_envi, file = "Output/07_PCA_RDA/Metabarcoding_PCA_NGRIP_Zooplankton_StratFamilies_FamScalingI_envi_notrans.pdf", width = 21, height = 16, units = "cm")


#########################################
###    (G) VARIATION PARTITIONING     ###
#########################################

names(envi_data)
envi_final <- envi_data[, -1] # remove age from final dataset

# all variables included
rda_all <- rda(seq_data_filtered ~., data = envi_final)
all_variance <- RsquareAdj (rda_all)$adj.r.squared

names(envi_final)
varp <- varpart (seq_data_filtered, ~NGRIP, ~Calanidae, ~het_pro, data = envi_final)
varp

plot (varp, digits = 2, Xnames = c('NGRIP', 'Calanidae', 'het_pro'), bg = c('navy', 'tomato', 'yellow'))


### TEST FOR SIGNIFICANCE ########################################

# (https://www.davidzeleny.net/anadat-r/doku.php/en:varpart_examples):
# Now, when we know both simple and conditional effect of each variables,
# we may want to know whether these variances are significant,
# and hence worth of interpreting. Results from varpart contain the column testable
# with logical values indicating whether given fraction is testable or not.
# To test each of them, we will need the models defined above, and the function anova,
# which (if applied on single object resulting from rda or cca method,
# returns Monte Carlo permutation test of the predictor effect).
# For this, we need to first define also partial ordination models with one
# variable as explanatory and the other as covariable (Condition): 

rda_all # all parameters included
rda_Cala <- rda (seq_data_filtered ~ Calanidae, data = envi_final) # Calanidae
rda_hetprot <- rda (seq_data_filtered ~ het_pro, data = envi_final) # heterotrophic protists
rda_ngrip <- rda (seq_data_filtered ~ NGRIP, data = envi_final) # NGRIP

# ANOVA
set.seed(123); anova(rda_all) # test significance of conditional effect
set.seed(123); anova(rda_Cala) # test simple (marginal) effect of Calanidae
set.seed(123); anova(rda_hetprot) # test simple (marginal) effect of heterotrophic protists
set.seed(123); anova(rda_ngrip) # test simple (marginal) effect of NGRIP

