setwd("")

library("tidyverse")
library("scales")
library("readxl")
library("openxlsx")

rm(list = ls())

dir.create("Output")
dir.create("Output/01_Subset_Phytoplankton")

########################################################################################################################################################
####################### 1 Creating Subset of data ######################################################################################################
########################################################################################################################################################


#############################################
###         CREATING THE DATASET          ###
#############################################

# read in data:
#KL-77-1_dataset10.2_lineageDB_age_depth_adj.xlsx and KL-77-2_dataset10.2_lineageDB_age_depth_adj.xlsx

# in these files, all rows with rank smaller than F (from F1 downwards) were deleted manually (order rank column in Excel and delete all rows with rank smaller F)
# otherwiese "reads.sums" will be counted several times and resulting read count per family will be artificially increased
# the modified Excel-files (name: "x_adj.xlsx") is then loaded as raw data set.
# Two datasets are loaded because samples were sequences in two runs:


# DATASET 1
dataset1 <- read.xlsx(("KL-77-1_nt0.2_lineageDB_age_depth_adj.xlsx"))
dataset1$age <- as.numeric(dataset1$age) #make numeric
dataset1 <- arrange(dataset1, age) #sort

# DATASET 2
dataset2 <- read.xlsx(("KL77-2_nt0.2_lineageDB_age_depth_adj.xlsx"))
dataset2$age <- as.numeric(dataset2$age)
dataset2 <- arrange(dataset2, age) 

# COMBINE DATASETS
alldata <- rbind(dataset1, dataset2)
alldata <- arrange(alldata, age) 

# remove "undetermined" from samples
data <- alldata[!grepl("Undetermined", alldata$sample),]

#############################################
###               OVERVIEW                ###
#############################################

# make rank a factor
class(data$rank)
data$rank <- as.factor(data$rank)
class(data$rank)

data_R <- subset(data, subset = rank %in% "R") # classified reads (root)
classified_reads <- sum(data_R$reads.sum)
data_U <- subset(data, subset = rank %in% "U") # unclassified reads
unclassified_reads <- sum(data_U$reads.sum)

reads_tot <- classified_reads + unclassified_reads # total reads
reads_tot

percent_classified_reads <- (classified_reads/reads_tot)*100
percent_classified_reads

# classified reads in samples
sample_classified_reads <- subset(data, subset = type %in% "sample")
sample_R <- subset(sample_classified_reads, subset = rank %in% "R")
sample_classified_reads <- sum(sample_R$reads.sum) # total classified reads in samples
percent_samples_classified_reads <- (sample_classified_reads/reads_tot)*100 # % of raw reads

# classified reads in library blanks
LB_classified_reads <- subset(data, subset = extract %in% "LB")
LB_R <- subset(LB_classified_reads, subset = rank %in% "R")
LB_classified_reads <- sum(LB_R$reads.sum)
percent_LB_classified_reads <- (LB_classified_reads/reads_tot)*100

# classified reads in extraction blanks
EB_classified_reads <- subset(data, subset = type %in% "blank")
EB_classified_reads <- subset(EB_classified_reads, subset = !extract %in% "LB")
EB_R <- subset(EB_classified_reads, subset = rank %in% "R")
EB_classified_reads <- sum(EB_R$reads.sum)
percent_EB_classified_reads <- (EB_classified_reads/reads_tot)*100

# total number of unique families
length(unique(alldata$family))

# counts on family level
family_read_count <- subset(alldata, subset = !is.na(alldata$family))
family_read_count <- subset(family_read_count, subset = rank %in% "F")
sum(family_read_count$reads.sum)


#############################################
###   SUBSETTING PHYTOPLANKTON FAMILIES   ###
#############################################

# import list of phytoplankton families that will be filtered for
fams <- read_delim("phyto_families.txt", delim = "\t", col_names = F)
fams_v <- unlist(fams) # vectorize
is.vector(fams_v)

phyto_fams <- subset(data, subset = family %in% fams_v) #subsetting

#summarize read counts over family even if lower taxa is known
phytoplankton_data <- phyto_fams %>%
  filter(!is.na(family)) %>%
  group_by(family, age) %>%
  summarise(sumcount_family = sum(reads.sum))

###turn to wide format
phytoplankton_wide <- phytoplankton_data %>%
  pivot_wider(names_from = family,
              values_from = sumcount_family)

phytos_by_age <- arrange(phytoplankton_wide, age)


# prepare dataset for export
phytos_by_age_final <- as.data.frame(t(phytos_by_age)) # transpose to age is in columns
phytos_by_age_final[is.na(phytos_by_age_final)] <- 0 # replace NAs by 0

ages <- phytos_by_age_final[1, ] # vector with ages
colnames(phytos_by_age_final) <- ages # change column names to age
phytos_by_age_final <- phytos_by_age_final[-1, ] # remove first column where age was before
ncol(phytos_by_age_final) # number of rows 
colnames(phytos_by_age_final) [43] <- "Blank" # name last row "Blank"

fam_names <- rownames(phytos_by_age_final) # vector with family names
phytos_by_age_final$family <- fam_names # create column with family names
rownames(phytos_by_age_final) <- NULL # change column name to "family"
ncol(phytos_by_age_final)
phytos_by_age_final <- phytos_by_age_final[, c(44, 1:43)] # rearrange so family names are first column

# export phytoplankton family counrs without taxonomic information
write.xlsx(phytos_by_age_final, file = "Output/01_Subset_Phytoplankton/Phytoplankton_notax.xlsx")

# total phytoplankton family counts
nrow(phytos_by_age_final)
all_phyto_fam_counts <- sum(phytos_by_age_final[, 2:(ncol(phytos_by_age_final)-1)]) # -1 for not counting blank reads
all_phyto_fam_counts # all phytoplankton counts in samples

### BLANKS
length(unique(phyto_fams$family))
LB_phyto <- subset(phyto_fams, subset = extract %in% "LB") # no reads in library blanks
EB_phyto <- subset(phyto_fams, subset = type %in% "blank")
EB_phyto <- subset(EB_phyto, subset = !extract %in% "LB")
unique(EB_phyto$extract) # check which samples are in the data frame
# StB091, StB107 and StB109 are not extraction blanks but samples (wrong assignment in Kraken)
EB_phyto <- subset(EB_phyto, subset = !extract %in% c("StB105", "StB107", "StB091"))
EB_phyto_reads <- sum(EB_phyto$reads.sum) # total reads in extraction blanks after phytoplankton filter


#############################################
###           GROUP INFORMATION           ###
#############################################

groups <- read.xlsx("Group_infos_phyto.xlsx")

phyto_groups <- merge(phytos_by_age_final, groups, by = "family")
ncol(phyto_groups)
names(phyto_groups)
phyto_groups <- phyto_groups[, c(1, 45, 2:44)] # rearrange dataframe and remove column 46 (colour information will be needed later)

write.xlsx(phyto_groups, file = "Output/01_Subset_Phytoplankton/Phytoplankton_groups_notax.xlsx")


#############################################
###             GROUP SUMMARY             ###
#############################################

library("readxl")
library("openxlsx")

rm(list = ls())

phyto <- read.xlsx("Output/01_Subset_Phytoplankton/Phytoplankton_groups_notax.xlsx")

class(phyto$group) # make group a factor
phyto$group <- as.factor(phyto$group)
class(phyto$group)

names(phyto)

# total reads
tot_reads <- sum(phyto[, 3:ncol(phyto)])
tot_reads

# reads per group
levels(phyto$group)
chloro <- subset(phyto, subset = group %in% "chlorophytes") # chlorophytes
bac <- subset(phyto, subset = group %in% "phototrophic bacteria") # phototrophic bacteria
pro <- subset(phyto, subset = group %in% "phototrophic protists") # phototrophic protists
alg <- subset(phyto, subset = group %in% "red algae") # red algae

chloro_reads <- sum(chloro[, 3:ncol(chloro)]) # reads chlorophytes
bac_reads <- sum(bac[, 3:ncol(bac)]) # reads phototrophic bacteria
pro_reads <- sum(pro[, 3:ncol(pro)]) # reads phototrophic protists
alg_reads <- sum(alg[, 3:ncol(alg)]) # reads red algae

chloro_reads
bac_reads
pro_reads
alg_reads

# relative to total reads
rel_chloro_reads <- (chloro_reads/tot_reads)*100
rel_bac_reads <- ((bac_reads/tot_reads)*100)
rel_pro_reads <- ((pro_reads/tot_reads)*100)
rel_algae_reads <- (alg_reads/tot_reads)*100

rel_chloro_reads # relative reads chlorophytes
rel_bac_reads # relative reads phototrophic bacteria
rel_pro_reads # relative reads phototrophic protists
rel_algae_reads  # relative reads red algae


########################################################################################################################################################
####################### 2 Rarefaction ##################################################################################################################
########################################################################################################################################################

rm(list = ls())

dir.create("Output/02_Rarefaction")

options(stringsAsFactors=FALSE)

library("openxlsx")


#############################################
###         PREPARE DATA FRAME            ###
#############################################

# read data and check visually if the df is as requested
# requirements: (1) blanks deleted
#               (2) "x" infront of the age, which are the column; otherwise they are not recognized as the column names

# import data
data <- read.xlsx("Output/01_Subset_Phytoplankton/Phytoplankton_notax.xlsx")

# add the "x" infront of the column names (age)
names(data)
ncol(data)
data_temp <- data[, 2:43] # remove columns where no "x" will be added (family and blank column)
names(data_temp)
ages <- names(data_temp) # extract ages
ages_mod <- paste("x", ages, sep = "") # add the "x"
ages_mod
colnames(data_temp) <- ages_mod # add modified ages to data frame
data_temp$family <- data$family # add family again to data frame
ncol(data_temp)
data_temp <- data_temp[, c(43, 1:42)] # rearage so family is first column


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
mincounts=sort(apply(specseq[2:43], 2, sum))[1]
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
save(genrare, file="Output/02_Rarefaction/rarefaction_pelagic.RDATA")


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
write.xlsx(totspec, file = "Output/02_Rarefaction/Rarefaction_Phytoplankton_FamiliesPerSample.xlsx")

str(famreads)	# count counts per families and sample/year

write.xlsx(famreads, file = "Output/02_Rarefaction/Rarefaction_Phytoplankton_FamilyCounts.xlsx")

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
write.xlsx(famdf, "Output/02_Rarefaction/Rarefaction_Phytoplankton_FamiliesCI95.xlsx")

#############################################
###       MEAN FAMILIES PER SAMPLE        ###
#############################################

rm(list = ls())

library("readxl")
library("openxlsx")

# load data coming from Rarefaction Script
phyto_rar <- read_excel("Output/02_Rarefaction/Rarefaction_Phytoplankton_FamiliesPerSample.xlsx")
colnames(phyto_rar) <- c("age", "SampleEff", "Nfamilies")

phyto_rar$age <- sub("x", "", phyto_rar$age) # removes the x infront of the age
phyto_rar$age <- factor(phyto_rar$age, levels = unique(phyto_rar$age)) # age as factor
is.factor(phyto_rar$age)


# calculate mean for each year (100 repeats for rarefaction)
nfam_mean <- data.frame(c(1:length(levels(phyto_rar$age))))
for(i in c(1:nrow(phyto_rar))) {
  nfam_mean[i, ] <- mean(phyto_rar$Nfamilies[phyto_rar$age == levels(phyto_rar$age)[i]])
}
nfam_mean <- na.omit(nfam_mean) # delete rows with NA

nfam_mean$age <- levels(phyto_rar$age)
nfam_mean <- nfam_mean[, c(2, 1)]
colnames(nfam_mean) <- c("age", "Nfamilies_mean")

write.xlsx(nfam_mean, file = "Output/02_Rarefaction/Rarefaction_Phytoplankton_FamiliesPerSample_mean.xlsx")


#############################################
###   MEAN COUNTS PER SAMPLE AND FAMILY   ###
#############################################

rm(list = ls())

library("readxl")
library("openxlsx")

# load data coming from Rarefaction Script
photos <- read_excel("Output/02_Rarefaction/Rarefaction_Phytoplankton_FamilyCounts.xlsx")

names(photos)
colnames(photos)[1] <- "age" # rename age column
photos <- photos[, -2] # removes sampleEff column

photos$age <- factor(photos$age, levels = unique(photos$age)) # make age a factor
is.factor(photos$age)
length(levels(photos$age)) # how many age levels (= samples)

# get subset of all resampling repeats per sample age (100 repeats for rarefaction)
a_photos <- subset(photos, subset = age %in% levels(photos$age)[1])
b_photos <- subset(photos, subset = age %in% levels(photos$age)[2])
c_photos <- subset(photos, subset = age %in% levels(photos$age)[3])
d_photos <- subset(photos, subset = age %in% levels(photos$age)[4])
e_photos <- subset(photos, subset = age %in% levels(photos$age)[5])
f_photos <- subset(photos, subset = age %in% levels(photos$age)[6])
g_photos <- subset(photos, subset = age %in% levels(photos$age)[7])
h_photos <- subset(photos, subset = age %in% levels(photos$age)[8])
i_photos <- subset(photos, subset = age %in% levels(photos$age)[9])
j_photos <- subset(photos, subset = age %in% levels(photos$age)[10])
k_photos <- subset(photos, subset = age %in% levels(photos$age)[11])
l_photos <- subset(photos, subset = age %in% levels(photos$age)[12])
m_photos <- subset(photos, subset = age %in% levels(photos$age)[13])
n_photos <- subset(photos, subset = age %in% levels(photos$age)[14])
o_photos <- subset(photos, subset = age %in% levels(photos$age)[15])
p_photos <- subset(photos, subset = age %in% levels(photos$age)[16])
q_photos <- subset(photos, subset = age %in% levels(photos$age)[17])
r_photos <- subset(photos, subset = age %in% levels(photos$age)[18])
s_photos <- subset(photos, subset = age %in% levels(photos$age)[19])
t_photos <- subset(photos, subset = age %in% levels(photos$age)[20])
u_photos <- subset(photos, subset = age %in% levels(photos$age)[21])
v_photos <- subset(photos, subset = age %in% levels(photos$age)[22])
w_photos <- subset(photos, subset = age %in% levels(photos$age)[23])
x_photos <- subset(photos, subset = age %in% levels(photos$age)[24])
y_photos <- subset(photos, subset = age %in% levels(photos$age)[25])
z_photos <- subset(photos, subset = age %in% levels(photos$age)[26])
aa_photos <- subset(photos, subset = age %in% levels(photos$age)[27])
ab_photos <- subset(photos, subset = age %in% levels(photos$age)[28])
ac_photos <- subset(photos, subset = age %in% levels(photos$age)[29])
ad_photos <- subset(photos, subset = age %in% levels(photos$age)[30])
ae_photos <- subset(photos, subset = age %in% levels(photos$age)[31])
af_photos <- subset(photos, subset = age %in% levels(photos$age)[32])
ag_photos <- subset(photos, subset = age %in% levels(photos$age)[33])
ah_photos <- subset(photos, subset = age %in% levels(photos$age)[34])
ai_photos <- subset(photos, subset = age %in% levels(photos$age)[35])
aj_photos <- subset(photos, subset = age %in% levels(photos$age)[36])
ak_photos <- subset(photos, subset = age %in% levels(photos$age)[37])
al_photos <- subset(photos, subset = age %in% levels(photos$age)[38])
am_photos <- subset(photos, subset = age %in% levels(photos$age)[39])
an_photos <- subset(photos, subset = age %in% levels(photos$age)[40])
ao_photos <- subset(photos, subset = age %in% levels(photos$age)[41])
ap_photos <- subset(photos, subset = age %in% levels(photos$age)[42])

# calculate mean for each sampling age
a <- colMeans(a_photos[, 2:ncol(a_photos)])
b <- colMeans(b_photos[, 2:ncol(b_photos)])
c <- colMeans(c_photos[, 2:ncol(c_photos)])
d <- colMeans(d_photos[, 2:ncol(d_photos)])
e <- colMeans(e_photos[, 2:ncol(e_photos)])
f <- colMeans(f_photos[, 2:ncol(f_photos)])
g <- colMeans(g_photos[, 2:ncol(g_photos)])
h <- colMeans(h_photos[, 2:ncol(h_photos)])
i <- colMeans(i_photos[, 2:ncol(i_photos)])
j <- colMeans(j_photos[, 2:ncol(j_photos)])
k <- colMeans(k_photos[, 2:ncol(k_photos)])
l <- colMeans(l_photos[, 2:ncol(l_photos)])
m <- colMeans(m_photos[, 2:ncol(m_photos)])
n <- colMeans(n_photos[, 2:ncol(n_photos)])
o <- colMeans(o_photos[, 2:ncol(o_photos)])
p <- colMeans(p_photos[, 2:ncol(p_photos)])
q <- colMeans(q_photos[, 2:ncol(q_photos)])
r <- colMeans(r_photos[, 2:ncol(r_photos)])
s <- colMeans(s_photos[, 2:ncol(s_photos)])
t <- colMeans(t_photos[, 2:ncol(t_photos)])
u <- colMeans(u_photos[, 2:ncol(u_photos)])
v <- colMeans(v_photos[, 2:ncol(v_photos)])
w <- colMeans(w_photos[, 2:ncol(w_photos)])
x <- colMeans(x_photos[, 2:ncol(x_photos)])
y <- colMeans(y_photos[, 2:ncol(y_photos)])
z <- colMeans(z_photos[, 2:ncol(z_photos)])
aa <- colMeans(aa_photos[, 2:ncol(aa_photos)])
ab <- colMeans(ab_photos[, 2:ncol(ab_photos)])
ac <- colMeans(ac_photos[, 2:ncol(ac_photos)])
ad <- colMeans(ad_photos[, 2:ncol(ad_photos)])
ae <- colMeans(ae_photos[, 2:ncol(ae_photos)])
af <- colMeans(af_photos[, 2:ncol(af_photos)])
ag <- colMeans(ag_photos[, 2:ncol(ag_photos)])
ah <- colMeans(ah_photos[, 2:ncol(ah_photos)])
ai <- colMeans(ai_photos[, 2:ncol(ai_photos)])
aj <- colMeans(aj_photos[, 2:ncol(aj_photos)])
ak <- colMeans(ak_photos[, 2:ncol(ak_photos)])
al <- colMeans(al_photos[, 2:ncol(al_photos)])
am <- colMeans(am_photos[, 2:ncol(am_photos)])
an <- colMeans(an_photos[, 2:ncol(an_photos)])
ao <- colMeans(ao_photos[, 2:ncol(ao_photos)])
ap <- colMeans(ap_photos[, 2:ncol(ap_photos)])

# combine means per age in one data frame
countfam_mean <- as.data.frame(rbind(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab, ac, ad, ae, af, ag, ah, ai, aj, ak, al, am, an, ao, ap))

family <- names(countfam_mean) # save vector with families for later

age <- levels(photos$age)
countfam_mean$age <- age
ncol(countfam_mean)
countfam_mean <- countfam_mean[, c(67, 1:66)]
rownames(countfam_mean) <- NULL

countfam_mean$age <- sub("x", "", countfam_mean$age) # removes the x infront of the age
countfam_mean <- as.data.frame(t(countfam_mean)) # transposing so age is in columns
age <- countfam_mean[1, ] # vector with ages
colnames(countfam_mean) <- age # ages to column names
countfam_mean <- countfam_mean[-1, ] # remove row with ages

countfam_mean <- as.data.frame(sapply(countfam_mean, as.numeric)) # counts numeric

countfam_mean$family <- family # add column with family
ncol(countfam_mean)
countfam_mean <- countfam_mean[, c(43, 1:42)] # make family the first column
rownames(countfam_mean) <- NULL # remove row names

write.xlsx(countfam_mean, file = "Output/02_Rarefaction/Rarefaction_Phytoplankton_FamilyCounts_means.xlsx")


########################################################################################################################################################
####################### 3 STRATIGRAMS ##################################################################################################################
########################################################################################################################################################

rm(list = ls())

dir.create("Output/03_Stratigram")

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
photos <- read_excel("Output/02_Rarefaction/Rarefaction_Phytoplankton_FamilyCounts_means.xlsx")

# create data frame with only columns with counts
ncol(photos)
names(photos)
photos_onlycounts <- photos[, -1] # remove columns with family information, so only counts are left

# create data frame with relative counts per column
photos_rel <- data.frame(c(1:nrow(photos_onlycounts)))
for( i in c(1:ncol(photos_onlycounts))) {
  photos_rel[, i] <- data.frame(photos_onlycounts[, i]) %>%
    transmute(rel_ab = (photos_onlycounts[, i]/sum((photos_onlycounts[, i])))) # calculation for relative abundance for each cell of the data frame
}
photos_rel <- photos_rel*100 # for making relative abundance in % (somehow it makes the column names weird if yoh do it right in the for-loop)

# add taxonomic information again
photos_rel <- cbind(photos_rel, photos$family) # add family column and class column
names(photos_rel)
ncol(photos_rel)
colnames(photos_rel)[43] <- "family" # rename family column correctly

# import group information
group_info <- read.xlsx("Group_infos_phyto.xlsx")

photos_rel <- merge(photos_rel, group_info, by = "family")

ncol(photos_rel)
names(photos_rel)
photos_rel <- photos_rel[, c(1, 44, 2:43)] # rearrange so family names are first column

write.xlsx(photos_rel, file = "Output/03_Stratigram/Phytoplankton_rel_abund.xlsx")


# summarize counts of groups
tot_chloro <- colSums(subset(photos_rel, subset = group %in% "chlorophytes")[, 3:ncol(photos_rel) ])
tot_photoprot <- colSums(subset(photos_rel, subset = group %in% "phototrophic protists")[, 3:ncol(photos_rel) ])
tot_photobac <- colSums(subset(photos_rel, subset = group %in% "phototrophic bacteria")[, 3:ncol(photos_rel) ])

tot_chloro <- c("total chlorophytes", "chlorophytes", tot_chloro)
tot_photoprot <- c("total phototrophic protists", "phototrophic protists", tot_photoprot)
tot_photobac <- c("total phototrophic bacteria", "phototrophic bacteria", tot_photobac)

# add rel. abundance of total counts of chlorophytes,
# phototrophoc protsits and phototrophic bacteria to the dataset
photos_rel <- rbind(photos_rel,
                    tot_chloro,
                    tot_photoprot,
                    tot_photobac)

# change table from wide onto long format (needed for plotting)
photos_long <- pivot_longer(photos_rel,
                            cols = !c(family, group),
                            names_to = c("age"),
                            values_to = "rel_abund")

# make age numeric
photos_long$age <- as.numeric(photos_long$age)

#############################################
###             (B) FILTERING             ###
#############################################

# For plotting I want to keep only families that occur in at least 20 samples and in at least 1 sample with a relative abundance >= 3%

# 1. remove all rows (in long format) of families occurring less than 20 times
long_biggerzero <-  filter(photos_long, rel_abund > 0) # subset of long format data with only relative abundance > 0
fam_biggerzero <- long_biggerzero$family # only family names
fam_counts <- as.data.frame(table(fam_biggerzero)) # counts how often families occur with relative abundance > 0
threshold_occurrence <- filter(fam_counts, Freq >= 20)$fam_biggerzero # subset of families that occur more than 20 times
threshold_occurrence <- as.character(threshold_occurrence) # make it a vector of characters

photos_filter1 <- subset(photos_long, subset = family %in% threshold_occurrence)

# 2. remove all rows (in long format) of families that don't occur in at least 1 sample with min. 3%
long_biggerpro <-  filter(photos_filter1, rel_abund >= 3) # subset of already filtered data set to select rows that have a relative >= 3% 
fam_biggerpro <- long_biggerpro$family # only family names 
fam_countspro <- as.data.frame(table(fam_biggerpro)) # counts how often families occur with relative abundance > 5
threshold_relabund <- filter(fam_countspro, Freq >= 1)$fam_biggerpro # subset of families that occur in at least 1 sample with relative abundance of 3%
threshold_relabund <- as.character(threshold_relabund) # make it a vector of characters

photos_filter2 <- subset(photos_filter1, subset = family %in% threshold_relabund)

# bring data into order in which they are plotted later
photos_filter2 <- photos_filter2[order(photos_filter2[, 3]),] # order first by column 3 (= age)

# for rel. abundance of main filtered families
photos_main <- pivot_wider(photos_filter2,
                           names_from = "age",
                           values_from = "rel_abund")

names(photos_main)
nrow(photos_main)

write.xlsx(photos_main, file = "Output/03_Stratigram/Phytoplankton_main_rel_abund.xlsx")


#############################################
###      (C) PLOTTING MAIN FAMILIES       ###
#############################################

# preparation of data
photos_filter2$age <- photos_filter2$age/1000 # make age in ka BP
names(photos_filter2)

photos_filter2$group <- factor(photos_filter2$group, levels = unique(photos_filter2$group)) # make group a factor with order it has right now


# add color information by subsetting groups, add color column and combine to one data frame again
levels(photos_filter2$group)
photos_filter2_chloro <- subset(photos_filter2, subset = group %in% "chlorophytes")
photos_filter2_chloro$co <- rep("green3")
photos_filter2_phopro <- subset(photos_filter2, subset = group %in% "phototrophic protists")
photos_filter2_phopro$co <- rep("olivedrab4")
photos_filter2_phobac <- subset(photos_filter2, subset = group %in% "phototrophic bacteria")
photos_filter2_phobac$co <- rep("cyan3")

photos_filter2_final <- rbind(photos_filter2_chloro,
                              photos_filter2_phopro,
                              photos_filter2_phobac)

# order
names(photos_filter2_final)

photos_filter2_final$rel_abund <- as.numeric(photos_filter2_final$rel_abund) # make re_abund numeric
photos_filter2_final$age <- as.numeric(photos_filter2_final$age) # male age numeric
photos_filter2_final$family <- factor(photos_filter2_final$family, levels = unique(photos_filter2_final$family)) # make family a factor with order it has right now
photos_filter2_final$group <- factor(photos_filter2_final$group, levels = unique(photos_filter2_final$group)) # make group a factor with order it has right now


# photos plot:
photos_plot <- ggplot(photos_filter2_final, aes(x = age, y = rel_abund, fill = group)) +
  geom_area(aes(fill = co)) +
   geom_segment(aes(xend = age, yend = 0), lwd = 0.5, color = "black") +
  expand_limits(x = 0, y = 0) + # axis start with zero
  scale_y_continuous(n.breaks = 3, # numbers mark y axis
                     expand = c(0,0),
                     name = "Relative abundance [%] \nbased on total reads of phytoplankton families") + 
  facet_abundanceh(vars(family), rotate_facet_labels = 70, space = "fixed", scales = "free_x") + # groups plots by family, makes every plot the same size (different scale)
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, color = c(NA, "black")), # rotates x axis by 90 degree, labels in middle of tick (vjust)
        axis.title.x = element_text(vjust = 0.5, size = 14),
        axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14, colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(colour = "black"),
        panel.grid = element_blank(),
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

photos_plot
ggsave(photos_plot, file = "Output/03_Stratigram/Phytoplankton_area.pdf", width = 33, height = 20, units = "cm")


#############################################
###               (D) GROUPS              ###
#############################################

rm(list = ls())

photos_rel <- read.xlsx("Output/03_Stratigram/Phytoplankton_rel_abund.xlsx")

names(photos_rel)
class(photos_rel$group)
photos_rel$group <- as.factor(photos_rel$group) # groups as factors
levels(photos_rel$group)

chloro <- subset(photos_rel, subset = group %in% "chlorophytes")
bac <- subset(photos_rel, subset = group %in% "phototrophic bacteria")
pro <- subset(photos_rel, subset = group %in% "phototrophic protists")
redal <- subset(photos_rel, subset = group %in% "red algae")

# create .xlsx with groups in different sheets
sheets <- list("Chlorophytes" = chloro,
               "Bacteria" = bac,
               "Protists" = pro,
               "Red algae" = redal)
write.xlsx(sheets, file = "Output/03_Stratigram/Phytoplankton_Groups.xlsx")


########################################################################################################################################################
####################### 4 ENVIRONMENTAL PARAMETERS #####################################################################################################
########################################################################################################################################################

rm(list = ls())

dir.create("Output/04_Environment")

library("openxlsx")

#############################################
###                  NGRIP                ###
#############################################

# read data
all_NGRIP <- read.xlsx("NGRIP_data.xlsx")

names(all_NGRIP)
all_NGRIP <- all_NGRIP[order(all_NGRIP$age), ] # order by age

# read in sample ages
samples <- read.xlsx("Output/02_Rarefaction/Rarefaction_Phytoplankton_FamiliesPerSample_mean.xlsx")
names(samples)
samples <- as.data.frame(samples[, "age"]) # get ages from samples
names(samples) <- "sample_age"
class(samples$sample_age)
samples$sample_age <- as.numeric(samples$sample_age) # age numeric
class(samples$sample_age)
samples$sample_age <- samples$sample_age/1000 # age in ka BP

# sample ages
samples$sample_age

# create empty column that will be filled with interpolated NGRIP data
samples$NGRIP_reg <- rep(NA, nrow(samples))

# LINEAR REGRESSION
for (index_SampleAge in 1:length(samples$sample_age)) { # goes through sample_age
  
  # difference from sample_age to all_NGRIP$age
  # get ages for interpolation
  diff <- samples$sample_age[index_SampleAge] - all_NGRIP$age
  diff_x2 <- max(diff[diff <= 0]) # gives age from all_NGRIP that is negative and closest to zero = older than sample
  diff_x1 <- min(diff[diff >= 0]) # gives age from all_NGRIP that is positive and closest to zero = younger than sample
  index_x2 <- which(diff == diff_x2) # index of all NGRIP$age for age2
  index_x1 <- which(diff == diff_x1) # index of all_NGRIP$age for age1
  #print(paste("index_x2:", index_x2, "index_x1:", index_x1, sep = " "))
  x2 <- all_NGRIP$age[index_x2] # age2
  x1 <- all_NGRIP$age[index_x1] #age1
  print(paste("age2:", x2, "age1:", x1, sep = " "))
  
  # get NGRIP-values for interpolation
  y2 <- all_NGRIP$d18O[index_x2]
  y1 <- all_NGRIP$d18O[index_x1]
  print(paste("d18O:", y2, "d18O:", y1, sep = " "))
  
  # slope m = (y2 - y1) / (x2 - x1)
  m <- (y2 - y1) / (x2 - x1)
  print(paste("slope:", m, sep = " "))
  
  # intercept n = y - mx
  n <- y2 - (m * x2)
  print(paste("intercept:", n, sep = " "))
  
  if (x2 != x1) {
    
    # linear interpolation: y = mx + n
    samples$NGRIP_reg[index_SampleAge] <- m * samples$sample_age[index_SampleAge] + n
  }
  
  else {samples$NGRIP_reg[index_SampleAge] <- all_NGRIP$d18O[index_x2]
  
  }
}

write.xlsx(samples, file = "Output/04_Environment/Shotgun_NGRIP.xlsx")


########################################################################################################################################################
####################### 5 ZOOPLANKTON ##################################################################################################################
########################################################################################################################################################

rm(list = ls())

dir.create("Output/05_Zooplankton")

library("tidyverse")
library("scales")
library("readxl")
library("openxlsx")

#############################################
###       (A) CREATING DATA SET           ###
#############################################

# DATASET 1
dataset1 <- read.xlsx(("KL-77-1_nt0.2_lineageDB_age_depth_adj.xlsx"))
dataset1$age <- as.numeric(dataset1$age) #make numeric
dataset1 <- arrange(dataset1, age) #sort

# DATASET 2
dataset2 <- read.xlsx(("KL77-2_nt0.2_lineageDB_age_depth_adj.xlsx"))
dataset2$age <- as.numeric(dataset2$age)
dataset2 <- arrange(dataset2, age)

# COMBINE DATASETS
alldata <- rbind(dataset1, dataset2)
alldata <- arrange(alldata, age) 

# remove "undetermined" from samples
data <- alldata[!grepl("Undetermined", alldata$sample),]

#############################################
###               OVERVIEW                ###
#############################################

# total number of unique families
length(unique(data$family))

# counts on family level
family_read_count <- subset(data, subset = !is.na(data$family))
family_read_count <- subset(family_read_count, subset = rank %in% "F")
sum(family_read_count$reads.sum)

# reads of zooplankton families
fams <- read_delim("zoo_families.txt", delim = "\t", col_names = F)
fams_v <- unlist(fams) # vectorize
is.vector(fams_v)

#############################################
###    SUBSETTING ZOOPLANKTON FAMILIES    ###
#############################################

zoo_fams <- subset(data, subset = family %in% fams_v)

#summarize read counts over family even if lower taxa is known
zoo_data <- zoo_fams %>%
  filter(!is.na(family)) %>%
  group_by(family, age) %>%
  summarise(sumcount_family = sum(reads.sum))

###turn to wide format if needed
zoo_wide <- zoo_data %>%
  pivot_wider(names_from = family,
              values_from = sumcount_family)

zoo_by_age <- arrange(zoo_wide, age)

# prepare dataset for export
zoo_by_age_final <- as.data.frame(t(zoo_by_age)) # transpose to age is in columns
zoo_by_age_final[is.na(zoo_by_age_final)] <- 0 # replace NAs by 0

ages <- zoo_by_age_final[1, ] # vector with ages
colnames(zoo_by_age_final) <- ages # change column names to age
zoo_by_age_final <- zoo_by_age_final[-1, ] # remove first column where age was before
ncol(zoo_by_age_final) # number of rows 

fam_names <- rownames(zoo_by_age_final) # vector with family names
zoo_by_age_final$family <- fam_names # create column with family names
rownames(zoo_by_age_final) <- NULL # change column name to "family"
ncol(zoo_by_age_final)
zoo_by_age_final <- zoo_by_age_final[, c(43, 1:42)] # rearrange so family names are first column

write.xlsx(zoo_by_age_final, file = "Output/05_Zooplankton/Zooplankton_notax.xlsx")

# total zooplankton family counts
nrow(zoo_by_age_final)
all_zoo_fam_counts <- sum(zoo_by_age_final[, 2:(ncol(zoo_by_age_final))])
all_zoo_fam_counts # all zooplankton counts in samples

### BLANKS
length(unique(zoo_fams$family))
LB_zoo <- subset(zoo_fams, subset = extract %in% "LB") # no reads in library blanks
EB_zoo <- subset(zoo_fams, subset = type %in% "blank")
EB_zoo <- subset(EB_zoo, subset = !extract %in% "LB")
unique(EB_zoo$extract) # check which samples are in the data frame
# StB091, StB107 and StB109 are not extraction blanks but samples (wrong assignment in Kraken)
EB_zoo <- subset(EB_zoo, subset = !extract %in% c("StB105", "StB107", "StB091"))
EB_zoo_reads <- sum(EB_zoo$reads.sum) # total reads in extraction blanks after phytoplankton filter


#############################################
###           GROUP INFORMATION           ###
#############################################

groups <- read.xlsx("Group_infos_zoo.xlsx")

zoo_groups <- merge(zoo_by_age_final, groups, by = "family")
ncol(zoo_groups)
zoo_groups <- zoo_groups[, c(1, 45, 2:44)]

write.xlsx(zoo_groups, file = "Output/05_Zooplankton/Zooplankton_groups_notax.xlsx")


#############################################
###           (B) RAREFACTION             ###
#############################################

rm(list = ls())

options(stringsAsFactors=FALSE)

library("openxlsx")

#############################################
###         PREPARE DATA FRAME            ###
#############################################

# read data and check visually if the df is as requested
# requirements: (1) blanks deleted
#               (2) "x" infront of the age, which are the column; otherwise they are not recognized as the column names

# import data
data <- read.xlsx("Output/05_Zooplankton/Zooplankton_notax.xlsx")

# add the "x" infront of the column names (age)
names(data)
ncol(data)
data_temp <- data[, 2:43] # remove columns where no "x" will be added (family and blank column)
names(data_temp)
ages <- names(data_temp) # extract ages
ages_mod <- paste("x", ages, sep = "") # add the "x"
ages_mod
colnames(data_temp) <- ages_mod # add modified ages to data frame
data_temp$family <- data$family # add family again to data frame
ncol(data_temp)
data_temp <- data_temp[, c(43, 1:42)] # rearage so family is first column


#####################################################################################
###  Rarefaction_process_data_family_level_Stefan's codes                         ###
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
mincounts=sort(apply(specseq[2:43], 2, sum))[1]
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
save(genrare, file="Output/05_Zooplankton/rarefaction_pelagic.RDATA")


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

write.xlsx(totspec, file = "Output/05_Zooplankton/Rarefaction_Zooplankton_FamiliesPerSample.xlsx")

str(famreads)	# count counts per families and sample/year

write.xlsx(famreads, file = "Output/05_Zooplankton/Rarefaction_Zooplankton_FamilyCounts.xlsx")

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
write.xlsx(famdf, "Output/05_Zooplankton/Rarefaction_Zooplankton_FamiliesCI95.xlsx")


#############################################
###     MEAN TOTAL COUNTS PER SAMPLE      ###
#############################################

rm(list = ls())

library("readxl")
library("openxlsx")

# load data coming from Rarefaction Script
zoo_rar <- read_excel("Output/05_Zooplankton/Rarefaction_Zooplankton_FamiliesPerSample.xlsx")
colnames(zoo_rar) <- c("age", "SampleEff", "Nfamilies")

zoo_rar$age <- sub("x", "", zoo_rar$age) # removes the x infront of the age
zoo_rar$age <- factor(zoo_rar$age, levels = unique(zoo_rar$age)) # make age a factor
is.factor(zoo_rar$age) # check


# calculate mean for each year (100 repeats for rarefaction)
nfam_mean <- data.frame(c(1:length(levels(zoo_rar$age))))
for(i in c(1:nrow(zoo_rar))) {
  nfam_mean[i, ] <- mean(zoo_rar$Nfamilies[zoo_rar$age == levels(zoo_rar$age)[i]])
}
nfam_mean <- na.omit(nfam_mean) # delete rows with NA

nfam_mean$age <- levels(zoo_rar$age)
nfam_mean <- nfam_mean[, c(2, 1)]
colnames(nfam_mean) <- c("age", "Nfamilies_mean")

write.xlsx(nfam_mean, file = "Output/05_Zooplankton/Rarefaction_Zooplankton_FamiliesPerSample_mean.xlsx")


#############################################
###   MEAN FAMILIES PER SAMPLE AND SAMPLE ###
#############################################

rm(list = ls())

library("readxl")
library("openxlsx")

# load data coming from Rarefaction Script
zoo <- read_excel("Output/05_Zooplankton/Rarefaction_Zooplankton_FamilyCounts.xlsx")

names(zoo)
colnames(zoo)[1] <- "age" # rename age column
zoo <- zoo[, -2] # removes sampleEff column

zoo$age <- factor(zoo$age, levels = unique(zoo$age)) # make age a factor
is.factor(zoo$age) # check
length(levels(zoo$age)) # hopw many age levels (= samples)

# get subset of all resampling repeats per sample age (100 repeats for rarefaction)
a_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[1])
b_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[2])
c_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[3])
d_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[4])
e_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[5])
f_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[6])
g_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[7])
h_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[8])
i_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[9])
j_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[10])
k_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[11])
l_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[12])
m_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[13])
n_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[14])
o_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[15])
p_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[16])
q_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[17])
r_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[18])
s_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[19])
t_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[20])
u_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[21])
v_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[22])
w_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[23])
x_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[24])
y_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[25])
z_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[26])
aa_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[27])
ab_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[28])
ac_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[29])
ad_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[30])
ae_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[31])
af_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[32])
ag_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[33])
ah_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[34])
ai_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[35])
aj_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[36])
ak_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[37])
al_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[38])
am_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[39])
an_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[40])
ao_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[41])
ap_zoo <- subset(zoo, subset = age %in% levels(zoo$age)[42])

# calculate mean for each sampling age
a <- colMeans(a_zoo[, 2:ncol(a_zoo)])
b <- colMeans(b_zoo[, 2:ncol(b_zoo)])
c <- colMeans(c_zoo[, 2:ncol(c_zoo)])
d <- colMeans(d_zoo[, 2:ncol(d_zoo)])
e <- colMeans(e_zoo[, 2:ncol(e_zoo)])
f <- colMeans(f_zoo[, 2:ncol(f_zoo)])
g <- colMeans(g_zoo[, 2:ncol(g_zoo)])
h <- colMeans(h_zoo[, 2:ncol(h_zoo)])
i <- colMeans(i_zoo[, 2:ncol(i_zoo)])
j <- colMeans(j_zoo[, 2:ncol(j_zoo)])
k <- colMeans(k_zoo[, 2:ncol(k_zoo)])
l <- colMeans(l_zoo[, 2:ncol(l_zoo)])
m <- colMeans(m_zoo[, 2:ncol(m_zoo)])
n <- colMeans(n_zoo[, 2:ncol(n_zoo)])
o <- colMeans(o_zoo[, 2:ncol(o_zoo)])
p <- colMeans(p_zoo[, 2:ncol(p_zoo)])
q <- colMeans(q_zoo[, 2:ncol(q_zoo)])
r <- colMeans(r_zoo[, 2:ncol(r_zoo)])
s <- colMeans(s_zoo[, 2:ncol(s_zoo)])
t <- colMeans(t_zoo[, 2:ncol(t_zoo)])
u <- colMeans(u_zoo[, 2:ncol(u_zoo)])
v <- colMeans(v_zoo[, 2:ncol(v_zoo)])
w <- colMeans(w_zoo[, 2:ncol(w_zoo)])
x <- colMeans(x_zoo[, 2:ncol(x_zoo)])
y <- colMeans(y_zoo[, 2:ncol(y_zoo)])
z <- colMeans(z_zoo[, 2:ncol(z_zoo)])
aa <- colMeans(aa_zoo[, 2:ncol(aa_zoo)])
ab <- colMeans(ab_zoo[, 2:ncol(ab_zoo)])
ac <- colMeans(ac_zoo[, 2:ncol(ac_zoo)])
ad <- colMeans(ad_zoo[, 2:ncol(ad_zoo)])
ae <- colMeans(ae_zoo[, 2:ncol(ae_zoo)])
af <- colMeans(af_zoo[, 2:ncol(af_zoo)])
ag <- colMeans(ag_zoo[, 2:ncol(ag_zoo)])
ah <- colMeans(ah_zoo[, 2:ncol(ah_zoo)])
ai <- colMeans(ai_zoo[, 2:ncol(ai_zoo)])
aj <- colMeans(aj_zoo[, 2:ncol(aj_zoo)])
ak <- colMeans(ak_zoo[, 2:ncol(ak_zoo)])
al <- colMeans(al_zoo[, 2:ncol(al_zoo)])
am <- colMeans(am_zoo[, 2:ncol(am_zoo)])
an <- colMeans(an_zoo[, 2:ncol(an_zoo)])
ao <- colMeans(ao_zoo[, 2:ncol(ao_zoo)])
ap <- colMeans(ap_zoo[, 2:ncol(ap_zoo)])

# combine means per age in one data frame
countfam_mean <- as.data.frame(rbind(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab, ac, ad, ae, af, ag, ah, ai, aj, ak, al, am, an, ao, ap))

family <- names(countfam_mean) # save vactor with families for later

age <- levels(zoo$age)
countfam_mean$age <- age
ncol(countfam_mean)
countfam_mean <- countfam_mean[, c(38, 1:37)]
rownames(countfam_mean) <- NULL

countfam_mean$age <- sub("x", "", countfam_mean$age) # removes the x infront of the age
countfam_mean <- as.data.frame(t(countfam_mean)) # transposing so age is in columns
age <- countfam_mean[1, ] # vector with ages
colnames(countfam_mean) <- age # ages to column names
countfam_mean <- countfam_mean[-1, ] # remove row with ages

countfam_mean <- as.data.frame(sapply(countfam_mean, as.numeric)) # counts numeric

countfam_mean$family <- family # add column with family
ncol(countfam_mean)
countfam_mean <- countfam_mean[, c(43, 1:42)] # make family the first column
rownames(countfam_mean) <- NULL # remove row names

write.xlsx(countfam_mean, file = "Output/05_Zooplankton/Rarefaction_Zooplankton_FamilyCounts_means.xlsx")


#############################################
###  (C) ZOOPLANKTON RELATIVE ABUNDANCE   ###
#############################################

rm(list = ls())

library("readxl")
library("openxlsx")
library("dplyr")

zoo <- read.xlsx("Output/05_Zooplankton/Rarefaction_Zooplankton_FamilyCounts_means.xlsx")

ncol(zoo)
names(zoo)
zoo_onlycounts <- zoo[, -1] # remove columns with taxonomic information, so only counts are left

# create data frame with relative counts per column
zoo_rel <- data.frame(c(1:nrow(zoo_onlycounts)))
for( i in c(1:ncol(zoo_onlycounts))) {
  zoo_rel[, i] <- data.frame(zoo_onlycounts[, i]) %>%
    transmute(rel_ab = (zoo_onlycounts[, i]/sum((zoo_onlycounts[, i])))) # calculation for relative abundance for each cell of the data frame
}
zoo_rel <- zoo_rel*100 # for making relative abundance in % (somehow it makes the column names weird if yoh do it right in the for-loop)

# add taxonomic information again
zoo_rel <- cbind(zoo_rel, zoo$family) # add family column and class column
names(zoo_rel)
ncol(zoo_rel)
colnames(zoo_rel)[43] <- "family" # rename family column correctly
zoo_rel <- zoo_rel[, c(43, 1:42)]
names(zoo_rel) <- names(zoo)

# import taxonomic information (created during subsetting)
group_info <- read.xlsx("Group_infos_zoo.xlsx")

zoo_rel <- merge(zoo_rel, group_info, by = "family")

ncol(zoo_rel)
names(zoo_rel)
zoo_rel <- zoo_rel[, c(1, 44, 2:43)] # rearrange so family names are first column

write.xlsx(zoo_rel, file = "Output/05_Zooplankton/Zooplankton_rel_abund_all.xlsx")

#MAIN FAMILIES
# the zooplankton dataset was screened manually for marine families,
# excluding parasites and freshwater families
# the main zooplankton families were defined:
zoo_rel_main <- subset(zoo_rel, subset = family %in% c("Apusomonadidae",
                                                       "Calanidae",
                                                       "Metridinidae",
                                                       "Salpingoecidae",
                                                       "Temoridae",
                                                       "Parameciidae",
                                                       "Sphaerozoidae",
                                                       "Vahlkampfiidae"))


write.xlsx(zoo_rel_main, file = "Output/05_Zooplankton/Zooplankton_rel_abund_main.xlsx")

#############################################
###      (D) ZOOPLANKTON STRATIGRAM       ###
#############################################

rm(list = ls())

library("tidypaleo")
library("readxl")
library("openxlsx")
library("tidyverse")
library("ggplot2")
library("ggh4x")

theme_set(theme_paleo(8))

zoo_rel <- read.xlsx("Output/05_Zooplankton/Zooplankton_rel_abund_all.xlsx")

zoo_rel_main <- read.xlsx("Output/05_Zooplankton/Zooplankton_rel_abund_main.xlsx")

# summarize remaining rare zooplankton families by groups
zoo_rel$family # all families
zoo_rel_rare <- zoo_rel[c(3:19, 21:22, 24:29, 32, 34:36),] # exclude families of main zooplankton
zoo_rel_rare$group <- as.factor(zoo_rel_rare$group) # group as a factor
levels(zoo_rel_rare$group)
rare_crust <- subset(zoo_rel_rare, subset = group %in% "crustaceous zooplankton")
rare_gel <- subset(zoo_rel_rare, subset = group %in% "gelatinous zooplankton")
rare_pro <- subset(zoo_rel_rare, subset = group %in% "heterotrophic protists")

# summarize counts of groups per age
sum_rare_crust <- as.data.frame(colSums(rare_crust[, 3:ncol(rare_crust)]))
names(sum_rare_crust) <- "rare crustaceous zooplankton"

sum_rare_gel <- as.data.frame(colSums(rare_gel[, 3:ncol(rare_gel)]))
names(sum_rare_gel) <- "rare gelantinous zooplankton"

sum_rare_pro <- as.data.frame(colSums(rare_pro[, 3:ncol(rare_pro)]))
names(sum_rare_pro) <- "rare heterotrophic protists"

rare_zooplankton <- cbind(sum_rare_crust,
                          sum_rare_pro,
                          sum_rare_gel)

# make it fit to the main zooplankton table
rare_zooplankton <- as.data.frame(t(rare_zooplankton)) # transpose
rare_fam <- rownames(rare_zooplankton)
groups <- c("rare", "rare", "rare")
rare_zooplankton$family <- rare_fam
rare_zooplankton$group <- groups
names(rare_zooplankton)
ncol(rare_zooplankton)
rare_zooplankton <- rare_zooplankton[, c(43, 44, 1:42)] # rearange
rownames(rare_zooplankton) <- NULL

# get all in one data frame with main zooplankton
zoo_final <- rbind(zoo_rel_main, rare_zooplankton)

# change table from wide onto long format (needed for plotting)
zoo_long <- pivot_longer(zoo_final,
                         cols = !c(family, group),
                         names_to = c("age"),
                         values_to = "rel_abund")

# make age numeric
zoo_long$age <- as.numeric(zoo_long$age)

# preparation of data
zoo_long$age <- zoo_long$age/1000 # make age in ka BP
names(zoo_long)

zoo_long$group <- factor(zoo_long$group, levels = unique(zoo_long$group)) # make group a factor with order it has right now

# add color information by subsetting groups, add color column and combine to one data frame again
levels(zoo_long$group)
zoo_long_pro <- subset(zoo_long, subset = group %in% "heterotrophic protists")
zoo_long_pro$co <- rep("coral")
zoo_long_crust <- subset(zoo_long, subset = group %in% "crustaceous zooplankton")
zoo_long_crust$co <- rep("darkgoldenrod")
zoo_long_rare <- subset(zoo_long, subset = group %in% "rare")
zoo_long_rare$co <- rep("black")

zoo_long_final <- rbind(zoo_long_pro,
                        zoo_long_crust,
                        zoo_long_rare)

# order
names(zoo_long_final)

zoo_long_final$rel_abund <- as.numeric(zoo_long_final$rel_abund) # make rel_abund numeric
zoo_long_final$age <- as.numeric(zoo_long_final$age) # male age numeric
zoo_long_final$family <- factor(zoo_long_final$family, levels = unique(zoo_long_final$family)) # make family a factor with order it has right now
zoo_long_final$group <- factor(zoo_long_final$group, levels = unique(zoo_long_final$group)) # make group a factor with order it has right now

# zoo plot:
zoo_plot <- ggplot(zoo_long_final, aes(x = age, y = rel_abund, fill = group)) +
  geom_area(aes(fill = co)) +
  geom_segment(aes(xend = age, yend = 0), lwd = 0.5, color = "black") +
  expand_limits(x = 0, y = 0) + # axis start with zero
  scale_y_continuous(n.breaks = 3, # numbers mark y axis
                     expand = c(0,0),
                     name = "Relative abundance [%] \nbased on total reads of zooplankton families") + 
  facet_abundanceh(vars(family), rotate_facet_labels = 70, space = "fixed", scales = "free_x") + # groups plots by family, makes every plot the same size (different scale)
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

zoo_plot
ggsave(zoo_plot, file = "Output/03_Stratigram/Zooplankton_area.pdf", width = 20, height = 20, units = "cm")


########################################################################################################################################################
####################### 6 PCA/RDA  #####################################################################################################################
########################################################################################################################################################

dir.create("Output/06_PCA_RDA")

rm(list = ls())

library("openxlsx")
library("vegan")
library("ggplot2")
library("tidypaleo")
library("openxlsx")
library("olsrr")
library("ggrepel")

theme_set(theme_paleo(8))

#############################################
###       (A) ZOOPLANKTON DATA            ###
#############################################

# read in zooplankton data
zoo_data <- read.xlsx("Output/05_Zooplankton/Zooplankton_rel_abund_main.xlsx")
names(zoo_data)

names(zoo_data)
zoo_data <- zoo_data[, -2] # remove group column
families <- zoo_data[, 1] # get family names
rownames(zoo_data) <- families # family as row name
zoo_data <- zoo_data[, -1] # remove family column
age <- as.numeric(names(zoo_data))/1000 # get sample ages

zoo_envi <- as.data.frame(t(zoo_data)) # transpose

# correlation zooplankton families
zoo_cor <- data.frame(cor(zoo_envi, method = "spearman"))
plot(zoo_cor)

# VIFs
# check collinearity (on untransformed data):
zoo_envi$age <- age
ncol(zoo_envi)
zoo_envi <- zoo_envi[, c(9, 1:8)]
model <- lm(zoo_envi)
ols_vif_tol(model)
# VIF should be <10, but better is <4


##############################
### WITHOUT SHPHAEROZOIDAE ###
##############################
names(zoo_envi)
zoo_envi2 <- zoo_envi[, -7]

# VIFs
model2 <- lm(zoo_envi2)
ols_vif_tol(model2)
# VIF should be <10, but better is <4

#############################
### WITHOUT VALKAMPFIIDAE ###
#############################

names(zoo_envi)
zoo_envi3 <- zoo_envi[, -9]

# VIFs
model3 <- lm(zoo_envi3)
ols_vif_tol(model3)
# VIF should be <10, but better is <4

############################
### WITHOUT PARAMECIIDAE ###
############################
names(zoo_envi)
zoo_envi4 <- zoo_envi[, -5]

# VIFs
model4 <- lm(zoo_envi4)
ols_vif_tol(model4)
# VIF should be <10, but better is <4


### get final dataset based on VIFs and filtering for main families:
zoo_data <- read.xlsx("Output/05_Zooplankton/Zooplankton_rel_abund_main.xlsx")
names(zoo_data)
zoo_red <- zoo_data[, c(1, 3:ncol(zoo_data))] # get data without group information
zoo_red <- as.data.frame(t(zoo_red)) # transpose
names(zoo_red) <- zoo_red[1, ] # family as column names
zoo_red <- zoo_red[-1, ] # delete family row
names(zoo_red)
zoo_red <- zoo_red[, c("Apusomonadidae",
                       "Calanidae",
                       "Salpingoecidae")] # get data frame with the final 3 zooplankton families

zoo_red <- as.data.frame(sapply(zoo_red, as.numeric)) # counts numeric

# summarize heterotrophic protists with Apusomonadidae + Salpingoecidae
zoo_red$het_prot <- zoo_red$Apusomonadidae + zoo_red$Salpingoecidae
names(zoo_red)
zoo_final <- zoo_red[, c(2, 4)]
zoo_final$age <- age
names(zoo_final)
zoo_final <- zoo_final[, c(3, 1:2)]


##################################################
### (B) HELLINGER TRANSFORMED ZOOPLANKTON DATA ###
##################################################

rm(list = ls())

# read in zooplankton data
zoo_data <- read.xlsx("Output/05_Zooplankton/Zooplankton_rel_abund_all.xlsx")
names(zoo_data)

zoo_data <- zoo_data[, -2] # remove group column
families <- zoo_data[, 1] # get family names

# add Apusomonadidae and Salpingoecidae (togther they will represent the heterotrophic zooplankton)
het_zoo <- subset(zoo_data, subset = family %in% c("Apusomonadidae", # subset with the two families
                                                   "Salpingoecidae"))
het_zoo_sum <- colSums(het_zoo[, 2:ncol(het_zoo)]) # sum of relative abundances
het_zoo_sum <- c("het_pro", het_zoo_sum)

zoo_data <- rbind(zoo_data, het_zoo_sum) # add sum of heterotrophic zooplankton to total zooplankton dataset
zoo_data$family

# remove Apusomonadidae and Salpingoecidae from zooplankton dataset (they are there with the sum now)
row_Apu <- which(zoo_data$family == "Apusomonadidae")
row_Salp <- which(zoo_data$family == "Salpingoecidae")

zoo_data <- zoo_data[-c(row_Apu, row_Salp), ]

# prep. for Hellinger transformation
family <- zoo_data$family
zoo_data_temp <- zoo_data[, 2:ncol(zoo_data)]
zoo_data_temp <- as.data.frame(sapply(zoo_data_temp, as.numeric)) # counts numeric

# Hellinger transformation = square root of relative abundance data
# here: count data already transformed to relative abundances,
# therefore only sqaure root transformation on seq_data needed for Hellinger transformed data
zoo_data_hell <- sqrt(zoo_data_temp)

# make data frame look nice
zoo_data_hell$family <- family
names(zoo_data_hell)
ncol(zoo_data_hell)
zoo_data_hell <- zoo_data_hell[, c(43, 1:42)]

# filter for the 3 zooplankton families
zoo_filtered <- subset(zoo_data_hell, subset = family %in% c("Calanidae",
                                                             "het_pro"))
names(zoo_filtered)
families_filtered <- zoo_filtered[, 1] # get family names
rownames(zoo_filtered) <- families_filtered # family as row name
zoo_filtered <- zoo_filtered[, -1] # remove family column
age <- as.numeric(names(zoo_filtered))/1000 # get sample ages

zoo_envi <- as.data.frame(t(zoo_filtered)) # transpose
zoo_envi$age <- age
names(zoo_envi)
zoo_envi <- zoo_envi[, c(3, 1:2)]

#############################################
###           (C) NGRIP DATA              ###
#############################################

# read in environmental data
ngrip_data <-read.xlsx("Output/04_Environment/Shotgun_NGRIP.xlsx") 
names(ngrip_data)
names(ngrip_data) <- c("age", "NGRIP")

# combine data frames of available environmental data
all_envi_data <- merge(ngrip_data, zoo_envi, by = "age") # ngrip_data with zoo families
envi_cor <- data.frame(cor(all_envi_data, method = "spearman"))
plot(all_envi_data)

# VIFs
model_all <- lm(all_envi_data)
ols_vif_tol(model_all)
# VIF should be <10, but better is <4

# export
write.xlsx(all_envi_data, file = "Output/06_PCA_RDA/Envi_for_RDA.xlsx")

#############################################
###         (D) SEQUENCING DATA           ###
#############################################

# import relative abundance of the main families (same as in stratigram)
# and create vector with families for PCA
seq_data_fams <-read.xlsx("Output/03_Stratigram/Phytoplankton_main_rel_abund.xlsx")
names(seq_data_fams)
seq_fams <- as.vector(seq_data_fams$family)

# import full relative abundance dataset for Hellinger transformation
seq_data <- read.xlsx("Output/03_Stratigram/Phytoplankton_rel_abund.xlsx")

seq_data$family
nrow(seq_data)

names(seq_data)
family <- seq_data$family
age <- as.numeric(names(seq_data[3:42]))/1000
seq_data_temp <- as.data.frame(sapply(seq_data[, 3:ncol(seq_data)], as.numeric)) # only numeric
rownames(seq_data_temp) <- family

# Hellinger transformation = square root of relative abundance data
# here: count data already transformed to relative abundances,
# therefore only sqaure root transformation on seq_data needed for Hellinger transformed data
seq_data_hell <- sqrt(seq_data_temp)

# filter for the main families
seq_data_hell <- subset(seq_data_hell, subset = family %in% seq_fams)

# transpose
seq_data_hell <- as.data.frame(t(seq_data_hell))

# Check the Axis lengths of DCA1, use CCA if >4
decorana(seq_data_hell)
# DCA1 axis length = 1.30407 --> OK

ncol(seq_data_hell)


###############################################
### (E) FORWARD MODEL WITH ZOOPLANKTON ONLY ###
###############################################

# dataset without age
names(all_envi_data)
envi_final <- all_envi_data[, -1]

# RDA
rda <- rda(seq_data_hell ~., data = envi_final)
plot(rda)
adjR2.all <- RsquareAdj (rda)$adj.r.squared
adjR2.all


data_envi_rda0 <- rda(seq_data_hell ~ 1, data = envi_final) # "empty" model with intercept only
rda <- rda(seq_data_hell ~., data = envi_final) # full model
sel.osR2 <- ordiR2step (data_envi_rda0, scope = formula (rda), R2scope = adjR2.all, direction = 'forward', permutations = 999)
sel.osR2
sel.osR2$anova


sel.osR2_adj <- sel.osR2
sel.osR2_adj$anova$`Pr(>F)` <- p.adjust (sel.osR2$anova$`Pr(>F)`, method = 'holm', n = ncol (envi_final))
sel.osR2_adj$anova


# check backward selection (gives the same result as forward selection)
sel.osR2_back <- ordiR2step (data_envi_rda0, scope = formula (rda), R2scope = adjR2.all, direction = 'backward', permutations = 999)
sel.osR2_back
sel.osR2_back$anova


sel.osR2_adj_back <- sel.osR2_back
sel.osR2_adj_back$anova$`Pr(>F)` <- p.adjust (sel.osR2_back$anova$`Pr(>F)`, method = 'holm', n = ncol (envi_final))
sel.osR2_adj_back$anova
# output: same result as with foreward selection


### SUMMARY RDA

# https://mb3is.megx.net/gustame/constrained-analyses/rda for interpretation

# final RDA
plot(rda)

adjR2.final <- RsquareAdj (rda)$adj.r.squared
adjR2.final

summary(rda)

screeplot(rda)

#################################################
### (G) PROJECTING RDA IN UNCONSTRAINED SPACE ###
#################################################

# run a  pca
pca <- rda(seq_data_hell, scale=TRUE)
summary <- summary(pca)
plot(pca)

# incorporate environmental factors as restricting variables
fit_all <- envfit(pca ~ . , envi_final, perm = 999)
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

# PLOTTING
# extract age information
age_position_pca <- as.data.frame(summary$sites[, 1:2])
age_position_pca$age <- as.numeric(rownames(age_position_pca))/1000 # make column with age

age_info <- read.xlsx("Age_Period_Group.xlsx") # import age_period table
age_position_pca <- merge(age_position_pca, age_info, by = "age")
age_position_pca$pch <- as.factor(age_position_pca$pch)


# extract family information
family_position_pca <- as.data.frame(summary$species[, 1:2])
# type I scaling for families:
family_scaling <- plot(pca, scaling = 1)
family_position_scaleII <- as.data.frame(family_scaling$species)

# add size information in colour of the data point
family_info <- read.xlsx("Group_infos_phyto.xlsx")
family_temp <- family_position_scaleII
family_temp$family <- rownames(family_temp)
family_temp <- merge(family_temp, family_info, by = "family")
names(family_temp)
#family_temp <- family_temp[, -4]

rownames(family_temp) <- family_temp$family

family_position_final <- family_temp
names(family_position_final)
family_position_final <- family_position_final[, -1] # remove family column

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

#PLOT
PCA_envi <- ggplot(family_position_final, aes(x = PC1, y = PC2)) +
  geom_point(data = family_position_final, aes(x = PC1, y = PC2), size = 1, stroke = 0.8)+#,  color = family_position_final$col) +
  
  geom_segment(aes(x = 0, y = 0, xend = NGRIP_PC1, yend = NGRIP_PC2), colour = "gray14", # NGRIP
               arrow = arrow(length = unit(0.35, "cm")), size = 0.5) +
  annotate("text", x = NGRIP_PC1+0.03, y = NGRIP_PC2+0.03, label = expression(paste(delta^"18", "O NGRIP")), colour = "gray14", angle = 0, size = 5, fontface = "bold") +
  
  geom_segment(aes(x = 0, y = 0, xend = Cala_PC1, yend = Cala_PC2), colour = "gray14", # Calanidae
               arrow = arrow(length = unit(0.35, "cm")),  size = 0.5) +
  annotate("text", x = Cala_PC1+0.2, y = Cala_PC2+0.03, label = "Calanidae",colour = "gray14", angle = 0, size = 5, fontface = "bold") +
  
  geom_segment(aes(x = 0, y = 0, xend = het_PC1, yend = het_PC2), colour = "gray14", # Salpingoecidae
               arrow = arrow(length = unit(0.35, "cm")),  size = 0.5) +
  annotate("text", x = het_PC1-0.6, y = het_PC2+0.03, label = "heterotrophic protists",colour = "gray14", angle = 0, size = 5, fontface = "bold") +
  
  geom_point(data = age_position_pca, aes(x = PC1, y = PC2, shape = group),size = 4) + # age
  geom_text_repel(data = family_position_final,
                  aes(label=rownames(family_position_final),
                      point.size = 10,
                      colour = col,
                      segment.color = "black"),
                  min.segment.length = 0.1, seed = 45, box.padding = 0.1)+
  scale_color_manual(values = c("cyan3", "green3", "olivedrab4", "coral2")) + # colors for bacteria, protists, chlorophytes, red algae (in this order!)
  scale_shape_manual(values=c(18, 1, 13)) + 
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_x_continuous(name = PC1_explain) +
  expand_limits(x = c(0.8, 0.6)) +
  scale_y_continuous(name = PC2_explain) +
  geom_vline(xintercept = 0, size = 0.4, linetype = "dashed", colour = adjustcolor("gray26", alpha.f = 0.7)) +
  geom_hline(yintercept = 0, size = 0.4, linetype = "dashed", colour = adjustcolor("gray26", alpha.f = 0.7))
PCA_envi

ggsave(PCA_envi, file = "Output/06_PCA_RDA/PCA_NGRIP_Zooplankton_StratFamilies_FamScalingI_envi_notrans.pdf", width = 21, height = 16, units = "cm")


#########################################################
### (H) VARIATION PARTITIONING ENVIONMENTAL VARIABLES ###
#########################################################

names(envi_final)

# all variables included
rda_all <- rda(seq_data_hell ~., data = envi_final)
all_variance <- RsquareAdj (rda_all)$adj.r.squared
all_variance

names(envi_final)
varp <- varpart (seq_data_hell, ~ NGRIP, ~ Calanidae, ~het_pro, data = envi_final)
varp

plot (varp, digits = 2, Xnames = c("NGRIP", "Calanidae", "heterotrophic protists"), bg = c('navy', 'tomato', 'yellow'))


#################################################################################
################## TEST FOR SIGNIFICANCE ########################################
#################################################################################

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
rda_Cala <- rda (seq_data_hell ~ Calanidae, data = envi_final) # Calanidae
rda_hetprot <- rda (seq_data_hell ~ het_pro, data = envi_final) # heterotrophic protists
rda_ngrip <- rda (seq_data_hell ~ NGRIP, data = envi_final) # NGRIP
rda_Cala_cond <- rda (seq_data_hell ~ Calanidae + Condition (het_pro) + Condition(NGRIP), data = envi_final)
rda_hetprot_cond <- rda (seq_data_hell ~ het_pro + Condition (Calanidae) + Condition(NGRIP), data = envi_final)
rda_ngrip_cond <- rda (seq_data_hell ~ NGRIP + Condition (Calanidae) + Condition(het_pro), data = envi_final)

# ANOVA
set.seed(123); anova(rda_all) # test significance of conditional effect
set.seed(123); anova(rda_Cala) # test simple (marginal) effect of Calanidae
set.seed(123); anova(rda_hetprot) # test simple (marginal) effect of heterotrophic protists
set.seed(123); anova(rda_ngrip) # test simple (marginal) effect of NGRIP

set.seed(123); anova(rda_hetprot_cond)
set.seed(123); anova(rda_Cala_cond)
set.seed(123); anova(rda_ngrip_cond)
