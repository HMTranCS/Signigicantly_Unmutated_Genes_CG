library(tidyverse)
library(readxl)
library(ggplot2)
library(likert)

setwd("~/Documents/School/BaileyResearch/Cancer_Genomics")
par(mar=c(1,1,1,1))

########pan-cancer cnv gene data#########
dat1 = read_tsv("Applications/CNVTree/TCGA_Data/7d64377f-2cea-4ee3-917f-8fcfbcd999e7")

# remove first three columns to only leave numeric data
dat1_data = select(dat1,!c("Gene Symbol", "Locus ID", "Cytoband"))

# create new column that records the total number of values that are either -1 or 0
totals1 = dat1_data %>%
  mutate(totals = rowSums(. == -1 | . == 0))
dat1 = bind_cols(select(dat1,`Gene Symbol`),totals1)

# plot and print
select(dat1, `Gene Symbol`, `totals`) -> totals_by_gene
plot(density(totals_by_gene$totals))
arrange(totals_by_gene,totals) %>% print(n=300)

raw_random = dat1[c(sample(25128,1000,replace = FALSE)),]

raw_random_data = select(raw_random,!c("Gene Symbol"))

# create new column that records the total number of values that are either -1 or 0
raw_random %>%
  mutate(totals = rowSums(. <= -1)) -> del_random_data
raw_random %>%
  mutate(totals = rowSums(. > 0)) -> amp_random_data
raw_random %>%
  mutate(totals = rowSums(. == 0)) -> neutral_random_data

random_counts = tibble(deletion=del_random_data$totals,amplification=amp_random_data$totals,neutral=neutral_random_data$totals)
average_random_counts = summarise(random_counts,avg_del=mean(random_counts$deletion),avg_amp=mean(random_counts$amplification),avg_neu=mean(random_counts$neutral))

########for oncogenes###########
dat2 = read_tsv("https://ongene.bioinfo-minzhao.org/ongene_human.txt")
oncogenes = pull(dat2,OncogeneName)

totals_by_gene = mutate(totals_by_gene,oncogene=`Gene Symbol`%in%oncogenes)

non_oncogene = sample(which(totals_by_gene$oncogene == FALSE),803,replace = FALSE)
is_oncogene = which(totals_by_gene$oncogene == TRUE)

sample_total_by_genes = totals_by_gene[c(non_oncogene,is_oncogene),]

ggplot() + 
  geom_density(data=sample_total_by_genes, aes(x=totals,fill=oncogene)) +
  scale_fill_manual(values=alpha(c('red','blue'),.3))

# test of significance using 
wilcox.test(sample_total_by_genes$totals~sample_total_by_genes$oncogene)

# test for normal distribution
shapiro.test(sample_total_by_genes$totals)

raw_oncogene = filter(dat1,`Gene Symbol`%in%oncogenes)

raw_oncogene_data = select(raw_oncogene,!c("Gene Symbol"))

raw_oncogene_data %>%
  mutate(totals = rowSums(. <= -1)) -> del_oncogene_data
raw_oncogene_data %>%
  mutate(totals = rowSums(. > 0)) -> amp_oncogene_data
raw_oncogene_data %>%
  mutate(totals = rowSums(. == 0)) -> neutral_oncogene_data

oncogene_counts = tibble(deletion=del_oncogene_data$totals,amplification=amp_oncogene_data$totals,neutral=neutral_oncogene_data$totals)
average_oncogene_counts = summarise(oncogene_counts,avg_del=mean(oncogene_counts$deletion),avg_amp=mean(oncogene_counts$amplification),avg_neu=mean(oncogene_counts$neutral))

#############################
#######for tumor suppressors#########
dat3 = read_tsv("https://bioinfo.uth.edu/TSGene/Human_TSGs.txt?csrt=4608055084028856876")
tsg = pull(dat3,GeneSymbol)
select(dat1, `Gene Symbol`, `totals`) -> totals_by_gene

totals_by_gene = mutate(totals_by_gene,tsg=`Gene Symbol`%in%tsg)

non_tsg = sample(which(totals_by_gene$tsg == FALSE),1217,replace = FALSE)
is_tsg = which(totals_by_gene$tsg == TRUE)

sample_total_by_genes = totals_by_gene[c(non_tsg,is_tsg),]

ggplot() + 
  geom_density(data=sample_total_by_genes, aes(x=totals,fill=tsg)) +
  scale_fill_manual(values=alpha(c('red','blue'),.3))

raw_tsg = filter(dat1,`Gene Symbol`%in%tsg)

# remove first three columns to only leave numeric data
raw_tsg_data = select(raw_tsg,!c("Gene Symbol"))

# create new column that records the total number of values that are either -1 or 0
raw_tsg_data %>%
  mutate(totals = rowSums(. <= -1)) -> del_tsg_data
raw_tsg_data %>%
  mutate(totals = rowSums(. > 0)) -> amp_tsg_data
raw_tsg_data %>%
  mutate(totals = rowSums(. == 0)) -> neutral_tsg_data

tsg_counts = tibble(deletion=del_tsg_data$totals,amplification=amp_tsg_data$totals,neutral=neutral_tsg_data$totals)
average_tsg_counts = summarise(tsg_counts,avg_del=mean(tsg_counts$deletion),avg_amp=mean(tsg_counts$amplification),avg_neu=mean(tsg_counts$neutral))

# test of significance using 
wilcox.test(sample_total_by_genes$totals~sample_total_by_genes$tsg)

# test for normal distribution
shapiro.test(sample_total_by_genes$totals)
####################
##########for significantly unmutated genes###################
sun_mg = read_tsv("/Users/hannahtran/Downloads/sunMG_top625.txt",col_names = FALSE)
select(dat1, `Gene Symbol`, `totals`) -> totals_by_gene

sun_mg_totals = mutate(totals_by_gene,sun_mg=`Gene Symbol`%in%sun_mg$X1)

non_sun_mg = sample(which(sun_mg_totals$sun_mg == FALSE),625,replace = FALSE)
is_sun_mg = which(sun_mg_totals$sun_mg == TRUE)

sample_sun_mg = sun_mg_totals[c(non_sun_mg,is_sun_mg),]

ggplot() + 
  geom_density(data=sample_sun_mg, aes(x=totals,fill=sun_mg)) +
  scale_fill_manual(values=alpha(c('red','blue'),.3))

sun_mg_restricted = filter(sun_mg_totals,totals > 8250 & totals < 9000 & sun_mg == TRUE)

raw_sun_mg = filter(dat1,`Gene Symbol`%in%sun_mg_restricted$`Gene Symbol`)

# remove first three columns to only leave numeric data
raw_sun_mg_data = select(raw_sun_mg,!c("Gene Symbol"))

# create new column that records the total number of values that are either -1 or 0
raw_sun_mg_data %>%
  mutate(totals = rowSums(. <= -1)) -> del_sun_mg_data
raw_sun_mg_data %>%
  mutate(totals = rowSums(. > 0)) -> amp_sun_mg_data
raw_sun_mg_data %>%
  mutate(totals = rowSums(. == 0)) -> neutral_sun_mg_data

sun_mg_counts = tibble(deletion=del_sun_mg_data$totals,amplification=amp_sun_mg_data$totals,neutral=neutral_sun_mg_data$totals)
average_sun_mg_counts = summarise(sun_mg_counts,avg_del=mean(sun_mg_counts$deletion),avg_amp=mean(sun_mg_counts$amplification),avg_neu=mean(sun_mg_counts$neutral))
#################

likert_levels = c("DELETION","NEUTRAL","AMPLIFICATION")

random_data = mutate_all(raw_random_data,~ replace(.,(. == -1 |. == -2),"DELETION")) %>%
  mutate_all(~ replace(., (. == 1 | . == 2),"AMPLIFICATION")) %>%
  mutate_all(~ replace(.,(. == 0),"NEUTRAL")) %>%
  select(!c("totals"))
oncogene_data = mutate_all(raw_oncogene_data,~ replace(.,(. == -1 |. == -2),"DELETION")) %>%
  mutate_all(~ replace(.,(. == 1 | . == 2),"AMPLIFICATION")) %>%
  mutate_all(~ replace(.,(. == 0),"NEUTRAL")) %>%
  select(!c("totals"))
tsg_data = mutate_all(raw_tsg_data,~ replace(.,(. == -1 |. == -2),"DELETION")) %>%
  mutate_all(~ replace(.,(. == 1 | . == 2),"AMPLIFICATION")) %>%
  mutate_all(~ replace(.,(. == 0),"NEUTRAL"))%>%
  select(!c("totals"))
sun_mg_data = mutate_all(raw_sun_mg_data,~ replace(.,(. == -1 |. == -2),"DELETION")) %>%
  mutate_all(~ replace(.,(. == 1 | . == 2),"AMPLIFICATION")) %>%
  mutate_all(~ replace(.,(. == 0),"NEUTRAL"))%>%
  select(!c("totals"))

random_all=unlist(random_data, use.names = FALSE)
oncogene_all=unlist(oncogene_data, use.names = FALSE)
tsg_all=unlist(tsg_data, use.names = FALSE)
sun_mg_all=unlist(sun_mg_data, use.names = FALSE)

length(random_all) = length(tsg_all)
length(oncogene_all) = length(tsg_all)
length(sun_mg_all) = length(tsg_all)

dat <- cbind(Random=random_all,Oncogene=oncogene_all,TumorSuppressor=tsg_all,Unmutated625=sun_mg_all)
df = data.frame(Random = factor(dat[,"Random"],likert_levels),
                Oncogene = factor(dat[,"Oncogene"],likert_levels),
                TumorSuppressor = factor(dat[,"TumorSuppressor"],likert_levels),
                Sun_Mg = factor(dat[,"Unmutated625"],likert_levels))
p <- likert(df)
plot(p,group.order = c("Random","Oncogene","TumorSuppressor","Sun_Mg"))
