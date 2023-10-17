library(tidyverse)
library("httr")
library("jsonlite")
library(fingerprint)
library(rcdk)
library(rprojroot)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChemmineR",force=TRUE)
library(ChemmineR)
setwd("C:/Users/atran/Documents/R Scripts/SCOPE/SCOPE_Lipid_Analysis--master/analyses/data/C18_IT_BTLE_vs_nESI_offline")

get_smiles_database

mdb <- read.SDFset("structures.sdf")
mdb

tibble(LM_ID=rep("test",length(mdb)),SMILES=rep("test",length(mdb)),
       lipid_name=rep("test",length(mdb)))->get_smiles_database
for (i in (seq_along(mdb))){
  if ("LM_ID" %in% names(mdb@SDF[[i]]@datablock) & 
      "SMILES" %in% names(mdb@SDF[[i]]@datablock)&
      "SYNONYMS" %in% names(mdb@SDF[[i]]@datablock)&
      "NAME" %in% names(mdb@SDF[[i]]@datablock)){
    
  mdb@SDF[[i]]@datablock[["LM_ID"]] ->get_smiles_database[i,1]
  mdb@SDF[[i]]@datablock[["SMILES"]] ->get_smiles_database[i,2]
  mdb@SDF[[i]]@datablock[["NAME"]] ->get_smiles_database[i,3]
  
  }
else {
  NA->get_smiles_database[i,1]
  NA->get_smiles_database[i,2]
  NA->get_smiles_database[i,3]
  }
}
get_smiles_database%>%filter(str_detect(get_smiles_database$LM_ID,"LMSP")&!is.na(LM_ID))->get_smiles_SPdb
get_smiles_SPdb

sum(is.na(get_smiles_database[["SMILES"]]))
#making database with lipid shorthand and SMILES from lipidmaps database
#some lipidmaps database does not contain SMILES input, so they become NA in get_smiles_database
"GlcCer(d18:2(4Z)/16:0(2OH[R]))"->test
str_split(test, "\\(")[[1]][1]
str_extract(test, "(?<=\\().*(?=/)")
str_extract(test, "(?<=/).*(?=\\))")

str_extract_lipidsubclass <- function(y) {map_chr(y,
 function(x){
 str_split (x, "\\(")[[1]][1] 
 })
}
str_extract_chain1<- function(y) {map_chr(y,
 function(x){
 str_split(x, "\\(|\\/")[[1]][2] 
 })
}

str_extract_chain2<- function(y) {map_chr(y,
  function(x){((str_split(x, "\\/")[[1]][2])%>%str_remove("\\)")%>%str_split("\\("))[[1]][1]
 })
}
str_extract_chain2_all<- function(y) {map_chr(y,
 function(x){(str_split(x, "\\/")[[1]][2])%>%str_remove("\\)")
 })
}
str_extract_lipidsubclass(test)
str_extract_chain1(test)
str_extract_chain2(test)
#string regular expressions for getting lipid class and acyl chain info... 
#extracts out all stuff 1. before "(" 2. after the first "(" and before the "/" 3. between "/" and ")" (gets chains only, no OH or EZ info)
#chain 2 all is everything after the "/"
#x is the sub_string 
str_extract_chain1(testdb)
get_smiles_SPdb[[3]]->testdb

read_csv("c18_hilic_test.csv")->tibble1
tibble1[1]->tibble1
tibble(abbrev=tibble1[[1]],smi="2",matches=1)->tibble1
tibble1
str_replace(tibble1[[1]], "Hex", "Glc")->tibble1[1]
str_replace(tibble1[[1]],"_","/")->tibble1[1]
tibble1

for (i in seq_along(tibble1[[1]])) {
  if(str_extract_chain2_all(tibble1[[1]][i]) %>%str_detect("OH")) {
tibble1[[1]][i]%>% str_extract_lipidsubclass() ==
testdb%>%str_extract_lipidsubclass()%>% str_match(str_extract_lipidsubclass(tibble1[[1]][i]))&
      
tibble1[[1]][i]%>% str_extract_lipidsubclass()%>% str_count("[:alpha:]") == 
testdb%>% str_extract_lipidsubclass()%>% str_count("[:alpha:]")&
      
tibble1[[1]][i]%>% str_extract_chain1() == testdb%>% str_extract_chain1()  &
tibble1[[1]][i]%>% str_extract_chain2() == testdb%>% str_extract_chain2() %>% str_match(str_extract_chain2(tibble1[[1]][i])) &
tibble1[[1]][i]%>% str_extract_chain2_all()%>%str_detect("OH")==testdb%>% str_extract_chain2_all() %>% str_detect("OH") ->lipid_str_detect 
  }else{
  
tibble1[[1]][i]%>% str_extract_lipidsubclass() ==
testdb%>%str_extract_lipidsubclass()%>% str_match(str_extract_lipidsubclass(tibble1[[1]][i]))&
    
tibble1[[1]][i]%>% str_extract_lipidsubclass()%>% str_count("[:alpha:]") == 
testdb%>%str_extract_lipidsubclass()%>% str_count("[:alpha:]")&

tibble1[[1]][i]%>% str_extract_chain1() == testdb%>% str_extract_chain1() &
      
tibble1[[1]][i]%>% str_extract_chain2() == testdb%>% str_extract_chain2()->lipid_str_detect
  }
  #lipid_str_detect is a string detection vector of true/false
  #made different for hydroxylated acylchain2 (first if)
    if (TRUE %in% lipid_str_detect) {
  get_smiles_SPdb[[2]][min(which(lipid_str_detect))]->tibble1[[2]][i]
  sum(lipid_str_detect, na.rm = TRUE) ->tibble1 [[3]][i]
  } else {
  NA->tibble1[[2]][i]
    NA->tibble1[[3]][i]
  }
}
tibble1
head(tibble1)
write.csv(tibble1,"test.csv")
read.csv("tibble1.csv")->tibble1
#matches for subclass: needs to have same number of letters and the subclass name must be in match somewhere
#e.g. "GlcCer" has "Cer" in there somewhere but does not have same number of letters 
#For chain2, exact match must be found (exactly as extracted from dataset) since there are no d or t to notate hydroxylation
#if there are multiple matches, picks lowest one
#counts number of matches


Case1_v_Ctrl_Negative <- read.csv("c18_hilic_test.csv", stringsAsFactors = F)
Case1_v_Ctrl_Positive <- read.csv("c18_hilic_test.csv", stringsAsFactors = F)
All_lipids <- rbind(Case1_v_Ctrl_Negative, Case1_v_Ctrl_Positive)

df <- as.data.frame(tibble1)
tibble1%>%filter(!is.na(tibble1[2]))->df
df
head(df)

sum(df$smi == "")

df.c <- df[df$smi != "", ]
df.na <- df[df$smi == "", ]

head(df.c)

fing <- parse.smiles(df.c$smi)
fing <- lapply(fing, get.fingerprint, type = "circular")
fp.sim <- fp.sim.matrix(fing)
row.names(fp.sim) <- df.c$abbrev
fp.dist <- as.dist(1 - fp.sim)
cls <- hclust(fp.dist)
plot(cls)

# Dendrogram of all identified lipids, structurally related using an ECFP_6 fingerprint, Euclidian distance and average linkage method.
save.image("fp_all_lipid.rdata")
library(httr)
library(jsonlite)
library(fingerprint)
library(rcdk)
library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)
library(limma)
library(digest)
library(factoextra)
library(NbClust)
library(reshape2)
library(scales)
library(tidyverse)


#Loading lipid clustering
tupgma <- upgma(fp.dist, method = "average")



## Building Log2FC heatmaps of statistically significant lipids
data_fn <- c("c18_hilic_test_01", "c18_hilic_test_02")
fn <- data_fn[1]
sig.df <- data.frame(matrix(nrow = 37, ncol = 2*length(data_fn)))
for( i in 1 : length(data_fn)) {
  fn <- data_fn[i]
  df <- read.csv(paste0(fn,".csv"))
  sig.df[,i] <- (df$Log2FC)
  sig.df[,i+2]<- (df$Case_v_Ctrl_Flag) 
  colnames(sig.df) <- c("Case1_v_Ctrl_Log2FC", "Case2_v_Ctrl_Log2FC","Case1_v_Ctrl_Flag", "Case2_v_Ctrl_Flag")
}
sig.df$Name <- df$Name
sig.df <- sig.df[c(5, 1:4)]
lipids <- sig.df$Name
sig.df <- unique(sig.df)
rownames(sig.df) <- sig.df$Name
sig.df$Name <- NULL
sig.df$Case1_v_Ctrl_Log2FC <- sig.df$Case1_v_Ctrl_Log2FC*abs(sig.df$Case1_v_Ctrl_Flag)
sig.df$Case2_v_Ctrl_Log2FC <- sig.df$Case2_v_Ctrl_Log2FC*abs(sig.df$Case2_v_Ctrl_Flag)
sig.df[sig.df == 0] <- NA
breaks <- seq(from=min(range(-5)), to=max(range(5)), length.out=10000)
midpoint <- which.min(abs(breaks - 0))
rampCol1 <- colorRampPalette(c("#a50026", "#a50026", "#a50026"))(midpoint)
rampCol2 <- colorRampPalette(c("#fee090","#e0f3f8","#abd9e9", "#74add1", "#4575b4"))(10000-(midpoint+1))
rampCols <- c(rampCol1,rampCol2)


##Lipid Dendogram -- Identifications annotated

p<- ggtree(tupgma, layout="circular", size=1, branch.length="none")
p <- p + geom_text(aes(label=label, angle=angle, fontface="bold"), hjust=-0.15, size=1.65)
p <- open_tree(p, angle = 3)
p <- gheatmap(p, sig.df [1], offset=20, width=.25, font.size=1.5, colnames = F)+
  scale_fill_gradient2(low = rampCol1, high = rampCol2, mid = "white", na.value = "grey70", midpoint = 0)
p <- gheatmap(p, sig.df [2], offset=30, width=.25, font.size=1.5, colnames = F)+
  scale_fill_gradient2(low = rampCol1, high = rampCol2, mid = "white", na.value = "grey70", midpoint = 0)



## Integrating node annotation by lipid class
sig.df_1 <- sig.df
cmpds <- read.csv("at_test_lipid_color_classification.csv", header= TRUE)
#color by the property (Using HG classification)
sortMat2 = cmpds
#colorby property (Property is desc. col. name like p-value)
HG = cmpds$HG
cols10=c("HexCer"="purple","OH_HexCer"="red","Cer" = "navy","SM" = "orange")
to_plot3= as.data.frame(cbind(sig.df_1[,c("Case1_v_Ctrl_Log2FC")], 
                              sig.df_1[,c("Case2_v_Ctrl_Log2FC")])) 
names(to_plot3)=c("1","2")
rownames(to_plot3) = rownames(sig.df_1)


## Case vs. Control dendrogram

#Generate Circular Dendrogram
t4 <- ggtree(tupgma, layout="circular", size=1.5) 
#%<+% is a pipe symbol to combine datasets more efficiently
#merging circ. dend. w/ pvalue color assignments
t4 <- t4 %<+% cmpds +
  geom_tippoint(aes(color=HG), size=3.0, alpha =0.7, shape = 16)+
  scale_color_brewer(palette = "Dark2")+
theme(legend.position="bottom",legend.text=element_text(size=10))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme(text = element_text(size = 0.01)) +
  geom_treescale(x = NULL, y = NULL, width = 0.5, offset = 30,
                 color = "white", linesize = 1E-100, fontsize = 1E-100)
t4

#plots circular dendrogram with layered heatmap
t4 <- gheatmap(t4, to_plot3[2], offset = 0, width = 0.12, colnames =T, colnames_angle = 50)  +
  scale_fill_gradient2(low = rampCol1, high = rampCol2,  mid = "white", na.value = "#fee090", midpoint = 0)
t4 <- gheatmap(t4, to_plot3[1], offset = 0.05, width = 0.12, colnames =T, colnames_angle = 50)  +
  scale_fill_gradient2(low = rampCol1, high = rampCol2,  mid = "white", na.value = "#fee090", midpoint = 0)
open_tree(t4, 50) %>% rotate_tree(50)
open_tree(t4, 30) 
t4
#coloring options
scale_fill_brewer(palette = "RdBu")
scale_color_manual(values =cols10) 
  

t4 <- gheatmap(t4, to_plot3[2], offset = 0, width = 0.12, colnames =T, colnames_angle = 50)  +
  scale_color_distiller(palette = "Set1")
t4 <- gheatmap(t4, to_plot3[1], offset = 0.05, width = 0.12, colnames =T, colnames_angle = 50)  +
  scale_color_distiller(palette = "Set1")
open_tree(t4, 50) %>% rotate_tree(50)

