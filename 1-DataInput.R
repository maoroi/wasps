### initial look at the data and tree of Epiponini
library(phytools)
library(corHMM)
library(hisse)
library(stringr)

setwd("C:/Users/Roi Maor/Desktop/Projects/Wasps/workspace")

tree <- read.nexus('Epiponini.nex')
# remove duplicate tips to leave one per species
dups <- c("Polybia33", "Agelaia14", "Pseudopolybia2") # Polybia rejecta_rejecta, Agelaia pallipes_morfo_sp2, Pseudopolybia compressa_morfo_sp2
tree <- drop.tip(tree, dups) 

taxa <- read.table('taxa_block.txt') 
dat <- read.csv('1.0 Correlates_SocEvo_trait_data.csv')
#dat <- read.csv('2.0 Correlates_SocEvo_Dimorphism_Metric.csv')
dat <- dat[order(dat$genus_species),]
dat$genus <- str_to_title(dat$genus) # capitalise genus names

# sort out nomenclature
taxa$binomial <- c(1:nrow(taxa))
for (i in 1:nrow(taxa)) {
    taxa[i,1] <- strsplit(taxa[i,1], "\\[")[[1]][1]
    taxa[i,2] <- gsub("names\"=\"", "", taxa[i,3])
    taxa[i,4] <- strsplit(taxa[i,4], "\"")[[1]][1]
    if (is.na(strsplit(taxa[i,4],"_")[[1]][2])) {  # if no subspecific denomination
        taxa[i,3] <- regmatches(taxa[i,4], regexpr("_", taxa[i,4]), invert = TRUE)[[1]][1]
        taxa[i,4] <- ""
    } else {    # split binomials from subspecific additions
        taxa[i,3] <- regmatches(taxa[i,4], regexpr("_", taxa[i,4]), invert = TRUE)[[1]][1]
        taxa[i,4] <- regmatches(taxa[i,4], regexpr("_", taxa[i,4]), invert = TRUE)[[1]][2]
    }
    taxa[i,5] <- paste0(taxa[i,2], "_", taxa[i,3])
}
colnames(taxa) <- c("tip","Genus","Species","notes", "binomial")

# swap tip labels for true binomials
for (j in 1:length(tree$tip.label)) {
    rr <- which(taxa$tip == tree$tip.label[j])
    tree$tip.label[j] <- taxa$binomial[rr]
}

# conform spelling to phylo
dat$genus_species[which(dat$genus_species == "Asteloeca_ujhelyii")] <- "Asteloeca_ujhelyi"  

## plotting data species (blue) on the tree
label <- rep("black", Ntip(tree))
names(label) <- tree$tip.label
label[which(!tree$tip.label %in% dat$genus_species)] <- "red"   # no data
plot(tree, type='fan', cex=0.7, tip.color = label[tree$tip.label])


# remove tips w/o data
ptree <- drop.tip(tree, setdiff(tree$tip.label, dat$genus_species))

# taxa w phylogenetic and trait data
dtaxa <- dat[which(dat$genus_species %in% tree$tip.label),]
length(unique(dtaxa$genus_species)) == 50 # species on tree w data - matching Phoebe's email from 26/05/2021


# plotting data on the pruned tree
label1 <- character(Ntip(ptree))
label2 <- numeric(Ntip(ptree))
names(label1) <- names(label2) <- ptree$tip.label
# map traits onto tips
for (i in 1:Ntip(ptree)){
    state <- dat$caste_differentiation[which(dat$genus_species == ptree$tip.label[i])]
    if (state == "Morphologically similar") {
        label1[ptree$tip.label[i]] <- "grey80"
    } else if (state == "Isometric") {
        label1[ptree$tip.label[i]] <- "skyblue"
    } else if (state == "Allometric") {
        label1[ptree$tip.label[i]] <- "blue"
    } else {
        label1[ptree$tip.label[i]] <- "red"
    }
    label2[ptree$tip.label[i]] <- log(dat$colony_size_max[which(dat$genus_species == ptree$tip.label[i])])
}

# plotting flat tree with coloured dots for trait state
#plot(ptree, cex=0.7, label.offset = 2.5)
pdf(file = "Epiponine_ColSize_CastDiff_names.pdf", height=6, width=9)
plotTree.wBars(ptree,label2, col=label1[ptree$tip.label], tip.labels = TRUE)
legend("bottomleft", legend = c("Allometric","Isometric","Morphologically similar"), 
       col = c("blue","lightblue","grey80"), pch=15, bty = "n", cex=1.5)
dev.off()

# plotting circular tree with coloured bars for trait state
plot(ptree, cex=0.7, type='fan', label.offset = 4) #open.angle=90
ring(2, ptree, style = "ring", offset = 1, col = label1[ptree$tip.label])

pdf(file = "Epiponine_ColSize_CastDiff_circplot.pdf", height=5, width=5)
plotTree.wBars(ptree,label2, type='fan', cex=0.5, col="grey40")
ring(2, ptree, style = "ring", offset = -3, col = label1[ptree$tip.label])
legend("topleft", legend = c("Morphologically similar","Isometric","Allometric","log(max_col_size)"), 
       col = c("grey90","lightblue","blue", "grey40"), pch=15, bty = "n")
dev.off()

pdf(file="ColSizeVs.Differentiation.pdf")
ggplot(dat) +
    aes(x=caste_differentiation, y=log(colony_size_max)) +
    geom_violin() +
    geom_jitter(width = 0.05) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="D") +
    theme(legend.position="none", plot.title = element_text(size=11)) +
    xlab("")
dev.off()

### evolutionary models

## from phytools: http://www.phytools.org/Cordoba2017/ex/4/PGLS.html
bm <- corBrownian(1, ptree)

library(nlme)
model1 <- gls(colony_size_max ~ caste_differentiation, data=dat, correlation=bm)

label <- label1
label <- gsub("grey80","Similar",label)
label <- gsub("skyblue","Isometric",label)
label <- gsub("blue","Allometric",label)
mkmod <- fitMk(ptree, label, model="ARD")
plot(mkmod, title="Mk model fit to differentitation data")

allo <- dat[which(dat$caste_differentiation == "Allometric"),]
iso <- dat[which(dat$caste_differentiation == "Isometric"),]


modpgls <- pgls(colony_size_max ~ caste_differentiation, data=dat, correlation=)