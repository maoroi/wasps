### initial look at the data and tree of Epiponini
library("phytools")
library("caper")
library("corHMM")
library("hisse")
library("diversitree")
library("tidyverse")
set.seed(88)

# 1. data and tree ------------------------------------------------------------
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
for (j in 1:Ntip(tree)) {
    rr <- which(taxa$tip == tree$tip.label[j])
    tree$tip.label[j] <- taxa$binomial[rr]
}
#write.nexus(tree, file="EpiponiniBinomials.nex")

# remove ourgroups (non-Epiponine)
Epis <- getMRCA(tree, c("Angiopolybia_pallens","Apoica_arborea"))
etree <- extract.clade(tree, Epis)
# spelling 
dat$genus_species[which(dat$genus_species == "Asteloeca_ujhelyii")] <- "Asteloeca_ujhelyi" 

# remove tips w/o data
dtree <- drop.tip(etree, setdiff(etree$tip.label, dat$genus_species))

# taxa w phylogenetic and trait data
dtaxa <- dat[which(dat$genus_species %in% dtree$tip.label),]
length(unique(dtaxa$genus_species)) == 50 # 50 species w data on this tree - agreeing with Phoebe's email from 26/05/2021
dtaxa$caste_differentiation[which(dtaxa$caste_differentiation == "Morphologically similar")] <- 1
dtaxa$caste_differentiation[which(dtaxa$caste_differentiation == "Isometric")] <- 2
dtaxa$caste_differentiation[which(dtaxa$caste_differentiation == "Allometric")] <- 3

# 2. data distribution plots --------------------------------------------------

## plotting data species (blue) on the full tree
label <- rep("black", Ntip(etree))
names(label) <- etree$tip.label
label[which(!etree$tip.label %in% dat$genus_species)] <- "red"   # no data
plot(etree, type='fan', cex=0.7, tip.color = label[etree$tip.label])

# plotting data on the pruned tree
label1 <- character(Ntip(dtree))
label2 <- numeric(Ntip(dtree))
names(label1) <- names(label2) <- dtree$tip.label
# map traits onto tips
for (i in 1:Ntip(dtree)){
    state <- dat$caste_differentiation[which(dat$genus_species == dtree$tip.label[i])]
    if (state == "Morphologically similar") {
        label1[dtree$tip.label[i]] <- "grey80"
    } else if (state == "Isometric") {
        label1[dtree$tip.label[i]] <- "skyblue"
    } else if (state == "Allometric") {
        label1[dtree$tip.label[i]] <- "blue"
    } else {
        label1[dtree$tip.label[i]] <- "red"
    }
    label2[dtree$tip.label[i]] <- log(dat$colony_size_max[which(dat$genus_species == dtree$tip.label[i])])
}

# plotting flat tree with coloured dots for trait state
#plot(dtree, cex=0.7, label.offset = 2.5)
pdf(file = "Epiponine_ColSize_CastDiff_names.pdf", height=6, width=9)
plotTree.wBars(dtree,label2, col=label1[dtree$tip.label], tip.labels = TRUE)
legend("bottomleft", legend = c("Allometric","Isometric","Morphologically similar"), 
       col = c("blue","lightblue","grey80"), pch=15, bty = "n", cex=1.5)
dev.off()

# plotting circular tree with coloured bars for trait state
plot(dtree, cex=0.7, type='fan', label.offset = 4) #open.angle=90
ring(2, dtree, style = "ring", offset = 1, col = label1[dtree$tip.label])

pdf(file = "Epiponine_ColSize_CastDiff_circplot.pdf", height=5, width=5)
plotTree.wBars(dtree,label2, type='fan', cex=0.5, col="grey40")
ring(2, dtree, style = "ring", offset = -3, col = label1[dtree$tip.label])
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


# 3. evolutionary models ------------------------------------------------------

## 3.1 phyogenetic signal, PGLS etc -------------------------------------------
# PGLS requires a continuous response and predictor to give reasonable results so I leave it for now

## from phytools: http://www.phytools.org/Cordoba2017/ex/4/PGLS.html
bm <- corBrownian(1, dtree)

library(nlme)
model1 <- gls(colony_size_max ~ caste_differentiation, data=dat, correlation=bm)

label <- label1
label <- gsub("grey80","Similar",label)
label <- gsub("skyblue","Isometric",label)
label <- gsub("blue","Allometric",label)
mkmod <- fitMk(dtree, label, model="ARD")
plot(mkmod, title="Mk model fit to differentitation data")

allo <- dat[which(dat$caste_differentiation == "Allometric"),]
iso <- dat[which(dat$caste_differentiation == "Isometric"),]


## 3.2 rates of evolution (HiSSE, MuHiSSE) ------------------------------------
### 3.2.1 Estimating sampling fraction  ---------------------------------------
## [Proportion sampled out of known species (incl. not in the tree), by differentiation mode, assuming unbiased sampling] 

# sampling proportion same as data distribution (assumption of unbiased sampling)
allEpi <- 246           # total Epiponini (Menezes et al. 2020)
f <- nrow(dat) / allEpi

### 3.2.2 Setting up (Mu)HiSSE analysis ---------------------------------------
# dtree has a polytomy, I resolved randomly
phy <- multi2di(dtree, random = TRUE, equiprob = TRUE)

for (i in 1:length(phy$tip.label)){
    phy$tip.state[i] <- as.numeric(dtaxa$caste_differentiation[which(dtaxa$genus_species == phy$tip.label[i])])
}

states <- data.frame(phy$tip.label, phy$tip.state, phy$tip.state)
states_trans <- states
for(i in 1:Ntip(phy)){
    if(states[i,2] == 1){
        states_trans[i,2] = 0
        states_trans[i,3] = 0
    }
    if(states[i,2] == 2){
        states_trans[i,2] = 0
        states_trans[i,3] = 1
    }
    if(states[i,2] == 3){
        states_trans[i,2] = 1
        states_trans[i,3] = 1
    }
}

freq = c(f, f, 0, f)

# null model - same diversification throughout; no direct transition from similar to allometric 
turnover <- c(1,1,0,1)  
extinction.fraction <- c(1,1,0,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))   # eliminate unused state
diag(trans.rate.mod) <- NA
dull.null <- MuHiSSE(phy=phy, data=states_trans, f=freq, 
                     turnover=turnover, eps=extinction.fraction, 
                     hidden.states=FALSE,
                     trans.rate=trans.rate.mod)
saveRDS(dull.null, file="wasp_null_constrained.RDS")

# null model - same diversification throughout; no direct transition from similar to allometric 
turnover <- c(1,1,0,1)  
extinction.fraction <- c(1,1,0,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0, include.diagonals = TRUE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,7,8,9,12))   # eliminate unused state
diag(trans.rate.mod) <- NA
free.null <- MuHiSSE(phy=phy, data=states_trans, f=freq, 
                     turnover=turnover, eps=extinction.fraction, 
                     hidden.states=FALSE,
                     trans.rate=trans.rate.mod)
saveRDS(free.null, file="wasp_null_free.RDS")


# MuSSE model - only differentiation mode affects diversification
turnover <- c(1,2,0,3)  
extinction.fraction <- c(1,1,0,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))   # eliminate unused state
diag(trans.rate.mod) <- NA
MuSSEmod <- MuHiSSE(phy=phy, data=states_trans, f=freq, 
                     turnover=turnover, eps=extinction.fraction, 
                     hidden.states=FALSE,
                     trans.rate=trans.rate.mod)
saveRDS(MuSSEmod, file="wasp_MuSSE_cons.RDS")

# MuSSE w free transitions
turnover <- c(1,2,0,3)  
extinction.fraction <- c(1,1,0,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0, include.diagonals = TRUE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,7,8,9,12))   # eliminate unused state
diag(trans.rate.mod) <- NA
MuSSEfreemod <- MuHiSSE(phy=phy, data=states_trans, f=freq, 
                    turnover=turnover, eps=extinction.fraction, 
                    hidden.states=FALSE,
                    trans.rate=trans.rate.mod)
saveRDS(MuSSEfreemod, file="wasp_MuSSE_free.RDS")


## ** is it worth testing a model with no direct transition isometric and allometric ?  ** 

# MuHiSSE models
TO_rates <- c(1,2,0,3, 4,5,0,6, 7,8,0,9, 10,11,0,12, 13,14,0,15, 16,17,0,18, 19,20,0,21, 22,23,0,24)
drops_cons <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64)
drops_free <- c(2,5,7,8,9,12, 14,17,19,20,21,24, 26,29,31,32,33,36, 38,41,43,44,45,48, 50,53,55,56,57,60, 62,65,67,68,69,72)

for (n in 2:4) {
    # constrained
    set_to_0 <- 4*(1:n)-1
    turnover <- TO_rates[1:(4*n)]
    extinction.fraction <- rep(1,4*n) 
    extinction.fraction[set_to_0] <- 0
    trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = FALSE)
    trans.rate.mod <- ParDrop(trans.rate, drops_cons[1:(4*n)])   
    trans.rate.mod[,set_to_0] <- 0
    diag(trans.rate.mod) <- NA
    tmp <- MuHiSSE(phy=phy, data=states_trans, f=freq, 
                        turnover=turnover, eps=extinction.fraction, 
                        hidden.states=TRUE,
                        trans.rate=trans.rate.mod)
    saveRDS(tmp, file=paste0("waspMuHiSSE",n,"cons.RDS"))

    # MuHiSSE w free transitions
    turnover <- TO_rates[1:(4*n)]
    extinction.fraction <- rep(1,4*n) 
    extinction.fraction[set_to_0] <- 0
    trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = TRUE)
    trans.rate.mod <- ParDrop(trans.rate, drops_free[1:(6*n)])   # eliminate unused state
    diag(trans.rate.mod) <- NA
    tmpfree <- MuHiSSE(phy=phy, data=states_trans, f=freq, 
                            turnover=turnover, eps=extinction.fraction, 
                            hidden.states=TRUE,
                            trans.rate=trans.rate.mod)
    saveRDS(tmpfree, file=paste0("waspMuHiSSE",n,"free.RDS"))
}


# 4. Post processing ----------------------------------------------------------
#* Function definitions -------------------------------------------------------

extRDS <- function(file){
    tmp <- readRDS(file)
    res <- tmp[c(1:3)]
}

extRDS_par <- function(file){
    tmp <- readRDS(file)
    res <- tmp$solution
}

allmods <- list.files(pattern = ".RDS")

## 4.1 summarising models -----------------------------------------------------
modat <- data.frame(t(sapply(allmods, extRDS)))
score <- as.data.frame(do.call(cbind, lapply(modat[,1:3], as.numeric)))
rownames(score) <- rownames(modat)
score$npar <- score$states <- score$type <- NA

score$type[which(is.na(str_extract(rownames(score), "free")))] <- "constrained"
score$type[which(!is.na(str_extract(rownames(score), "free")))] <- "free"

# extracting model parameters from file name
for (i in 1:nrow(score)){
    score$states[i] <- as.numeric(str_extract(rownames(score)[i],"[0-9]"))    # no. of hidden states
    score$states[which(is.na(score$states))] <- 1
    if (score$type[i] == "constrained") { 
        score$npar[i] <- 7*score$states[i]+2
    } else if (score$type[i] == "free") { 
        score$npar[i] <- 9*score$states[i]+2
    } 
    score["wasp_null_constrained.RDS",6] <- 6 # surely theres a better way to do this
    score["wasp_null_free.RDS",6] <- 8
    score["wasp_MuSSE_cons.RDS",6] <- 8
    score["wasp_MuSSE_free.RDS",6] <- 10
} 

#model selection
ord <- score[order(score$AICc, decreasing = FALSE),]
write.csv(ord,file='best-fitting models Noll tree.csv')


## 4.2 parameter estimates ----------------------------------------------------
pardat <- data.frame(t(sapply(allmods, extRDS_par)))
params <- pardat[which(is.na(str_extract(names(pardat), ".10.")))]      # remove state "10" (disabled in all models)
params <- params[which(is.na(str_extract(names(params), ".00..11.|.11..00.")))]  # remove 'diagonal' change (disabled in all)
