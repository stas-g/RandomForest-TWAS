library(randomForest)
library(magrittr)
library(parallel)
library(plyr)
library(annotSnpStats)

setwd('/scratch/ng414/twas')
DIST <- 1e7

args_ <- commandArgs(trailingOnly = TRUE)
if(length(args_) < 2) stop('must provide chromosome number and data name')
k <- as.numeric(args_[1])
a <- as.numeric(args_[2]) #a = {1, 2} = {knight, raj}

a <- c('knight', 'raj')[a]

message('this is chromosome: ', k)
message('this is dataset: ', a)

FILE.M <- sprintf("model_output/chrm%s/%s_gwas_mods/ichip", k, a)
FILE.P <- sprintf("model_output/chrm%s/%s_gwas_preds/ichip", k, a)

#--------------------------------------------------------------------------------
#FUNCTIONS for grouping traits
#--------------------------------------------------------------------------------

#z = correlation matrix for various traits
sieve <- function(z){
    z <- abs(z)
    cor.m <- colMeans(z)
    print(cor.m)
    #does the column/trait have at least one cor_pw < 0.3?
    ind <- colSums(z < 0.3) > 0
    #while there are still traits with any cor_pw < 0.3 and there are more than two traits in total, delete the cor_pw < 0.3 trait with the smallest average cor_pw
    while(any(ind)){
      #if z is a 2x2 matrix and there is still some cor_pw < 0.3, then there are no subgroups of interest, so stop
      if(nrow(z) < 3) {print('empty'); return(NULL)}
      i <- names(which.min(cor.m[ind]))
      i <- which(colnames(z) == i)
      z <- z[-i, -i]
      cor.m <- colMeans(z)
      print(cor.m)
      ind <- colSums(z < 0.3) > 0
    }
    #names of the traits forming the largest group of variables with every cor_pw > 0.3
    names(ind)
}

#iterative application of the sieve algorithm to split all variables in z into "best" correlated groups
#if sieve returns NA, then there are no suitabe groups and we perform STL
#if sieve returns a vector of length 5, then there is one group comprising all the tissue types and we perform full MTL on all the tissues
#if sieve returns a vector of length 2 or 3, then the seived out 3, resp. 2, traits can potentially be formed into another cluster, so apply sieve again
#z = correlation matrix as above

group <- function(z){
  N <- ncol(z)
  gr.list <- list()
  while(sum(lengths(gr.list)) < N & !is.null(z)){
      k <- length(gr.list)
      tr <- colnames(z)
      #best group from traits tr in full z cor matrix:
      tmp <- sieve(z)
      if(is.null(tmp)) {
        #if can't find a cluster in z, then each trait in z forms its own singleton group
        gr.list <- c(gr.list, as.list(tr))
        } else {
          if(length(tmp) == ncol(z) - 1){
            gr.list <- c(gr.list, list(tmp), list(setdiff(tr, tmp)))
          } else {
            gr.list <- c(gr.list, list(tmp))
            #if all traits in z form a cluster, then z is now empty and we can stop; else delete traits forming a new cluster from z and apply sieve again
            if(length(tmp) == ncol(z)) z <- NULL else z <- z[setdiff(tr, tmp), setdiff(tr, tmp)]
          }
    }
    message('')
  }
  gr.list
}

#--------------------------------------------------------------------------------
#run MTL model
#fit MTL-RF given a gene j and a group of traits tr; save the model & return GWAS predictions
run.mtl <- function(j, tr, a){
  mod.f <- file.path(FILE.M, sprintf("%s-%s.rds", j, paste(tr, collapse = '-')))
  pred.f <- file.path(FILE.P, sprintf("%s-%s.rds", j, paste(tr, collapse = '-')))

  ##status == TRUE: fit the model, status == FALSE: proceed to next gene
  if(!file.exists(mod.f)) status <- TRUE
  if(file.exists(mod.f)) {
    #if file exists but cannot be read set global status to TRUE
    tryCatch(mod <- readRDS(mod.f), error = function(e) status <<- TRUE)
    if(exists('mod')) { if(class(mod) != 'randomForest') {
       status <- TRUE
       pred.sp <- readRDS(pred.f)
        } else { status <- FALSE }
     }
  }

  if(status) {
    y <- D.mult[trait.id %in% tr, j]
    X <- X.mult[trait.id %in% tr, ]
    id <- droplevels(trait.id[trait.id %in% tr])
    X.new <- new.geno[tr]
    #create tissue id for X.new
    new.id <- factor(rep(tr, sapply(X.new, nrow)), levels = tr)
    ##
    if(a == 'raj'){
      pp <- intersect(colnames(new.geno[[1]]), colnames(new.geno[[2]]))
      for(i in 1 : 2) X.new[[i]] <- new.geno[[i]][, pp]
    }
    X.new <- do.call(rbind, X.new)

    pos.y <- t(p22[j, c("start_position", "end_position")])
    if(a == 'knight'){
      pos.x <- mysnps$position
    } else {
      pos.x <- mysnps[[1]][colnames(X), ]$position
    }
    #selecting SNPs in the DIST envelope within the probe
    sel.x <- pmin(abs(pos.x - pos.y[1]), abs(pos.x - pos.y[2])) < DIST
    X <- X[, sel.x, drop = FALSE]
    X.new <- X.new[, sel.x, drop = FALSE]

    suppressWarnings(X <- data.frame(id = id, X))
    suppressWarnings(X.new <- data.frame(id = new.id, X.new))

    set.seed(123)
    mod <- randomForest(X, y, ntree = 500, importance = TRUE)
    pred <- predict(mod, X.new)
    #split pred values by trait
    pred.sp <- split(pred, X.new$id) %>% do.call(cbind, .)
    colnames(pred.sp) <- paste0(tr, ".", j, "_", k)

    saveRDS(mod, file = mod.f)
    saveRDS(pred.sp, file = pred.f)
    message(sprintf("MTL: chrm %s gene %s: %s ok!", k, j, paste(tr, collapse = ' ')))
    }

  message(sprintf("MTL: chrm %s gene %s: %s done!", k, j, paste(tr, collapse = ' ')))

  pred.sp
}

#run STL-RF model for probe j, trace tr and dataset a
run.stl <- function(j, tr, a){
  mod.f <- file.path(FILE.M, sprintf("%s-%s.rds", j, tr))
  pred.f <- file.path(FILE.P, sprintf("%s-%s.rds", j, tr))

  ##status == TRUE: fit the model, status == FALSE: proceed to next gene
  if(!file.exists(mod.f)) status <- TRUE
  if(file.exists(mod.f)) {
    #if file exists but cannot be read set global status to TRUE
    tryCatch(mod <- readRDS(mod.f), error = function(e) status <<- TRUE)
    if(exists('mod')) { if(class(mod) != 'randomForest') {
       status <- TRUE
       pred <- readRDS(pred.f)
        } else { status <- FALSE }
     }
  }

  if(status){
    y <- D[[tr]][, j]
    if(a == 'raj') X22 <- X22[[tr]]
    pos.y <- t(p22[j, c("start_position", "end_position")])
    if(a == 'knight'){
      pos.x <- mysnps$position
    } else {
      pos.x <- mysnps[[tr]]$position
    }
    sel.x <- pmin(abs(pos.x - pos.y[1]), abs(pos.x - pos.y[2])) < DIST
    X <- X[[tr]][, sel.x, drop = FALSE]
    X.new <- new.geno[[tr]][, sel.x, drop = FALSE]

    set.seed(123)
    mod <- randomForest(X, y, ntree = 500, importance = TRUE)
    pred <- predict(mod, X.new)

    pred <- data.frame(x = pred)
    colnames(pred) <- paste0(tr, ".", j, "_", k)

    saveRDS(pred, file = pred.f)
    message(sprintf("STL: chrm %s gene %s: %s ok!", k, j, tr))
  }

  message(sprintf("STL: chrm %s gene %s: %s done!", k, j, tr))
  pred
}

#-----------------------------------------------------------------------------------------
#loading data and initial pre-processing
#-----------------------------------------------------------------------------------------

##KNIGHT
if(a == 'knight'){
  pheno.names <- c('IFN', 'LPS24', 'LPS2', 'CD14', 'BCELL')

  X22 <- readRDS(sprintf('knight-genotypes/geno-chrm%s.rds', k))
  (load("knight-expression.RData")) # D
  D <- D[-6]; names(D) <- pheno.names
  (load("probes.RData")) # probes

  p22 <- probes[which(probes$chromosome_name == k), ]

  pvals <- readRDS(file = sprintf('model_output/chrm%s/pvals_min.rds', k))
  #which genes have sufficient signal in at least one tissue type?
  pvals <- pvals[apply(pvals, 1, FUN = function(x) any(x <= 1e-07)), ]

  tmp <- get(load(sprintf("gwas-data/ichip/knight-%s.RData", k)))
  if(k == 16){
    tmp <- tmp[, colnames(tmp) != "rs908989"] #missing SNP
  }

  ii <- intersect(colnames(X22), colnames(tmp))
  X22 <- X22[, ii]
  tmp <- tmp[, ii]
  mysnps <- snps(tmp)

  new.geno <- vector('list', 5)
  names(new.geno) <- pheno.names
  tmp <- as(tmp, 'numeric')
  for(i in 1 : 5) new.geno[[i]] <- tmp
}


##RAJ
if(a == 'raj'){
  pheno.names <- c('raj-cd4', 'raj-cd14')

  D <- X22 <- new.geno <- mysnps <-vector('list', 2)
  names(D) <- names(X22) <- names(new.geno) <- names(mysnps) <- pheno.names
  for(p in pheno.names){
    trans <- readRDS(sprintf("model_output/trans_%s.rds", p))
    (load(sprintf("%s-expression.RData", p))) # expr
    colnames(expr) <- trans[colnames(expr)]
    D[[p]] <- expr
    rm(expr)

    X22[[p]] <- readRDS(sprintf("%s-genotypes/ichip/geno-chrm%s.rds", p, k))
    tmp <- get(load(file = sprintf("gwas-data/ichip/%s-%s.RData", p, k)))
    if(k == 1 & p == 'raj-cd4'){
      tmp <- tmp[, colnames(tmp) != 'rs8177971'] #missing SNP
    }
    if(k == 16){
      tmp <- tmp[, colnames(tmp) != "rs908989"] #missing SNP
    }
    pp <- intersect(colnames(X22[[p]]), colnames(tmp))

    X22[[p]] <- X22[[p]][, pp]
    mysnps[[p]] <- snps(tmp[, pp])
    new.geno[[p]] <- as(tmp[, pp], 'numeric')
    }
  rm(A)

  (load("raj-cd4-probes.RData")) # probes (same as raj-cd14-probes.RData)
  p22 <- probes[which(probes$seqname == k), c('ID', 'seqname', 'RANGE_START', 'RANGE_STOP')]
  colnames(p22) <- c("Name", "chromosome_name", "start_position", "end_position")

  pvals <- readRDS(sprintf('model_output/chrm%s/pvals_min_raj.rds', k))
  #which genes have sufficient signal in at least one tissue type?
  pvals <- pvals[apply(pvals, 1, FUN = function(x) any(x <= 1e-07)), ]
}

#-----------------------------------------------------------------------------------------
# constructing .mult objects in preparation for analysis
#-----------------------------------------------------------------------------------------

if(a == 'knight'){
  X.mult <- D.mult <- vector("list", length(pheno.names))
  names(D.mult) <- names(X.mult) <- pheno.names
  for(i in 1 : length(D)){

  		ii <- sort(intersect(rownames(D[[i]]), p22$Name))
  		E22 <- D[[i]][ii, ]

      ## order samples
      m <- intersect(colnames(E22), rownames(X22))
      N22 <- X22[m, ]
      E22 <- t(E22[ ,m])

      E22 <- scale(E22)
      N22 <- as(N22, 'numeric')
      set.seed(100)

      X.mult[[i]] <- N22
      D.mult[[i]] <- E22
      message(i)
    }
    #getting rid of the nuisance-some extra gene in chrm 1 Bcell
    if(k == 1){
      D.mult[['BCELL']] <- D.mult[['BCELL']][, colnames(D.mult[['BCELL']]) %in% colnames(D.mult[[1]])]
    }
}


if(a == 'raj'){
  X.mult <- D.mult <- vector("list", length(pheno.names))
  names(D.mult) <- names(X.mult) <- pheno.names
  for(i in 1 : length(D)){

      ii <- sort(intersect(rownames(D[[i]]), p22$Name))
      E22 <- D[[i]][ii, ]

      ## order samples
      m <- intersect(colnames(E22), rownames(X22[[i]]))
      N22 <- X22[[i]][m, ]
      E22 <- t(E22[ ,m])

      E22 <- scale(E22)
      N22 <- as(N22, 'numeric')

      X.mult[[i]] <- N22
      D.mult[[i]] <- E22
      message(i)
    }
}

#trait id vector
z <- sapply(X.mult, nrow) #same as sapply(D.mult, nrow)
trait.id <- factor(rep(names(z), z), levels = pheno.names)

#creating stacked X-matrix
X <- X.mult
if(a == 'raj'){
  #only include probes featuring for both cell types
  pp <- intersect(colnames(X.mult[[1]]), colnames(X.mult[[2]]))
  for(i in 1 : 2) X.mult[[i]] <- X.mult[[i]][, pp]
}
#stacking for different cell types
X.mult <- do.call(rbind, X.mult)
##--------_>>>>>>>>>>>>>>>>>

##weeding out probes with nothing in the SNP catchment area
ww <- sapply(rownames(pvals), FUN = function(j){
  pos.y <- t(p22[j, c("start_position", "end_position")])
  if(a == 'knight'){
    pos.x <- mysnps$position
    # pos.x <- snps(X22)$bp
  } else {
    pp <- unique(rownames(mysnps[[1]]), rownames(mysnps[[2]]))
    pos.x <- mysnps[[1]][pp, ]$position
  }
  sel.x <- pmin(abs(pos.x - pos.y[1]), abs(pos.x - pos.y[2])) < DIST
  sum(sel.x)
  }
)

iind <- names(which(ww > 0))
##--------_>>>>>>>>>>>>>>>>>

#-----------------------------------------------------------------------------------------
# running grouping algorithm
#-----------------------------------------------------------------------------------------
#for each gene j (with at least one trait pval_min < 1e-7) return best groups
trait.groups <- lapply(iind, FUN = function(j){
  message()
  message(j)
  #accumulating expression data across cell types and binding it in a data.frame
  zz <- lapply(D.mult, FUN = function(x) as.data.frame(t(x[, j, drop = FALSE])))
  zz <- t(rbind.fill(zz))
  colnames(zz) <- names(D.mult)
  z <- cor(zz, method = "pearson", use = "complete.obs")

  group(z)
  })
  names(trait.groups) <- iind

#check that for each probe groups contain all cell types
all(sapply(trait.groups, FUN = function(x) sum(lengths(x))) == length(pheno.names))

#sieve out groups for which the pval_min requirement is not satisfied
groups.fin <- lapply(names(trait.groups), FUN = function(j){
  x <- trait.groups[[j]]
  out <- lapply(x, FUN = function(a){
          if(all(pvals[j, a] >  1e-07)) {
            return(NA) } else {
              a
            }
    })
  out[!is.na(out)]
  })
names(groups.fin) <- names(trait.groups)

D <- D.mult
D.mult <- do.call(rbind, D.mult)
#--------------------------------------------------------------------------------
#RUN procedure
#--------------------------------------------------------------------------------

res <- mclapply(names(groups.fin), FUN = function(j){
    out <- lapply(1 : length(groups.fin[[j]]), FUN = function(p){
        tr <- groups.fin[[j]][[p]]
        message(' ')
        message(j, ' ', tr)
        if(length(tr) == 1){
          message('STL')
          run.stl(j, tr, a = a)
        } else {
            message('MTL')
            run.mtl(j, tr, a = a)
        }
  }) %>% do.call(cbind, .)
}, mc.cores = 3) %>% do.call(cbind, .)
rownames(res) <- rownames(new.geno[[1]])

saveRDS(res, file = sprintf("/scratch/ng414/twas/model_output/chrm%s/gwas_mtl_%s_ichip_res.rds", k, a))

#--------------------------------------------------------------------------------------------
##update vector containing info on which chrm has already been completed

mm <-  readRDS(sprintf('ind-%s-ichip_preds.rds', a))
mm[k] <- TRUE
saveRDS(mm, file = sprintf('ind-%s-ichip_preds.rds', a))


q('no')

plot(x, 100 * x, col = 'red', type = 'l', ylim = c(-100, 100))
lines(x,  x/100, col = 'firebrick1')
lines(x, 2 * x, col = 'blue')
lines(x,  x/2, col = 'dodgerblue')



#------------------------------------------------------------------------

#create progress indicators
for(a in c('knight', 'raj')){
    file_ <- readRDS(sprintf('/scratch/wallace/twas/ind-%s-ichip.rds', a))
    print(file_)
}

for(a in c('knight', 'raj')){
  for(k in 1 : 22){
    val <- file.exists(sprintf("/scratch/wallace/twas/model_output/chrm%s/gwas_mtl_%s_ichip_res.rds", k, a))
    message('dat: ', a, '; chrm: ', k, ' -> ', val)
      }
    message('')
}


lapply(c('knight', 'raj'), FUN = function(a) {
  aa <- ifelse(a == 'knight', 1, 2)
  out <- sapply(1 : 22, FUN = function(k) {
  file.exists(sprintf("/scratch/ng414/twas/model_output/chrm%s/gwas_mtl_%s_ichip_res.rds", k, a))
    })
       ind <- which(!out)
       if(length(ind) > 1){
         message(sprintf('export dat=%s', aa))
         message(sprintf('  for k in {%s};', paste(ind, collapse = ',')))
         message('    do')
         message('      export chr=$k')
         message('      sbatch mtl_gwas_ichip-26072017.sh')
         message('done')
         message('')
     } else {
       message(sprintf('export dat=%s k=%s; sbatch mtl_gwas_ichip-26072017.sh', aa, ind))
       message('')
     }
})


for a in {1,2}
  do export dat=$a
    for k in {1..22};
      do
        export chr=$k
          sbatch mtl_gwas_ichip-26072017.sh
  done
done


export dat=2
for k in {1..3};
    do
      export chr=$k
      sbatch mtl_gwas_ichip-26072017.sh
done



#checking all colnames are there
zz <- lapply(c('knight', 'raj'), FUN = function(a) {
  lapply(1 : 22, FUN = function(k) {
    message(a, ' chrm ', k)
    x <- readRDS(sprintf("/scratch/wallace/twas/model_output/chrm%s/gwas_mtl_%s_ichip_res.rds", k, a))
    colnames(x)
  })
})











##
