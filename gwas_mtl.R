library(randomForest)
# devtools::install_github("stas-g/funfun")
library(funfun)
library(magrittr)
library(parallel)
library(plyr)
library(annotSnpStats)

setwd('/scratch/wallace/twas')
#progress ind
DIST <- 1e7

args_ <- commandArgs(trailingOnly = TRUE)
if(length(args_) == 0) stop('must provide chromosome number and tissue type')
print(args_)
message('class(args_) ', class(args_))
k <- as.numeric(args_[1])
message('this is chromosome ', k)
a <- as.numeric(args_[2]) #a = {1, 2} = {knight, raj}
d <- as.numeric(args_[3]) #d = {1, 2, 3} = {t1dgc+wtccc, jia, ms}

a <- c('knight', 'raj')[a]
d <- c('t1dgc+wtccc', 'jia', 'ms')[d]

message('this is chromosome: ', k)
message('this is dataset: ', a)
message('this is disease: ', d)

FILE <- sprintf("model_output/chrm%s/%s_gwas_preds/%s", k, a, d)

# for(a in c('t1dgc+wtccc', 'jia', 'ms')){
#   for(k in 1 : 22){
#       FILE <- sprintf("model_output/chrm%s/raj_gwas_preds/%s", k, d)
#       x <- list.files(FILE)
#       # file.remove(file.path(FILE, x[grep('-raj-cd4-raj-cd14', x)]))
#       print(x[grep('-raj-cd4-raj-cd14', x)])
#   }
# }

#
# for(k in 1 : 22){
# for(a in c('knight', 'raj', 'ichip')){
#   for(d in c('t1dgc+wtccc', 'jia', 'ms')){
#     # FILE <- sprintf("rm /scratch/wallace/twas/model_output/chrm%s/%s_gwas_preds/%s/*", k, a, d)
#     FILE  <- sprintf("rm /scratch/wallace/twas/model_output/chrm%s/gwas_mlt_%s_%s_res.rds", k, a, d)
#     message(FILE)
#   }
# }
# }
#
#
ind <- rep(FALSE, 22)
for(a in c('knight', 'raj')){
  for(d in c('t1dgc+wtccc', 'jia', 'ms')){
    file_ <- sprintf('ind-%s-%s_preds.rds', a, d)
    message(file_)
    saveRDS(ind, file = file_)
  }
}

##check whether chrm k has already been completed, if so stop
mm <- readRDS(sprintf('ind-%s-%s_preds.rds', a, d))
if(mm[k]) stop(sprintf('chrm %s already done!', k))

#--------------------------------------------------------------------------------
#z = correlation matrix for various traits
sieve <- function(z){
    z <- abs(z)
    cor.m <- colMeans(z)
    #does the column/trait have at least one cor_pw < 0.3
    ind <- colSums(z < 0.3) > 0
    print(cor.m)
    while(any(ind)){
      #if z is a 2x2 matrix and there is still some cor_pw < 0.3, then there are no subgroups of interest, so stop
      if(nrow(z) < 3) {print('empty'); return(NULL)}
      i <- names(which.min(cor.m[ind]))
      i <- which(colnames(z) == i)
      z <- z[-i, -i]
      cor.m <- colMeans(z)
      ind <- colSums(z < 0.3) > 0
      print(cor.m)
    }
    names(ind)
}

#iterative application of the sieve algorithm to split all variables in z into "best" correlated groups
#if sieve returns NA, then perform STL
#if sieve returns a vector of length 5, then perform full MTL on all the tissues
#if sieve returns a vector of length 2 or 3, then sieve for another cluster
#z = correlation matrix as above

group <- function(z){
  N <- ncol(z)
  gr.list <- list()
  while(sum(lengths(gr.list)) < N & !is.null(z)){
      k <- length(gr.list)
      tr <- colnames(z)
      #best group from traits tr of z
      tmp <- sieve(z)
      if(is.null(tmp)) {
        gr.list <- c(gr.list, as.list(tr))
        } else {
          if(length(tmp) == ncol(z) - 1){
            gr.list <- c(gr.list, list(tmp), list(setdiff(tr, tmp)))
          } else {
            gr.list <- c(gr.list, list(tmp))
            if(length(tmp) == ncol(z)) z <- NULL else z <- z[setdiff(tr, tmp), setdiff(tr, tmp)]
            #if(length(tr) != length(tmp)) z <- z[setdiff(tr, tmp), setdiff(tr, tmp)] #!!! what happens if we get a whole group (after the split)
          }
    }
    message('')
  }
  gr.list
}

#--------------------------------------------------------------------------------
#run MTL model
#fit mlt rf given a gene j and a group of traits tr; save the model & return test predictions
run.mtl <- function(j, tr, a){
  pred.f <- file.path(FILE, sprintf("%s-%s.rds", j, paste(tr, collapse = '-')))

  ##status == TRUE: fit the model, status == FALSE: proceed to next gene
  if(!file.exists(pred.f)) status <- TRUE
  if(file.exists(pred.f)) {
    #if file exists but cannot be read set global status to TRUE
    tryCatch(pred.sp <- readRDS(pred.f), error = function(e) status <<- TRUE)
    if(exists('pred.sp')) { if(!(class(pred.sp) %in% c('matrix', 'numeric', 'data.frame'))) status <- TRUE else status <- FALSE }
  }

  if(status) {
    y <- D.mult[trait.id %in% tr, j]
    X <- X.mult[trait.id %in% tr, ]
    id <- droplevels(trait.id[trait.id %in% tr])
    X.new <- new.geno[tr]
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
    sel.x <- pmin(abs(pos.x - pos.y[1]), abs(pos.x - pos.y[2])) < DIST
    X <- X[, sel.x, drop = FALSE]
    X.new <- X.new[, sel.x, drop = FALSE]

    # ind.g <- intersect(colnames(X), colnames(X.new))
    # X.new <- X.new[, ind.g]
    # X <- X[, ind.g]

    suppressWarnings(X <- data.frame(id = id, X))
    suppressWarnings(X.new <- data.frame(id = new.id, X.new))

    set.seed(123)
    mod <- randomForest(X, y, ntree = 500, importance = TRUE)
    pred <- predict(mod, X.new)
    #split pred values by trait
    pred.sp <- split(pred, X.new$id) %>% do.call(cbind, .)
    colnames(pred.sp) <- paste0(tr, ".", j, "_", k)

    saveRDS(pred.sp, file = pred.f)
    message(sprintf("MTL: chrm %s gene %s: %s ok!", k, j, paste(tr, collapse = ' ')))
    }

  message(sprintf("MTL: chrm %s gene %s: %s done!", k, j, paste(tr, collapse = ' ')))

  pred.sp
}

#run STL model
run.stl <- function(j, tr, a){

  pred.f <- file.path(FILE, sprintf("%s-%s.rds", j, tr))

  ##status == TRUE: fit the model, status == FALSE: proceed to next gene
  if(!file.exists(pred.f)) status <- TRUE
  if(file.exists(pred.f)) {
    #if file exists but cannot be read set global status to TRUE
    tryCatch(pred <- readRDS(pred.f), error = function(e) status <<- TRUE)
    if(exists('pred')) { if(!(class(pred) %in% c('numeric', 'data.frame'))) status <- TRUE else status <- FALSE }
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

    # ind.g <- intersect(colnames(X), colnames(X.new))
    # X.new <- X.new[, ind.g]
    # X <- X[, ind.g]

    set.seed(123)
    mod <- randomForest(X, y, ntree = 500, importance = TRUE)
    pred <- predict(mod, X.new)

    pred <- data.frame(x = pred)
    colnames(pred) <- paste0(tr, ".", j, "_", k)

    saveRDS(pred, file = pred.f)
    message(sprintf("STL: chrm %s gene %s: %s ok!", k, j, paste(tr, collapse = ' ')))
  }

  message(sprintf("STL: chrm %s gene %s: %s done!", k, j, paste(tr, collapse = ' ')))
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

  #loading GWAS data
  if(d == 'jia' | d == "t1dgc+wtccc"){
    tmp <- readRDS(sprintf("gwas-data/%s/pro-geno/knight-%s.rds", d, k))
    } else {
      tmp <- readRDS(sprintf("gwas-data/%s/pro-geno/gwas-%s.rds", d, k))
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

    X22[[p]] <- readRDS(sprintf("%s-genotypes/%s/geno-chrm%s.rds", p, d, k))
    tmp <- readRDS(file = sprintf("gwas-data/%s/pro-geno/%s-%s.rds", d, p, k))
    pp <- intersect(colnames(X22[[p]]), colnames(tmp))

    X22[[p]] <- X22[[p]][, pp]
    mysnps[[p]] <- snps(tmp[, pp])
    new.geno[[p]] <- as(tmp[, pp], 'numeric')
    }

  (load("raj-cd4-probes.RData")) # probes (same as raj-cd14-probes.RData)
  p22 <- probes[which(probes$seqname == k), c('ID', 'seqname', 'RANGE_START', 'RANGE_STOP')]
  colnames(p22) <- c("Name", "chromosome_name", "start_position", "end_position")

  pvals <- readRDS(sprintf('model_output/chrm%s/pvals_min_raj.rds', k))
  #which genes have sufficient signal in at least one tissue type?
  pvals <- pvals[apply(pvals, 1, FUN = function(x) any(x <= 1e-07)), ]
}
#-----------------------------------------------------------------------------------------

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
      set.seed(100)

      X.mult[[i]] <- N22
      D.mult[[i]] <- E22
      message(i)
    }
}

#trait id vector
z <- sapply(X.mult, nrow) #same as sapply(D.mult, nrow)
trait.id <- factor(rep(names(z), z), levels = pheno.names)

X <- X.mult
if(a == 'raj'){
  pp <- intersect(colnames(X.mult[[1]]), colnames(X.mult[[2]]))
  for(i in 1 : 2) X.mult[[i]] <- X.mult[[i]][, pp]
}
X.mult <- do.call(rbind, X.mult)

#-----------------------------------------------------------------------------------------
# running grouping algorithm
#-----------------------------------------------------------------------------------------
#for each gene j (with at least one trait pval_min < 1e-7) return best groups
trait.groups <- lapply(rownames(pvals), FUN = function(j){
  message()
  message(j)
  zz <- lapply(D.mult, FUN = function(x) as.data.frame(t(x[, j, drop = FALSE])))
  zz <- t(rbind.fill(zz))
  colnames(zz) <- names(D.mult)
  z <- cor(zz, method = "pearson", use = "complete.obs")

  group(z)
  })

all(sapply(trait.groups, FUN = function(x) sum(lengths(x))) == length(pheno.names))
names(trait.groups) <- rownames(pvals)

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
}, mc.cores = 1) %>% do.call(cbind, .)

rownames(res) <- rownames(new.geno[[1]])

saveRDS(res, file = sprintf("/scratch/ng414/twas/model_output/chrm%s/gwas_mtl_%s_%s_res.rds", k, a, d))

#--------------------------------------------------------------------------------------------
##update vector containing info on which chrm has already been completed

mm <- readRDS(file = sprintf('ind-%s-%s_preds.rds', a, d))
mm[k] <- TRUE
saveRDS(mm, file = sprintf('ind-%s-%s_preds.rds', a, d))


q('no')

##dealing with hla --> not finished (28.09.17)
load('/scratch/ng414/twas/model_output/hla_probes.RData') #hla.r,hla.k
k <- 6
for(a in c('knight', 'raj')){
  if(a == 'knight') {hla <- hla.k} else {hla <- hla.r}
  for d in c('t1dgc+wtccc', 'jia', 'ms'){
      res <- readRDS(sprintf("/scratch/ng414/twas/model_output/chrm%s/gwas_mtl_%s_%s_res.rds", k, a, d))
      res <- res[, !(colnames(res) %in% hla)]
      ###!!! need to extract probe names from colnames, then subset by hla (28.09.17)
      saveRDS(res, sprintf("/scratch/ng414/twas/model_output/chrm%s/gwas_mtl_%s_%s_res.rds", k, a, d))
  }
}


for k in {1..22}
  do -r rm model_output/chrm$k/knight_gwas_preds/ms/*
  done


a <- c('knight', 'raj')[a]
d <- c('t1dgc+wtccc', 'jia', 'ms')[d]

export dat=1
for d in {1,2,3}
  do export dd=$d
    for k in {1..22};
      do
        export chr=$k
          sbatch mtl_gwas-19072017.sh
  done
done


export dat=2
for d in {1,2}
  do export dd=$d
    for k in {1..22};
      do
        export chr=$k
          sbatch mtl_gwas-19072017.sh
  done
done

export dat=2
for d in {1,2}
  do export dd=$d
    for k in {1..3};
      do
        export chr=$k
          sbatch mtl_gwas-19072017.sh
  done
done


#PRINT ind vectors
for(a in c('knight', 'raj')){
  for(d in c('t1dgc+wtccc', 'jia', 'ms')){
    file_ <- readRDS(sprintf('/scratch/wallace/twas/ind-%s-%s.rds', a, d))
    message(a, ' ', d)
    print(file_)

  }
}

for(a in c('knight', 'raj')){
  for(d in c('t1dgc+wtccc', 'jia', 'ms')){
    for(k in 1 : 22){
      val <- file.exists(sprintf("/scratch/wallace/twas/model_output/chrm%s/gwas_mtl_%s_%s_res.rds", k, a, d))
      message('dat: ', a, '; disease: ', d, '; chrm: ', k, ' -> ', val)
    }
  message('')
  }
message('')
}


lapply(c('knight', 'raj'), FUN = function(a) {
  aa <- ifelse(a == 'knight', 1, 2)
  lapply(c('t1dgc+wtccc', 'jia', 'ms'), FUN = function(d) {
    dd <- which(c('t1dgc+wtccc', 'jia', 'ms') == d)
    out <- sapply(1 : 22, FUN = function(k) {
   file.exists(sprintf("/scratch/wallace/twas/model_output/chrm%s/gwas_mtl_%s_%s_res.rds", k, a, d))
    })
       ind <- which(!out)
       if(length(ind) > 1){
         message(sprintf('export dat=%s d=%s', aa, dd))
         message(sprintf('  for k in {%s};', paste(ind, collapse = ',')))
         message('    do')
         message('      export chr=$k')
         message('      sbatch mtl_gwas-19072017.sh')
         message('done')
         message('')
     } else {
       message(sprintf('export dat=%s d=%s k=%s; sbatch mtl_gwas-19072017.sh', aa, dd, ind))
       message('')
     }
  })
})






xx <- best.iind[[22]]
zz <- lapply(1 : 5, FUN = function(i) paste0(names(xx)[i], ".", xx[[i]], "_", 22)) %>% unlist(.)

zz2 <- lapply(1 : length(groups.fin), FUN = function(i){
  z <- unlist(groups.fin[[i]])
  paste0(z, ".", names(groups.fin)[i], "_", 22)
  }) %>% unlist(.)










##
