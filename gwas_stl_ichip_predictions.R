library(randomForest)
library(glmnet)
library(magrittr)
library(parallel)
library(plyr)
library(annotSnpStats)

setwd('/rds/user/ng414/hpc-work/twas')

DIST <- 1e7

args_ <- commandArgs(trailingOnly = TRUE)
if(length(args_) < 4) stop('must provide chromosome number, data type and method')
print(args_)
message('class(args_) ', class(args_))
k <- as.numeric(args_[1])
a <- as.numeric(args_[2]) #a = {1, 2} = {knight, raj}
m <- as.numeric(args_[3]) #m = {1, 2} = {lasso, stl_rf}

a <- c('knight', 'raj')[a]
m <- c('lasso', 'rf')[m]

message('this is chromosome: ', k)
message('this is dataset: ', a)
message('this is method: ', m)

FILE <- sprintf("model_output/chrm%s/%s_gwas_%s_mods/ichip", k, a, m)

#--------------------------------------------------------------------------------
#STL routine (lasso, RF)
run.stl <- function(j, m){
  p <- pheno.names[i]
  mod.f <- file.path(FILE.F, sprintf("%s-%s.rds", j, p))

  ##status == TRUE: fit the model, status == FALSE: proceed to next gene
  if(!file.exists(mod.f)) status <- TRUE
  if(file.exists(mod.f)) {
    #if file exists but cannot be read set global status to TRUE
    tryCatch(mod <- readRDS(mod.f), error = function(e) status <<- TRUE)
    if(exists('mod')) { if(!(class(mod) %in% c('cv.glmnet', 'randomForest'))) status <- TRUE else status <- FALSE }
  }

    y <- E22[, j]
    pos.y <- t(p22[j, c("start_position", "end_position")])
    pos.x <- mysnps$position
    sel.x <- pmin(abs(pos.x - pos.y[1]), abs(pos.x - pos.y[2])) < DIST
    X <- N22[, sel.x, drop = FALSE]
    X.new <- new.geno[, sel.x, drop = FALSE]

  if(status){
    set.seed(123)
    if(m == 'lasso'){
      mod <- cv.glmnet(X, y, family = "gaussian", alpha = 1, standardize = FALSE)
      } else {
          mod <- randomForest(X, y, ntree = 500, importance = TRUE)
      }
    saveRDS(mod, file = mod.f)

    message(sprintf("STL: chrm %s gene %s: %s ok!", k, p, j))
  }

  message(sprintf("STL: chrm %s gene %s: %s done!", k, p, j))
  pred <- predict(mod, X.new)
  pred
}

#-----------------------------------------------------------------------------------------
#loading data and initial pre-processing
#-----------------------------------------------------------------------------------------

##KNIGHT
if(a == 'knight'){
  #cell types`
  pheno.names <- c('IFN', 'LPS24', 'LPS2', 'CD14', 'BCELL')
  #probes that passed initial univariate test weeding
  best.iind <- readRDS('model_output/best_iind.rds')

  X22 <- readRDS(sprintf('knight-genotypes/geno-chrm%s.rds', k))
  (load("probes.RData")) # probes
  p22 <- probes[which(probes$chromosome_name == k), ]

  (load("knight-expression.RData")) # D
  D <- D[-6] #get rid of neutrophils

  new.geno <- get(load(sprintf("gwas-data/ichip/knight-%s.RData", k)))
  if(k == 16){
    new.geno <- new.geno[, colnames(new.geno) != "rs908989"] #missing SNP
  }

  ii <- intersect(colnames(X22), colnames(new.geno))
  X22 <- X22[, ii]
  new.geno <- new.geno[, ii]
  mysnps <- snps(new.geno)
  new.geno <- as(new.geno, 'numeric')
}

##RAJ
if(a == 'raj'){
  #cell types
  pheno.names <- c('raj-cd4', 'raj-cd14')
  #probes that passed initial univariate test weeding
  best.iind <- readRDS('model_output/best_iind_raj.rds')

  D <- X22 <- new.geno0 <- mysnps0 <- vector('list', 2)
  names(D) <- names(X22) <- names(new.geno0) <- names(mysnps0) <- pheno.names
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
    mysnps0[[p]] <- snps(tmp[, pp])
    new.geno0[[p]] <- as(tmp[, pp], 'numeric')
    }
  rm(A) #memory management

  (load("raj-cd4-probes.RData")) # probes (same as raj-cd14-probes.RData)
  p22 <- probes[which(probes$seqname == k), c('ID', 'seqname', 'RANGE_START', 'RANGE_STOP')]
  colnames(p22) <- c("Name", "chromosome_name", "start_position", "end_position")
}

#-----------------------------------------------------------------------------------------
# more pre-processing, running stl
#-----------------------------------------------------------------------------------------

pred <- vector('list', length(pheno.names)) #creating empty lists for results
names(pred) <- pheno.names
##iterating over cell types ---->>>>
for(i in 1 : length(pheno.names)) {
    message('CELL ', pheno.names[i])
    message('')

  if(a == 'knight'){
        iind <- best.iind[[k]][[i]]

    		ii <- sort(intersect(rownames(D[[i]]), p22$Name))
    		E22 <- D[[i]][ii, ]

        ## order samples
        m <- intersect(colnames(E22), rownames(X22))
        N22 <- X22[m, ]
        E22 <- t(E22[, m])

        E22 <- scale(E22)
        N22 <- as(N22, 'numeric')
    }

  if(a == 'raj'){
        iind <- best.iind[[k]][[i]]

        ii <- sort(intersect(rownames(D[[i]]), p22$Name))
        E22 <- D[[i]][ii, ]

        ## order samples
        m <- intersect(colnames(E22), rownames(X22[[i]]))
        N22 <- X22[[i]][m, ]
        E22 <- t(E22[ ,m])

        E22 <- scale(E22)
        N22 <- as(N22, 'numeric')

        mysnps <- mysnps0[[i]]
        new.geno <- new.geno0[[i]]
      }

  ##-------->>>>>>>>>>>>>>>>>
  ##weeding out probes with nothing in the SNP catchment area (which is not uncommon due to clustering of SNPs in the iCHIP data)
  ww <- sapply(iind, FUN = function(j){
    pos.y <- t(p22[j, c("start_position", "end_position")])
    pos.x <- mysnps$position
    sel.x <- pmin(abs(pos.x - pos.y[1]), abs(pos.x - pos.y[2])) < DIST
    sum(sel.x)
    }
  )

  #final list of probes for trainig and predictions
  iind <- names(which(ww > 1))
  ##-------->>>>>>>>>>>>>>>>>

  #RUN procedure
  res <- mclapply(iind, FUN = function(j){
          message(sprintf('%s: %s - %s' , toupper(m), pheno.names[i], j))

          run.stl(j, m)
              }, mc.cores = 10) %>% do.call(cbind, .)

    dimnames(res) <- list(rownames(new.geno), iind)
    pred[[i]] <- res
}

saveRDS(pred, file = sprintf("model_output/chrm%s/%s_gwas_%s_ichip_preds.rds", k, a, m))

q('no')

#-----------------------------------------------------------------------------------------
# aggregating results (Rsq)
#-----------------------------------------------------------------------------------------
# devtools::install_github('stas-g/MLhelper') #lasso.betas here
library(MLhelper)
library(glmnet)
library(randomForest)

#obtain variable importances from randomForest model mod and order from most to least important
rf.varimp <- function(mod){
  imp <- importance(mod)
  imp[order(imp[, 1], decreasing = TRUE), ]
}

get.rsq <- function(m){
  mclapply(1 : 22 , FUN = function(k) {
    out <- lapply(c('knight', 'raj'), FUN = function(a) {
        path <- sprintf("model_output/chrm%s/%s_gwas_%s_mods/ichip", k, a, m)
        files_ <- list.files(path)
        out <- lapply(files_, FUN = function(z){
              mod <- readRDS(file.path(path, z))
              j <- strsplit(z, "-|\\.")[[1]][1]
              cell <- strsplit(z, "-|\\.")[[1]][2]
              if(m == 'lasso'){
                rsq <- mod$glmnet.fit$dev.ratio[mod$glmnet.fit$lambda == mod$lambda.1se]
                if(rsq == 0) {
                  best.snp <- best.snp.coef <- NA
                } else {
                    best.snp.coef <- lasso.betas(mod)[1]
                    best.snp <- names(best.snp.coef)
                  }
              } else {
                rsq <- mod$rsq[500]
                best.snp.coef <- rf.varimp(mod)[, 1][1]
                best.snp <- names(best.snp.coef)
              }
              message(sprintf('chr %s dataset %s disease %s: %s', k, a, d, z))
              data.frame(probe = j, cell = cell, rsq = rsq, best.snp = best.snp, best.snp.coef = best.snp.coef)
              }) %>% do.call(rbind, .)
        out$dat <- a
        out
        }) %>% do.call(rbind, .)
      out$chr <- k
      out
      }, mc.cores = 16
  ) %>% do.call(rbind, .)
}

rsq.lasso <- get.rsq('lasso')
rsq.rf <- get.rsq('stl_rf')

rsq.lasso$ID <- paste(rsq.lasso$probe, rsq.lasso$cell, sep = '.')
ID <- sapply(best.iind, FUN = function(z) sapply(1 : 5, FUN = function(i) paste(z, names(z)[i], sep = '.')))





##
