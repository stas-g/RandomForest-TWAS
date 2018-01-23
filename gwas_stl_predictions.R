library(randomForest)
library(glmnet)
library(magrittr)
library(parallel)
library(plyr)
library(annotSnpStats)

setwd('/rds/user/ng414/hpc-work/twas')

DIST <- 1e7

args_ <- commandArgs(trailingOnly = TRUE)
if(length(args_) < 4) stop('must provide chromosome number, data type, disease and method')
print(args_)
message('class(args_) ', class(args_))
k <- as.numeric(args_[1])
a <- as.numeric(args_[2]) #a = {1, 2} = {knight, raj}
d <- as.numeric(args_[3]) #d = {1, 2, 3} = {t1dgc+wtccc, jia, ms}
m <- as.numeric(args_[4]) #m = {1, 2} = {lasso, stl_rf}

a <- c('knight', 'raj')[a]
d <- c('t1dgc+wtccc', 'jia', 'ms')[d]
m <- c('lasso', 'rf')[m]

message('this is chromosome: ', k)
message('this is dataset: ', a)
message('this is disease: ', d)
message('this is method: ', m)

FILE <- sprintf("model_output/chrm%s/%s_gwas_%s_mods/%s", k, a, m, d)

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

  #loading GWAS data
  if(d == 'jia' | d == "t1dgc+wtccc"){
    new.geno <- readRDS(sprintf("gwas-data/%s/pro-geno/knight-%s.rds", d, k))
    } else {
      new.geno <- readRDS(sprintf("gwas-data/%s/pro-geno/gwas-%s.rds", d, k))
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

    X22[[p]] <- readRDS(sprintf("%s-genotypes/%s/geno-chrm%s.rds", p, d, k))
    tmp <- readRDS(file = sprintf("gwas-data/%s/pro-geno/%s-%s.rds", d, p, k))
    pp <- intersect(colnames(X22[[p]]), colnames(tmp))

    X22[[p]] <- X22[[p]][, pp]
    mysnps0[[p]] <- snps(tmp[, pp])
    new.geno0[[p]] <- as(tmp[, pp], 'numeric')
    }

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

  iind <- best.iind[[k]][[i]]
  ii <- sort(intersect(rownames(D[[i]]), p22$Name))
  E22 <- D[[i]][ii, ]

  if(a == 'knight'){
        ## order samples
        m <- intersect(colnames(E22), rownames(X22))
        N22 <- X22[m, ]
        E22 <- t(E22[, m])

        E22 <- scale(E22)
        N22 <- as(N22, 'numeric')
    }

  if(a == 'raj'){
        ## order samples
        m <- intersect(colnames(E22), rownames(X22[[i]]))
        N22 <- X22[[i]][m, ]
        E22 <- t(E22[ ,m])

        E22 <- scale(E22)
        N22 <- as(N22, 'numeric')

        mysnps <- mysnps0[[i]]
        new.geno <- new.geno0[[i]]
      }

  #RUN procedure
  res <- mclapply(iind, FUN = function(j){
          message(sprintf('%s: %s - %s' , toupper(m), pheno.names[i], j))

          run.stl(j, m)
              }, mc.cores = 10) %>% do.call(cbind, .)

    dimnames(res) <- list(rownames(new.geno), iind)
    pred[[i]] <- res
}

saveRDS(pred, file = sprintf("model_output/chrm%s/%s_gwas_%s_%s_preds.rds", k, a, m, d))

q('no')

#-----------------------------------------------------------------------------------------
# aggregating results (Rsq)
#-----------------------------------------------------------------------------------------
devtools::install_github('stas-g/MLhelper') #lasso.betas here
library(MLhelper)

#obtain variable importances from randomForest model mod and order from most to least important
rf.varimp <- function(mod){
  imp <- importance(mod)
  imp[order(imp[, 1], decreasing = TRUE), ]
}

get.rsq <- function(m){
  mclapply(1 : 22 , FUN = function(k) {
    out <- lapply(c('knight', 'raj'), FUN = function(a) {
      out <- lapply(c('t1dgc+wtccc', 'jia', 'ms'), FUN = function(d) {
        if(a == 'raj' & d == 'ms') return(NA)
        path <- sprintf("model_output/chrm%s/%s_gwas_%s_mods/%s", k, a, m, d)
        files_ <- list.files(path)
        out <- lapply(files_, FUN = function(z){
              mod <- readRDS(file.path(path, z))
              j <- strsplit(z, "-|\\.")[[1]][1]
              cell <- strsplit(z, "-|\\.")[[1]][2]
              if(m == 'lasso'){
                rsq <- mod$glmnet.fit$dev.ratio[mod$glmnet.fit$lambda == mod$lambda.1se]
                if(rsq == 0) best.snp <- NA else best.snp <- lasso.betas(mod)[1]
              } else {
                rsq <- mod$rsq[500]
                best.snp <- rf.varimp(mod)[, 1][1]
              }
              message(sprintf('chr %s dataset %s disease %s: %s', k, a, d, z))
              data.frame(probe = j, cell = cell, rsq = rsq, best.snp = best.snp)
              }) %>% do.call(rbind, .)
         out$disease <- d
         out
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


##
