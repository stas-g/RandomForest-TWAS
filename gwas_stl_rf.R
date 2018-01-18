library(randomForest)
# devtools::install_github("stas-g/funfun")
library(funfun)
library(magrittr)
library(parallel)
library(plyr)
library(annotSnpStats)

# setwd('/scratch/ng414/twas')
setwd('/rds/user/ng414/hpc-work/twas')
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

FILE <- sprintf("/scratch/ng414/twas/model_output/chrm%s/%s_gwas_stl_rf_preds/%s", k, a, d)

# #create progress indicators
# ind <- rep(FALSE, 22)
# for(a in c('knight', 'raj')){
#   for(d in c('t1dgc+wtccc', 'jia', 'ms')){
#     file_ <- sprintf('ind-stl-rf-%s-%s_preds.rds', a, d)
#     message(file_)
#     saveRDS(ind, file = file_)
#   }
# }
#
#
# ##check whether chrm k has already been completed, if so stop
# mm <- readRDS(sprintf('ind-stl-rf-%s-%s_preds.rds', a, d))
# if(mm[k]) stop(sprintf('chrm %s already done!', k))

#--------------------------------------------------------------------------------
#run STL model
run.stl <- function(j){
  p <- pheno.names[i]
  pred.f <- file.path(FILE, sprintf("%s-%s.rds", j, p))

  ##status == TRUE: fit the model, status == FALSE: proceed to next gene
  if(!file.exists(pred.f)) status <- TRUE
  if(file.exists(pred.f)) {
    #if file exists but cannot be read set global status to TRUE
    tryCatch(pred <- readRDS(pred.f), error = function(e) status <<- TRUE)
    if(exists('pred')) { if(!(class(pred) %in% c('matrix', 'numeric', 'data.frame'))) status <- TRUE else status <- FALSE }
  }

  if(status){
    y <- E22[, j]
    pos.y <- t(p22[j, c("start_position", "end_position")])
    pos.x <- mysnps$position
    sel.x <- pmin(abs(pos.x - pos.y[1]), abs(pos.x - pos.y[2])) < DIST
    X <- N22[, sel.x, drop = FALSE]
    X.new <- new.geno[, sel.x, drop = FALSE]

    set.seed(123)
    mod <- randomForest(X, y, ntree = 500, importance = TRUE)
    pred <- predict(mod, X.new)

    saveRDS(pred, file = pred.f)
    message(sprintf("STL: chrm %s gene %s: %s ok!", k, p, j))
  }

  message(sprintf("STL: chrm %s gene %s: %s done!", k, p, j))
  pred
}

#-----------------------------------------------------------------------------------------
#loading data and initial pre-processing
#-----------------------------------------------------------------------------------------

##KNIGHT
if(a == 'knight'){
  pheno.names <- c('IFN', 'LPS24', 'LPS2', 'CD14', 'BCELL')
  best.iind <- readRDS('model_output/best_iind.rds')

  X22 <- readRDS(sprintf('knight-genotypes/geno-chrm%s.rds', k))
  (load("probes.RData")) # probes
  p22 <- probes[which(probes$chromosome_name == k), ]

  (load("knight-expression.RData")) # D
  D <- D[-6]
  # D <- D[[i]]

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
  pheno.names <- c('raj-cd4', 'raj-cd14')
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
#-----------------------------------------------------------------------------------------
# constructing .mult objects in preparation for analysis
#-----------------------------------------------------------------------------------------

##iterating over cell types ---->>>>
pred <- vector('list', length(pheno.names))
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

  ##--------_>>>>>>>>>>>>>>>>>

  #--------------------------------------------------------------------------------
  #RUN procedure
  #--------------------------------------------------------------------------------

  res <- mclapply(iind, FUN = function(j){
          message(j)

          run.stl(j)
              }, mc.cores = 1) %>% do.call(cbind, .)

    rownames(res) <- rownames(new.geno)
    colnames(res) <- iind

    pred[[i]] <- res
}

names(pred) <- pheno.names

saveRDS(pred, file = sprintf("/scratch/ng414/twas/model_output/chrm%s/gwas_stl_rf_%s_%s_res.rds", k, a, d))

#--------------------------------------------------------------------------------------------
##update vector containing info on which chrm has already been completed

mm <- readRDS(sprintf('ind-stl-rf-%s-%s_preds.rds', a, d))
mm[k] <- TRUE
saveRDS(mm, file = sprintf('ind-stl-rf-%s-%s_preds.rds', a, d))


q('no')

#------------------------------------------------------------------------

#create progress indicators
for(a in c('knight', 'raj')){
  for(d in c('t1dgc+wtccc', 'jia', 'ms')){
    file_ <- readRDS(sprintf('/scratch/wallace/twas/ind-lasso-%s-%s.rds', a, d))
    message(a, ' ', d)
    print(file_)
  }
}



lapply(c('knight', 'raj'), FUN = function(a) {
  aa <- ifelse(a == 'knight', 1, 2)
  lapply(c('t1dgc+wtccc', 'jia', 'ms'), FUN = function(d) {
    dd <- which(c('t1dgc+wtccc', 'jia', 'ms') == d)
    out <- sapply(1 : 22, FUN = function(k) {
   file.exists(sprintf("/scratch/ng414/twas/model_output/chrm%s/gwas_stl_rf_%s_%s_res.rds", k, a, d))
    })
       ind <- which(!out)
       if(length(ind) > 1){
         message(sprintf('export dat=%s d=%s', aa, dd))
         message(sprintf('  for k in {%s};', paste(ind, collapse = ',')))
         message('    do')
         message('      export chr=$k')
         message('      sbatch gwas_lasso-08082017.sh')
         message('done')
         message('')
     } else {
       message(sprintf('export dat=%s d=%s chr=%s; sbatch gwas_stl_rf-20092017.sh', aa, dd, ind))
       message('')
     }
  })
})


export dat=1
for dd in {1,2,3}
  do export d=$dd
    for k in {1..22};
      do
        export chr=$k
          sbatch gwas_stl_rf-20092017.sh
  done
done


export dat=2
for dd in {1,2}
  do export d=$dd
    for k in {1..22};
      do
        export chr=$k
          sbatch gwas_stl_rf-20092017.sh
  done
done

export dat=2
for dd in {1,2}
  do export d=$dd
    for k in {1..3};
      do
        export chr=$k
          sbatch gwas_stl_rf-20092017.sh
  done
done













##
