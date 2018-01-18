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

a <- c('knight', 'raj')[a]

message('this is chromosome: ', k)
message('this is dataset: ', a)

FILE <- sprintf("model_output/chrm%s/%s_gwas_stl_rf_mods/ichip", k, a)
# FILE <- sprintf("model_output/chrm%s/%s_ichip_stl_rf_mods", k, a)

# #create progress indicators
# ind <- rep(FALSE, 22)
# for(a in c('knight', 'raj')){
#     file_ <- sprintf('ind-%s-slt-rf-ichip.rds', a)
#     message(file_)
#     saveRDS(ind, file = file_)
# }

#
# k <- 6
# a <- 'raj'
# message(sprintf("rm /scratch/ng414/twas/model_output/chrm%s/%s_ichip_stl_rf_mods/*", k, a))
#
#
# ##check whether chrm k has already been completed, if so stop
# mm <- readRDS(sprintf('ind-%s-slt-rf-ichip.rds', a))
# if(mm[k]) stop(sprintf('chrm %s already done!', k))

#--------------------------------------------------------------------------------
#run STL model
run.stl <- function(j){
  p <- pheno.names[i]
  mod.f <- file.path(FILE, sprintf("%s-%s.rds", j, p))

  ##status == TRUE: fit the model, status == FALSE: proceed to next gene
  if(!file.exists(mod.f)) status <- TRUE
  if(file.exists(mod.f)) {
    #if file exists but cannot be read set global status to TRUE
    tryCatch(mod <- readRDS(mod.f), error = function(e) status <<- TRUE)
    if(exists('mod')) { if(class(mod) != 'randomForest') status <- TRUE else status <- FALSE }
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
    saveRDS(mod, mod.f)
    message(sprintf("STL: chrm %s gene %s: %s ok!", k, p, j))
  }

  message(sprintf("STL: chrm %s gene %s: %s done!", k, p, j))
  rsq <- mod$rsq[500]
  names(rsq) <- paste(p, j, sep = ".")
  rsq
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
  rm(A)

  (load("raj-cd4-probes.RData")) # probes (same as raj-cd14-probes.RData)
  p22 <- probes[which(probes$seqname == k), c('ID', 'seqname', 'RANGE_START', 'RANGE_STOP')]
  colnames(p22) <- c("Name", "chromosome_name", "start_position", "end_position")
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# constructing .mult objects in preparation for analysis
#-----------------------------------------------------------------------------------------
##iterating over cell types ---->>>>
rsq <- vector('list', length(pheno.names))
names(rsq) <- pheno.names

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
  ##weeding out probes with nothing in the SNP catchment area

  ww <- sapply(iind, FUN = function(j){
    pos.y <- t(p22[j, c("start_position", "end_position")])
    pos.x <- mysnps$position
    sel.x <- pmin(abs(pos.x - pos.y[1]), abs(pos.x - pos.y[2])) < DIST
    sum(sel.x)
    }
  )

  iind <- names(which(ww > 1))
  ##--------_>>>>>>>>>>>>>>>>>

  #--------------------------------------------------------------------------------
  #RUN procedure
  #--------------------------------------------------------------------------------

  res <- mclapply(iind, FUN = function(j){
          message(pheno.names[i])

          message(j)
          run.stl(j)
              }, mc.cores = 7) %>% do.call(c, .)

          rsq[[i]] <- res
}

saveRDS(rsq, file = sprintf("model_output/chrm%s/gwas_stl_rf_%s_ichip_rsq.rds", k, a))

#--------------------------------------------------------------------------------------------
##update vector containing info on which chrm has already been completed

mm <- readRDS(sprintf('ind-%s-slt-rf-ichip.rds', a))
mm[k] <- TRUE
saveRDS(mm, file = sprintf('ind-%s-slt-rf-ichip.rds', a))


q('no')

#------------------------------------------------------------------------

#create progress indicators
for(a in c('knight', 'raj')){
    file_ <- readRDS(sprintf('/scratch/ng414/twas/ind-%s-slt-rf-ichip.rds', a))
    print(file_)
}

for(a in c('knight', 'raj')){
  for(k in 1 : 22){
    val <- file.exists(sprintf("/scratch/ng414/twas/model_output/chrm%s/gwas_stl_rf_%s_ichip_rsq.rds", k, a))
    message('dat: ', a, '; chrm: ', k, ' -> ', val)
      }
    message('')
}


lapply(c('knight', 'raj'), FUN = function(a) {
  aa <- ifelse(a == 'knight', 1, 2)
  out <- sapply(1 : 22, FUN = function(k) {
  file.exists(sprintf("/scratch/ng414/twas/model_output/chrm%s/gwas_stl_rf_%s_ichip_rsq.rds", k, a))
    })
       ind <- which(!out)
       if(length(ind) > 1){
         message(sprintf('export dat=%s', aa))
         message(sprintf('  for k in {%s};', paste(ind, collapse = ',')))
         message('    do')
         message('      export chr=$k')
         message('      sbatch gwas_ichip_lasso-07082017.sh')
         message('done')
         message('')
     } else {
       message(sprintf('export dat=%s k=%s; sbatch gwas_stl_ichip_mods-09112017.sh', aa, ind))
       message('')
     }
})


for a in {1,2}
  do export dat=$a
    for k in {1..22};
      do
        export chr=$k
          sbatch gwas_stl_ichip_mods-09112017.sh
  done
done

export dat=2
for k in {1..3};
    do
      export chr=$k
        sbatch gwas_stl_ichip_mods-09112017.sh
done



for(k in 1 : 22){
  for(a in c('knight', 'raj')){
    pred.rf <- readRDS(sprintf("/scratch/wallace/twas/model_output/chrm%s/gwas_stl_rf_%s_ichip_res.rds", k, a))
    pred.l <- readRDS(sprintf("/scratch/wallace/twas/model_output/chrm%s/gwas_lasso_%s_ichip_res.rds", k, a))
    message(all.equal(sapply(pred.l, dim), sapply(pred.rf, dim)))
    message(a, ' - ', k)
  }
}








##
