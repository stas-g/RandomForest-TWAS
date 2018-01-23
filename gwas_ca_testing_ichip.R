#Cochran-Armitage 31.07.2017
library(annotSnpStats)
# devtools::install_github("stas-g/funfun")
library(randomFunctions)
library(magrittr)
library(parallel)

setwd('/rds/user/ng414/hpc-work/twas')
#---------------------------------------------------------
#----------------------------------------------------------

a <- c('raj', 'knight')

vw <- function(x, cc, stratum){
  N <- tapply(x, list(cc, stratum), length)
  M <- tapply(x, list(cc, stratum), mean)
  V <- tapply(x, list(cc, stratum), var)

  m <- M[2, ] - M[1, ]
  v <- V[1 ,]/N[1, ] + V[2, ]/N[2, ]
  D <- sum(m/v)/sum(1/v)
  se.D <- sqrt(1/sum(1/v))
  p.val <- 2 * pnorm(abs(D/se.D), lower.tail = FALSE)
  return(list(D = D, se.D = se.D, p.val = p.val))
}

do.ca <- function(x, cc, stratum = rep(1, length(cc)), ctrl){
    require(randomFunctions)

    if(var(x) == 0) return(data.frame(control = NA, case = NA, diff = NA, p.val = NA, V = 0))

    means <- tapply(x, list(cc, stratum), mean)
    d <- setdiff(rownames(means), ctrl)
    m.expr <- rowMeans(means)
    if(length(unique(stratum)) == 2){
      diff <- vw(x, cc, stratum)[['D']]
    } else {
      diff <- mean(means[d, ] - means[ctrl, ])
    }
    ca.test.aff <- Cochran.Armitage.test(exposure = x, cc = cc, stratum = stratum)
    res <- data.frame(control = m.expr[ctrl], case = m.expr[d], diff = diff, p.val = ca.test.aff[['p.value']], V = var(x), stringsAsFactors = FALSE)
    return(res)
}


#lasso = 1, rf-mtl = 2, rf-stl = 3
gwas.ca <- function(a, experiment){

  res <- mclapply(1 : 22, FUN = function(k){
    if(experiment == 1) {
      pred <- readRDS(sprintf("model_output/chrm%s/gwas_lasso_%s_ichip_res.rds", k, a))
      } else {
          if(experiment == 2){
            pred <- readRDS(sprintf("model_output/chrm%s/gwas_mtl_%s_ichip_res.rds", k, a))
            } else {
              pred <- readRDS(sprintf("model_output/chrm%s/gwas_stl_rf_%s_ichip_res.rds", k, a))
          }
    }
    if(experiment %in% c(1, 3)){
      if(a == 'raj') names(pred) <- c('raj-CD4', 'raj-CD14')
      names_ <- lapply(1 : length(pred), FUN = function(i) paste0(names(pred)[i], ".", colnames(pred[[i]]), "_", k)) %>% unlist(.)
      pred <- do.call(cbind, pred)
      colnames(pred) <- names_
    }

    if(a == 'knight') {
      new.geno <- get(load(sprintf("/rds/user/ng414/hpc-work/twas/gwas-data/ichip/knight-%s.RData", k)))
      } else {
        new.geno <- get(load(file = sprintf("/rds/user/ng414/hpc-work/twas/gwas-data/ichip/raj-cd14-%s.RData", k)))
      }

    mysamples <- samples(new.geno)

    # if(!all.equal(rownames(mysamples), rownames(pred))) stop("mismatch between samples in preds and cc!")

    phenotype <- mysamples[, 'phenotype']
    aff <- mysamples[, 'affected']
    country <- mysamples[, 'country']

    DISEASES <- setdiff(unique(phenotype), c(NA,"CONTROL"))

    #diff = case - control
    out <- lapply(1 : ncol(pred), FUN = function(i){
        x <- pred[, i]
        res.aff <- do.ca(x, aff, country, ctrl = "1")
        # res.aff <- do.ca(x, sample(aff), country, ctrl = "1")

        res <- lapply(DISEASES, FUN = function(d) {
          message(d)
          if(d == 'RA'){
            use <- !is.na(x) & phenotype %in% c("CONTROL", d)
            res <- do.ca(x[use], phenotype[use], country[use], ctrl = "CONTROL")
            # res <- do.ca(x[use], sample(phenotype[use]), country[use], ctrl = "CONTROL")
          } else {
            use <- !is.na(x) & phenotype %in% c("CONTROL", d) & country == 'UK'
            res <- do.ca(x[use], phenotype[use], ctrl = "CONTROL")
            #res <- do.ca(x[use], sample(phenotype[use]), ctrl = "CONTROL")
          }
          res
        }
      ) %>% do.call(rbind, .)
        ans <- rbind(res.aff, res)
        ans$disease <- c('AUTO', DISEASES)

        z <- strsplit(colnames(pred)[i], "\\.|_")[[1]]
        ans$cell <- z[1]
        if(a == 'knight'){
          ans$probe <- paste(z[2], z[3], sep = '_')
          ans$chr <- z[4]
        } else {
            ans$cell[res$cell == 'raj-cd14'] <- 'raj-CD14'
            ans$cell[res$cell == 'raj-cd4'] <- 'raj-CD4'
            ans$probe <- z[2]
            ans$chr <- z[3]
      }

        message(colnames(pred)[i])
        ans
          }) %>% do.call(rbind, .)

    rownames(out) <- NULL

    message("chrm ", k, " done!")
    out
      }, mc.cores = 16)
    res <- do.call(rbind, res)
    return(res)
}


#--------------------------------
# reading in and assembling results
# #//RF-MTL
# raj.res <- gwas.ca('raj', 2)
# saveRDS(raj.res, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_ca0.rds')
# saveRDS(raj.res, file = '/scratch/ng414/twas/model_output/gwas_ichip_raj_ca0.rds')

knight.res <- gwas.ca('knight', 2)
knight.res <- do.call(cbind, knight.res)
knight.res <- as.data.frame(knight.res, stringsAsFactors = FALSE)
knight.res$control <- as.numeric(knight.res$control)
knight.res$case <- as.numeric(knight.res$case)
knight.res$diff <- as.numeric(knight.res$diff)
knight.res$p.val <- as.numeric(knight.res$p.val)
knight.res$V <- as.numeric(knight.res$V)
rownames(knight.res) <- with(knight.res, paste(disease, cell, probe, sep = '.'))
saveRDS(knight.res, file = '/scratch/ng414/twas/model_output/gwas_ichip_knight_ca.rds')
# saveRDS(knight.res, file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_ca.rds') #old version here

# #//LASSO
# raj.res.l <- gwas.ca('raj', 1)
# saveRDS(raj.res.l, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_lasso_ca0.rds')
# saveRDS(raj.res.l, file = '/scratch/ng414/twas/model_output/gwas_ichip_raj_lasso_ca0.rds')

knight.res.l <- gwas.ca('knight', 1)
saveRDS(knight.res.l, file = '/scratch/ng414/twas/model_output/gwas_ichip_knight_lasso_ca.rds')
saveRDS(knight.res.l, file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_lasso_ca.rds')

# #//RF-STL
# raj.res.stl <- gwas.ca('raj', 3)
# saveRDS(raj.res.stl, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_rf_stl_ca0.rds')
# saveRDS(raj.res.stl, file = '/scratch/ng414/twas/model_output/gwas_ichip_raj_rf_stl_ca0.rds')

knight.res.stl <- gwas.ca('knight', 3)
knight.res.stl <- subset(knight.res.stl,!(probe %in% hla.k))
saveRDS(knight.res.stl, file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_rf_stl_ca.rds')
saveRDS(knight.res.stl, file = '/scratch/ng414/twas/model_output/gwas_ichip_knight_rf_stl_ca.rds')

##permutation results
# raj.res.x <- gwas.ca('raj')
# saveRDS(raj.res.x, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_ca_x.rds')
# knight.res.x <- gwas.ca('knight')
# saveRDS(knight.res.x, file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_ca_x.rds')

#

q('no')

####RAJ#### -->>
##adding genes
load("raj-probes.RData")
probes.r <- probes; rm(probes)
load("knight-probes.RData")
raj.genes <- readRDS("gene_details.rds")
colnames(raj.genes) <- c("probe", "ensid", "symbol", "chr")
load('model_output/hla_probes.RData') #hla.r,hla.k

# raj.res <- readRDS('/scratch/wallace/twas/model_output/gwas_ichip_raj_ca0.rds')
raj.res <- gwas.ca('raj', 2)
raj.res <- subset(raj.res, !(raj.res$probe %in% hla.r))
#which probes have mpre than one gene associated to them?
probes.bad <- names(table(raj.genes$probe))[which(table(raj.genes$probe) > 1)]
probes.bad.raj <- unique(raj.res$probe[raj.res$probe %in% probes.bad])
raj.res <- subset(raj.res, !(probe %in% probes.bad.raj))

raj.res <- split(raj.res, raj.res$disease) %>% lapply(., FUN = function(x){
  tmp <- merge(x, raj.genes, by = 'probe', all.x = TRUE)
  out <- subset(tmp, !is.na(symbol))
  print(all(out$chr.x == out$chr.y))
  out <- out[, c(1 : 9, 11)]
  colnames(out)[c(9, 10)] <- c('chr', 'gene')
  rownames(out) <- paste0(out$cell, ".", out$probe, "_", out$chr)
  out
  })
raj.res <- lapply(raj.res, FUN = function(z) subset(z, !(z$probe %in% hla.r)))
saveRDS(raj.res, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_ca.rds')

# raj.res.l <- readRDS('/scratch/wallace/twas/model_output/gwas_ichip_raj_lasso_ca.rds')
raj.res <- gwas.ca('raj', 1)
raj.res.l <- split(raj.res.l, raj.res.l$disease) %>% lapply(., FUN = function(x){
  tmp <- merge(x, raj.genes, by = 'probe', all.x = TRUE)
  ind <- paste0(tmp$cell, ".", tmp$probe, "_", tmp$chr.x)
  tb <- table(ind)
  keep <- names(tb)[tb == 1]
  out <- tmp[ind %in% keep, ]
  out <- out[!is.na(out$symbol), ]
  print(all(out$chr.x == out$chr.y))
  out <- out[, c(1 : 9, 11)]
  colnames(out)[c(9, 10)] <- c('chr', 'gene')
  rownames(out) <- paste0(out$cell, ".", out$probe, "_", out$chr)
  out
  })
raj.res.l <- lapply(raj.res.l, FUN = function(z) subset(z, !(z$probe %in% hla.r)))
saveRDS(raj.res.l, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_lasso_ca.rds')


#rf-stl
raj.res.rf <- readRDS('/scratch/wallace/twas/model_output/gwas_ichip_raj_rf_stl_ca0.rds')
raj.res.rf <- split(raj.res.rf, raj.res.rf$disease) %>% lapply(., FUN = function(x){
  tmp <- merge(x, raj.genes, by = 'probe', all.x = TRUE)
  ind <- paste0(tmp$cell, ".", tmp$probe, "_", tmp$chr.x)
  tb <- table(ind)
  keep <- names(tb)[tb == 1]
  out <- tmp[ind %in% keep, ]
  out <- out[!is.na(out$symbol), ]
  print(all(out$chr.x == out$chr.y))
  out <- out[, c(1 : 9, 11)]
  colnames(out)[c(9, 10)] <- c('chr', 'gene')
  rownames(out) <- paste0(out$cell, ".", out$probe, "_", out$chr)
  out
  })
raj.res.rf <- lapply(raj.res.rf, FUN = function(z) subset(z, !(z$probe %in% hla.r)))
saveRDS(raj.res.rf, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_rf_stl_ca.rds')


#--------------
knight.res <- readRDS('/scratch/ng414/twas/model_output/gwas_ichip_knight_ca.rds')
knight.res <- split(knight.res, knight.res$disease) %>% lapply(., FUN = function(x){
      x$gene <- probes[x$probe, "hgnc_symbol"]
      rownames(x) <- paste0(x$cell, ".", x$probe, "_", x$chr)
      x
  })
knight.res <- lapply(knight.res, FUN = function(z) subset(z, !(z$probe %in% hla.k)))
saveRDS(knight.res, file = '/scratch/ng414/twas/model_output/gwas_ichip_knight_ca.rds')

knight.res.l <- readRDS('/scratch/wallace/twas/model_output/gwas_ichip_knight_lasso_ca.rds')
knight.res.l <- split(knight.res.l, knight.res.l$disease) %>% lapply(., FUN = function(x){
      x$gene <- probes[x$probe, "hgnc_symbol"]
      rownames(x) <- paste0(x$cell, ".", x$probe, "_", x$chr)
      x
  })
knight.res.l <- lapply(knight.res.l, FUN = function(z) subset(z, !(z$probe %in% hla.k)))
saveRDS(knight.res.l, file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_lasso_ca.rds')

knight.res.rf <- readRDS('/scratch/wallace/twas/model_output/gwas_ichip_knight_rf_stl_ca.rds')
knight.res.rf <- split(knight.res.rf, knight.res.rf$disease) %>% lapply(., FUN = function(x){
      x$gene <- probes[x$probe, "hgnc_symbol"]
      rownames(x) <- paste0(x$cell, ".", x$probe, "_", x$chr)
      x
  })
knight.res.rf <- lapply(knight.res.rf, FUN = function(z) subset(z, !(z$probe %in% hla.k)))
saveRDS(knight.res.rf, file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_rf_stl_ca.rds')


##HLA ADJUSTMENT
load('/scratch/ng414/twas/model_output/hla_probes.RData')

setwd('/scratch/wallace/twas')
raj.res <- readRDS('model_output/gwas_ichip_raj_ca.rds')
raj.res.l <- readRDS('model_output/gwas_ichip_raj_lasso_ca.rds')
raj.res.rf <- readRDS('model_output/gwas_ichip_raj_rf_stl_ca.rds')
knight.res.l <- readRDS('model_output/gwas_ichip_knight_lasso_ca.rds')
knight.res <- readRDS('model_output/gwas_ichip_knight_ca.rds')
knight.res.rf <- readRDS('model_output/gwas_ichip_knight_rf_stl_ca.rds')

knight.res <- lapply(knight.res, FUN = function(z) subset(z, !(probe %in% hla.k)))
knight.res.l <- lapply(knight.res.l, FUN = function(z) subset(z, !(probe %in% hla.k)))
knight.res.rf <- lapply(knight.res.rf, FUN = function(z) subset(z, !(probe %in% hla.k)))
raj.res <- lapply(raj.res, FUN = function(z) subset(z, !(probe %in% hla.r)))
raj.res.l <- lapply(raj.res.l, FUN = function(z) subset(z, !(probe %in% hla.r)))
raj.res.rf <- lapply(raj.res.rf, FUN = function(z) subset(z, !(probe %in% hla.r)))

saveRDS(raj.res, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_ca.rds')
saveRDS(raj.res.l, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_lasso_ca.rds')
saveRDS(raj.res.rf, file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_rf_stl_ca.rds')
saveRDS(knight.res.l, file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_lasso_ca.rds')
saveRDS(knight.res, file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_ca.rds')
saveRDS(knight.res.rf, file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_rf_stl_ca.rds')

saveRDS(raj.res, file = '/scratch/ng414/twas/model_output/gwas_ichip_raj_ca.rds')
saveRDS(raj.res.l, file = '/scratch/ng414/twas/model_output/gwas_ichip_raj_lasso_ca.rds')
saveRDS(raj.res.rf, file = '/scratch/ng414/twas/model_output/gwas_ichip_raj_rf_stl_ca.rds')
saveRDS(knight.res.l, file = '/scratch/ng414/twas/model_output/gwas_ichip_knight_lasso_ca.rds')
saveRDS(knight.res, file = '/scratch/ng414/twas/model_output/gwas_ichip_knight_ca.rds')
saveRDS(knight.res.rf, file = '/scratch/ng414/twas/model_output/gwas_ichip_knight_rf_stl_ca.rds')


scp ng414@darwin:/scratch/wallace/twas/model_output/gwas_ichip_{knight,raj}_ca.rds model_output/
scp ng414@darwin:/scratch/wallace/twas/model_output/gwas_ichip_{knight,raj}_ca_x.rds model_output/

###analysis (08.08.2017)
#james' function 04.08.2017
library('magrittr')
QQ=function(x,l=0.99,add=FALSE,minx=TRUE,ab=T,...) {
  if (max(abs(x),na.rm=T)<1.1) x=-log10(x)
  if (missing(l) & length(x)<10000 & length(x)>500) l=(length(x)-500)/length(x)
  if (missing(l) & length(x)<501) l=0
  n=length(x); q=-log10((n:1)/(n+1)); x=sort(x)
  n1=round(l*n)
  if (minx) {
    if (add) points(c(0,q[n1:n]),c(0,x[n1:n]),...) else plot(c(0,q[n1:n]),c(0,x[n1:n]),...)
    lines(c(0,q[n1]),c(0,x[n1]),...)
  } else if (add) points(q,x,...) else plot(q,x,...)
  if (ab) abline(0,1,col="red")
}



raj.res <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_ca.rds')
knight.res <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_ca.rds')

##RAJ
par(mfrow = c(4, 4))
invisible(u <- split(raj.res, raj.res$disease) %>% lapply(., FUN =
  function(x) {
    split(x, x$cell) %>% lapply(., FUN = function(z) {
      # p <- p.adjust(z$p.val, method = 'BH')
      p <- z$p.val
      QQ(p, l = 0, main = paste(z$cell[1], z$disease[1]))
    }
      )
    }))

dev.print(pdf, file = 'pictures/ichip_ca_qq.pdf')

#-------------------
#permutation plots
raj.res.x <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_ichip_raj_ca_x.rds')
knight.res.x <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_ichip_knight_ca_x.rds')

par(mfrow = c(4, 4))
invisible(u <- split(raj.res.x, raj.res.x$disease) %>% lapply(., FUN =
  function(x) {
    split(x, x$cell) %>% lapply(., FUN = function(z) {
      # p <- p.adjust(z$p.val, method = 'BH')
      p <- z$p.val
      QQ(p, l = 0, main = paste(z$cell[1], z$disease[1]))
    }
      )
    }))

dev.print(pdf, file = 'pictures/ichip_ca_perm_qq.pdf')
#-------------------
##KNIGHT
DISEASES <- unique(knight.res$disease)
par(mfrow = c(3, 5))
use <- knight.res$disease %in% DISEASES[1 : 3]
invisible(u <- split(knight.res[use, ], knight.res$disease[use]) %>% lapply(., FUN =
  function(x) {
    split(x, x$cell) %>% lapply(., FUN = function(z) {
      # p <- p.adjust(z$p.val, method = 'BH')
      p <- z$p.val
      QQ(p, l = 0, main = paste(z$cell[1], z$disease[1]))
    }
      )
    }))
dev.print(pdf, file = 'pictures/ichip_knight_ca_qq.pdf')

par(mfrow = c(4, 5))
use <- knight.res$disease %in% DISEASES[4 : 7]
invisible(u <- split(knight.res[use, ], knight.res$disease[use]) %>% lapply(., FUN =
  function(x) {
    split(x, x$cell) %>% lapply(., FUN = function(z) {
      # p <- p.adjust(z$p.val, method = 'BH')
      p <- z$p.val
      QQ(p, l = 0, main = paste(z$cell[1], z$disease[1]))
    }
      )
    }))
dev.print(pdf, file = 'pictures/ichip_knight_ca_qq2.pdf')



#-------------------------------
##investigation yo
a <- 'raj'

pred <- readRDS(sprintf("/scratch/wallace/twas/model_output/chrm%s/gwas_mtl_%s_ichip_res.rds", k, a))

if(a == 'knight') {
  new.geno <- get(load(sprintf("/scratch/wallace/twas/gwas-data/ichip/knight-%s.RData", k)))
  } else {
    new.geno <- get(load(file = sprintf("/scratch/wallace/twas/gwas-data/ichip/raj-cd14-%s.RData", k)))
  }

mysamples <- samples(new.geno)

# if(!all.equal(rownames(mysamples), rownames(pred))) stop("mismatch between samples in preds and cc!")

phenotype <- mysamples[, 'phenotype']
aff <- mysamples[, 'affected']
country <- mysamples[, 'country']

DISEASES <- setdiff(unique(phenotype), c(NA,"CONTROL"))

#diff = case - control
out <- lapply(1 : ncol(pred), FUN = function(i){
    x <- pred[, i]
    res.aff <- do.ca(x, aff, country, ctrl = "1")

    res <- lapply(DISEASES, FUN = function(d) {
      message(d)
      if(d == 'RA'){
        use <- !is.na(x) & phenotype %in% c("CONTROL", d)
        var(x)
      } else {
        use <- !is.na(x) & phenotype %in% c("CONTROL", d) & country == 'UK'
        var(x)
      }
      res
    }
  ) %>% do.call(rbind, .)
    ans <- rbind(res.aff, res)
    ans$disease <- c('AUTO', DISEASES)

    z <- strsplit(colnames(pred)[i], "\\.|_")[[1]]
    ans$cell <- z[1]
    if(a == 'knight'){
      ans$probe <- paste(z[2], z[3], sep = '_')
      ans$chr <- z[4]
    } else {
        ans$probe <- z[2]
        ans$chr <- z[3]
  }

    message(colnames(pred)[i])
    ans
      }) %>% do.call(rbind, .)

rownames(out) <- NULL

message("chrm ", k, " done!")
out
  }, mc.cores = 11)
res <- do.call(rbind, res)

















##
