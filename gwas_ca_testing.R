#-----------------------------------------------------------------------------------
#Cochran-Armitage 27.07.2017
library(annotSnpStats)
library(randomFunctions)
library(magrittr)
library(parallel)

setwd('/rds/user/ng414/hpc-work/twas')
#---------------------------------------------------------

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


do.ca <- function(x, cc, stratum = rep(1, length(cc)), ctrl = "0"){
    require(randomFunctions)

    if(var(x) == 0) return(data.frame(control = NA, case = NA, diff = NA, p.vals = NA, V = 0))

    means <- tapply(x, list(cc, stratum), mean)
    d <- setdiff(rownames(means), ctrl)
    m.expr <- rowMeans(means)
    if(length(unique(stratum)) == 2){
      diff <- vw(x, cc, stratum)[['D']]
    } else {
      diff <- mean(means[d, ] - means[ctrl, ])
    }
    ca.test.aff <- Cochran.Armitage.test(exposure = x, cc = cc, stratum = stratum)
    res <- data.frame(control = m.expr[ctrl], case = m.expr[d], diff = diff, p.vals = ca.test.aff[['p.value']], V = var(x), stringsAsFactors = FALSE)
    return(res)
}

#----------------------------------------------------------

a <- c('raj', 'knight')
dd <- c('jia', 'ms', 't1dgc+wtccc', 'ichip')

#exp_ = {mlt, lasso, stl_rf}
gwas.ca <- function(a, dd, exp_ ){
  cc.table <- c(phenotype = 'jia', case = 'ms', cc = 't1dgc+wtccc')
  cc.name <- names(cc.table)[match(dd, cc.table)]

  dd.table <- c(JIA = 'jia', MS = 'ms', T1D = 't1dgc+wtccc')
  d <- names(dd.table)[match(dd, dd.table)]

  res <- mclapply(1 : 22, FUN = function(k){
    if(exp_ == 'mtl'){
      pred <- readRDS(sprintf("model_output/chrm%s/gwas_mtl_%s_%s_res.rds", k, a, dd))
    } else {
      pred <- readRDS(sprintf("model_output/chrm%s/gwas_%s_%s_%s_res.rds", k, exp_, a, dd))
      if(a == 'raj') names(pred) <- c('raj-CD4', 'raj-CD14')
      names_ <- lapply(1 : length(pred), FUN = function(i) paste0(names(pred)[i], ".", colnames(pred[[i]]), "_", k)) %>% unlist(.)
      pred <- do.call(cbind, pred)
      colnames(pred) <- names_
    }

    if(a == 'knight') {
      if(dd == 'ms'){
      new.geno <- readRDS(file = sprintf("gwas-data/%s/pro-geno/gwas-%s.rds", dd, k))
      } else {
        new.geno <- readRDS(file = sprintf("gwas-data/%s/pro-geno/%s-%s.rds", dd, a, k))
      }
        }

    if(a == 'raj'){
      #wlog using raj-cd14 genomatrix to extract sample info (same as for cd4)
      new.geno <- readRDS(file = sprintf("gwas-data/%s/pro-geno/raj-cd14-%s.rds", dd, k))
    }

    mysamples <- samples(new.geno)

    # if(!all.equal(rownames(mysamples), rownames(pred))) stop("mismatch between samples in preds and cc!")

    cc <- mysamples[, cc.name]

    #diff = case - control
    res <- lapply(1 : ncol(pred), FUN = function(i){
        x <- pred[, i]
        if(dd == 't1dgc+wtccc'){
          stratum <- mysamples$cohort
        } else {
          stratum <- rep(1, nrow(mysamples))
        }

        res <- do.ca(x, cc, stratum, ctrl = "0")
        message(sprintf("%s %s: %s", a, d, colnames(pred)[i]))
        res
          }) %>% do.call(rbind, .)
    res$disease <- d
    z <- strsplit(colnames(pred), "\\.|_") %>% do.call(rbind, .)
    res$cell <- z[, 1]
    if(a == 'knight'){
      res$probe <- paste(z[, 2], z[, 3], sep = '_')
      res$chr <- z[, 4]
    } else {
        res$cell[res$cell == 'raj-cd14'] <- 'raj-CD14'
        res$cell[res$cell == 'raj-cd4'] <- 'raj-CD4'
        res$probe <- z[, 2]
        res$chr <- z[, 3]
    }

    message(sprintf("dataset %s disease %s: chrm %s done!", a, dd, k))
    message('')
    res
      }, mc.cores = 10)
    res <- do.call(rbind, res)
    return(res)
}


###----------------------------------------------------------
#------------------------------------------------------------
##ASSEMBLING

load("raj-probes.RData")
probes.r <- probes; rm(probes)
load("/scratch/wallace/twas/knight-probes.RData")
raj.genes <- readRDS("/scratch/wallace/twas/gene_details.rds")
colnames(raj.genes) <- c("probe", "ensid", "symbol", "chr")
load('model_output/hla_probes.RData') #hla.r,hla.k

knight.res <- lapply(c('t1dgc+wtccc', 'jia', 'ms'), FUN = function(d){
  message(d)
  out <- gwas.ca('knight', d, exp_ = 'mtl')
  message(d, ' done!')
  out
  })
names(knight.res) <- c('T1D', 'JIA', 'MS')
knight.res <- lapply(knight.res, FUN = function(x){
    x$gene <- probes[x$probe, "hgnc_symbol"]
    rownames(x) <- paste0(x$cell, ".", x$probe, "_", x$chr)
    x
  })

#excluding HLA probes
knight.res <- lapply(knight.res, FUN = function(z) subset(z, !(z$probe %in% hla.k)))
saveRDS(knight.res, file = 'model_output/gwas_knight_ca.rds')
#--------

raj.res <- lapply(c('t1dgc+wtccc', 'jia'), FUN = function(d){
  out <- gwas.ca('raj', d, exp_ = 'mtl')
  message(d, ' done!')
  out
  })
names(raj.res) <- c('T1D', 'JIA')
raj.res <- lapply(raj.res, FUN = function(z) subset(z, !(z$probe %in% hla.r)))
saveRDS(raj.res, file = '/scratch/wallace/twas/model_output/gwas_raj_ca0.rds')

raj.res <- lapply(raj.res, FUN = function(x){
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
raj.res <- lapply(raj.res, FUN = function(z) subset(z, !(z$probe %in% hla.r)))
saveRDS(raj.res, file = '/scratch/wallace/twas/model_output/gwas_raj_ca.rds')
#--------

##LASSO
knight.res.l <- lapply(c('t1dgc+wtccc', 'jia', 'ms'), FUN = function(d){
  message(d)
  out <- gwas.ca('knight', d, exp_ = 'lasso')
  message(d, ' done!')
  out
  })
names(knight.res.l) <- c('T1D', 'JIA', 'MS')
knight.res.l <- lapply(knight.res.l, FUN = function(x){
    x$gene <- probes[x$probe, "hgnc_symbol"]
    rownames(x) <- paste0(x$cell, ".", x$probe, "_", x$chr)
    x
  })
knight.res.l <- lapply(knight.res.l, FUN = function(z) subset(z, !(z$probe %in% hla.k)))
saveRDS(knight.res.l, file = '/scratch/wallace/twas/model_output/gwas_knight_lasso_ca.rds')

raj.res.l <- lapply(c('t1dgc+wtccc', 'jia'), FUN = function(d){
  out <- gwas.ca('raj', d, exp_ = 'lasso')
  message(d, ' done!')
  out
  })
names(raj.res.l) <- c('T1D', 'JIA')
raj.res.l <- lapply(raj.res.l, FUN = function(z) subset(z, !(z$probe %in% hla.r)))
saveRDS(raj.res.l, file = '/scratch/wallace/twas/model_output/gwas_raj_lasso_ca0.rds')
saveRDS(raj.res.l, file = '/scratch/ng414/twas/model_output/gwas_raj_lasso_ca0.rds')

raj.res.l <- lapply(raj.res.l, FUN = function(x){
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
saveRDS(raj.res.l, file = '/scratch/wallace/twas/model_output/gwas_raj_lasso_ca.rds')
saveRDS(raj.res.l, file = '/scratch/ng414/twas/model_output/gwas_raj_lasso_ca.rds')
#--------

##STL-RF
knight.res.rf <- lapply(c('t1dgc+wtccc', 'jia', 'ms'), FUN = function(d){
  message(d)
  out <- gwas.ca('knight', d, exp_ = 'stl_rf')
  message(d, ' done!')
  out
  })
names(knight.res.rf) <- c('T1D', 'JIA', 'MS')
knight.res.rf <- lapply(knight.res.rf, FUN = function(x){
    x$gene <- probes[x$probe, "hgnc_symbol"]
    rownames(x) <- paste0(x$cell, ".", x$probe, "_", x$chr)
    x
  })
knight.res.rf <- lapply(knight.res.rf, FUN = function(z) subset(z, !(z$probe %in% hla.k)))
saveRDS(knight.res.rf, file = '/scratch/wallace/twas/model_output/gwas_knight_stl_rf_ca.rds')

raj.res.rf <- lapply(c('t1dgc+wtccc', 'jia'), FUN = function(d){
  out <- gwas.ca('raj', d, exp_ = 'stl_rf')
  message(d, ' done!')
  out
  })
names(raj.res.rf) <- c('T1D', 'JIA')
raj.res.rf <- lapply(raj.res.rf, FUN = function(z) subset(z, !(z$probe %in% hla.r)))
saveRDS(raj.res.rf, file = '/scratch/ng414/twas/model_output/gwas_raj_stl_rf_ca0.rds')

raj.res.rf <- lapply(raj.res.rf, FUN = function(x){
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
saveRDS(raj.res.rf, file = '/scratch/ng414/twas/model_output/gwas_raj_stl_rf_ca.rds')


q('no')

#--------------------------------
#reading in and assembling results

ksum <- function(x, thr = 0.05) {
    y <- split(x, x$cell) %>% lapply(., function(z) { z$padj <- p.adjust(z$p.vals, method = "BH"); subset(z, padj < thr })
    y
}

##HLA ADJUSTMENT
load('/scratch/ng414/twas/model_output/hla_probes.RData')

knight.res <- readRDS('/scratch/wallace/twas/model_output/gwas_knight_ca.rds')
raj.res <- readRDS('/scratch/wallace/twas/model_output/gwas_raj_ca.rds')
knight.res.l <- readRDS('/scratch/wallace/twas/model_output/gwas_knight_lasso_ca.rds')
raj.res.l <- readRDS('/scratch/wallace/twas/model_output/gwas_raj_lasso_ca.rds')
knight.res.rf <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_knight_stl_rf_ca.rds')
raj.res.rf <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_raj_stl_rf_ca.rds')


knight.res <- lapply(knight.res, FUN = function(z) subset(z, !(probe %in% hla.k)))
knight.res.l <- lapply(knight.res.l, FUN = function(z) subset(z, !(probe %in% hla.k)))
knight.res.rf <- lapply(knight.res.rf, FUN = function(z) subset(z, !(probe %in% hla.k)))
raj.res <- lapply(raj.res, FUN = function(z) subset(z, !(probe %in% hla.r)))
raj.res.l <- lapply(raj.res.l, FUN = function(z) subset(z, !(probe %in% hla.r)))
raj.res.rf <- lapply(raj.res.rf, FUN = function(z) subset(z, !(probe %in% hla.r)))


saveRDS(knight.res, file = '/scratch/wallace/twas/model_output/gwas_knight_ca.rds')
saveRDS(raj.res, file = '/scratch/wallace/twas/model_output/gwas_raj_ca.rds')
saveRDS(knight.res.l, file = '/scratch/wallace/twas/model_output/gwas_knight_lasso_ca.rds')
saveRDS(raj.res.l, file = '/scratch/wallace/twas/model_output/gwas_raj_lasso_ca.rds')
saveRDS(knight.res.rf, file = '/scratch/wallace/twas/model_output/gwas_knight_stl_rf_ca.rds')
saveRDS(raj.res.rf, file = '/scratch/wallace/twas/model_output/gwas_raj_stl_rf_ca.rds')

saveRDS(knight.res, file = '/scratch/ng414/twas/model_output/gwas_knight_ca.rds')
saveRDS(raj.res, file = '/scratch/ng414/twas/model_output/gwas_raj_ca.rds')
saveRDS(knight.res.l, file = '/scratch/ng414/twas/model_output/gwas_knight_lasso_ca.rds')
saveRDS(raj.res.l, file = '/scratch/ng414/twas/model_output/gwas_raj_lasso_ca.rds')
saveRDS(knight.res.rf, file = '/scratch/ng414/twas/model_output/gwas_knight_stl_rf_ca.rds')
saveRDS(raj.res.rf, file = '/scratch/ng414/twas/model_output/gwas_raj_stl_rf_ca.rds')



ksum <- function(x, thr = 0.01) {
    y <- split(x,x$cell) %>% lapply(., function(z) { z$padj <- p.adjust(z$p.vals, method="BH", n = nrow(z)); z }) %>% do.call(rbind, .)
    subset(y, padj < thr)
}

knight.res <- lapply(knight.res, ksum, thr = 1) %>% do.call(rbind, .)
rownames(knight.res) <- with(knight.res, paste(probe, cell, disease, sep = "."))
#old (STL) results
kdata <- readRDS(file = '/scratch/wallace/twas/model_output/rf_knight_res.rds')
sum(kdata$padj < 0.05)
sum(knight.res$padj < 0.05)
sum(rownames(kdata)[kdata$padj < 0.05] %in% rownames(knight.res)[knight.res$padj < 0.05])
sum(rownames(knight.res)[knight.res$padj < 0.05] %in% rownames(kdata)[kdata$padj < 0.05])
ii <- rownames(kdata)[kdata$padj < 0.05][which(!(rownames(kdata)[kdata$padj < 0.05] %in% rownames(knight.res)[knight.res$padj < 0.05]))]

##
scp ng414@darwin:/scratch/wallace/twas/model_output/gwas_{knight,raj}_ca.rds model_output/
##

###analysis (09.08.2017)
#james' function 04.08.2017
library(magrittr)
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


# knight.res <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_knight_ca.rds')
# raj.res <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_raj_ca.rds')

raj.res <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_raj_lasso_ca.rds')
knight.res <- readRDS(file = '/scratch/wallace/twas/model_output/gwas_knight_lasso_ca.rds')

par(mfrow = c(3, 5))
invisible(u <- lapply(1 : 3, FUN = function(i){
  split(knight.res[[i]], knight.res[[i]]$cell) %>% lapply(., FUN = function(z){
      p <- z$p.vals
      QQ(p, l = 0, main = paste(z$cell[1], z$disease[1]))
    })
  }))

# dev.print(pdf, file = 'pictures/gwas_knight_ca_qq.pdf')

par(mfrow = c(2, 2))
invisible(u <- lapply(1 : 2, FUN = function(i){
      split(raj.res[[i]], raj.res[[i]]$cell) %>% lapply(., FUN = function(z){
          p <- z$p.vals
          QQ(p, l = 0, main = paste(z$cell[1], z$disease[1]))
        })
      }))
# dev.print(pdf, file = 'pictures/gwas_raj_ca_qq.pdf')













##
