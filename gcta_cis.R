
help_message <- "This script constructs prediction models for gene expression based on SNPs in the vicinity of genes, depending on GCTA.
Usage: Rscript gcta_cis.R [options]

Options:
  --gffs_file=value   An RData file containing gene annotation information, including a data frame named 'gffs' with columns Gene, chr, start, end, strand, Qsymbols, and Note.
  --exp_file=value   An RData file containing a gene expression matrix with a matrix named 'exp_m', where row names represent gene IDs and column names represent individual IDs.
  --genodir=value   Directory containing genotype data in plink's binary format, including .bed, .bim, and .fam files, saved per chromosome.
  --gfile_prefix=value   Pattern for the file names of genotype data, where chromosome number is represented by %s, exclusive of suffix.
  --plinkdir=value   Path to the plink program.
  --gctadir=value   Path to the GCTA program.
  --out_dir=value   Output directory.
  --extend=value   The range of cis-regulatory regions around genes. SNPs within a distance less than 'extend' from genes are used to build gene expression prediction models. e.g., extend=100000.
  --ncor=value   The number of threads to use.
  --help   Display this help message.
"

args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args) {
  cat(help_message)
  quit("no", status = 0)
}

for (arg in args) {
  if (startsWith(arg, "--gffs_file=")) {
    gffs_file <- gsub("--gffs_file=", "", arg)
  }
  else if (startsWith(arg, "--exp_file=")) {
    exp_file <- gsub("--exp_file=", "", arg)
  }
  else if (startsWith(arg, "--genodir=")) {
    genodir <- gsub("--genodir=", "", arg)
  }
  else if (startsWith(arg, "--gfile_prefix=")) {
    gfile_prefix <- gsub("--gfile_prefix=", "", arg)
  }
  else if (startsWith(arg, "--plinkdir=")) {
    plinkdir <- gsub("--plinkdir=", "", arg)
  }
  else if (startsWith(arg, "--gctadir=")) {
    gctadir <- gsub("--gctadir=", "", arg)
  }
  else if (startsWith(arg, "--out_dir=")) {
    out_dir <- gsub("--out_dir=", "", arg)
  }
  else if (startsWith(arg, "--extend=")) {
    extend <- as.numeric(gsub("--extend=", "", arg))
  }
  else if (startsWith(arg, "--ncor=")) {
    ncor <- as.numeric(gsub("--ncor=", "", arg))
  }
}


gcta_gene_cis_pvalue <- function(counts_matrix, gffs, chr, extend=1e5, ncor=2, geno_dir, gfile_name, tmp_dir, plink_dir, gcta_dir, filter.maf=0.05){
    #library(parallel)
    library(pbmcapply)
    require(snpStats)
    fam_file <- sprintf("%s/%s.fam",geno_dir,gfile_name)
    vars.all <- read.table(fam_file,stringsAsFactors=F,header=F)[,1]
    #vars.all <- setdiff(vars.all, c("C126","W190","W196","W232"))
    vars.u <- intersect(colnames(counts_matrix), vars.all)
    gene_all <- intersect(rownames(counts_matrix), gffs$Gene[gffs$chr==chr])
    counts_m <- counts_matrix[gene_all,vars.u]
    bim_file <- sprintf("%s/%s.bim",geno_dir,gfile_name)
    snp_all <- read.table(bim_file,stringsAsFactors=F,header=F)
    bfile <- sprintf("%s/%s",geno_dir,gfile_name)
    
    func <-  function(gene){
       #print(gene)
       left <- gffs[match(gene,gffs$Gene),"start"] - extend; right <- gffs[match(gene,gffs$Gene),"end"] + extend
       snp_s <- snp_all[snp_all[,4] >= left & snp_all[,4] <= right,2]
       
       #f.plink <- tmp_f <- sprintf("%s/%s%s",tmp_dir,basename(tempfile()),Sys.getpid())
       f.plink <- tmp_f <- sprintf("%s/%s_%d_%sk", tmp_dir, gene, length(vars.all), extend/1000)
       snp_file <- sprintf("%s.snp",tmp_f)
       ind_file <- sprintf("%s.ind",tmp_f)
       #write(paste(vars.u,vars.u,sep="\t"),file=ind_file)
       write(paste(vars.all,vars.all,sep="\t"),file=ind_file)  ##############
       write(snp_s, file=snp_file)
       system(sprintf("%s --bfile %s --extract %s --keep %s --noweb --make-bed --maf %f --out %s > /dev/null 2>&1", plink_dir,bfile,snp_file,ind_file,filter.maf,f.plink))
       
       grm_file <- sprintf("%s.grm",tmp_f)
       system(sprintf("%s --bfile %s --autosome --maf 0.05 --make-grm-inbred --out %s > /dev/null 2>&1", gcta_dir, f.plink,grm_file))
       phe_file <- sprintf("%s.phen",tmp_f)
       write.table(cbind(vars.u,vars.u,counts_m[gene,vars.u]), file=phe_file, quote=F, row.names=F, col.names=F, sep="\t")
       reml_file <- sprintf("%s.reml",tmp_f)
       system(sprintf("%s --reml --reml-pred-rand --grm %s --pheno %s --out %s > /dev/null 2>&1", gcta_dir, grm_file, phe_file, reml_file))
       
       summ_file <- sprintf("%s.hsq", reml_file)
       blp_file <- sprintf("%s.indi.blp", reml_file)
       if(file.exists(summ_file) & file.exists(blp_file)){
          reml_res <- readLines(summ_file)
          Vg_line <- grep("V\\(G\\)/Vp", reml_res)
          Pval_line <- grep("Pval", reml_res)
          Vg <- as.numeric(strsplit(reml_res[Vg_line], split="\t")[[1]][2])
          pvalue <- as.numeric(strsplit(reml_res[Pval_line], split="\t")[[1]][2])
          vc <- data.frame(Vg=Vg, pvalue=pvalue, stringsAsFactors=F)
          
          blp_res <- read.table(blp_file, stringsAsFactors=F, header=F, row.names=1)
          cis_fit <- blp_res[vars.u,3]; names(cis_fit) <- vars.u
          trans_fit <- blp_res[vars.u,5]; names(trans_fit) <- vars.u  ############
          
          system(sprintf("%s --bfile %s --blup-snp %s --out %s > /dev/null 2>&1", gcta_dir,f.plink,blp_file, tmp_f))
          pred_file <- sprintf("%s.snp.blp", tmp_f)
          
          plink.data <- snpStats::read.plink(f.plink)
          geno <- plink.data$genotype@.Data; mode(geno) <- 'integer' ## codes:1,3 NA:0
          geno.o <- geno; codes=c('AA'=0,'AB'=NA,'BB'=1,'NA'=NA)
          geno[geno.o==0]<- codes['NA']; geno[geno.o==1]<- codes['AA']
          geno[geno.o==2]<- codes['AB']; geno[geno.o==3]<- codes['BB']
          g.maf <- colMeans(geno,na.rm=T); f.na <- is.na(geno); geno[f.na] <- matrix(g.maf,nrow=nrow(geno),ncol=ncol(geno),byrow=T)[f.na]
          #geno <- scale(geno)
          
          betas <- read.table(pred_file, stringsAsFactors=F, header=F, row.names=1)
          snps <- intersect(rownames(betas), colnames(geno))
          geno <- geno[ ,match(snps,colnames(geno))]
          betas <- betas[match(snps,rownames(betas)), ]
          maps <- plink.data$map[match(snps,plink.data$map$snp.name), ]
          ids <- which(betas[,1] != maps[,"allele.2"])
          betas[ids,2] <- -betas[ids,2]
          cis_pred <- rowSums(t(t(geno) * betas[,2]))[vars.all]
       }else{
          vc <- data.frame(Vg=NA, pvalue=NA, stringsAsFactors=F)
          cis_fit <- rep(NA, length(vars.u)); names(cis_fit) <- vars.u
          trans_fit <- rep(NA, length(vars.u)); names(trans_fit) <- vars.u
          cis_pred <- rep(NA, length(vars.all)); names(cis_pred) <- vars.all
       }
       
       system(sprintf("rm %s*",tmp_f))
       return(list(res_vc=vc, res_cis=cis_fit, res_trans=trans_fit, cis_pred=cis_pred))
    }
    
    #system.time(res <- mclapply(gene_all, func, mc.cores=getOption("mc.cores",ncor)))
    res <- pbmclapply(gene_all, func, mc.cores=ncor, ignore.interactive=TRUE)
    res_vc <- do.call(rbind, lapply(res,function(x) x$res_vc)); rownames(res_vc) <- gene_all
    res_cis <- t(do.call(cbind,lapply(res,function(x) x$res_cis))); rownames(res_cis) <- gene_all
    res_trans <- t(do.call(cbind,lapply(res,function(x) x$res_trans))); rownames(res_trans) <- gene_all
    res_pred <- t(do.call(cbind,lapply(res,function(x) x$cis_pred))); rownames(res_pred) <- gene_all
    return(list(res_vc=res_vc, res_cis=res_cis, res_trans=res_trans, res_pred=res_pred))
}


(load(gffs_file))
(load(exp_file))
chrs <- unique(gffs$chr[match(rownames(exp_m),gffs$Gene)])

for(chr in chrs){
   log_file <- sprintf("%s/gcta_cis_pvalue_%s.log", out_dir, chr)
   if(file.exists(log_file)){ print(sprintf("%s already exists!", log_file)); next }
   system(sprintf("touch %s", log_file))
   print(sprintf("Processing %s",chr))
   gfilename <- sprintf(gfile_prefix, chr)
   res <- gcta_gene_cis_pvalue(counts_matrix=exp_m, gffs=gffs, chr=chr, extend=extend, ncor=ncor, geno_dir=genodir, gfile_name=gfilename, tmp_dir=out_dir, plink_dir=plinkdir, gcta_dir=gctadir, filter.maf=0.05)
   save(res, file=sprintf("%s/gcta_cis_pred_%s.RData", out_dir, chr))
}



