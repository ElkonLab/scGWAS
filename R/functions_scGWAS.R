
conv_hs <- function(exp)
orth.path <- paste0("/specific/elkon/eldadshulman/", 
                    "GWAS_dif/annots/ortholog_one2one.RDS")
orth <- readRDS(orth.path)

exp <- merge(orth[,c(2,4)], exp, by.x = "Mouse.gene.name", by.y = "g")

exp <- exp[,-1]
colnames(exp)[1] <- "g"
exp <- exp[!duplicated(exp$g),]

# This function generate covrate files for MAGMA's gene property analysis.
# exp, is a data.fram of the count matrix, where the first column, "Gene" is the gene name.
# cors is the number of cores to use.
generate_covs_cells <- function(exp, cors = 25){
  cov.dir <- paste0("/cov_files")
  if(!dir.exists(cov.dir)) dir.create(cov.dir)
  genes.means =  rowMeans(exp[,-1])
  write_cov <- function(ItemNumber){
    df <- data.frame(GENE =  exp$gene, E = exp[,ItemNumber], A = genes.means)
    cell.name <- colnames(exp)[ItemNumber]
    cover_name <- paste0(cov.dir, "/", cell.name, ".cov")
    write.table(x = df, file = cover_name, quote = F, col.names = T, row.names = F, sep = "\t")
  }
  
  parallel::mclapply(X = 2:ncol(exp), FUN = write_cov, mc.cores = cors)
  
}


######


gene_prop_cells <- function(cor = 10, raw.file.path, cover.files.path, magma.path){
  
  cell.n <- gsub(pattern = "[.]cov$", replacement = "", x = cover.files.path)
  cover.df <- data.frame(cell = cell.n, path = cover)
  mag_prop <- function(IterNumber){
    mag.com <- paste0(magma.path, "magma",
                      " --gene-results ", raw.file.path, 
                      " --gene-covar ",  cover.df$path[IterNumber], 
                      " --model condition-hide=A direction=greater --out ", 
                      "/", cover.df$cell[IterNumber])
    system(command =  mag.com, wait = T)
  }
  
  parallel::mclapply(X = 1:length(cover.files), FUN = mag_prop, mc.cores = cor)
  
}

