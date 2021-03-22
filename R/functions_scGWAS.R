

generate_covs_cells <- function(exp, out.dir, cors = 25){
  if(!dir.exists(out.dir)) dir.create(out.dir)
  genes.means =  rowMeans(exp[,-1])
  write_cov <- function(ItemNumber){
    df <- data.frame(GENE =  exp$gene, E = exp[,ItemNumber], A = genes.means)
    cell.name <- colnames(exp)[ItemNumber]
    cover_name <- paste0(out.dir, "/", cell.name, ".cov")
    write.table(x = df, file = cover_name, quote = F, col.names = T,
                row.names = F, sep = "\t")
  }
  
  parallel::mclapply(X = 2:ncol(exp), FUN = write_cov, mc.cores = cors)
  
}


######


gene_prop_cells <- function(cor = 10, raw.file.path, cover.files.dir.path,
                            magma.path, out.dir){
  cover.files.path = list.files(path = cover.files.dir.path, 
                                pattern = "[.]cov$")
  cell.n <- gsub(pattern = "[.]cov$", replacement = "", x = cover.files.path)
  cover.df <- data.frame(cell = cell.n, path = cover)
  mag_prop <- function(IterNumber){
    mag.com <- paste0(magma.path, "magma",
                      " --gene-results ", raw.file.path, 
                      " --gene-covar ",  cover.df$path[IterNumber], 
                      " --model condition-hide=A direction=greater --out ",
                      out.dir, "/", cover.df$cell[IterNumber])
    system(command =  mag.com, wait = T)
  }
  
  parallel::mclapply(X = 1:length(cover.files), FUN = mag_prop, 
                     mc.cores = cor)
  
}

#######
get_scores <- function(gsa.files.path, cors = 20){
    prop.dir <- list.dirs(path =gsa.files.path, recursive = F)
    
    tstat_calc <- function(x){
      val <- fread(x, data.table = F, skip = "VARIABLE", nrows = 1,
                   select = c(4,6))
      tstat <- val$BETA/val$SE
      tstat
    }
    cells <- list.files(path = prop.dir[i], pattern = "gsa[.]out$",
                        full.names = T)
    cells.n <- list.files(path = prop.dir[i], pattern = "gsa[.]out$",
                          full.names = F)
    cells.n <- gsub(pattern = "[.]gsa[.]out", replacement = "", x = cells.n )
      
    ts <- mclapply(X = cells, FUN = tstat_calc, mc.cores = cors)
    ts <- unlist(ts)
    ts <- data.frame(cells = cells.n, ts = ts)
    ts
    }
  
###### 

