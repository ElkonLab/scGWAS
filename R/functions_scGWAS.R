

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


df.genes <- function(exp = exp, meta = meta, model, branch = NULL, cors = 10){
  
  if(!is.null(branch)){
    ref <- levels(meta$branch)[levels(meta$branch) != branch]
    meta$branch <- relevel(meta$branch, ref = "")
  }
  
  var <- gsub(model, pattern = "~", replacement = "")
  var <- gsub(var, pattern = ":", replacement = " ")
  var <- gsub(var, pattern = "[+]", replacement = " ")
  var <- gsub(var, pattern = "  ", replacement = " ")
  var <- gsub(var, pattern = "  ", replacement = " ")
  var <- strsplit(var, split = " ")[[1]]
  var <- c("cells", var)
  met <- meta[,which(colnames(meta) %in% var)]
  
  model <- paste0("exp ", model)    
  model <- as.formula(model)
  df.gene <- function(genename, .exp = exp, .meta = met, .model = model){
    exp.gene <- exp[genename,]
    exp.gene <- data.frame(cells = names(exp.gene), exp = as.numeric(exp.gene))
    exp.gene <- merge(exp.gene, met, by = "cells")
    exp.gene <- es[complete.cases(exp.gene),]
    full_model_fit <- speedglm::speedglm(data = exp.gene , formula = model, 
                                         family = stats::gaussian(), acc = 0.001, model = FALSE, 
                                         y = FALSE)
    
    s <- summary(full_model_fit)
    df <- s$coefficients[4,]
    df$gene <- genename
    df <- df[,c(5,1:4)]
    colnames(df)[4] <- "t"
    df
    
  }
  
  plyr::rbind.fill(mclapply(FUN = df_branch, X = rownames(exp), 
                            mc.cores = 10, mc.preschedule = T))
  
}


##########


branch_assign <- function(traj, progenitors, branch1, branch2, b.names){
  traj <- dplyr::arrange(traj, Pseudotime) %>% 
    dplyr::mutate(rank = 1:nrow(traj))
  traj <- rbind.data.frame(traj[1,], traj)
  traj$rank[1] <- 0
  traj$cell.type <- ifelse(traj$State %in% progenitors, "Progenitor",
                           b.names[1])
  traj$cell.type <- ifelse(traj$State %in% branch2, b.names[2], 
                           traj$cell.type)
  
  traj$Branch <- traj$cell.type
  traj$Branch <- ifelse(test = ((traj$cell.type ==  "Progenitor") & 
                                  (traj$rank%%2 == 0)),
                        yes = b.names[1], no = traj$Branch)
  traj$Branch <- ifelse(test = (traj$cell.type ==  "Progenitor") & 
                          (traj$rank%%2 != 0), 
                        yes = b.names[2], 
                        no =  traj$Branch)
  traj$Branch <- factor(traj$Branch, levels = c(b.names[1],
                                                b.names[2]))
  traj$cell.type <- factor(traj$cell.type, levels = c("Progenitor", 
                                                      b.names[1],
                                                      b.names[2]))
  traj
}
