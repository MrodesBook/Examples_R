# APY Function

APY_inverse = function(GRM, corelist) {
  # Index core and non-core
  core = corelist$core_status == 1
  noncore = corelist$core_status == 0
  
  ncore = sum(core)
  nnoncore = sum(noncore)
  
  # Partition G to core and non-core
  Gcc = GRM[core, core]
  Gcn = GRM[core, noncore]
  Gnc = GRM[noncore, core]
  Gnn = GRM[noncore, noncore]
  
  # Gpart = rbind(cbind(Gcc, Gcn), cbind(Gnc, Gnn)) 
  
  # APY inverse based on above Gpart:
  Gcc_inv = solve(Gcc)
  
  Mnn_inv = matrix(data = 0, nrow = nnoncore, ncol = nnoncore)
  for(i in 1:nnoncore) 
  { 
    Mnn_inv[i,i] = 1 / (Gnn[i,i] - Gnc[i,] %*% Gcc_inv %*% Gcn[,i]) 
  } 
  
  APY11 = Gcc_inv + Gcc_inv %*% Gcn %*% Mnn_inv %*% Gnc %*% Gcc_inv 
  APY12 = -1 * Gcc_inv %*% Gcn %*% Mnn_inv 
  APY21 = -1 * Mnn_inv %*% Gnc %*% Gcc_inv 
  # APY21 = t(APY12)
  APY22 = Mnn_inv 
  
  Ginv_APY = rbind(cbind(APY11, APY12), 
                   cbind(APY21, APY22))
  
  # Create index of ID's
  # Note - order is not the same as in the original G
  apy_ref = tibble(geno_id = c(rownames(Gcc), rownames(Gnn)), apy_order = seq(1:nrow(Ginv_APY)))
  orig_ref = full_join(corelist, apy_ref, by = "geno_id")
  
  list(APY_Ginv = Ginv_APY, index_file = orig_ref)
}

