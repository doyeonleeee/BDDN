library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(tidyr)
library(dplyr)

# data load
load("data/fitting.rda")
load("data/dn_full.rda")

get_offDiag = function(matrix) {
  matrix[upper.tri(matrix)]
}

get_protein_pairs = function(protein_names) {
  combn(protein_names, 2, FUN = function(x) paste(x[1], "-", x[2]))
}

protein_names = colnames(Y) 
protein_pairs = get_protein_pairs(protein_names)

time_info_hf = c("0hr", "0.083hr", "0.167hr", "0.25hr", "0.33hr", "0.5hr", "1hr", "1.5hr", "2hr",
                 "4hr", "8hr", "12hr", "24hr", "25hr", "48hr", "72hr", "96hr")
t = gsub("[hr,min]","", time_info_hf)
t = sort(as.numeric(t))

off_list = list()
for (i in 1:dim(dn)[1]) {
  mat = dn[i, , ]
  colnames(mat) = colnames(Y)
  
  off_diag = get_offDiag(mat)
  off_list[[i]] = data.frame(
    Time = t[i],
    Pair = protein_pairs,
    Value = off_diag
  )
}

off_diag_data = do.call(rbind, off_list)


############## Functional clustering ############## 
library(fda)
library(reshape2)
library(data.table)

dn_wide = reshape2::dcast(off_diag_data, Pair ~ Time, value.var = "Value")
dn_matrix = as.matrix(dn_wide[,-1])
rownames(dn_matrix) = dn_wide$Pair

### 1. Basis
knots = c(seq(0,max(t),1)) # Location of knots
n_knots   = length(knots) # Number of knots
n_order   = 3 # order of basis functions
n_basis   = length(knots) + n_order - 2;
basis = create.bspline.basis(rangeval = range(t), nbasis = n_basis)

# generate functional data object
fdn = Data2fd(argvals = t, y = t(dn_matrix), basisobj = basis) 


### 2. FPCA
pca_fd = pca.fd(fdn, nharm = 3)
pca_fd$harmonics
pca_fd$values
pca_fd$varprop
fpca_scores  = pca_fd$scores 


### 3. K-means clustering
derivatives = t(apply(dn_matrix, 1, diff)) 
early_slope = rowMeans(derivatives[, 1:10]) 
mid_slope = rowMeans(derivatives[, 11:13])
late_slope = rowMeans(derivatives[, 14:16])

fpca_scores_extended = cbind(fpca_scores, early_slope, mid_slope, late_slope)

set.seed(4)
kmeans_result2 = kmeans(fpca_scores_extended, centers = 5)

cluster_protein_list2 = list()
for (cc in sort(unique(kmeans_result2$cluster))) {
  cluster_data = which(kmeans_result2$cluster == cc)
  
  protein_in_cluster = rownames(dn_matrix)[cluster_data]
  cluster_protein_list2[[paste("Cluster", cc)]] = protein_in_cluster
}

fpca_scores_df2 = data.frame(PC1 = fpca_scores[, 1], 
                             PC2 = fpca_scores[, 2], 
                             cluster = as.factor(kmeans_result2$cluster))


cluster_color = list(
  Cluster = c(
    "1" = "#369fce",
    "2" = "#E78AC3",  
    "3" = "#6abf6a",
    "4" = "#9786fc",
    "5" = "#ADADAD"
  ),
  Cluster_Name = c(
    "Negative Chage (Late Peak)" = "#369fce",
    "Positive Chage (Late Peak)" = "#E78AC3",  
    "Negative Chage (Early Peak)" = "#6abf6a",
    "Positive Chage (Early Peak)" = "#9786fc",
    "Small Change" = "#ADADAD"
  )
)

cluster_names = c("Negative Chage (Late Peak)", 
                  "Positive Chage (Late Peak)", 
                  "Negative Chage (Early Peak)", 
                  "Positive Chage (Early Peak)", 
                  "Small Change")
names(cluster_names) = 1:5

# Adding cluster names to the data
fpca_scores_df2$Cluster_Name = factor(fpca_scores_df2$cluster, levels = 1:5, labels = cluster_names)

# cluster 1
cluster1_data = which(kmeans_result2$cluster == 1)
dn_long1 = as.data.frame(t(dn_matrix[cluster1_data, ])) %>%
  mutate(Time = time_info_hf) %>%
  pivot_longer(cols = -Time, names_to = "Protein", values_to = "DN_Value")
dn_long1$Time = as.numeric(gsub("[hr,min]","", dn_long1$Time))


# cluster 2
cluster2_data = which(kmeans_result2$cluster == 2)
dn_long2 = as.data.frame(t(dn_matrix[cluster2_data, ])) %>%
  mutate(Time = time_info_hf) %>%
  pivot_longer(cols = -Time, names_to = "Protein", values_to = "DN_Value")
dn_long2$Time = as.numeric(gsub("[hr,min]","", dn_long2$Time))

# cluster 3
cluster3_data = which(kmeans_result2$cluster == 3)
dn_long3 = as.data.frame(t(dn_matrix[cluster3_data, ])) %>%
  mutate(Time = time_info_hf) %>%
  pivot_longer(cols = -Time, names_to = "Protein", values_to = "DN_Value")
dn_long3$Time = as.numeric(gsub("[hr,min]","", dn_long3$Time))

# cluster 4
cluster4_data = which(kmeans_result2$cluster == 4)
dn_long4 = as.data.frame(t(dn_matrix[cluster4_data, ])) %>%
  mutate(Time = time_info_hf) %>%
  pivot_longer(cols = -Time, names_to = "Protein", values_to = "DN_Value")
dn_long4$Time = as.numeric(gsub("[hr,min]","", dn_long4$Time))




#################### Dynamic differential network #################### 
B.psamp2=fit_mcmc_cpp2$B.psamp
A.psamp2=fit_mcmc_cpp2$A.psamp
B.psamp2 = B.psamp2[,,101:2000] 
A.psamp2 = A.psamp2[,,101:2000]

p = dim(A.psamp2)[1]
q = dim(B.psamp2)[2]
niter = dim(A.psamp2)[3]

### igraph 
library(igraph)

cluster_graph_list1 = list()
edge_count_list1 = list()
vertex_count_list1 = list()
protein_names = colnames(Y)

for (i in 1:dim(dn)[1]) {
  A = dn[i, , ]
  colnames(A) = protein_names
  rownames(A) = protein_names
  
  cluster_graph_list1[[i]] = list()
  edge_count_list1[[i]] = list()
  vertex_count_list1[[i]] = list()
  
  cluster_graph_list2[[i]] = list()
  edge_count_list2[[i]] = list()
  vertex_count_list2[[i]] = list()
  
  for (cluster_idx in seq_along(cluster_protein_list2)) {
    
    protein_pairs = cluster_protein_list2[[cluster_idx]]
    proteins_in_cluster = unique(unlist(strsplit(protein_pairs, " - ")))
    
    sub_A1 = matrix(0, nrow = nrow(A), ncol = ncol(A))
    rownames(sub_A1) = rownames(A)
    colnames(sub_A1) = colnames(A)
    
    for (pair in protein_pairs) {
      proteinA = strsplit(pair, " - ")[[1]][1]
      proteinB = strsplit(pair, " - ")[[1]][2]
      if (A[proteinA, proteinB] != 0) {
        sub_A1[proteinA, proteinB] = A[proteinA, proteinB]
        sub_A1[proteinB, proteinA] = A[proteinB, proteinA]
      }
    }
    
    G1 = graph_from_adjacency_matrix(sub_A1, mode = "undirected", weighted = TRUE, diag = FALSE)
    V(G1)$name = colnames(sub_A1)
    #E(G1)$color = ifelse(E(G1)$weight > 0, "blue", "red")
    E(G1)$color = cluster_color$Cluster[as.character(cluster_idx)]
    E(G1)$abs_weight = abs(E(G1)$weight)
    E(G1)$linetype = ifelse(E(G1)$weight > 0, "solid", "dashed")
    
    G1_layout = G1
    E(G1_layout)$weight = abs(E(G1_layout)$weight)
    
    cluster_graph_list1[[i]][[cluster_idx]] = G1
    edge_count_list1[[i]][[cluster_idx]] = ecount(G1)
    vertex_count_list1[[i]][[cluster_idx]] = vcount(G1)
    
  }
}

### ggraph
library(ggraph)
library(tidygraph)

cluster_color = list(
  Cluster = c(
    "1" = "#369fce",  
    "2" = "#E78AC3",  
    "3" = "#6abf6a",  
    "4" = "#9786fc",  
    "5" = "white"
  )
)

ggraph_time_list1 = list() 

for (i in 1:dim(dn)[1]) { 
  A = dn[i, , ]
  colnames(A) = protein_names
  rownames(A) = protein_names
  
  G_time = graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  V(G_time)$name = colnames(A)
  
  E(G_time)$color = "white"
  for (cluster_idx in seq_along(cluster_protein_list2)) {
    protein_pairs = cluster_protein_list2[[cluster_idx]]
    
    for (pair in protein_pairs) {
      proteinA = strsplit(pair, " - ")[[1]][1]
      proteinB = strsplit(pair, " - ")[[1]][2]
      
      edge_id = get.edge.ids(G_time, c(proteinA, proteinB))
      if (edge_id > 0) {
        E(G_time)$color[edge_id] = cluster_color$Cluster[as.character(cluster_idx)]
      }
    }
  }
  
  edges_to_remove = which(E(G_time)$color == cluster_color$Cluster["5"])  
  if (length(edges_to_remove) > 0) {
    G_time = delete_edges(G_time, edges_to_remove)
  }
  
  
  G_time_tbl = as_tbl_graph(G_time) %>%
    activate(edges) %>%
    mutate(edge_color = E(G_time)$color)
  
  G_time_tbl = G_time_tbl %>%
    activate(nodes) %>%
    mutate(row_num = row_number())
  
  
  g_plot1 = ggraph(G_time_tbl, layout = "linear", circular = TRUE) +
    geom_edge_link(aes(color = edge_color, width = abs(weight/2), linetype = ifelse(weight > 0, "solid", "dashed")), show.legend = FALSE) +
    scale_edge_color_identity() + # use edge_color
    scale_edge_width(range = c(1.3, 1.7)) +  # width range
    geom_node_point(size = 3, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    theme_void() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  

  ggraph_time_list1[[i]] = g_plot1
}


for (i in 1:length(ggraph_time_list1)) {  
  print(ggraph_time_list1[[i]])
}