
library(KEGGgraph)
library(igraph)
library(ggraph)


download.file("https://rest.kegg.jp/get/ko04612/kgml", "ko04612.xml") # Antigen Processing
download.file("https://rest.kegg.jp/get/ko04660/kgml", "ko04660.xml") # T-Cell Signaling
download.file("https://rest.kegg.jp/get/ko05164/kgml", "ko05164.xml") # Influenza A (Contains MX1)
download.file("https://rest.kegg.jp/get/ko04620/kgml", "ko04620.xml") # Toll-like Receptor (Innate Immune Response)

# 2 & 3 Combined: Parse and Clean in one step
clean_kgml <- function(file) {
  g <- parseKGML2Graph(file, expandGenes = TRUE, genesOnly = FALSE)
  ig <- graph_from_graphnel(g)
  # THE FIX: Remove nodes with empty names
  ig <- delete_vertices(ig, V(ig)[name == ""])
  return(ig)
}


ig1 <- clean_kgml("ko04612.xml")
ig2 <- clean_kgml("ko04660.xml")
ig3 <- clean_kgml("ko05164.xml")
ig4 <- clean_kgml("ko04620.xml")

# 4. Weld the 4 pathways
stitched_network <- igraph::union(ig1, ig2)
stitched_network <- igraph::union(stitched_network, ig3)
stitched_network <- igraph::union(stitched_network, ig4)

# 4. Weld the 4 pathways into your focused immune network
stitched_network <- igraph::union(ig1, ig2)
stitched_network <- igraph::union(stitched_network, ig3)
stitched_network <- igraph::union(stitched_network, ig4)

# 5. Clean up any resulting empty nodes from the union operation
stitched_network <- delete_vertices(stitched_network, V(stitched_network)[name == ""])

# ==================== THE FIX: CORRECT MX1 ID ====================

# FIXED: MX1 is K14754. 
target_KOs <- c("ko:K14754", "ko:K06751", "ko:K06752", "ko:K06451", "ko:K06458")

gene_dictionary <- c(
  "ko:K14754" = "MX1",
  "ko:K06751" = "MHC-I",
  "ko:K06752" = "MHC-II",
  "ko:K06451" = "CD3e",
  "ko:K06458" = "CD8a"
)

# ==================== THE "FUZZY SEARCH" LOGIC ====================

# Set default styling for the background nodes
V(stitched_network)$node_color <- "lightblue3"
V(stitched_network)$node_size  <- 3.5
V(stitched_network)$label_text <- gsub("ko:", "", V(stitched_network)$name)
V(stitched_network)$label_size <- 2.3
V(stitched_network)$label_color <- "gray40"
V(stitched_network)$is_target  <- FALSE 

# Loop through the targets and fish them out of any protein complexes
cat("\n--- SEARCHING FOR YOUR GENES ---\n")
for (ko in target_KOs) {
  matches <- grepl(ko, V(stitched_network)$name) 
  
  if (any(matches)) {
    cat("✅ Found:", gene_dictionary[ko], "\n")
    V(stitched_network)$node_color[matches] <- "firebrick1"
    V(stitched_network)$node_size[matches] <- 8
    V(stitched_network)$label_text[matches] <- gene_dictionary[ko]
    V(stitched_network)$label_size[matches] <- 4.5
    V(stitched_network)$label_color[matches] <- "firebrick4"
    V(stitched_network)$is_target[matches] <- TRUE
  } else {
    cat("❌ WARNING: Could not find", gene_dictionary[ko], "\n")
  }
}
cat("--------------------------------\n\n")

# ==================== EDGE HIGHLIGHTING ====================

target_vertices <- V(stitched_network)[V(stitched_network)$is_target]

# Identify lines connecting to or from the targets
incident_edges <- E(stitched_network)[.from(target_vertices) | .to(target_vertices)]

# Default all lines to faint gray
E(stitched_network)$edge_color <- "gray45"
E(stitched_network)$edge_alpha <- 0.3
E(stitched_network)$edge_width <- 0.7

# Force the links directly touching your 5 targets to glow bright red
E(stitched_network)$edge_color[incident_edges] <- "firebrick1"
E(stitched_network)$edge_alpha[incident_edges] <- 1.0
E(stitched_network)$edge_width[incident_edges] <- 1.6

# ==================== PLOTTING ====================

ggraph(stitched_network, layout = "nicely") + 
  geom_edge_link(aes(color = edge_color, alpha = edge_alpha, width = edge_width), 
                 arrow = arrow(length = unit(1.8, 'mm')), 
                 end_cap = circle(3, 'mm'), show.legend = FALSE) + 
  scale_edge_color_identity() +
  scale_edge_alpha_identity() +
  scale_edge_width_identity() +
  geom_node_point(aes(color = node_color, size = node_size)) +
  geom_node_text(aes(label = label_text, size = label_size, color = label_color), 
                 repel = TRUE, fontface = "bold", 
                 bg.color = "white", bg.r = 0.15, max.overlaps = Inf) +
  scale_color_identity() + 
  scale_size_identity() +
  theme_graph() +
  labs(title = "Targeted Anas platyrhynchos Immune Response Network",
       subtitle = "Different pathway network mappings with selected Innate and Adaptive biomarkers")

```