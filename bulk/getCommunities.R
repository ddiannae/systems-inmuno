library(igraph)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Three arguments must be provided", call.=FALSE)
  
}

VERT <- args[1]
INTER <-  args[2]
OUT <- args[3]

getComInfo <- function(cmembership, network){
  comp.info <- lapply(unique(cmembership), function(idc){
    mem <- names(cmembership[cmembership == idc])
    com <- induced.subgraph(network, mem)
    intras <- igraph::as_data_frame(com, what = "edges") %>% 
      filter(interaction_type == "Intra") %>% nrow()
    prs <- page.rank(com)
    chrs <- table(V(com)$chr)
    edges <- length(E(com))
    return(data.frame(com_id = idc,
                      pg_gene = names(which.max(prs$vector))[1], 
                      chr = names(which.max(chrs))[1], order = length(V(com)), 
                      size = edges,
                      intra_fraction = intras/edges))
  })
  comp.info <- plyr::compact(comp.info)
  comp.info <- plyr::ldply(comp.info)
  comp.info <- comp.info %>% arrange(desc(size))
  return(comp.info)
}

interactions <- read_tsv(INTER)  %>% 
  dplyr::rename("from" = "source", "to" = "target",
                "weight" = "mi")
vertices <- read_tsv(VERT) %>% 
  dplyr::rename("name" = "ensembl")
vertices <- vertices[!duplicated(vertices$name), ]

net <- graph_from_data_frame(interactions, 
                             directed=F, vertices = vertices)
comm <- cluster_louvain(graph = net)
cat("Modularity: ", modularity(comm), "\n")
cat("Communities: ", length(comm))

names(comm$membership) <- comm$names
df_comm <- tibble(comm$names, comm$membership)
colnames(df_comm) <- c("ensembl", "community")

comm_info <- getComInfo(comm$membership, net)

cat("Saving files\n")
write_tsv(df_comm, paste0(OUT, "_membership.tsv"))
write_tsv(comm_info, paste0(OUT, "_comminfo.tsv"))
