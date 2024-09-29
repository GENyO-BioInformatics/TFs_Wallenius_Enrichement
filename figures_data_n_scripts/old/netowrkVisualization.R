# library
library(igraph)

mynetwork <- read.delim("../data/network.tsv")

# Create data
network <- graph_from_data_frame(mynetwork)
network$color <- ifelse(V(g)$gender == "Male", "lightblue", "pink")

# When ploting, we can use different layouts:
par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(network, layout=layout.fruchterman.reingold, main="fruchterman.reingold",
     edge.arrow.size=.5, vertex.color="grey", vertex.size=15, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.2)
dev.off()

# See the complete list with
# help(layout)

