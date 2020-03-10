# connectivity matrix

library(stats)
library(spatstat)
library(natverse)
library(knitr)
library(plyr)
library(ggplot2)
library(igraph)
library(RColorBrewer)


setwd("./") # set to this code's directory

# clean everythign up.  
rm(list=ls())

#close any open rgl windows
while (rgl.cur() > 0) { rgl.close() }

# set up for 3d plots based on rgl package
rgl::setupKnitr()

load("LC6_Analysis_connMatrix_data.RData")


# Table S1 ,conn matrix -------------------------------------------------------------------------------------------

library(reshape2)
longData <- melt(conn_all)
longData<-longData[longData$value!=0,]
# longData<-melt(conn_all[1:69,1:69])
mat_names <- c(paste("LC6_", anno_rightLC6$skid, sep = ""), 
               paste("BiL_", biL_skid, sep = ""), 
               paste("BiR_", biR_skid, sep = ""), 
               paste("Ipsi_", ipsi_skid, sep = ""),
               "LC6_all", "Bi_LHS", "Bi_RHS", "Ipsi_all")
dev.new()
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  # geom_point(aes(size=value)) +
  scale_fill_gradient(low="grey90", high="red", trans = "log", breaks = c(1,10,100,1000), labels = c(1,10,100,1000)) +
  guides(fill = guide_colourbar(title = "No. connection")) +
  # guides(fill = guide_colourbar(label = c("1","10","100","1000"))) +
  # geom_text(aes(label = value)) +
  labs(x="postsynaptic neurons", y="presynaptic neurons", title="Connectivity Matrix") +
  scale_x_continuous(breaks = seq(1,conn_all_dim),labels = mat_names, position = "top", expand = c(0,0)) +
  scale_y_reverse(breaks = seq(1,conn_all_dim), labels = mat_names, expand = c(0,0)) +
  # theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="#993333", size=10, angle=90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(face="bold", color="#993333", size=10, angle=0),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black"))




# Figure 6B,conn matrix compact -----------------------------------------------------------------------------------

conn_cpt <- conn_all[(length(neu_all)+1):conn_all_dim, (length(neu_all)+1):conn_all_dim]
longData <- melt(conn_cpt)
# longData<-longData[longData$value!=0,]
mat_names <- c("LC6_all", "Bi_LHS", "Bi_RHS", "Ipsi_all")
dev.new()
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  # geom_point(aes(size=value)) +
  scale_fill_gradient(low="grey90", high="red", trans = "log", breaks = c(1,10,100,1000), labels = c(1,10,100,1000)) +
  guides(fill = guide_colourbar(title = "No. connection")) +
  geom_text(aes(label = value)) +
  labs(x="postsynaptic neurons", y="presynaptic neurons", title="Connectivity Matrix") +
  scale_x_continuous(breaks = seq(1,4),labels = mat_names, position = "top", expand = c(0,0)) +
  scale_y_reverse(breaks = seq(1,4), labels = mat_names, expand = c(0,0)) +
  # theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="#993333", size=10, angle=90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(face="bold", color="#993333", size=10, angle=0),
        panel.background = element_blank())


