# codes used for analyze LC6 and downstream neurons
# !!! open with encoding UTF-8
# All data is saved in the .RData file
# For a paritcular figure, search for, eg. "5A" (main figures) or "S5A" (supplimentary)


# load libraries 
library(natverse)
library(rjson)
library(alphashape3d)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(sf)
library(ggExtra)
library(cowplot)
library(alphahull)
library(reshape2)

setwd("C:/Users/zhaoa/Dropbox (HHMI)/LC6 downstream paper new/Code")
setwd("./") # set to this code's directory
# source("someFunc.R") # load some useful functions

# clean everythign up.  
rm(list=ls())

#close any open rgl windows
while (rgl.cur() > 0) { rgl.close() }

# set up for 3d plots based on rgl package
rgl::setupKnitr()



# Buchner 71 eye map from Andrew Straw----------------------------------------------------------------------------------------------

Npt <- 699
buchner <- read.csv("buchner71_tp.csv", header = FALSE)
buchner <- buchner[1:Npt,]
buchner <- buchner / pi * 180

dev.new()
plot(buchner)

range(buchner[buchner[,2] > 70 & buchner[,2] < 110, 1]) # [-7, 160] as horizontal range
# range(buchner[buchner,2]) 

buchner_phi <- c(-10, 160) # for longitude


# func to generate polygon from set of points ---------------------------------------------------------------------

mkpoly <- function(xy) {
  xy <- as.data.frame(xy)[,1:6]
  xyset <- list() # vertices for a list of polygons
  if (dim(xy)[1] > 2) {
    N <- 1
    xyset[[N]] <- xy[1,]
    xy <- xy[-1, ]
    ii <- c()
    while (dim(xy)[1] >= 1) {
      ii[1] <- match(tail(xyset[[N]], 1)[2], xy[, 1])
      ii[2] <- match(tail(xyset[[N]], 1)[2], xy[, 2])
      if (!is.na(ii[1])) {
        xyset[[N]] <- rbind(xyset[[N]], xy[ii[1], ])
        xy <- xy[-ii[1], ]
      } else if (!is.na(ii[2])){
        xytmp <- xy[ii[2], c(2,1,5,6,3,4)]
        colnames(xytmp) <- colnames(xyset[[N]])
        xyset[[N]] <- rbind(xyset[[N]], xytmp)
        xy <- xy[-ii[2], ]
      } else {
        N <- N + 1
        xyset[[N]] <- xy[1, ]
        xy <- xy[-1, ]
      }
    }
  }
  return(xyset)
}

# load neurons ----------------------------------------------------------------------------------------------------

# LC6 neurons
anno_LC6 <- catmaid_query_by_annotation("^LC6 neuron$")
anno_RHS <- catmaid_query_by_annotation("^LO-R")
anno_rightLC6 <- anno_LC6[anno_LC6[,1] %in% anno_RHS[,1],]
LC6_skid <- anno_rightLC6[,"skid"]
neu <-  read.neurons.catmaid(LC6_skid, .progress='text')
LC6 <- neu

# target neurons
anno_ipsi <- catmaid_query_by_annotation("^putative ss2036$")
ipsi_skid <- anno_ipsi[,"skid"]
neu_ipsi <-  read.neurons.catmaid(ipsi_skid, .progress='text') 

# anno_biL <- catmaid_query_by_annotation("^LC6 target - bilateral$")
# biL_skid <- anno_biL$skid # LHS x4 , 
# biL_skid <- biL_skid[1:2] #select 2
biL_skid <- c(3149978, 3262999)
neu_biL <- read.neurons.catmaid(biL_skid, .progress='text') 

biR_skid <- c(3154013, 3155729) #RHS x2
neu_biR <- read.neurons.catmaid(biR_skid, .progress='text')

neu_target <- c(neu_biL, neu_biR, neu_ipsi)

# TM5
neu_JSON <- fromJSON(file = "Tm5_LC6 mapping.json")
neu_skid <- c()
for (j in 1:length(neu_JSON)) {
  neu_skid[j] <- neu_JSON[[j]]$skeleton_id
}
TM5 = read.neurons.catmaid(neu_skid, .progress='text')

# load glomerulus volumne mesh 
glo_vol <- catmaid_get_volume("v14.LC6_glomerulus_preSyn_R")

# whole brain mesh
v14 <- catmaid_get_volume(439, rval = 'mesh3d')  

# define a plane separating lobula portion of the LC6
a = -0.6; b = 1; c = 1.3; d = -230000

# make a layer using LC6 ----------------------------------------------------------------------------------------

# get all end points and branch points
for (j in 1:length(neu)) {
  tar <- neu[[j]]
  xyz_ep <-  tar$d[tar$EndPoints, ] %>% xyzmatrix()
  xyz_bp = tar$d[tar$BranchPoints, ] %>% xyzmatrix()
  xyz_LO <- rbind(xyz_ep, xyz_bp) %>% 
    as_tibble() %>% 
    mutate(LO = a*X + b*Y + c*Z + d) %>%
    filter(LO > 0) %>%
    select(X,Y,Z)
  if (j == 1) {
    xyz_node <- xyz_LO
  } else {
    xyz_node <- bind_rows(xyz_node, xyz_LO)
  }
}

# fit curved surface
polyfitorder <- 2 # 2nd order surface
gridMargin <- 20000 #add margin on the edge
xyz_node <- data.matrix(xyz_node)
X <- xyz_node[, 1];  Y <- xyz_node[, 2]; Z <- xyz_node[, 3]
fitlm <- lm(Z ~ poly(X, Y, degree = polyfitorder, raw = TRUE)) #linear model fit

# make an underlying grid to interpret the fit as a point set
dx2 <- 500 
dy2 <- 500
xx <- seq(range(X)[1] - gridMargin, range(X)[2] + gridMargin, by = dx2)
yy <- seq(range(Y)[1] - gridMargin, range(Y)[2] + gridMargin, by = dy2)
xygrid <- expand.grid(xx, yy)
xygrid <- setNames(data.frame(xygrid), c('X', 'Y'));
valfit <- predict(fitlm, xygrid) #generate values from the fit
xyz_lm <- cbind(xygrid, valfit)
dist_min <- apply(xyz_lm, 1, function(pt) {min(rowSums(sweep(xyz_node, 2, pt) ^ 2))}) #calculate min distance
ii <- dist_min > 20e6 # old, for visulization
xyz_layer_25e6 <- xyz_lm[!ii,] # pts
valfit_sup <- valfit
valfit_sup[ii] <- NA
m_surf <-  matrix(valfit_sup, nrow = length(xx), ncol = length(yy)) # surface

ii <- dist_min > 5e6
xyz_layer <- xyz_lm[!ii,] # pts

# coarse version for plotting
dx2c <- 1000
dy2c <- 1000
xxc <- seq(range(X)[1] - gridMargin, range(X)[2] + gridMargin, by = dx2c)
yyc <- seq(range(Y)[1] - gridMargin, range(Y)[2] + gridMargin, by = dy2c)
xygridc <- expand.grid(xxc, yyc)
xygridc <- setNames(data.frame(xygridc), c('X', 'Y'));
# valfitc <- predict(fitlm, xygridc) #generate values from the fit
# xyz_lmc <- cbind(xygridc, valfitc)
dist_min <- apply(xygridc, 1, function(pt) {min(rowSums(sweep(xyz_node[,c(1,2)], 2, pt) ^ 2))}) #calculate min distance
ii <- dist_min > 20e6
xy_layer_coarse <- xygridc[!ii,] # pts

# make alpha mesh
msh.a <- ashape3d(xyz_node, alpha = 20000) # 20000 look ok
msh <- as.mesh3d(msh.a)
neu_lo <- nlapply(neu, subset, function(x) pointsinside(x, msh))

# Figure 5B
nopen3d()
par3d('windowRect' = c(100,100,1700,1700))
plot3d(neu[[10]],  col= "#d7191c", lwd = 5, soma=T, WithNodes = F)
plot3d(neu[[51]],  col= "#2c7bb6", lwd = 5, soma=T, WithNodes = F)
plot3d(neu,  col= 'grey90', soma=T, WithNodes = F)
shade3d(glo_vol, col= "cyan", alpha = 1)
surface3d(xx,yy,m_surf, color = "#fdae61", alpha = 1)
plot3d(TM5[[1]], col = 'gold4', lwd = 5) #TM5
plot3d(TM5[[2]], col = 'brown', lwd = 5)
rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0))
rgl.pop(type = "light")
rgl.light(theta = 30, phi = 30)

# # save
# rgl.snapshot(filename = "LC6_3d.png",fmt = "png")


# 2d projections --------------------------------------------------------------------------------------------------
# projection dendrites onto the layer, both contour and center-of-mass

# xyz_layer <- xyz_layer
row.names(xyz_layer) <-  seq(1, dim(xyz_layer)[1])

ind_pj <- list() #index of projected grid points
xyz_com <- list() # center-of-mass
xyz_pj_com <- list() # xyz of com projecting on grid
for (j in 1:length(neu)){
  tar <- neu[[j]]
  xyz_ep <-  tar$d[tar$EndPoints, ] %>% 
    xyzmatrix() %>%
    as_tibble() %>% 
    mutate(LO = a*X + b*Y + c*Z + d) %>%
    filter(LO > 0) %>%
    select(X,Y,Z) %>%
    data.matrix()
  xyz_bp <-  tar$d[tar$BranchPoints, ] %>% 
    xyzmatrix() %>%
    as_tibble() %>% 
    mutate(LO = a*X + b*Y + c*Z + d) %>%
    filter(LO > 0) %>%
    select(X,Y,Z) %>%
    data.matrix()
  
  # center-of-mass
  xyz_eb <- rbind(xyz_ep,xyz_bp)
  xyz_com[[j]] <- colSums(xyz_eb)/dim(xyz_eb)[1]
  
  # project dendrite end points and com to the fitted grid by shortest distance
  xyz_dend_pj <- rbind(xyz_com[[j]], xyz_ep) #append com at the beginning
  Nsigma <- 5 #exclude > 5 sigma
  com_dist <- as.matrix(dist(xyz_dend_pj))
  thrhd_dist <- sd(com_dist[,1])*Nsigma
  xyz_dend_pj <- cbind(xyz_dend_pj, com_dist[,1]) %>%
    as_tibble() %>%
    filter(V4 < thrhd_dist)%>%
    select(X,Y,Z) %>%
    data.matrix()
  
  ind_min <- c()
  ind_min <-apply(xyz_dend_pj, 1, function(pt) {which.min(rowSums(sweep(xyz_layer, 2, pt) ^ 2))}) # index in xyz_layer with min distance
  ind_com <- ind_min[1] #index of grid point that's closest to com
  ind_min <- unique(ind_min[-1]) #index of grid points
  xyz_pj_com[[j]] <- xyz_layer[ind_com,]
  ind_pj[[j]] <- row.names(xyz_layer[ind_min,])
}

# projection down onto (x,y) plane
xy_pj <- list()
for (j in 1:length(ind_pj)){
  xy_tmp <- list()
  xy_ashape <- list()
  ii <- sort(as.integer(ind_pj[[j]]))
  xy_tmp <- xyz_layer[ii, 1:2]
  xy_ashape <- ashape(xy_tmp + matrix(runif(dim(xy_tmp)[1]*2, 1e-9, 2e-9), ncol = 2), alpha = 6000)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_pj[[j]] <- list(xy=xy_tmp, ashape=xy_ashape, edge=xy_edge)
}


# ahull
xy_ashape_grid <- ashape(xyz_layer[,1:2] + matrix(runif(dim(xyz_layer)[1]*2, 1e-9, 2e-9), ncol = 2), alpha = 6000)


# # use chull for the grid projection
# xy_grid <- xy_ashape_grid$edges[,c("x1","y1")]
# hpts <- chull(xy_grid)
# hpts <- c(hpts, hpts[1])
# xy_edge_grid <- xy_grid[hpts,] # hull edge points

# ### ### use alpha-hull for the grid projection
xy_grid_ahull <- mkpoly(xy_ashape_grid$edges)[[1]][,3:4]
xy_edge_grid <- xy_grid_ahull # hull edge points

# get edge points for each LC6 projection
xy_poly <- list()
for (j in 1:length(ind_pj)) {
  ls_poly <- mkpoly(xy_pj[[j]]$ashape$edges)
  xy_poly[[j]] <- ls_poly[[1]][,3:4]
}


# use TM5 to determine the center and meridian  -------------------------------------------------------------------

tar <- TM5[[1]]
ng <- as.ngraph(tar)
distal_points <- igraph::graph.dfs(ng, root=2137, unreachable=FALSE, neimode='out')$order
proximal_points <- setdiff(igraph::V(ng), distal_points)
# points3d(xyzmatrix(tar$d[proximal_points,]), col = 'red', size = 5)
TM5_u <- colMeans(xyzmatrix(tar$d[proximal_points,])) #upper pt

tar <- TM5[[2]]
ng <- as.ngraph(tar)
distal_points <- igraph::graph.dfs(ng, root=873, unreachable=FALSE, neimode='out')$order
proximal_points <- setdiff(igraph::V(ng), distal_points)
TM5_c <- colMeans(xyzmatrix(tar$d[proximal_points,])) # central pt

grid_c <- which.min(rowSums((sweep(xyz_layer, 2, c(TM5_c), "-"))^2))
grid_u <- which.min(rowSums((sweep(xyz_layer, 2, c(TM5_u), "-"))^2))


# --- map equator and coord system 
ymax <- max(xy_edge_grid[,2])
ymin <- min(xy_edge_grid[,2])
center_new <- xyz_layer[grid_c,1:2]
x_med_new <- as.numeric(center_new[1])
y_eq_new <- as.numeric(center_new[2])
angR_new <- acos((x_med_new - xyz_layer[grid_u,1])/sqrt(sum((xyz_layer[grid_c,1:2] - xyz_layer[grid_u,1:2])^2)))

# PLOT,PAPER, 2d projection with dendrites
ang_2 <- pi/2 - angR_new
rot_2 <- matrix(c(cos(ang_2), sin(ang_2), -sin(ang_2), cos(ang_2)), ncol = 2)
xy_layer_coarse_rot <- sweep(xy_layer_coarse, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
xy_layer_coarse_rot <- sweep(t(rot_2 %*% t(xy_layer_coarse_rot)), 2, c(x_med_new, y_eq_new), '+')

# boundary
xy_edge_grid_rot <- sweep(xy_edge_grid, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
xy_edge_grid_rot <- sweep(t(rot_2 %*% t(xy_edge_grid_rot)), 2, c(x_med_new, y_eq_new), '+')

# Figure 5C
windows(record = F, width = 8, height = 8)
# pdf('LC6_2d_lo.pdf', width = 8, height = 8, family = "Courier")
plot(xy_layer_coarse_rot, col="#fdae61", cex = 1, pch = ".", ylim = rev(range(xy_layer_coarse_rot[,2])), asp = 1)
for (j in 1:length(neu_lo)) {
  neu_lo_xy <- neu_lo[[j]]$d[,c("X","Y")]
  pp <- sweep(neu_lo_xy, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
  pp <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
  # points(neu_lo_xy, col = "grey80", pch = ".", cex = 1)
  points(pp, col = "grey80", pch = ".", cex = 1)
  points(colMeans(pp)[1], colMeans(pp)[2], pch = 20, col = "blue", cex = 1.5)
}
ng=as.ngraph(neu[[10]])
distal_points=igraph::graph.dfs(ng, root=449, unreachable=FALSE, neimode='out')$order
distal_tree=subset(neu[[10]], distal_points)
pp <- sweep(distal_tree$d[,c("X","Y")], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
distal_tree$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
plot(distal_tree,  col= "#d7191c", lwd = 2, soma=T, WithNodes = F, add = T)
ng=as.ngraph(neu[[51]])
distal_points=igraph::graph.dfs(ng, root=636, unreachable=FALSE, neimode='out')$order
distal_tree=subset(neu[[51]], distal_points)
pp <- sweep(distal_tree$d[,c("X","Y")], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
distal_tree$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
plot(distal_tree,  col= "#2c7bb6", lwd = 2, soma=T, WithNodes = F, add = T)
# lines(rbind(xyz_layer[grid_c,1:2], xyz_layer[grid_c,1:2]+(xyz_layer[grid_u,1:2]-xyz_layer[grid_c,1:2])*2), lwd = 3, col = 'cyan')
pp <- sweep(xyz_layer[grid_u,1:2], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
pp <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
lines(rbind(c(x_med_new,y_eq_new), c(x_med_new,y_eq_new)+(pp-c(x_med_new,y_eq_new))*2), lwd = 3, col = 'cyan')
points(matrix(c(x_med_new,y_eq_new), ncol =2), pch = 18, col = 'brown', cex = 2)
points(c(x_med_new,y_eq_new)+(pp-c(x_med_new,y_eq_new)), pch = 18, col = 'gold4', cex = 2)
lines(c(300000,310000), c(300000, 300000), col = 'black', lwd = 3)
text(x = 305000, 295000, labels = "10 µm")
# polygon(xy_edge_grid_rot)


# wrap the lobula layer onto a hemisphere, use line-polygon boundary to calculate radial distance --------------

# get edge points for each LC6 projection
xy_poly <- list()
for (j in 1:length(ind_pj)) {
  ls_poly <- mkpoly(xy_pj[[j]]$ashape$edges)
  xy_poly[[j]] <- ls_poly[[1]][,3:4]
}

grid_bdpt <- xy_edge_grid
grid_bdpt <- rbind(grid_bdpt, grid_bdpt[1,])

library(sf)

poly_st <- st_polygon(list(data.matrix(grid_bdpt)))
xy_ori <- c(x_med_new, y_eq_new)
R <- (ymax-ymin)*2
line_st <- st_linestring(t(matrix(c(x_med_new, y_eq_new, x_med_new+100000, y_eq_new+100000), nrow = 2)))
int_st = st_intersection(line_st, poly_st) # the 2nd point is the intersection
xy_com <- list()
for (j in 1:length(xy_poly)) {
  xy_pj_com <- xyz_pj_com[[j]][c('X','Y')]
  colnames(xy_pj_com) <- c("x1", "y1")
  xy_poly[[j]] <- rbind(xy_pj_com, xy_poly[[j]])
  xy_poly[[j]] %<>% 
    as_tibble() %>%
    mutate(phiC = acos((x1-x_med_new)/sqrt((x1-x_med_new)^2+(y1-y_eq_new)^2))*(-1)^(y1<y_eq_new) + 2*pi*(y1<y_eq_new)) %>%  #angle
    mutate(thetaC = NA) %>% #radius
    transmute(x1, y1, phiC, thetaC) %>%
    data.matrix()
  for (k in 1:dim(xy_poly[[j]])[1]) {
    alpha <- xy_poly[[j]][k,'phiC']
    line_st <- st_linestring(rbind(xy_ori, c(xy_ori[1] + R*cos(alpha), xy_ori[2] + R*sin(alpha))))
    int <- data.matrix(st_intersection(line_st, poly_st))[2,]
    xy_poly[[j]][k,'thetaC'] <- pi/2 * dist(rbind(xy_poly[[j]][k,1:2], xy_ori)) / dist(rbind(int,xy_ori))
  }
  # now turn inside out and a 90-rotation about x-axis to have the front edge on the x-z plane, 
  # angle wrt looking from behiind the eye
  xy_poly[[j]] %<>%
    as_tibble() %>%
    mutate(thetaC = pi - thetaC) %>% # turn inside out
    mutate(phiC = phiC - pi/2 - angR_new) %>% # align front edge to x-z plane
    mutate(x = sin(thetaC)*cos(phiC), y = sin(thetaC)*sin(phiC), z = cos(thetaC)) %>% #(x,y,z) coord
    mutate(xr = x, yr = -z, zr = y) %>% # +90 rotation about x-axis
    mutate(xrr = xr, yrr = yr, zrr = zr) %>% # do nothing
    mutate(theta = acos(zrr/sqrt(xrr^2+yrr^2+zrr^2)), phi = acos(xrr/sqrt(xrr^2+yrr^2))) %>%
    mutate(theta_deg = theta/pi*180, phi_deg = phi/pi*180/180*diff(buchner_phi)+buchner_phi[1]) %>% # longitude use buchner_phi
    select(x1, y1, theta_deg, phi_deg) %>%
    data.matrix()
  xy_com[[j]] <- xy_poly[[j]][1,]
  xy_poly[[j]] <- xy_poly[[j]][-1,]
}


# Figure 5D
windows(record = F, width = 8, height = 8)
bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (j == 10) {
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "#d7191c", density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 3, pch = 3)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
  else if (j == 51) {
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "#2c7bb6", density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 3, pch = 3)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
  else {
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "cyan", border = 'black', density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
}
xy_bd <- xy_bd[-1,]
hpts_2 <- chull(xy_bd)
hpts_2 <- c(hpts_2, hpts_2[1])
xy_bd_chull <- xy_bd[hpts_2,] # hull edge points
polygon(xy_bd_chull)
lines(rbind(c(-11,180), c(-2,180)), lwd = 3)
text(-5, 180, labels = expression(paste("9",degree)), pos = 1, offset = 0.3)
lines(rbind(c(-11,180), c(-11,171)), lwd = 3)
text(-11, 175, labels = expression(paste("9",degree)), pos = 2, offset = 0.2)
lines(rbind(c(0,0), c(0,180)), lwd = 3) 
text(0, -5, labels = "front", pos = 1, offset = 0)
lines(rbind(c(90,0), c(90,180)), lwd = 3)
text(90, -5, labels = expression(paste("side 90",degree)), pos = 1, offset = 0)
lines(rbind(c(-12,90), c(162,90)), lwd = 3)
text(-17, 90, labels = "equator", pos = 1, offset = 0, srt = 90)


# Figure S6A
bd_phi2 <- seq(buchner_phi[1], buchner_phi[2], by = 2)
bd_theta2 <- seq(1, 180, by = 2)
pt_grid <- expand.grid(bd_phi2, bd_theta2) 
pt_grid <- cbind(pt_grid, rep(0, dim(pt_grid)[1]))
colnames(pt_grid) <- c('x','y','paint')
for (j in 1:length(xy_poly)) {
  ii <- sp::point.in.polygon(pt_grid[,"x"], pt_grid[,"y"], xy_poly[[j]][,"phi_deg"], xy_poly[[j]][,"theta_deg"])
  pt_grid[,"paint"] <- pt_grid[,"paint"] + ii
}
pt_grid <- pt_grid[pt_grid[,"paint"] != 0, ] 
pt_grid[,"paint"] <- factor(pt_grid[,"paint"], labels = seq(1,5, by = 1))

# PLOT
windows(record = F, width = 8, height = 8)
xy_bd_chull_df <- as.data.frame(xy_bd_chull)
x_tick <- c(0, 45, 90, 135, 180) + buchner_phi[1]
y_tick <- c(0, 45, 90, 135, 180)
dev.new()
ggplot() + 
  geom_point(data = pt_grid, aes(x = x, y = y, colour = paint)) +
  scale_color_brewer(palette =  "Blues") +
  scale_y_reverse() +
  guides(col = guide_legend(title = "groupings")) +
  geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg)) +
  geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
  annotate("text", x = -5, y = 185, label = "9°") +
  geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
  annotate("text", x = -15, y = 175.5, label = "9°") +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
  geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
  geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  theme(axis.title = element_blank(), axis.text = element_blank()) +
  theme_void() +
  labs(title = "2d projection on LO")




# LC6 glomerulus ------------------------------------------------------------------------------------------------

# all the connectors on LC6
LC6_pre <- data.frame()
LC6_post <- data.frame()
conn_glu <- list()
conn_pre <- list()
conn_post <- list()
for (j in 1:length(neu)) {
  tar <- neu[[j]]
  conn_glu[[j]] <- tar$connectors %>% 
    as_tibble() %>%
    mutate(glu = a*x + b*y + c*z + d) %>%
    filter(glu < 0)
  conn_pre[[j]] <- conn_glu[[j]] %>% 
    filter(prepost == 0) %>%
    select(treenode_id, connector_id, x, y, z) %>%
    as.data.frame()
  conn_post[[j]] <- conn_glu[[j]] %>% 
    filter(prepost == 1) %>%
    select(treenode_id, connector_id, x, y, z) %>%
    as.data.frame()
  LC6_pre <- rbind(LC6_pre, cbind(rep(j, dim(conn_pre[[j]])[1]),conn_pre[[j]]))
  LC6_post <- rbind(LC6_post, cbind(rep(j, dim(conn_post[[j]])[1]),conn_post[[j]]))
}
colnames(LC6_pre)[1] <- "ID"
colnames(LC6_post)[1] <- "ID"

if (exists("LC6LC6")) {
  rm(LC6LC6)
}

for (j in 1:length(neu)) {
  # tar_pre <- neu[[j]]
  for (k in 1:length(neu)) {
    # tar_post <- neu[[k]]
    cft <- catmaid_get_connectors_between(pre_skids = neu[[j]]$skid, post_skids = neu[[k]]$skid)
    if (!is.null(cft)) {
      cft_nrow <- dim(cft)[1]
      cft_comb <- cbind(rep(j,cft_nrow), rep(k,cft_nrow), cft[,1:8])
      if (exists("LC6LC6")) {
        LC6LC6 <- rbind(LC6LC6,cft_comb)
      } else {
        LC6LC6 <- cft_comb
      }
    }
  }
}
LC6LC6 <- as.data.frame(LC6LC6)
colnames(LC6LC6) <- c('pre_ind','post_ind','pre_skid','post_skid','conn_id','pre_node_id','post_node_id','x','y','z')
rownames(LC6LC6) <- seq(1,dim(LC6LC6)[1])


# -- look at pairwise distance vs pairwise connection
a1 = -0.6; b1 = 1; c1 = 1.3; d1 = -200000 #separate LO and glomerulus
dist_Nconn <- matrix(ncol = 7) #pairwise dist and connections
# colnames(dist_Nconn) <- c("from","to", "dist", "fromto_LO", "fromto_glu", "tofrom_LO", "tofrom_glu")
for (j in 1:(length(conn_pre)-1)) {
  for (k in (j+1):length(conn_pre)) {
    dist_d <- dist(rbind(xyz_com[[j]], xyz_com[[k]]))
    conn_fromto <- catmaid_get_connectors_between(pre_skids = neu[[j]]$skid, post_skids = neu[[k]]$skid)
    if (is.null(conn_fromto)) {
      fromto_LO <- 0
      fromto_glu <- 0
    } else {
      mat_conn <- data.matrix(conn_fromto[,c("connector_x", "connector_y", "connector_z")])
      fromto_LO <- sum(mat_conn %*% c(a1,b1,c1) + d1 > 0)
      fromto_glu <- dim(conn_fromto)[1] - fromto_LO
    }
    conn_tofrom <- catmaid_get_connectors_between(pre_skids = neu[[k]]$skid, post_skids = neu[[j]]$skid)
    if (is.null(conn_tofrom)) {
      tofrom_LO <- 0
      tofrom_glu <- 0
    } else {
      mat_conn <- data.matrix(conn_tofrom[,c("connector_x", "connector_y", "connector_z")])
      tofrom_LO <- sum(mat_conn %*% c(a1,b1,c1) + d1 > 0)
      tofrom_glu <- dim(conn_tofrom)[1] - tofrom_LO
    }
    dist_Nconn_tmp <- matrix(c(j,k,dist_d, fromto_LO, fromto_glu, tofrom_LO, tofrom_glu), nrow = 1)
    dist_Nconn <- rbind(dist_Nconn, dist_Nconn_tmp)
  }
}
dist_Nconn <- dist_Nconn[-1,]
colnames(dist_Nconn) <- c("from","to","dist_com","fromto_LO","fromto_glu","tofrom_LO","tofrom_glu")
dist_Nconn %<>% 
  as_tibble() %>%
  mutate(Nconn_glu = tofrom_glu + fromto_glu) %>%
  as.data.frame()
colSums(dist_Nconn)
dist_com_mean <- mean(dist_Nconn$dist_com)
Nconn_glu_tot <- sum(dist_Nconn$Nconn_glu)


# dist to a neuron avg over all partners weighted by num of conn
dist_Nconn %<>% as_tibble() %>%
  mutate(distxN = dist_com*Nconn_glu) %>%
  as.data.frame()
LC6_avgdist_wt <- matrix(ncol = 2, nrow = length(conn_pre))
LC6_avgdist <- c()
LC6_nnbdist <- c()
# N_nb <- c()
for (j in 1:length(conn_pre)) {
  tib_tmp <- dist_Nconn %>%
    as_tibble() %>%
    mutate(bool = from == j | to == j) %>%
    filter(bool == TRUE) %>%
    # transmute(Nconn_glu, distxN) %>%
    data.matrix()
  # N_nb <- c(N_nb, sum(tib_tmp[,"distxN"] != 0))
  LC6_avgdist[j] <- sum(tib_tmp[,"dist_com"])/dim(tib_tmp)[1]
  LC6_nnbdist[j] <- min(tib_tmp[,"dist_com"])
  LC6_avgdist_wt[j,] <- c(j, sum(tib_tmp[,"distxN"])/sum(tib_tmp[,"Nconn_glu"]))
}
# avg_nb <- mean(N_nb)
# randomize connection partners
dist_nonzero <- dist_Nconn %>%
  filter(distxN != 0) %>%
  transmute(Nconn_glu, distxN) %>%
  data.matrix()
mat_avgdist_rand <- matrix(ncol = 2, nrow = length(conn_pre))
for (j in 1:length(conn_pre)) {
  # -- maintain connection num for given LC6
  ii_rand <- sample.int(length(conn_pre)-1)
  tib_tmp <- dist_Nconn %>%
    as_tibble() %>%
    mutate(bool = from == j | to == j) %>%
    filter(bool == TRUE) %>%
    select(dist_com, Nconn_glu) %>%
    mutate(dist_com = dist_com[ii_rand]) %>%
    mutate(distxN = dist_com*Nconn_glu) %>%
    as.data.frame()
    mat_avgdist_rand[j,] <- c(j, sum(tib_tmp[,"distxN"])/sum(tib_tmp[,"Nconn_glu"]))
}


# Figure 5E
dat_ggplot <- data.frame(rbind(LC6_avgdist_wt, mat_avgdist_rand))
gp  <- factor(c(rep(1,length(conn_pre)),rep(2,length(conn_pre))), labels = c("EM data","Randomized"))
dat_ggplot <- cbind(dat_ggplot, gp)
colnames(dat_ggplot) <- c("neu", "dist", "group")
dat_ggplot$dist <- dat_ggplot$dist/1000

windows(record = F, width = 8, height = 8)
p <- ggplot(dat_ggplot, aes(x = neu, y = dist, colour = gp)) + 
  geom_point(shape = 16, size = 3) +
  theme_void() +
  theme(legend.position = "top", plot.margin = unit(c(3,-5.5,4,3), "mm")) + 
  labs(title = "Avg distance over connection num") + 
  xlim(0, 65) + 
  ylim(30, 80) +
  xlab("neurons index") + 
  ylab("Avg dist [um]") +
  coord_fixed(ratio = 4) +
  geom_hline(yintercept = mean(LC6_nnbdist)/1000, linetype = 2) +
  geom_hline(yintercept = mean(LC6_avgdist)/1000, linetype = 2) +
  scale_colour_discrete(name="Groups")
ggMarginal(p, margins = "y", size = 2, type = "boxplot", outlier.size =3, groupColour = TRUE, groupFill = TRUE)

# KS test
ks.test(LC6_avgdist_wt, mat_avgdist_rand)

# Mann-Whitney
wilcox.test(LC6_avgdist_wt, mat_avgdist_rand, alternative = 'less')


# --- glomerulus dolphin compartment

N_gp <- 11

# for making ahull
conn_LC6 <- LC6_post # only post-synapses to divide the dolphin
glu_div <- quantile(conn_LC6$x, probs = seq(0,1,length.out = N_gp)) #separate into divisions of equal pts based on x
glu_div[1] <- glu_div[1]-1
conn_LC6 <- conn_LC6 %>%
  as_tibble() %>%
  mutate(gp_x = 0)
for (j in 1:length(glu_div)) {
  conn_LC6 %<>% as_tibble() %>%
    mutate(gp_x = gp_x + (x>glu_div[j]))
}

# assign group
conn_LC6LC6 <- LC6LC6
conn_LC6LC6 <- conn_LC6LC6 %>%
  as_tibble() %>%
  mutate(gp_x = 0)
for (j in 1:length(glu_div)) {
  conn_LC6LC6 %<>% as_tibble() %>%
    mutate(gp_x = gp_x + (x>glu_div[j]))
}

ID_wt <- list() #neuron skid with synapse weight in each division
for (j in 1:length(glu_div)) {
  df_tmp <- conn_LC6LC6 %>%
    filter(gp_x == j) %>%
    # transmute(ID) %>%
    data.frame()
  mat_tmp <- matrix(ncol = 2, nrow = length(neu))
  for (k in 1:length(neu)) {
    mat_tmp[k,] <- c(k, sum(df_tmp$pre_ind == k | df_tmp$post_ind == k))  
  }
  ID_wt[[j]] <- mat_tmp
}

# avg radius of projection
avg_size_xy <- c()
for (j in 1:length(neu)) {
  tmp <- as.data.frame(xy_poly[[j]])
  avg_size_xy <- c(avg_size_xy, (max(tmp$theta_deg)-min(tmp$theta_deg)+max(tmp$phi_deg)-min(tmp$phi_deg))/4)
}
r_xy <- mean(avg_size_xy)
# half-width sqrt(-log(0.5))*r_xy*2
# plot Gaussian around each com with synap# as height
ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])

n_lvl <- 11 # 10 compartments
plvl <- list() # ggplot
grid_Gaussian_ls <- list()
grid_Gaussian_sum <- matrix(ncol = 4)
colnames(grid_Gaussian_sum) <- c("X","Y","Z","gp")
plvl_sum <- ggplot()
getPalette <- colorRampPalette(brewer.pal(9, "RdYlBu"))
dolphin_col <- getPalette(n_lvl - 1)
pal_dolphin <- c()
gprange <- list()
for (j in 1:(length(glu_div)-1)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(conn_pre)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- ID_wt[[j]][k,2]
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  }
  
  breaks <- seq(0,180,length.out = 9)
  grid_Gaussian$equalSpace <- cut(grid_Gaussian$Z, breaks) #cut into 18 levels
  grid_Gaussian <- cbind(grid_Gaussian, rep(j,dim(grid_Gaussian)[1]))
  colnames(grid_Gaussian) <- c("X","Y","Z","equalSpace","gp")

  plvl[[j]] <- ggplot(grid_Gaussian, aes(X, Y, z= Z)) + 
    geom_raster(aes(fill = equalSpace), interpolate = T) +
    geom_contour(color = "black", alpha = 0.5) +
    theme_void() +
    theme(axis.title = element_blank()) +
    scale_fill_manual(values = (brewer.pal(8, 'Blues')), breaks = seq(0, 190, length.out = 8)) +
    scale_y_reverse() +
    xlim(0, 180) +
    ylim(0, 180) +
    labs(title = paste("Group", j))
  
  cutoffZ <- max(grid_Gaussian$Z) 
  cutoff_per <- 0.5
  grid_Gaussian <- grid_Gaussian[grid_Gaussian$Z > cutoffZ*cutoff_per, ]
  grid_Gaussian_sum <- rbind(grid_Gaussian_sum, grid_Gaussian[,c("X","Y","Z","gp")])
  grid_Gaussian_ls[[j]] <- grid_Gaussian[,c("X","Y","Z","gp")]
  gprange[[j]] <- range(grid_Gaussian_ls[[j]]$Z)
}
grid_Gaussian_sum <- grid_Gaussian_sum[-1,]
grid_Gaussian_sum$gp <- factor(grid_Gaussian_sum$gp)


# plot, all panels, 
dev.new()
plot_grid(plvl[[1]],
          plvl[[2]],
          plvl[[3]],
          plvl[[4]],
          plvl[[5]],
          plvl[[6]],
          plvl[[7]],
          plvl[[8]],
          plvl[[9]],
          plvl[[10]])


# Figure 5F
# color contours
xy_bd_chull_df <- as.data.frame(xy_bd_chull)
gg_cont <- ggplot(grid_Gaussian_sum, aes(X, Y, z = Z)) +
  coord_fixed(ratio = 1) +
  scale_y_reverse()
for (j in 1:(length(glu_div)-1)) {
  gg_cont <- gg_cont + 
    geom_contour(data = grid_Gaussian_ls[[j]], aes(X,Y,z=Z), breaks = seq(gprange[[j]][1], gprange[[j]][2], length.out = 3), color = dolphin_col[j], alpha = 0.9, lwd = 1)
}
gg_cont <- gg_cont + 
  geom_polygon(data = xy_bd_chull_df, aes(x=phi_deg, y=theta_deg, z=0),colour = "black", alpha = 0) +
  theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
  geom_segment(aes(x = -8, y = 90, xend = 160, yend = 90), size = 2) +
  geom_segment(aes(x = 10, y = 180, xend = 19, yend = 180), size=2, lineend = "round") +
  annotate("text", x = 14.5, y = 185, label = "9°") +
  geom_segment(aes(x = 10, y = 180, xend = 10, yend = 171), size=2, lineend = "round") +
  annotate("text", x = 5, y = 175.5, label = "9°") +
  annotate("text", x = 50, y = 0, label = paste(sprintf("%.1f %%", 100*cutoff_per)))

# make alpha hull
conn_LC6_hull <- as.data.frame(conn_LC6)[c(-5646,-5679,-5640,-5631,-5654,-5652),]
gp_x <- factor(conn_LC6_hull[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6_hull[,"gp_x"] <- gp_x
ii_ahull <- ashape(unique(conn_LC6_hull[, c("x","y")]), alpha = 7000)$edges[, 3:6]
ii_ahull <- as.data.frame(ii_ahull)
# make separate alpha hulls
ii_ahull_ls <- list()
for (j in 1:(length(glu_div)-1)) {
  ii_ahull_ls[[j]] <- ashape(unique(conn_LC6_hull[conn_LC6_hull$gp_x == j, c("x","y")]), alpha = 7000)$edges[, 3:6]
  ii_ahull_ls[[j]] <- as.data.frame(ii_ahull_ls[[j]])
}

# com of each division
range_syn <- matrix(nrow = length(glu_div)-1, ncol = 4)
for (j in 1:(length(glu_div)-1)) {
  tmp <- conn_LC6 %>%
    as_tibble() %>%
    filter(gp_x == j) %>%
    transmute(x, y) %>%
    as.data.frame() 
  range_syn[j,1:2] <- range(tmp$x)
  range_syn[j,3:4] <- range(tmp$y)
}

# num of syn in each div
N_syn_LC6 <- c()
for (j in 1:(length(glu_div)-1)) {
  N_syn_LC6[j] <- sum(conn_LC6LC6$gp_x == j)
}

# dolphin in colors
conn_LC6LC6 <- as.data.frame(conn_LC6LC6)
gp_x <- factor(conn_LC6LC6[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6LC6[,"gp_x"] <- gp_x

pglu <- ggplot(conn_LC6LC6) + 
  geom_point(aes(x = x, y = y, colour = gp_x), shape = 16)
for (j in 1:(length(glu_div)-1)) { #add ahulls
  pglu <- pglu + 
    geom_segment(data = ii_ahull_ls[[j]], aes(x = x1, y = y1, xend = x2, yend = y2)) +
    annotate("text", x = mean(range_syn[j,c(1,2)]), y = mean(range_syn[j,c(3,4)]), label = paste(sprintf("%.1f %%", 100*N_syn_LC6[j]/sum(N_syn_LC6))))
}
pglu <- pglu +
  scale_colour_manual(values = dolphin_col) +
  guides(col = guide_legend(title = "groupings")) +
  theme(legend.position = "bottom") +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  scale_y_reverse() +
  ylim(205000, 180000)+
  coord_fixed(ratio = 1) +
  geom_segment(aes(x = 380000, y = 205000, xend = 390000, yend = 205000), size=2, lineend = "round") +
  annotate("text", x = 385000, y = 204000, label = "10 µm") + 
  labs(title = "LC6 Synapses distribution in glomerulus")

windows(record = F, width = 16, height = 10)
ggdraw() +
  draw_plot(pglu + theme(legend.justification = "bottom"), 0, 0.51, 1, 0.5) +
  draw_plot(gg_cont + theme(axis.title = element_text()), 0, 0, 0.5, 0.5)



# Figure S6B, 10 individual compartments
ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
plvl <- list()
simdata <- list()
simdata_df <- list()
mat_names <- paste("compartment", seq(1,10),sep = " ")

for (j in 1:(n_lvl-1)) {
  getPalette <- colorRampPalette(c("white", dolphin_col[j]))
  pal_tar_red <- getPalette(n_lvl - 1)
  
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(conn_pre)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- ID_wt[[j]][k,2]
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  }
  simdata[[j]] <- grid_Gaussian
  simdata_df[[j]] <- simdata[[j]]
  colnames(simdata_df[[j]]) <- c("x","y","z")
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  plvl[[j]] <- ggplot(simdata_df[[j]], aes(x, y, z = z)) +
    geom_raster(aes(fill = equalSpace), interpolate = F) +
    scale_fill_manual(values = pal_tar_red, guide=FALSE) +
    theme_void() +
    geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE, lwd = 2) +
    # geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
    # annotate("text", x = -5, y = 185, label = "9°") +
    # geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
    # annotate("text", x = -15, y = 175.5, label = "9°") +
    geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
    geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
    geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
    scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
    scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
    # labs(title = mat_names[j]) +
    geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = seq(0,1,by = 0.1), color = dolphin_col[j], alpha = 1, lwd = 1)+
    coord_fixed(ratio = 1)
}

dev.new()
plvl[[1]]
dev.new()
plvl[[2]]
dev.new()
plvl[[3]]
dev.new()
plvl[[4]]
dev.new()
plvl[[5]]
dev.new()
plvl[[6]]
dev.new()
plvl[[7]]
dev.new()
plvl[[8]]
dev.new()
plvl[[9]]
dev.new()
plvl[[10]]



# target neuron RF wt by synapses ---------------------------------------------------------------------------------

a1 = -0.6; b1 = 1; c1 = 1.3; d1 = -200000 #separate LO and glomerulus
conn_target <- list()
for (j in 1:length(neu_target)) {
  tb_conn <- matrix(ncol = 7)
  for (k in 1:length(neu)) {
    conn_fromto <- catmaid_get_connectors_between(pre_skids = neu_target[[j]]$skid, post_skids = neu[[k]]$skid)
    if (is.null(conn_fromto)) {
      fromto_LO <- 0
      fromto_glu <- 0
    } else {
      mat_conn <- data.matrix(conn_fromto[, c("connector_x", "connector_y", "connector_z")])
      fromto_LO <- sum(mat_conn %*% c(a1, b1, c1) + d1 > 0)
      fromto_glu <- dim(conn_fromto)[1] - fromto_LO
    }
    conn_tofrom <- catmaid_get_connectors_between(pre_skids = neu[[k]]$skid, post_skids = neu_target[[j]]$skid)
    if (is.null(conn_tofrom)) {
      tofrom_LO <- 0
      tofrom_glu <- 0
    } else {
      mat_conn <- data.matrix(conn_tofrom[, c("connector_x", "connector_y", "connector_z")])
      tofrom_LO <- sum(mat_conn %*% c(a1, b1, c1) + d1 > 0)
      tofrom_glu <- dim(conn_tofrom)[1] - tofrom_LO
    }
    Nconn_glu <- fromto_glu + tofrom_glu
    tb_conn_tmp <- matrix(c(j, k, fromto_LO, fromto_glu, tofrom_LO, tofrom_glu, Nconn_glu), nrow = 1)
    tb_conn <- rbind(tb_conn, tb_conn_tmp)
  }
  conn_target[[j]] <- as.data.frame(tb_conn[-1,])
  colnames(conn_target[[j]]) <- c("target","LC6","fromto_LO","fromto_glu","tofrom_LO","tofrom_glu", "Nconn_glu")
}

conn_tt <- matrix(ncol = length(neu_target), nrow = length(neu_target)) #target to target
for (j in 1:length(neu_target)) {
  for (k in 1:length(neu_target)) {
    conn_fromto <- catmaid_get_connectors_between(pre_skids = neu_target[[j]]$skid, post_skids = neu_target[[k]]$skid)
    if (!is.null(conn_fromto)) {
      mat_conn <- data.matrix(conn_fromto[, c("connector_x", "connector_y", "connector_z")])
      fromto <- dim(conn_fromto)[1]
      conn_tt[j,k] <- fromto
    }
  }
}


# Figure S7
# PLOT,  Gaussian around each com with synap# as height,  cp data via binning
ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
plvl <- list()
simdata <- list()
simdata_df <- list()
LC6_tar_median <- list()
for (j in 1:length(conn_target)) {
  LC6_tar_median[[j]] <- (quantile(conn_target[[j]]$tofrom_glu, c(0.0)))
}
mat_names <- c(paste("biL_", biL_skid, sep = ""), 
               paste("biR_", biR_skid, sep = ""), 
               paste("ipsi_", ipsi_skid, sep = ""))

n_lvl <- 11
breaks_3 <- seq(0,70,length.out = n_lvl)
getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
pal_tar <- getPalette(n_lvl - 1)

x_tick <- c(16, 16+45, 16+90, 16+135, 16+180)
y_tick <- c(0, 45, 90, 135, 180)

for (j in 1:length(neu_target)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(neu)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- conn_target[[j]][k,"tofrom_glu"]
    if (A >= LC6_tar_median[[j]]) { # selected neuron
      grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
    }
  }
  
  grid_Gaussian_cut <- grid_Gaussian %>%
    as.data.frame()
  
  # -- whole range no binning
  simdata[[j]] <- grid_Gaussian_cut
  simdata_df[[j]] <- simdata[[j]]
  colnames(simdata_df[[j]]) <- c("x","y","z")
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))

  simdata_df[[j]]$z <- simdata_df[[j]]$z / mean(head(sort(simdata_df[[j]]$z, decreasing = T)))
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  plvl[[j]] <-
    ggplot(simdata_df[[j]], aes(x, y, z = z)) +
    geom_raster(aes(fill = equalSpace), interpolate = F) +
    scale_fill_manual(values = pal_tar) +
    geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
    geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
    annotate("text", x = -5, y = 185, label = "9°") +
    geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
    annotate("text", x = -15, y = 175.5, label = "9°") +
    geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
    geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
    geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
    scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
    scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
    labs(title = mat_names[j]) +
    theme_void()+
    coord_fixed(ratio = 1)

  plvl[[j]] <- plvl[[j]] +
    geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.7), color = "blue", alpha = 1, lwd = 1)
}

windows(record = F, width = 8, height = 8)
plvl[[1]]
windows(record = F, width = 8, height = 8)
plvl[[2]]
windows(record = F, width = 8, height = 8)
plvl[[3]]
windows(record = F, width = 8, height = 8)
plvl[[4]]
windows(record = F, width = 8, height = 8)
plvl[[5]]
windows(record = F, width = 8, height = 8)
plvl[[6]]
windows(record = F, width = 8, height = 8)
plvl[[7]]
windows(record = F, width = 8, height = 8)
plvl[[8]]
windows(record = F, width = 8, height = 8)
plvl[[9]]


# Figure 6C
ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
plvl <- list()
simdata <- list()
simdata_df <- list()
mat_names <- c("bi_all", "ipsi_all")

conn_target_agglo_sum <- list()
conn_target_agglo_sum[[1]] <- conn_target[[1]] + conn_target[[2]] + conn_target[[3]] + conn_target[[4]]
conn_target_agglo_sum[[2]] <- conn_target[[5]]+conn_target[[6]]+conn_target[[7]]+conn_target[[8]]+conn_target[[9]]

for (j in 1:length(conn_target_agglo_sum)) {  
    grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
    grid_Gaussian$Z = 0
    colnames(grid_Gaussian) <- c("X","Y","Z")
    for (k in 1:length(neu)) {
      x0 <- (xy_com[[k]]["phi_deg"])
      y0 <- (xy_com[[k]]["theta_deg"])
      A <- conn_target_agglo_sum[[j]][k,"tofrom_glu"]
      if (A > 0) { # selected neuron
        grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
      }
    }
    simdata[[j]] <- grid_Gaussian
    simdata_df[[j]] <- simdata[[j]]
    colnames(simdata_df[[j]]) <- c("x","y","z")
    # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
    # simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
    # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
    # simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
    
    simdata_df[[j]]$z <- simdata_df[[j]]$z / mean(head(sort(simdata_df[[j]]$z, decreasing = T)))
    simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
    
    plvl[[j]] <-
      ggplot(simdata_df[[j]], aes(x, y, z = z)) +
      geom_raster(aes(fill = equalSpace), interpolate = F) +
      scale_fill_manual(values = pal_tar) +
      geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
      geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
      annotate("text", x = -5, y = 185, label = "9°") +
      geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
      annotate("text", x = -15, y = 175.5, label = "9°") +
      geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
      geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
      geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
      scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
      scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
      labs(title = mat_names[j]) +
      theme_void() +
      coord_fixed(ratio = 1)

    plvl[[j]] <- plvl[[j]] +
      geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 1)
}

windows(record = F, width = 8, height = 8)
plvl[[1]] 
windows(record = F, width = 8, height = 8)
plvl[[2]] 




# # physiology data --------------------------------------------------------------------------------------------
# 
# library(R.matlab)
# expBi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss825_mean_RFmap_n=4.mat") #bi
# expBi <- as.matrix(expBi[[1]])
# expIpsi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss2036_mean_RFmap_n=7.mat") #ipsi
# expIpsi <- as.matrix(expIpsi[[1]])
# 
# slope_ex <- (90-83.28)/(117-16)
# (90-83.28)/(117-16)*(9+16)
# 
# n_lvl <- 11
# breaks_3 <- seq(0,1,length.out = n_lvl)
# getPalette <- colorRampPalette(brewer.pal(9, "Oranges"), space = "Lab", bias = 1.8)
# pal_tar_ex <- getPalette(n_lvl - 1)
# 
# x2 <- seq(0-4.25,117-4.25,by = 9) # -18 to 99, with sim data
# y2 <- seq(54, 108,by = 9)
# xygrid2 <- expand.grid(x2, y2)
# 
# expBi_df <- data.frame(xygrid2, as.vector(t(expBi)))
# colnames(expBi_df) <- c("x","y","z")
# expBi_df$equalSpace <- cut(expBi_df$z, breaks_3)
# expBi_df$z <- expBi_df$z*1
# dev.new()
# ggplot(expBi_df, aes(x, y, z = z)) +
#   geom_raster(data = expBi_df, aes(x,y,fill = equalSpace), interpolate = F) +
#   scale_fill_manual(values = pal_tar_ex, guide_legend("dF/F"), labels = paste(seq(0.1,1,length.out = 10))) +
#   geom_segment(aes(x = -9, y = 112, xend = 0, yend = 112), size=2, lineend = "round") +
#   annotate("text", x = 0, y = 116, label = "9°") +
#   geom_segment(aes(x = -9, y = 112, xend = -9, yend = 103), size=2, lineend = "round") +
#   annotate("text", x = -14.5, y = 102, label = "9°") +
#   # geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size=2, lineend = "round") +
#   geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
#   geom_segment(aes(x = 16, y = 90, xend = 117, yend = 90-righty), size=2, lineend = "round") +
#   # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
#   scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
#   scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
#   geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size = 2, lineend = "round") +
#   coord_fixed(ratio = 1)
# 
# expIpsi_df <- data.frame(xygrid2, as.vector(t(expIpsi)))
# colnames(expIpsi_df) <- c("x","y","z")
# expIpsi_df[14, "z"] <- expIpsi_df[14, "z"] - 0.01
# expIpsi_df$equalSpace <- cut(expIpsi_df$z, breaks_3)
# expIpsi_df$z <- expIpsi_df$z*1
# dev.new()
# ggplot(expIpsi_df, aes(x, y, z = z)) +
#   # geom_raster(aes(fill = equalSpace), interpolate = F) +
#   # scale_fill_manual(values = pal_tar_ex, guide_legend("synp den"), labels = paste(seq(0.1,1,length.out = 10))) +
#   geom_raster(data = expIpsi_df, aes(x,y,fill = equalSpace), interpolate = F) +
#   scale_fill_manual(values = pal_tar_ex, guide_legend("dF/F"), labels = paste(seq(0.1,1,length.out = 10))) +
#   geom_segment(aes(x = -9, y = 112, xend = 0, yend = 112), size=2, lineend = "round") +
#   annotate("text", x = 0, y = 116, label = "9°") +
#   geom_segment(aes(x = -9, y = 112, xend = -9, yend = 103), size=2, lineend = "round") +
#   annotate("text", x = -14.5, y = 102, label = "9°") +
#   # geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size=2, lineend = "round") +
#   geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
#   geom_segment(aes(x = 16, y = 90, xend = 117, yend = 90-righty), size=2, lineend = "round") +
#   # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
#   scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
#   scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
#   geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size = 2, lineend = "round") +
#   labs(title = paste("Target ipsi exp")) +
#   coord_fixed(ratio = 1)
# 
# 
# 
# 
# # Figure 6D
# # -- add up the same type
# mat_names <- c("bi_all", "ipsi_all")
# conn_target_agglo_sum <- list()
# conn_target_agglo_sum[[1]] <- conn_target[[1]] + conn_target[[2]] + conn_target[[3]] + conn_target[[4]]
# conn_target_agglo_sum[[2]] <- conn_target[[5]]+conn_target[[6]]+conn_target[[7]]+conn_target[[8]]+conn_target[[9]]
# 
# 
# ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
# plvl <- list()
# simdata <- list()
# simdata_df <- list()
# LC6_tar_median <- list()
# 
# for (j in 1:length(conn_target_agglo_sum)) {  
#   grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
#   ii <- grid_Gaussian[,1] > 16
#   grid_Gaussian[ii,2] <- slope_ex*(grid_Gaussian[ii,1] - 16) + grid_Gaussian[ii,2]
#   grid_Gaussian$Z = 0
#   colnames(grid_Gaussian) <- c("X","Y","Z")
#   for (k in 1:length(neu)) {
#     x0 <- (xy_com[[k]]["phi_deg"])
#     y0 <- (xy_com[[k]]["theta_deg"])
#     A <- conn_target_agglo_sum[[j]][k,"tofrom_glu"]
#     # if (A >= LC6_tar_median[[j]][j2]) { # selected neuron
#     grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
#     # }
#   }
#   grid_Gaussian[ii,2] <- round(- slope_ex*(grid_Gaussian[ii,1] - 16) + grid_Gaussian[ii,2])
#   simdata[[j]] <- grid_Gaussian
#   simdata_df[[j]] <- simdata[[j]]
#   colnames(simdata_df[[j]]) <- c("x","y","z")
#   simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
#   simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
#   simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
#   simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
#   
#   if (j <= 1) {
#     plvl[[j]] <- ggplot(expBi_df, aes(x, y, z = z))+
#       geom_raster(data = expBi_df, aes(x,y,fill = equalSpace), interpolate = F) +
#       scale_fill_manual(values = pal_tar_ex, guide_legend("dF/F"), labels = paste(seq(0.1,1,length.out = 10))) +
#       labs(title = paste("Target bi exp vs ", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glu"]),sep = " ")) +
#       geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)
#   } else {
#     plvl[[j]] <- ggplot(expIpsi_df, aes(x, y, z = z))+
#       geom_raster(data = expIpsi_df, aes(x,y,fill = equalSpace), interpolate = F) +
#       scale_fill_manual(values = pal_tar_ex, guide_legend("dF/F"), labels = paste(seq(0.1,1,length.out = 10))) +
#       labs(title = paste("Target ipsi exp vs", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glu"]),sep = " ")) +
#       geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)
#   }
#   
#   plvl[[j]] <- plvl[[j]] +
#     geom_segment(aes(x = -9, y = 112, xend = 0, yend = 112), size=2, lineend = "round") +
#     annotate("text", x = 0, y = 116, label = "9°") +
#     geom_segment(aes(x = -9, y = 112, xend = -9, yend = 103), size=2, lineend = "round") +
#     annotate("text", x = -14.5, y = 102, label = "9°") +
#     # geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size=2, lineend = "round") +
#     geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
#     geom_segment(aes(x = 16, y = 90, xend = 117, yend = 90-righty), size=2, lineend = "round") +
#     # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
#     scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
#     scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
#     geom_segment(aes(x = 16, y = 20, xend = 16, yend = 160), size = 2, lineend = "round") +
#     # geom_segment(aes(x = 0, y = 90, xend = 180, yend = 90), size = 2) +
#     coord_fixed(ratio = 1)
# }
# 
# windows(record = F, width = 8, height = 8)
# plvl[[1]] 
# windows(record = F, width = 8, height = 8)
# plvl[[2]] 



# target neurons in dolphin glomerulus ----------------------------------------------------------------------------------------------

N_gp <- 11

conn_LC6_tar <- list()
LC6_wt <- list() #neuron skid with synapse weight in each division for each target
for (j in 1:length(neu_target)) {
  tmp <- matrix(ncol = 4)
  for (k in 1:length(neu)) {
    tmp2 <- matrix(ncol = 3)
    conn_tofrom <- catmaid_get_connectors_between(pre_skids = neu[[k]]$skid, post_skids = neu_target[[j]]$skid)
    if (!is.null(conn_tofrom)) {
      tmp2 <- rbind(tmp2, data.matrix(conn_tofrom[, c("connector_x", "connector_y", "connector_z")]))
    }
    if (dim(tmp2)[1] == 2) {
      tmp <- rbind(tmp, c(rep(k,dim(tmp2)[1]-1), tmp2[-1,]))
    } else {
      tmp <- rbind(tmp, cbind(rep(k,dim(tmp2)[1]-1), tmp2[-1,]))
    }
  }
  tmp <- tmp[-1,] %>%
    as_tibble() %>%
    mutate(gp_x = 0)
  colnames(tmp) <- c("ID", "x","y","z","gp_x")
  
  for (k in 1:(length(glu_div)-1)) {
    tmp %<>% as_tibble() %>%
      mutate(gp_x = gp_x + (x>glu_div[k]))
  }
  
  ID_wt <- list() #neuron skid with synapse weight in each division
  for (k in 1:length(glu_div)) {
    df_tmp <- tmp %>%
      filter(gp_x == k) %>%
      transmute(ID) %>%
      data.frame()
    mat_tmp <- matrix(ncol = 2, nrow = length(neu))
    for (m in 1:length(neu)) {
      mat_tmp[m,] <- c(m, sum(df_tmp$ID %in% m))  
    }
    ID_wt[[k]] <- mat_tmp
  }
  
  conn_LC6_tar[[j]] <- as.data.frame(tmp)
  LC6_wt[[j]] <- ID_wt
}
# combine conn data
conn_LC6_tar_bi <- rbind(conn_LC6_tar[[1]], conn_LC6_tar[[2]], conn_LC6_tar[[3]], conn_LC6_tar[[4]])
gp_x <- factor(conn_LC6_tar_bi[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6_tar_bi[,"gp_x"] <- gp_x
conn_LC6_tar_ipsi <- rbind(conn_LC6_tar[[5]], conn_LC6_tar[[6]], conn_LC6_tar[[7]], conn_LC6_tar[[8]], conn_LC6_tar[[9]])
gp_x <- factor(conn_LC6_tar_ipsi[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6_tar_ipsi[,"gp_x"] <- gp_x


# PLOT, target with synapse
nopen3d()
plot3d(neu_target[[5]], col = 'lightblue', alpha = 0.7)
points3d(conn_LC6_tar[[5]][,c('x','y','z')], col = 'blue')
plot3d(neu_target[[6]], col = 'pink', alpha = 0.7)
points3d(conn_LC6_tar[[6]][,c('x','y','z')], col = 'red')


# histogram of num of synp in each compartment
LC6_ipsi_dol <- matrix(ncol = 5, nrow = 10)
for (j in 1:5) {
  for (k in 1:10) {
    LC6_ipsi_dol[k,j] <- sum(LC6_wt[[4+j]][[k]][,2])
  }
}


# Figure S8B, bar plot
num_tot <- colSums(LC6_ipsi_dol)
for (j in 1:5) {
  LC6_ipsi_dol[,j] <- LC6_ipsi_dol[,j]/num_tot[j]
}
dev.new()
barplot(LC6_ipsi_dol[,c(3,2,4,5,1)], main = 'LC6 to ipsi in dolphin compartments', ylim = c(-0.1, 1.2), col = dolphin_col )
# legend(x = 5, y = 800, legend = c(paste(seq(1,10))), cex = 1.3, fill = dolphin_col, horiz = F)
text(x = seq(0.7,5.5, length.out = 5), y = rep(1.1,5), labels = num_tot[c(3,2,4,5,1)])
text(x = seq(0.7,5.5, length.out = 5), y = rep(-0.05,5), labels = ipsi_skid[c(3,2,4,5,1)])


# separate targets
LC6_wt_biL <- list() # first 2 are biL, 2 biR, last 5 are ipsi
LC6_wt_biR <- list()
LC6_wt_bi <- list()
LC6_wt_ipsi <- list()
for (j in 1:(length(glu_div)-1)) {
  LC6_wt_biL[[j]] <- LC6_wt[[1]][[j]] + LC6_wt[[2]][[j]]
  LC6_wt_biR[[j]] <- LC6_wt[[3]][[j]] + LC6_wt[[4]][[j]]
  LC6_wt_bi[[j]] <- LC6_wt[[1]][[j]] + LC6_wt[[2]][[j]] + LC6_wt[[3]][[j]] + LC6_wt[[4]][[j]]
  LC6_wt_ipsi[[j]] <- LC6_wt[[5]][[j]] + LC6_wt[[6]][[j]] + LC6_wt[[7]][[j]] + LC6_wt[[8]][[j]] + LC6_wt[[9]][[j]]
}

# num of synapse per division
N_syn <- matrix(nrow = length(glu_div)-1, ncol = 3)
for (j in 1:(length(glu_div)-1)) {
  N_syn[j,1] <- colSums(LC6_wt_biL[[j]])[2]
  N_syn[j,2] <- colSums(LC6_wt_biR[[j]])[2]
  N_syn[j,3] <- colSums(LC6_wt_ipsi[[j]])[2]
}

# com of each division
com_syn <- cbind(rowMeans(range_syn[,1:2]), rowMeans(range_syn[,3:4]))


# Figure S8A
# -- bi all
plvl <- list()
grid_Gaussian_ls <- list()
grid_Gaussian_sum <- matrix(ncol = 4)
colnames(grid_Gaussian_sum) <- c("X","Y","Z","gp")
plvl_sum <- ggplot()
pal_dolphin <- c()
gprange <- list()
for (j in 1:(length(glu_div)-1)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(conn_pre)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- LC6_wt_biR[[j]][k,2] + LC6_wt_biL[[j]][k,2]
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  }
  breaks <- seq(0,180,length.out = 9)
  grid_Gaussian$equalSpace <- cut(grid_Gaussian$Z, breaks) #cut into 18 levels
  grid_Gaussian <- cbind(grid_Gaussian, rep(j,dim(grid_Gaussian)[1]))
  colnames(grid_Gaussian) <- c("X","Y","Z","equalSpace","gp")
  
  cutoffZ <- max(grid_Gaussian$Z) 
  cutofflevel <- 0.9
  grid_Gaussian <- grid_Gaussian[grid_Gaussian$Z > cutoffZ*cutofflevel, ]
  grid_Gaussian_sum <- rbind(grid_Gaussian_sum, grid_Gaussian[,c("X","Y","Z","gp")])
  grid_Gaussian_ls[[j]] <- grid_Gaussian[,c("X","Y","Z","gp")]
  gprange[[j]] <- range(grid_Gaussian_ls[[j]]$Z)
}
grid_Gaussian_sum <- grid_Gaussian_sum[-1,]
grid_Gaussian_sum$gp <- factor(grid_Gaussian_sum$gp)


xy_bd_chull_df <- as.data.frame(xy_bd_chull)
gg_cont <- ggplot(grid_Gaussian_sum, aes(X, Y, z = Z)) +
  geom_raster(aes(fill = gp)) +
  scale_fill_manual(values = dolphin_col) +
  coord_fixed(ratio = 1) +
  scale_y_reverse()

gg_cont <- gg_cont + 
  labs(title = paste("bi,", cutofflevel*100, "% cutoff"))+
  geom_polygon(data = xy_bd_chull_df, aes(x=phi_deg, y=theta_deg, z=0),colour = "black", alpha = 0) +
  theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
  annotate("text", x = -5, y = 185, label = "9°") +
  geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
  annotate("text", x = -15, y = 175.5, label = "9°") +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
  geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
  geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
  annotate("text", x = 5, y = 175.5, label = "9°") 

# dolphin in colors
pglu <- ggplot(conn_LC6_tar_bi) +
  geom_point(aes(x = x, y = y, colour = gp_x), shape = 16) 
  
for (j in 1:(length(glu_div)-1)) { #add ahulls
  pglu <- pglu + 
    geom_segment(data = ii_ahull_ls[[j]], aes(x = x1, y = y1, xend = x2, yend = y2))
}
pglu <- pglu +
  # geom_point(aes(x = x, y = y, colour = gp_x), shape = ".") +
  scale_colour_manual(values = dolphin_col) +
  guides(col = guide_legend(title = "groupings")) +
  theme(legend.position = "bottom") +
  theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  # scale_y_reverse() +
  ylim(205000, 180000)+
  # scale_x_reverse() +
  # xlab("x") + 
  # ylab("y") +
  coord_fixed(ratio = 1) +
  geom_segment(aes(x = 380000, y = 205000, xend = 390000, yend = 205000), size=2, lineend = "round") +
  annotate("text", x = 385000, y = 204000, label = "10 µm")+
  labs(title = "Synapses distribution in glumerulus")
for (j in 1:(length(glu_div)-2)) {
  glu_bd <- data.frame(x1 = range_syn[j,2], x2 = range_syn[j,2], y1 = range_syn[j,3], y2 = range_syn[j,4])
  # pglu <- pglu + geom_segment(data = glu_bd, aes(x = x1, y = y1, xend = x2, yend = y2))
}
for (j in 1:(length(glu_div)-1)) {
  pglu <- pglu + 
    # annotate("text", x = com_syn[j,1], y = com_syn[j,2], label = paste(sprintf("%.1f %%", 100*N_syn[j,2]/sum(N_syn[,2]))))
    annotate("text", x = com_syn[j,1], y = com_syn[j,2], label = paste(sprintf("%.1f ", 100*N_syn[j,2]/sum(N_syn[,2]))))
}

# PLOT, glu with 3 groups
windows(record = F, width = 16, height = 10)
ggdraw() +
  draw_plot(pglu + theme(legend.justification = "bottom"), 0, 0.51, 1, 0.5) +
  draw_plot(gg_cont + theme(axis.title = element_text()), 0, 0, 0.5, 0.5) 



# -- ipsi
ID_wt <- LC6_wt_ipsi
plvl <- list()
grid_Gaussian_ls <- list()
grid_Gaussian_sum <- matrix(ncol = 4)
colnames(grid_Gaussian_sum) <- c("X","Y","Z","gp")
plvl_sum <- ggplot()
# dolphin_col <- rainbow(N_gp + 1, alpha = 1)[1:(N_gp-1)]
pal_dolphin <- c()
gprange <- list()
for (j in 1:(length(glu_div)-1)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(conn_pre)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- ID_wt[[j]][k,2]
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  }
  
  breaks <- seq(0,180,length.out = 9)
  grid_Gaussian$equalSpace <- cut(grid_Gaussian$Z, breaks) 
  grid_Gaussian <- cbind(grid_Gaussian, rep(j,dim(grid_Gaussian)[1]))
  colnames(grid_Gaussian) <- c("X","Y","Z","equalSpace","gp")
  
  cutoffZ <- max(grid_Gaussian$Z) # j=1 -- 179, j=4 -- 108
  cutofflevel <- 0.9
  grid_Gaussian <- grid_Gaussian[grid_Gaussian$Z > cutoffZ*cutofflevel, ]
  grid_Gaussian_sum <- rbind(grid_Gaussian_sum, grid_Gaussian[,c("X","Y","Z","gp")])
  grid_Gaussian_ls[[j]] <- grid_Gaussian[,c("X","Y","Z","gp")]
  gprange[[j]] <- range(grid_Gaussian_ls[[j]]$Z)
}
grid_Gaussian_sum <- grid_Gaussian_sum[-1,]
grid_Gaussian_sum$gp <- factor(grid_Gaussian_sum$gp)


# color contours
xy_bd_chull_df <- as.data.frame(xy_bd_chull)
gg_cont <- ggplot(grid_Gaussian_sum, aes(X, Y, z = Z)) +
  geom_raster(aes(fill = gp)) +
  scale_fill_manual(values = dolphin_col) +
  coord_fixed(ratio = 1) +
  scale_y_reverse()
# for (j in 1:(length(glu_div)-1)) {
#   gg_cont <- gg_cont + 
#     geom_contour(data = grid_Gaussian_ls[[j]], aes(X,Y,z=Z), breaks = seq(gprange[[j]][1], gprange[[j]][2], length.out = 3), color = dolphin_col[j], alpha = 0.9, lwd = 3)
# }
titletext <- paste("ipsi,", cutofflevel*100, "% cutoff")
gg_cont <- gg_cont + 
  labs(title = titletext)+
  geom_polygon(data = xy_bd_chull_df, aes(x=phi_deg, y=theta_deg, z=0),colour = "black", alpha = 0) +
  theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
  annotate("text", x = -5, y = 185, label = "9°") +
  geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
  annotate("text", x = -15, y = 175.5, label = "9°") +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
  geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
  geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
  annotate("text", x = 5, y = 175.5, label = "9°") 

# dolphin in colors
pglu <- ggplot(conn_LC6_tar_ipsi) +
  geom_point(aes(x = x, y = y, colour = gp_x), shape = 16) 
  # geom_segment(data = ii_ahull, aes(x = x1, y = y1, xend = x2, yend = y2)) +
for (j in 1:(length(glu_div)-1)) { #add ahulls
  pglu <- pglu + 
    geom_segment(data = ii_ahull_ls[[j]], aes(x = x1, y = y1, xend = x2, yend = y2))
}
pglu <- pglu +
  scale_colour_manual(values = dolphin_col) +
  guides(col = guide_legend(title = "groupings", override.aes = list(shape = rep(".", 10)))) +
  theme_void() +
  theme(legend.position = "bottom") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  # scale_y_reverse() +
  ylim(205000, 180000)+
  coord_fixed(ratio = 1) +
  geom_segment(aes(x = 380000, y = 205000, xend = 390000, yend = 205000), size=2, lineend = "round") +
  annotate("text", x = 385000, y = 204000, label = "10 µm")+
  labs(title = "Synapses distribution in glumerulus")
for (j in 1:(length(glu_div)-2)) {
  glu_bd <- data.frame(x1 = range_syn[j,2], x2 = range_syn[j,2], y1 = range_syn[j,3], y2 = range_syn[j,4])
}
for (j in 1:(length(glu_div)-1)) {
  pglu <- pglu + 
    annotate("text", x = com_syn[j,1], y = com_syn[j,2], label = paste(sprintf("%.1f ", 100*N_syn[j,3]/sum(N_syn[,3]))))
}

# PLOT, glo with 3 groups
windows(record = F, width = 16, height = 10)
ggdraw() +
  draw_plot(pglu + theme(legend.justification = "bottom"), 0, 0.51, 1, 0.5) +
  draw_plot(gg_cont + theme(axis.title = element_text()), 0, 0, 0.5, 0.5) 



# Figure 6A, all targets with volume ref

nopen3d()
par3d('windowRect' = c(100,100,1100,1100))
shade3d(v14, alpha=0.05)
shade3d(glo_vol, col= "gold", alpha = 0.5)
rgl.viewpoint(userMatrix = rotationMatrix(1*90/180*pi+pi/2,1,0,0) %*% rotationMatrix(0,0,1,0), zoom = 0.7)
tar <- neu_biL[[1]]
plot3d(tar, lwd = 2, col = "blue")
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = "blue", size = 20)
tar <- neu_biR[[1]]
plot3d(tar, lwd = 2, col = "red")
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = "red", size = 20)
tar <- neu_ipsi[[1]]
plot3d(tar, lwd = 2, col = "cyan")
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = "cyan", size = 20)

# rgl.snapshot(filename = "LC6_targest.png",fmt = "png")
# rgl.snapshot(filename = "glo.png",fmt = "png")


# -- cp LM
nopen3d()
par3d('windowRect' = c(100,100,1100,1100))
# shade3d(v14, alpha=0.05)
shade3d(glo_vol, col = 'grey90', alpha = 0.3)
ipsi_pal <- brewer.pal(6, "RdYlBu")[c(1,2,3,5,6)]
tar <- neu_ipsi[[2]]
plot3d(tar, col = ipsi_pal[1])
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = ipsi_pal[1], size = 20)
tar <- neu_ipsi[[3]]
plot3d(tar, col = ipsi_pal[5])
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = ipsi_pal[5], size = 20)
rgl.viewpoint(userMatrix = rotationMatrix(1*90/180*pi+pi/2,1,0,0) %*% rotationMatrix(0,0,1,0), zoom = 1)

# Figure S7, individual targets RF
for (j in 1:2) {
  nopen3d()
  par3d('windowRect' = c(100,100,1100,1100))
  rgl.viewpoint(userMatrix = rotationMatrix(1*90/180*pi+pi/2,1,0,0) %*% rotationMatrix(0,0,1,0), zoom = 1)
  # shade3d(v14, alpha=0.05)
  shade3d(glo_vol, col = 'grey90', alpha = 0.3)
  tar <- neu_biR[[j]]
  # tar <- neu_biL[[j]]
  # tar <- neu_ipsi[[j]]
  plot3d(tar, col = ipsi_pal[2*j], lwd = 3)
  points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = ipsi_pal[2*j], size = 20)
  title3d(paste("biR_", biR_skid[j]))
  rgl.snapshot(filename = paste("biR_",  j, ".png", sep = ''),fmt = "png")
}

# rgl.snapshot(filename = "ipsi_x2.png",fmt = "png")


# data from Mai ---------------------------------------------------------------------------------------------------

library(cowplot)
library(R.matlab)
expBi3 <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss825_bi_RF_indiv_tseries.mat") #bi
expBi3 <- expBi3[[1]]
expBi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss825_mean_RFmap_n=4.mat") #bi
# expBi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss825_mean_AUC_RFmap_n=4.mat") #bi
expBi <- as.matrix(expBi[[1]])


expIpsi3 <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss2036_ipsi_RF_indiv_tseries.mat") #ipsi
expIpsi3 <- expIpsi3[[1]]
expIpsi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss2036_mean_RFmap_n=7.mat") #ipsi
# expIpsi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss2036_mean_AUC_RFmap_n=7.mat") #ipsi
expIpsi <- as.matrix(expIpsi[[1]])

# indexing
ind_mai <- c(t(matrix(seq(1:98), byrow = F, ncol = 14)))
tar_pal <- brewer.pal(4,"RdYlBu")[c(1,3,4,2)]


# x2 <- seq(-18 - 2.25, 99 - 2.25, by = 9) # -18 to 99, with sim data
# y2 <- seq(54, 108,by = 9)
# xygrid2 <- expand.grid(x2, y2) # looming center position 

# loom position from Matt
loom_theta_mat <- read.csv('loom_center_theta.txt', sep = ',', header = F) / pi * 180
loom_theta_mat <- 180 - loom_theta_mat
loom_theta_mat <- loom_theta_mat[seq(7,1),]
loom_phi_mat <- read.csv('loom_center_phi.txt', sep = ',', header = F) / pi * 180
loom_phi_mat <- loom_phi_mat[seq(7,1),]

loom_theta <- melt(t(as.matrix(loom_theta_mat)))$value 
loom_phi <- melt(t(as.matrix(loom_phi_mat)))$value 

xygrid2 <- cbind(loom_phi, loom_theta) 

# shim
shim_xy <- read.csv('shim.csv', sep = ',', header = F)
shim_xy <- as.matrix(shim_xy)
shim_xy <- t(shim_xy) / pi * 180
shim_xy[,2] <- - shim_xy[,2] + 90
shim_xy <- shim_xy[shim_xy[,1] > min(loom_phi)-4.5 & shim_xy[,1] < max(loom_phi)+4.5 & shim_xy[,2] > min(loom_theta)-4.5 & shim_xy[,2] < max(loom_theta)+4.5, ]
# head(shim_xy,1); tail(shim_xy,1); max(shim_xy[,1]); min(shim_xy[,2])
ptadd <- c(max(shim_xy[,1]), min(shim_xy[,2]))
shim_xy <- rbind(shim_xy, ptadd)
shim_xy <- as.data.frame(shim_xy)
colnames(shim_xy) <- c('x','y')


# construct data frame
expBi_df <- data.frame(xygrid2, as.vector(t(expBi)))
colnames(expBi_df) <- c("x","y","z")
expIpsi_df <- data.frame(xygrid2, as.vector(t(expIpsi)))
colnames(expIpsi_df) <- c("x","y","z")

# # slope of LED panel, cf. S5B
# slope_ex <- (90-83.28)/(117-16)
# lefty <- slope_ex*abs(x2[1]) 
# righty <- slope_ex*abs(tail(x2,1))


# not used --------------------------------------------------------------------------------------------------------

# # plot indiv traces
# mai <- expBi3
# plt_exp <- list()
# for (j in 1:dim(mai)[3]) {
#   amp_max <- quantile(na.omit(c(mai[,,j])), probs = c(0.98)) # normalize to 98%
#   for (k in 1:dim(mai)[1]) {
#     # amp_max <- quantile(na.omit(c(exp_raw[[j]][,,k])), probs = c(0.98)) # normalize to 98% 
#     # exp_raw[[j]][,,k] <- exp_raw[[j]][,,k] / amp_max
#     df <- data.frame(x = seq(1,dim(mai)[2]), y = mai[ind_mai[k],,j]/amp_max)
#     plt_exp[[(j-1)*dim(mai)[1] + k]] <- ggplot(df, aes(x,y)) + 
#       geom_point() 
#       # ylim(0,1)
#       # theme_void() 
#   }
#   
# }
# dev.new()
# plot_grid(plotlist = plt_exp, ncol = 14, labels = seq(1,98))
# 
# # plot mean
# mai_mean <- rowMeans(mai, dims = 2, na.rm = T)
# amp_max <- quantile(na.omit(c(mai_mean)), probs = c(0.98)) # normalize to 98% 
# for (k in 1:dim(mai_mean)[1]) {
#   df <- data.frame(x = seq(1,dim(mai_mean)[2]), y = mai_mean[ind_mai[k],]/amp_max)
#   plt_exp[[k]] <- ggplot(df, aes(x,y)) + 
#     geom_point() +
#     ylim(0,1)
#   # theme_void() 
# }
# dev.new()
# plot_grid(plotlist = plt_exp, ncol = 14, labels = seq(1,98))

# stim_end <- 200 # loom during the first 150 points, 4 sec
# 
# # normalize, mean aross flies, then normalize again
# exp_raw <- list(expBi3, expIpsi3)
# pval <- list()
# j <- 1
#   # normalize within fly
#   for (k in 1:dim(exp_raw[[j]])[3]) {
#     amp_max <- quantile(na.omit(c(exp_raw[[j]][,,k])), probs = c(0.98)) # normalize to 98% 
#     exp_raw[[j]][,,k] <- exp_raw[[j]][,,k] / amp_max
#   }
#   
# # t-test
# mu_test <- sum(exp_raw[[j]][ind_mai[c(12,13,14,27,28,42,56)], 1:stim_end,]) / 7 / dim(exp_raw[[j]])[3]
# pval_tmp <- c()
# for (k in 1:dim(exp_raw[[j]])[1]) {
#   pval_tmp[k] <- t.test(colSums(exp_raw[[j]][k, 1:stim_end,]), mu = mu_test)$p.value
# }
# pval[[j]] <- pval_tmp
# 
#   # mean and normalize across flies
#   mai_mean <- rowMeans(mai, dims = 2, na.rm = T)
#   amp_max <- quantile(na.omit(c(mai_mean)), probs = c(0.98)) # normalize to 98% 
#   
#   max(na.omit(c(mai_mean))) / amp_max
#   min(na.omit(c(mai_mean))) / amp_max


# # -- coutour from data
# # x2 <- seq(0-4.25,117-4.25,by = 9) # -18 to 99, with sim data
# # y2 <- seq(54, 108,by = 9)
# x3 <- seq(-18,99,by = 9) # -18 to 99, with sim data
# y3 <- seq(54, 108,by = 9)
# xygrid3 <- expand.grid(x3, y3)
# 
# exp_df <- data.frame(xygrid3, as.vector(t(expBi)))
# exp_df <- data.frame(xygrid3, as.vector(t(expIpsi)))
# colnames(exp_df) <- c("x","y","z")
# 
# exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
# x4 <- seq(-18, 99, by = 1)
# y4 <- seq(54, 108, by = 1)
# xygrid4 <- expand.grid(x = x4, y = y4)
# loess_fit <- predict (exp_loess, xygrid4, se = T)
# exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
# colnames(exp_interp) <- c('x','y','z')
# 
# dev.new()
# image(x= x4, y= y4, z = loess_fit$fit,  asp = 1)
# points(exp_df)
# 
# cut <- quantile(exp_interp$z, probs = c(0.7))
# # cut <- 0.7 * max(exp_interp$z)
# exp_interp_cut <- exp_interp[exp_interp$z > cut, ]
# 
# dim(exp_interp_cut)[1]



#  Figure 6E -- time series mean + EM contour  --------------------------------------------------------------------

# -- put in contour's coord
# exp data range [-0.1, 1.4] vs [1, 280]/200
# need to reverse y-axis

# # starting locations
# tseries_x0 <- sort(unique(expBi_df[,1])) - 4.5
# tseries_y0 <- sort(unique(expBi_df[,2])) + 4.5


# -- add up the same type
mat_names <- c("bi", "ipsi")
conn_target_agglo_sum <- list()
conn_target_agglo_sum[[1]] <- conn_target[[1]] + conn_target[[2]] + conn_target[[3]] + conn_target[[4]]
conn_target_agglo_sum[[2]] <- conn_target[[5]]+conn_target[[6]]+conn_target[[7]]+conn_target[[8]]+conn_target[[9]]


exp_raw <- list(expBi3, expIpsi3)
exp_raw_mean <- list(expBi, expIpsi)
stim_end <- 250 # loom during the first 150 points, 4 sec
stim_start <- 0

# for exp data contour 
# x4 <- seq(0-4.25,117-4.25,by = 9)
# y4 <- seq(54, 108, by = 1)
# xygrid4 <- expand.grid(x = x4, y = y4)
x4 <- seq(-25, 112,by = 1) # range(loom_phi)
y4 <- seq(58, 120, by = 1) # range(loom_theta)
xygrid4 <- expand.grid(x = x4, y = y4)


ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
plvl <- list()
simdata <- list()
simdata_df <- list()
LC6_tar_median <- list()
fdr <- list()
ct_cut_pt <- list()

for (j in 1:length(conn_target_agglo_sum)) { 
  for (k in 1:dim(exp_raw[[j]])[3]) {
    amp_max <- quantile(na.omit(c(exp_raw[[j]][,,k])), probs = c(0.98)) # normalize to 98% 
    exp_raw[[j]][,,k] <- exp_raw[[j]][,,k] / amp_max
  }
  mai_mean <- rowMeans(exp_raw[[j]], dims = 2, na.rm = T)
  amp_max <- quantile(na.omit(c(mai_mean)), probs = c(0.98)) # normalize to 98% 
  mai_mean <- mai_mean / amp_max
  
  max(na.omit(c(mai_mean))) / amp_max
  min(na.omit(c(mai_mean))) / amp_max
  
  # t-test
  # mu_test <- sum(exp_raw[[j]][ind_mai[c(12,13,14,27,28,42,56)], stim_start:stim_end, ]) / 7 / dim(exp_raw[[j]])[3]
  mu_test <- sum(exp_raw[[j]][ind_mai[c(13,14,28,42)], stim_start:stim_end, ]) / 4 / dim(exp_raw[[j]])[3]
  pval_tmp <- c()
  for (k in 1:dim(exp_raw[[j]])[1]) {
    pval_tmp[k] <- t.test(colSums(exp_raw[[j]][k, stim_start:stim_end,]), mu = mu_test)$p.value
  }
  fdr[[j]] <- p.adjust(pval_tmp, method = 'BH')
  
  exp_df <- data.frame(xygrid2, as.vector(t(exp_raw_mean[[j]])))
  colnames(exp_df) <- c("x","y","z")
  exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
  loess_fit <- predict (exp_loess, xygrid4, se = T)
  exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
  colnames(exp_interp) <- c('x','y','z')
  exp_interp <- exp_interp[!is.na(exp_interp[,'z']), ]
  
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  # ii <- grid_Gaussian[,1] > 16
  # grid_Gaussian[ii,2] <- slope_ex*(grid_Gaussian[ii,1] - 16) + grid_Gaussian[ii,2]
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(neu)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- conn_target_agglo_sum[[j]][k,"tofrom_glu"]
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
    # }
  }
  # grid_Gaussian[ii,2] <- round(- slope_ex*(grid_Gaussian[ii,1] - 16) + grid_Gaussian[ii,2])
  simdata[[j]] <- grid_Gaussian
  simdata_df[[j]] <- simdata[[j]]
  colnames(simdata_df[[j]]) <- c("x","y","z")
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / mean(head(sort(simdata_df[[j]]$z, decreasing = T)))
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  plvl[[j]] <- ggplot()
  for (k in 1:dim(mai_mean)[1]) {
    df <- data.frame(x = seq(1,dim(mai_mean)[2])/200*6 + loom_phi[k]-4.5, y = mai_mean[ind_mai[k],]*(-6) + loom_theta[k]+4.5)
    if (fdr[[j]][ind_mai[k]] >= 0.05 ) {
      plvl[[j]] <- plvl[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', alpha = 0.3, lwd = 1)
    } else {
      plvl[[j]] <- plvl[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
    }
  }
  plvl[[j]] <- plvl[[j]] + labs(title = paste("target", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glu"]),sep = " ")) +
    geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
    geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2) +
    geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61), color = "red", alpha = 1, lwd = 2) +
    geom_segment(aes(x = -12, y = 135, xend = -3, yend = 135), size=2, lineend = "round") +
    annotate("text", x = -3, y = 139, label = "9°") +
    geom_segment(aes(x = -12, y = 135, xend = -12, yend = 124), size=2, lineend = "round") +
    annotate("text", x = -17.5, y = 124, label = "9°") +
    # geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
    # geom_segment(aes(x = 16, y = 90, xend = 117, yend = 90-righty), size=2, lineend = "round") +
    geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
    geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
    geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
    # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
    scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
    scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
    theme_void() +
    coord_fixed(ratio = 1)
  
  # area
  ct_cut_pt[[j]] <- simdata_df[[j]][simdata_df[[j]]$z > 0.7, c('x','y')]
}

windows(record = F, width = 8, height = 8)
plvl[[1]] 
windows(record = F, width = 8, height = 8)
plvl[[2]] 


# # areas of RF
# dim(ct_cut_pt[[1]])[1] / sum(ii_inpoly)
# 
# dim(ct_cut_pt[[2]])[1] / sum(ii_inpoly)


# individual t series ---------------------------------------------------------------------------------------------

# Figure S5C

# color
# getPalette <- colorRampPalette(brewer.pal(9, "RdYlBu"))
# ct_col <- getPalette(7)
ct_col <- brewer.pal(9, "RdYlBu")[c(1,2,3,4,7,8,9)]
  
# -- add up the same type
mat_names <- c("bi_all", "ipsi_all")
conn_target_agglo_sum <- list()
conn_target_agglo_sum[[1]] <- conn_target[[1]] + conn_target[[2]] + conn_target[[3]] + conn_target[[4]]
conn_target_agglo_sum[[2]] <- conn_target[[5]]+conn_target[[6]]+conn_target[[7]]+conn_target[[8]]+conn_target[[9]]

# -- difference in EM contour
simdata_df_diff <- simdata_df[[1]][,1:3]
simdata_df_diff$z <- simdata_df_diff$z - simdata_df[[2]]$z
simdata_df_diff$z <- simdata_df_diff$z / max(abs(simdata_df_diff$z))
simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(0,max(abs(simdata_df_diff$z)),length.out = n_lvl))

exp_raw <- list(expBi3, expIpsi3)
exp_raw_mean <- list(expBi, expIpsi)
stim_end <- 250 # loom during the first 150 points, 4 sec
stim_start <- 0


plvl_ol <- list() # all exp contour
plvl_ol_em <- list() # all exp contour with EM bi - ipsi heat map
plvl_bi <- list()
plvl_ipsi <- list()
bi_ct <- list()
ipsi_ct <- list()
indi_max <- list() #peak value for indivial fly
fdr <- list() # False discovery rate

for (j in 1:length(conn_target_agglo_sum)) { 
  # -- fly mean
  max_tmp <- c()
  for (k in 1:dim(exp_raw[[j]])[3]) {
    amp_max <- quantile(na.omit(c(exp_raw[[j]][,,k])), probs = c(0.98)) # normalize to 98% 
    exp_raw[[j]][,,k] <- exp_raw[[j]][,,k] / amp_max
    max_tmp <- c(max_tmp, round(amp_max,1))
  }
  indi_max[[j]] <- max_tmp
  mai_mean <- rowMeans(exp_raw[[j]], dims = 2, na.rm = T)
  # amp_max <- quantile(na.omit(c(mai_mean)), probs = c(0.98)) # normalize to 98% 
  # mai_mean <- mai_mean / amp_max
  
  # t-test
  # mu_test <- sum(exp_raw[[j]][ind_mai[c(12,13,14,27,28,42,56)], stim_start:stim_end, ]) / 7 / dim(exp_raw[[j]])[3]
  mu_test <- sum(exp_raw[[j]][ind_mai[c(13,14,28,42)], stim_start:stim_end, ]) / 4 / dim(exp_raw[[j]])[3]
  pval_tmp <- c()
  for (k in 1:dim(exp_raw[[j]])[1]) {
    pval_tmp[k] <- t.test(colSums(exp_raw[[j]][k, stim_start:stim_end,]), mu = mu_test)$p.value
  }
  fdr[[j]] <- p.adjust(pval_tmp, method = 'BH')
  
  exp_df <- data.frame(xygrid2, as.vector(t(exp_raw_mean[[j]])))
  colnames(exp_df) <- c("x","y","z")
  exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
  loess_fit <- predict (exp_loess, xygrid4, se = T)
  exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
  colnames(exp_interp) <- c('x','y','z')
  exp_interp <- exp_interp[!is.na(exp_interp[,'z']), ]
  
  plvl_ol[[j]] <- ggplot()
  plvl_ol_em[[j]] <- ggplot() +
    geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
    geom_raster(data = simdata_df_diff, aes(x, y, z = z, fill = equalSpace), interpolate = F) +
    scale_fill_manual(values = pal_tar) +
    geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE)
  for (k in 1:dim(mai_mean)[1]) {
    df <- data.frame(x = seq(1,dim(mai_mean)[2])/200*6 + loom_phi[k]-4.5, y = mai_mean[ind_mai[k],]*(-6) + loom_theta[k]+4.5)
    if (fdr[[j]][ind_mai[k]] >= 0.05 ) {
      plvl_ol[[j]] <- plvl_ol[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', alpha = 0.3, lwd = 1)
      # plvl_ol_em[[j]] <- plvl_ol_em[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', alpha = 0.3, lwd = 1)
    } else {
      plvl_ol[[j]] <- plvl_ol[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
      # plvl_ol_em[[j]] <- plvl_ol_em[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
    }
  }
  plvl_ol[[j]] <- plvl_ol[[j]] + labs(title = paste("Target ipsi exp vs", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glu"]),sep = " ")) +
    geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
    # geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2) +
    # geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61), color = "red", alpha = 1, lwd = 2) +
    # geom_segment(aes(x = -12, y = 115, xend = -3, yend = 115), size=2, lineend = "round") +
    # annotate("text", x = -3, y = 119, label = "9°") +
    # geom_segment(aes(x = -12, y = 115, xend = -12, yend = 106), size=2, lineend = "round") +
    # annotate("text", x = -17.5, y = 105, label = "9°") +
    geom_segment(aes(x = -12, y = 135, xend = -3, yend = 135), size=2, lineend = "round") +
    annotate("text", x = -3, y = 139, label = "9°") +
    geom_segment(aes(x = -12, y = 135, xend = -12, yend = 124), size=2, lineend = "round") +
    annotate("text", x = -17.5, y = 124, label = "9°") +
    geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
    geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
    geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
    # geom_segment(aes(x = 0, y = y2[1]-4.5, xend = 0, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
    # geom_segment(aes(x = 90, y = y2[1]-4.5, xend = 90, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
    # geom_segment(aes(x = x2[1], y = 90-lefty, xend = 0, yend = 90), size=1, lineend = "round") +
    # geom_segment(aes(x = 0, y = 90, xend = tail(x2,1), yend = 90-righty), size=1, lineend = "round") +
    # geom_segment(aes(x = min(loom_phi), y = 90, xend = max(loom_phi), yend = 90), size=1, lineend = "round") +
    scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
    scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
    theme_void() +
    coord_fixed(ratio = 1)
  
  # for (j in 1:length(conn_target_agglo_sum)) { 
    
  # -- fly indiv
  for (k in 1:dim(exp_raw[[j]])[3]) {
    exp_df <- data.frame(xygrid2, rowSums(exp_raw[[j]][ind_mai, stim_start:stim_end, k]))
    colnames(exp_df) <- c("x","y","z")
    exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
    loess_fit <- predict (exp_loess, xygrid4, se = T)
    exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
    colnames(exp_interp) <- c('x','y','z')
    exp_interp <- exp_interp[!is.na(exp_interp[,'z']), ]
    
    if (j == 1) {
      plvl_bi[[k]] <- ggplot()
      for (m in 1:dim(exp_raw[[j]])[1]) {
        df <- data.frame(x = seq(1,dim(exp_raw[[j]])[2])/200*6 + loom_phi[m]-4.5, y = exp_raw[[j]][ind_mai[m],,k]*(-6) + loom_theta[m]+4.5)
        plvl_bi[[k]] <- plvl_bi[[k]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
      }
      bi_ct[[k]] <- exp_interp
      
      plvl_bi[[k]] <- plvl_bi[[k]] + labs(title = paste(k, ", bi, 60%" ,sep = " ")) +
        geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
        # geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2) +
        geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[2*(k-1)+1], alpha = 1, lwd = 2) +
        # geom_segment(aes(x = -12, y = 115, xend = -3, yend = 115), size=2, lineend = "round") +
        # annotate("text", x = -3, y = 119, label = "9°") +
        # geom_segment(aes(x = -12, y = 115, xend = -12, yend = 106), size=2, lineend = "round") +
        # annotate("text", x = -17.5, y = 105, label = "9°") +
        geom_segment(aes(x = loom_phi[1]-4.5, y = loom_theta[1]+4.5, xend = loom_phi[1]-4.5, yend = loom_theta[1]-9+4.5), size=2, lineend = "round") +
        annotate("text", x = loom_phi[1]-4.5-4.5, y = loom_theta[1]-4.5+4.5, label = "ΔF/F") +
        annotate("text", x = loom_phi[1]-4.5-4.5, y = loom_theta[1]+4.5, label = paste(indi_max[[j]][k])) +
        geom_segment(aes(x = loom_phi[1]-4.5, y = loom_theta[1]+4.5, xend = loom_phi[1]+9-4.5, yend = loom_theta[1]+4.5), size=2, lineend = "round") +
        annotate("text", x = loom_phi[1]+4.5-4.5, y = loom_theta[1]+4+4.5, label = "5s") +
        geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
        geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
        geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
        # geom_segment(aes(x = 0, y = y2[1]-4.5, xend = 0, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
        # geom_segment(aes(x = 90, y = y2[1]-4.5, xend = 90, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
        # geom_segment(aes(x = x2[1], y = 90-lefty, xend = 0, yend = 90), size=1, lineend = "round") +
        # geom_segment(aes(x = 0, y = 90, xend = tail(x2,1), yend = 90-righty), size=1, lineend = "round") +
        # geom_segment(aes(x = min(loom_phi), y = 90, xend = max(loom_phi), yend = 90), size=1, lineend = "round") +
        scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
        scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
        theme_void() +
        coord_fixed(ratio = 1)
      
      # if (k == 1) {
      #   plvl_bi[[k]] <- plvl_bi[[k]] + 
      #     geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[1,14], yend = loom_theta_mat[1,14]), size = 1, lineend = "round", col = 'blue') + # exp
      #     geom_segment(aes(x = loom_phi_mat[7,1], y = loom_theta_mat[7,1], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'blue') +
      #     geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[7,1], yend = loom_theta_mat[7,1]), size = 1, lineend = "round", col = 'blue') + 
      #     geom_segment(aes(x = loom_phi_mat[1,14], y = loom_theta_mat[1,14], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'blue')  
      # }
      
    } else {
      plvl_ipsi[[k]] <- ggplot()
      for (m in 1:dim(exp_raw[[j]])[1]) {
        df <- data.frame(x = seq(1,dim(exp_raw[[j]])[2])/200*6 + loom_phi[m]-4.5, y = exp_raw[[j]][ind_mai[m],,k]*(-6) + loom_theta[m]+4.5)
        plvl_ipsi[[k]] <- plvl_ipsi[[k]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
      }
      ipsi_ct[[k]] <- exp_interp
      
      plvl_ipsi[[k]] <- plvl_ipsi[[k]] + labs(title = paste(k, ", ipsi, 60%" ,sep = " ")) +
        geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
        # geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = ct_col[k], alpha = 0, lwd = 2) +
        geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[k], alpha = 1, lwd = 2) +
        geom_segment(aes(x = loom_phi[1]-4.5, y = loom_theta[1]+4.5, xend = loom_phi[1]-4.5, yend = loom_theta[1]-9+4.5), size=2, lineend = "round") +
        annotate("text", x = loom_phi[1]-4.5-4.5, y = loom_theta[1]-4.5+4.5, label = "ΔF/F") +
        annotate("text", x = loom_phi[1]-4.5-4.5, y = loom_theta[1]+4.5, label = paste(indi_max[[j]][k])) +
        geom_segment(aes(x = loom_phi[1]-4.5, y = loom_theta[1]+4.5, xend = loom_phi[1]+9-4.5, yend = loom_theta[1]+4.5), size=2, lineend = "round") +
        annotate("text", x = loom_phi[1]+4.5-4.5, y = loom_theta[1]+4+4.5, label = "5s") +
        geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
        geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
        geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
        # geom_segment(aes(x = 0, y = y2[1]-4.5, xend = 0, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
        # geom_segment(aes(x = 90, y = y2[1]-4.5, xend = 90, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
        # geom_segment(aes(x = x2[1], y = 90-lefty, xend = 0, yend = 90), size=1, lineend = "round") +
        # geom_segment(aes(x = 0, y = 90, xend = tail(x2,1), yend = 90-righty), size=1, lineend = "round") +
        # geom_segment(aes(x = min(loom_phi), y = 90, xend = max(loom_phi), yend = 90), size=1, lineend = "round") +
        scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
        scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
        theme_void() +
        coord_fixed(ratio = 1)
    }
    if (j == 1) {
      plvl_ol[[j]] <- plvl_ol[[j]] + geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[2*(k-1)+1], alpha = 1, lwd = 2)
      plvl_ol_em[[j]] <- plvl_ol_em[[j]] + geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[2*(k-1)+1], alpha = 1, lwd = 2)
    } else {
      plvl_ol[[j]] <- plvl_ol[[j]] + geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[k], alpha = 1, lwd = 2)
      plvl_ol_em[[j]] <- plvl_ol_em[[j]] + geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[k], alpha = 1, lwd = 2)
    }
  }
}
  
  
  
# Figure 3B and S5C
windows(record = F, width = 8, height = 8)
plvl_bi[[1]] 
windows(record = F, width = 8, height = 8)
plvl_bi[[2]] 
windows(record = F, width = 8, height = 8)
plvl_bi[[3]] 
windows(record = F, width = 8, height = 8)
plvl_bi[[4]] 

windows(record = F, width = 8, height = 8)
plvl_ipsi[[1]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[2]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[3]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[4]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[5]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[6]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[7]]

# Figure 3C
windows(record = F, width = 8, height = 8)
plvl_ol[[1]] 
windows(record = F, width = 8, height = 8)
plvl_ol[[2]] 

# Figure SX
windows(record = F, width = 8, height = 8)
plvl_ol_em[[1]] 
windows(record = F, width = 8, height = 8)
plvl_ol_em[[2]] 

# Figure 6D
# contour map + exp contour
windows(record = F, width = 8, height = 8)
plvl_ol[[1]] + geom_contour(data = simdata_df[[1]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)

windows(record = F, width = 8, height = 8)
plvl_ol[[2]] + geom_contour(data = simdata_df[[2]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)



# -- exp and EM overlap
# polygon from EM
simdata_ct <- list()
simdata_poly <- list()
for (j in 1:2) {
  simdata_ct[[j]] <- simdata_df[[j]][simdata_df[[j]]$z > 0.71, 1:2]
  ash <- ashape(simdata_ct[[j]]+matrix(runif(dim(simdata_ct[[j]])[1]*2, 1e-9, 2e-9), ncol = 2), alpha = 5)
  simdata_poly[[j]] <- mkpoly(ash$edges)[[1]][,3:4]
  # dev.new()
  # polygon(simdata_poly[[j]])
}


# polygon from exp
exp_raw <- list(expBi3, expIpsi3)
stim_end <- 200 # loom during the first 150 points, 4 sec
stim_start <- 50
ol_ratio <- matrix(ncol = 2, nrow = 2*(4+7)) # percentage of overlap = intersection / union, bi-biEM, ipsi-biEM, etc

NN <- 1
for (jem in 1:2) {
  for (j in 1:2) {
    for (k in 1:dim(exp_raw[[j]])[3]) {
      exp_df <- data.frame(xygrid2, rowSums(exp_raw[[j]][ind_mai, stim_start:stim_end, k]))
      colnames(exp_df) <- c("x","y","z")
      exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
      loess_fit <- predict (exp_loess, xygrid4, se = T)
      exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
      colnames(exp_interp) <- c('x','y','z')
      exp_interp <- exp_interp[!is.na(exp_interp[,'z']), ]
      
      exp_interp_ct <- exp_interp[(exp_interp$z > 0.61*max(exp_interp$z)),]
      
      N_in <- sp::point.in.polygon(exp_interp_ct[,1], exp_interp_ct[,2], simdata_poly[[jem]][,1], simdata_poly[[jem]][,2])
      N_in <- sum(N_in)
      N_out <- dim(exp_interp_ct)[1] - N_in
      N_em <- sum(simdata_ct[[jem]]$y >= 54 | simdata_ct[[jem]]$y <= 108)
      ol_ratio[NN,] <- c(j + 2*(jem %/% 2), N_in / (N_out + N_em))
      NN <- NN + 1
    }
  }
}

ol_ratio <- as.data.frame(ol_ratio)
colnames(ol_ratio) <- c('pair', 'ratio')

library(Hmisc)
windows(record = F, width = 8, height = 8)
ggplot(ol_ratio, aes(x = pair, group = pair, y = ratio)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 1) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4), labels=c("bi-biEM", "ipsi-biEM", "bi-ipsiEM", "ipsi-ipsiEM"))+
  theme_bw()

# ggplot(ol_ratio, aes(x = pair, group = pair, y = ratio)) + geom_boxplot()+
#   scale_x_continuous(breaks=c(1,2,3,4), labels=c("bi-biEM", "ipsi-biEM", "bi-ipsiEM", "ipsi-ipsiEM")) +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) 
  


# -- difference in EM contour
simdata_df_diff <- simdata_df[[1]][,1:3]
simdata_df_diff$z <- simdata_df_diff$z - simdata_df[[2]]$z
simdata_df_diff$z <- simdata_df_diff$z / max(abs(simdata_df_diff$z))
# simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(0,max(abs(simdata_df_diff$z)),length.out = n_lvl))
# simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(min(simdata_df_diff$z),max(simdata_df_diff$z),length.out = n_lvl))
simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(-1,1,length.out = n_lvl))

# contour map
windows(record = F, width = 8, height = 8)
ggplot(simdata_df_diff, aes(x, y, z = z)) +
  geom_raster(aes(fill = equalSpace), interpolate = F) +
  scale_fill_manual(values = pal_tar) +
  geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +  
  plvl_ol[[1]]

ggplot(simdata_df_diff, aes(x, y, z = z)) +
  geom_raster(aes(fill = equalSpace), interpolate = F) +
  scale_fill_manual(values = pal_tar) +
  geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
  geom_segment(aes(x = 10, y = 180, xend = 19, yend = 180), size=2, lineend = "round") +
  annotate("text", x = 14.5, y = 185, label = "9°") +
  geom_segment(aes(x = 10, y = 180, xend = 10, yend = 171), size=2, lineend = "round") +
  annotate("text", x = 5, y = 175.5, label = "9°") +
  # geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
  # geom_segment(aes(x = 16, y = 90, xend = 180, yend = 79.1), size=2, lineend = "round") + #90-slope_ex*(180-16)
  geom_segment(aes(x = min(loom_phi), y = 90, xend = 160, yend = 90), size=1, lineend = "round") +
  # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 180), size = 2, lineend = "round") +
  geom_segment(aes(x = 90, y = 0, xend = 90, yend = 180), size = 2, lineend = "round") +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  theme_void() +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[1,14], yend = loom_theta_mat[1,14]), size = 1, lineend = "round", col = 'red') + # exp
  geom_segment(aes(x = loom_phi_mat[7,1], y = loom_theta_mat[7,1], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[7,1], yend = loom_theta_mat[7,1]), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,14], y = loom_theta_mat[1,14], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'red')
  labs(title = "bi - ipsi") +
  coord_fixed(ratio = 1)


windows(record = F, width = 8, height = 8)
ggplot(simdata_df_diff, aes(x, y, z = z)) +
  geom_raster(aes(fill = equalSpace), interpolate = F) +
  scale_fill_manual(values = pal_tar) +
  geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
  geom_segment(aes(x = 10, y = 190, xend = 19, yend = 190), size=2, lineend = "round") +
  annotate("text", x = 14.5, y = 195, label = "9°") +
  geom_segment(aes(x = 10, y = 190, xend = 10, yend = 181), size=2, lineend = "round") +
  annotate("text", x = 5, y = 185.5, label = "9°") +
  geom_segment(aes(x = min(loom_phi), y = 90, xend = 160, yend = 90), size=1, lineend = "round") +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 180), size = 2, lineend = "round") +
  geom_segment(aes(x = 90, y = 0, xend = 90, yend = 180), size = 2, lineend = "round") +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  theme_void() +
  labs(title = "bi - ipsi and bi contour") +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[1,14], yend = loom_theta_mat[1,14]), size = 1, lineend = "round", col = 'red') + # exp
  geom_segment(aes(x = loom_phi_mat[7,1], y = loom_theta_mat[7,1], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[7,1], yend = loom_theta_mat[7,1]), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,14], y = loom_theta_mat[1,14], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'red') +
  geom_contour(data = bi_ct[[1]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[1], alpha = 1, lwd = 2) +
  geom_contour(data = bi_ct[[2]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[3], alpha = 1, lwd = 2) +
  geom_contour(data = bi_ct[[3]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[5], alpha = 1, lwd = 2) +
  geom_contour(data = bi_ct[[4]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[7], alpha = 1, lwd = 2) +
  coord_fixed(ratio = 1)

windows(record = F, width = 8, height = 8)
ggplot(simdata_df_diff, aes(x, y, z = z)) +
  geom_raster(aes(fill = equalSpace), interpolate = F) +
  scale_fill_manual(values = pal_tar) +
  geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
  geom_segment(aes(x = 10, y = 190, xend = 19, yend = 190), size=2, lineend = "round") +
  annotate("text", x = 14.5, y = 195, label = "9°") +
  geom_segment(aes(x = 10, y = 190, xend = 10, yend = 181), size=2, lineend = "round") +
  annotate("text", x = 5, y = 185.5, label = "9°") +
  geom_segment(aes(x = min(loom_phi), y = 90, xend = 160, yend = 90), size=1, lineend = "round") +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 180), size = 2, lineend = "round") +
  geom_segment(aes(x = 90, y = 0, xend = 90, yend = 180), size = 2, lineend = "round") +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  theme_void() +
  labs(title = "bi - ipsi and ipsi contour") +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[1,14], yend = loom_theta_mat[1,14]), size = 1, lineend = "round", col = 'red') + # exp
  geom_segment(aes(x = loom_phi_mat[7,1], y = loom_theta_mat[7,1], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[7,1], yend = loom_theta_mat[7,1]), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,14], y = loom_theta_mat[1,14], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'red') +
  geom_contour(data = ipsi_ct[[1]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[1], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[2]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[2], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[3]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[3], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[4]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[4], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[5]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[5], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[6]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[6], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[7]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[7], alpha = 1, lwd = 2) +
  coord_fixed(ratio = 1)




# user ------------------------------------------------------------------------------------------------------------

user_stat <- catmaid_get_contributor_stats(c(LC6_skid, ipsi_skid, biL_skid, biR_skid))

user_stat$node_contributors[order(user_stat$node_contributors$n),]
user_stat$pre_contributors[order(user_stat$pre_contributors$n), ]
user_stat$post_contributors[order(user_stat$post_contributors$n), ]
user_stat$review_contributors[order(user_stat$review_contributors$n), ]

udf <- catmaid_get_user_list()

# make a df of contributions 
jd <- full_join(user_stat$post_contributors, user_stat$pre_contributors, by = "id")
jd <- full_join(jd, user_stat$review_contributors, by = 'id')
jd <- full_join(jd, user_stat$node_contributors, by = 'id')
jd[is.na(jd)] <- 0
jd <- left_join(jd, udf) %>%
  # mutate(n = n.x + n.y + n.x.x + n.y.y) %>%
  mutate(n = n.x + n.y ) %>%
  select(id, n, full_name) %>%
  arrange(desc(n)) 

