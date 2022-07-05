
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(vegan)

#ggplot theme
mytheme <- theme(axis.line = element_line(linetype = "solid"), 
                 axis.line.x = element_line(color = "black", size = 0.7), 
                 axis.line.y = element_line(color = "black", size = 0.7), 
                 axis.ticks = element_line(colour = "black"), 
                 panel.grid.major = element_line(linetype = "blank"), 
                 panel.background = element_rect(fill = NA), 
                 legend.key = element_rect(fill = NA), 
                 legend.background = element_rect(fill = NA), 
                 legend.text=element_text(size=13),
                 legend.title=element_text(size=13),
                 #legend.position = "bottom", 
                 #legend.direction = "horizontal",
                 axis.title=element_text(size=17,face="bold"),
                 axis.text.x = element_text( size = 15),
                 axis.text.y = element_text( size = 15))

path_results <- "path/to/your/result/folder"

data(dune)
data("dune.env")


# NMDS
dune.mds <- metaMDS(dune,
                distance = "bray",
                trymax = 1000)
#extract site NMDS scores (x and y coordinates)
site.scores <- as.data.frame(scores(nmds))

#fit species and environmental varables
dune.envfit <- envfit(dune.mds, dune.env, permutations = 9999)
dune.spp.fit <- envfit(dune.mds, dune, permutations = 9999)

#adding environmental variables to dataframe for ploting
site.scores$Management <- dune.env$Management
site.scores$Land_use <- dune.env$Use
site.scores$site <- paste("site",
                          rownames(site.scores))
head(site.scores)

#creating species scores dataframe
spp.scrs <- as.data.frame(scores(dune.spp.fit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs <- cbind(spp.scrs, pval = dune.spp.fit$vectors$pvals)
#subset data to show species significant at 0.05
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) 

head(spp.scrs)

#
# function for ellipses
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the management factor
df_ell <- data.frame() #sets up a data frame before running the function.
for(g in levels(as.factor(site.scores$Land_use))){
  df_ell <- rbind(df_ell,
                  cbind(as.data.frame(with(data.scores[site.scores$Land_use==g, ],
                        veganCovEllipse(cov.wt(cbind(NMDS1, NMDS2),
                                                wt=rep(1/length(NMDS1),
                                                length(NMDS1)))$cov,
                                                center=c(mean(NMDS1),
                                                          mean(NMDS2))))),
                                                Land_use=g))
}


# data for labelling the ellipse
NMDS.mean <- aggregate(site.scores[,c("NMDS1", "NMDS2")],
                    list(group = site.scores$Land_use), mean)

pal <- wes_palette("Darjeeling2",
                    length(site.scores$Land_use),
                    type = "continuous")

pdf(file = file.path(path_results, "NMDS.pdf"),
    width = 12, # The width of the plot in inches
    height = 8) # The height of the plot in inches
ggplot(site.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 4, aes(shape = Land_use, colour = site)) +
  geom_segment(data = spp.scrs,
                aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                arrow = arrow(length = unit(0.25, "cm")),
                colour = "grey10", lwd = 0.3) +
  ggrepel::geom_text_repel(data = spp.scrs,
                            aes(x = NMDS1, y = NMDS2, label = Species),
                            cex = 3,
                            direction = "both",
                            segment.size = 0.25) +
  geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, group = Land_use)) +
  labs(x = "NMDS1", colour = "Sites", y = "NMDS2", shape = "Land use") +
  scale_color_manual(values = pal) +
  mytheme
dev.off()
