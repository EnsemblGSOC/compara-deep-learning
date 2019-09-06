###Mateus plots

##Libraries
library(data.table)
library(ggplot2)
library(plyr)
library(UpSetR)

##Read data
wd <- "./"
#uncompress predicions.tar.gz first
files_path <- file.path(wd, "predictions_pfam_model_composition_1")
files <- list.files(files_path, full.names = T, pattern = ".txt")

all <- do.call(rbind, lapply(files, function(f){
  goc <- strsplit(basename(f), "_")[[1]][6]
  identity <- strsplit(basename(f), "_")[[1]][7]
  homology_type <- strsplit(basename(f), "_")[[1]][8]
  
  df <- fread(f)
  
  data.table("goc" = goc
             , "identity" = identity
             , "homology_type" = homology_type
             , "accuracy" = nrow(subset(df, V6 == 1))/10000
             , "wrong" = nrow(subset(df, V6 == 0))
             , "right" = nrow(subset(df, V6 == 1))
             , "error" = nrow(subset(df, V6 == "NaN")
             ))
  
}))

##Subset
all2 <- subset(all, select = c(goc, identity, homology_type, accuracy))

##Dataframe for plots
dfp <- all2

##Plot 1 - Overall accuracy
dfp$coll <- paste(dfp$goc, dfp$identity, dfp$homology_type, sep = "-")
dfp$coll <- factor(dfp$coll, levels = unique(dfp[order(dfp$accuracy),]$coll))

p1 <- ggplot(dfp, aes(x = coll, y = accuracy)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

##Plot 2 - Violin accuracy all
p2 <- ggplot(dfp, aes(x = "", y = accuracy)) +
  geom_violin() +
  ylim(0, 1) +
  geom_point() +
  xlab("Observations") +
  ylab("Accuracy\n") +
  theme(axis.title = element_text(size = 16)
        , axis.text = element_text(size = 14))
p2

##Plot 3 - Violin accuracy homology type
p3 <- ggplot(subset(dfp, homology_type != "samples"), aes(x = homology_type, y = accuracy)) +
  geom_violin(width = 1.2) +
  geom_point() +
  ylim(0, 1) +
  xlab("\nHomology type") +
  ylab("Accuracy\n") +
  theme(axis.title = element_text(size = 16)
        , axis.text = element_text(size = 14))
p3

##Plot 4 - Violin accuracy homology type facet by goc
# p4 <- ggplot(dfp, aes(x = homology_type, y = accuracy)) +
#   geom_boxplot() +
#   geom_point() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   facet_wrap(~goc, scales = "free_x")
# p4

##Plot 5 - Facet by homology type coloured by goc
# dfp$identity <- as.numeric(dfp$identity)
# 
# p5 <- ggplot(subset(dfp, goc != "models"), aes(x = identity, y = accuracy)) +
#   geom_point(aes(color = goc)) +
#   scale_x_continuous(breaks = c(25,50,75,100)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   facet_wrap(~homology_type, scales = "free_x")
# p5

##Plot 6 - Test
dfp$identity <- factor(dfp$identity, levels = c(25,50,75,100))
dfp$goc <- factor(dfp$goc, levels = c(100,75,50,25,0,"nan"))

p6 <- ggplot(subset(dfp, homology_type != "model" & identity != "samples")
             # ggplot(subset(dfp, homology_type == "many2many")
             , aes(x = identity, y = accuracy)) +
  geom_boxplot() +
  geom_point(aes(color = goc)
             # , width = .1
             , size = 3
             # , stroke = 1.5
             # , shape = 21
             # , alpha = 0.8
  ) +
  geom_point(aes(color = goc)
             , color = "black"
             , size = 3
             , shape = 21
  ) +
  scale_color_manual(name = "GOC"
                     , values = c("#990000"
                                  ,"#CC0000"
                                  ,"#FF0000"
                                  ,"#FF8000"
                                  ,"#FFB266"
                                  ,"#FFFF66")) +
  ylim(c(0,1)) +
  facet_wrap(~homology_type, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
        , axis.title = element_text(size = 16)
        , axis.text = element_text(size = 14)
        , strip.text.x = element_text(size = 12)) +
  xlab("Identity") +
  ylab("Accuracy\n")

p6

##Plot 7 - heatmap
dfp$goc <- factor(dfp$goc, levels = rev(c(100,75,50,25,0,"nan")))
p7 <- ggplot(subset(dfp, homology_type != "model" & identity != "samples")
             , aes(x = identity, y = goc)) +
  geom_tile(aes(fill = accuracy)) +
  facet_wrap(~homology_type, scales = "free_x") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)
        , axis.text = element_text(size = 14)
        , strip.text.x = element_text(size = 12)
        # , panel.border = element_blank()
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
  ) +
  xlab("Identity") +
  ylab("GOC\n") +
  scale_fill_gradient(name = "Accuracy", low = "yellow", high = "red", limits = c(0,1))
p7

##Saving plots
ggsave(plot = p1
       , filename = "distribution.png"
       , path = file.path(wd, "plots")
       , width = 8
       , height = 6
       )

ggsave(plot = p2
       , filename = "violin.png"
       , path = file.path(wd, "plots")
       , width = 8
       , height = 6
)

ggsave(plot = p3
       , filename = "violin_homology.png"
       , path = file.path(wd, "plots")
       , width = 8
       , height = 6
)

ggsave(plot = p6
       , filename = "boxplot.png"
       , path = file.path(wd, "plots")
       , width = 8
       , height = 6
)

ggsave(plot = p7
       , filename = "heatmap.png"
       , path = file.path(wd, "plots")
       , width = 8
       , height = 6
)

