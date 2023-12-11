### Morphological diversity in true and false crabs reveals a common middle 
### ground – the megalopa phase

### Florian Braig, Carolin Haug, Joachim T. Haug.

## Start workflow --------------------------------------------------------------

library(Momocs)
library(dispRity)
library(ggplot2)
library(RColorBrewer)
library(ape)

set.seed(1234)

## MOMOCS OUTLINE ANALYSIS -----------------------------------------------------

lf <- list.files("D:/Documents/AnovsBra_datset_december23", 
                 full.names = TRUE)
coo <- import_jpg(lf, auto.notcentered = TRUE, fun.notcentered = NULL, 
                  threshold = 0.8) 
avb_all <- "App_2_Datalist"
# Import data list from supplemental material however you like into the 
# R session
df <- data.frame(dev = as.factor(avb_all$developmental_phase),
                 major = as.factor(avb_all$general_group),
                 no = as.factor(avb_all$no),
                 group = as.factor(avb_all$major_out),
                 fine_group = as.factor(avb_all$major_group))
df$indi <- as.factor(paste0(df$dev, df$major))

anobra <- Out(coo, fac = df)

panel(anobra, names = TRUE, cex.names = 0.5)


# ELLIPTIC FOURIER ANALYSIS AND PCA OF COMPLETE DATA SET------------------------

anobra <- coo_center(anobra)
anobra <- coo_scale(anobra)
anobra.l <- coo_slidedirection(anobra, direction = "left")
stack(anobra.l, first.point = TRUE)
inspect(anobra.l)

calibrate_harmonicpower_efourier(anobra.l, id = 1:1567, 100, drop = 1, 
                                 thresh = c(95, 99, 99.9))
calibrate_deviations_efourier(anobra.l, range = c(3, 6, 11, 24), 
                              norm.centsize = TRUE, dist.method = edm_nearest, 
                              interpolate.factor = 1, dist.nbpts = 60)

anobra.e <- efourier(anobra.l, nb.h = 11, norm = FALSE, start = TRUE)
boxplot(anobra.e, drop = 1)
anobra.pca <- PCA(anobra.e)
scree(anobra.pca, nax = c(1:20))
anobra.scores <- anobra.pca %>% as_df(17)
PCcontrib(anobra.pca, nax = c(1:10))

## PLOTS OF THE MORPHOSPACE ----------------------------------------------------

#Plot with shapes in background and polygons with alpha
plot(anobra.pca, fac = "dev", col = c("orange", "cyan", "purple", "red"), 
     fill = c("orange", "cyan", "purple", "red"))

#Same plot but with differnt colours
plot(anobra.pca, fac = "dev")

# Plot with every data point having its respective shape
plot_PCA(anobra.pca, ~ dev, points = TRUE, chull = TRUE, chullfilled = FALSE,
         axes = c(1, 2), zoom = 1, morphospace_position = "xy",
         palette = pal_manual(c('#1f78b4', '#33a02c', '#e31a1c', '#ff7f00')))

#Plot with shapes in background and polygons
plot_PCA(anobra.pca, ~ major, points = TRUE, chull = TRUE, chullfilled = FALSE,
         axes = c(1, 2), zoom = 0.9, 
         palette = pal_manual(c('#1f78b4', '#33a02c', '#e31a1c')))

#Plot with labels
ggplot(anobra.scores, aes(x = PC1, y = PC2, color = dev, label = no)) + 
  geom_point() +
  geom_text() +
  coord_fixed (ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
  theme_gray()

#Plot development by colour, outgroups by shape, FOR PUBLICATION
ggplot(anobra.scores, aes(x = PC1, y = PC2, fill = dev, shape = major)) + 
  geom_point(size = 15, color = "black") +
  scale_fill_manual(values = c("#d7191c", "#ffd700", "#abd9e9", "#2c7bb6"),
                    guide = guide_legend(override.aes = list(color =
                                                               c("#d7191c", "#ffd700", "#abd9e9", "#2c7bb6")))) +
  scale_shape_manual(values = c(23, 24, 21)) +
  coord_fixed (ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
  theme_void(base_size = 60)

# Plot for supplemental material
ggplot(anobra.scores, aes(x = PC1, y = PC2, fill = group, shape = dev)) + 
  geom_point(size = 15, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#d7191c", "#abd9e9", "#1b9e77", "#d95f02", "#7570b3", "#d7191c", "#abd9e9", "#1b9e77", "#d95f02", "#7570b3", "#d7191c", "#abd9e9")) +
  scale_shape_manual(values = c(23, 24, 21, 22)) +
  facet_wrap(~ group) +
  coord_fixed (ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
  theme_light(base_size = 60)


## Values for Figure 2----------------------------------------------------------
# gives occupation ranges of PC1

#ALL IN ONE PLOT
df_box <- data.frame(PC1 = anobra.scores$PC1, group = anobra.scores$group, 
                     dev = anobra.scores$dev)
ggplot(df_box, aes(x = PC1, y = group, fill = dev)) + 
  geom_boxplot(color = "black", outlier.shape = NA) + 
  scale_fill_manual(values = c("#d7191c", "#ffa500", "#ffd700", 
                               "#abd9e9", "#2c7bb6")) +
  theme_classic(base_size = 30)


# Carcinidae

max(anobra.scores$PC1[anobra.scores$group == "Carcinidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$group == "Carcinidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$group == "Carcinidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$group == "Carcinidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$group == "Carcinidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$group == "Carcinidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$group == "Carcinidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$group == "Carcinidae" & anobra.scores$dev == "zoea"])

# Dromiidae

max(anobra.scores$PC1[anobra.scores$group == "Dromiidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$group == "Dromiidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$group == "Dromiidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$group == "Dromiidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$group == "Dromiidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$group == "Dromiidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$group == "Dromiidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$group == "Dromiidae" & anobra.scores$dev == "zoea"])


# Diogenidae

max(anobra.scores$PC1[anobra.scores$fine_group == "Diogenidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Diogenidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Diogenidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Diogenidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Diogenidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Diogenidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Diogenidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Diogenidae" & anobra.scores$dev == "zoea"])


# Earliest

max(anobra.scores$PC1[anobra.scores$group == "Earliest" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$group == "Earliest" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$group == "Earliest" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$group == "Earliest" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$group == "Earliest" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$group == "Earliest" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$group == "Earliest" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$group == "Earliest" & anobra.scores$dev == "zoea"])


# Galatheidae

max(anobra.scores$PC1[anobra.scores$group == "Galatheidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$group == "Galatheidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$group == "Galatheidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$group == "Galatheidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$group == "Galatheidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$group == "Galatheidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$group == "Galatheidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$group == "Galatheidae" & anobra.scores$dev == "zoea"])


# Grapsidae

max(anobra.scores$PC1[anobra.scores$group == "Grapsidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$group == "Grapsidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$group == "Grapsidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$group == "Grapsidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$group == "Grapsidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$group == "Grapsidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$group == "Grapsidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$group == "Grapsidae" & anobra.scores$dev == "zoea"])


# Hippidae

max(anobra.scores$PC1[anobra.scores$fine_group == "Hippidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Hippidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Hippidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Hippidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Hippidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Hippidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Hippidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Hippidae" & anobra.scores$dev == "zoea"])


# Lithodidae

max(anobra.scores$PC1[anobra.scores$group == "Lithodidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$group == "Lithodidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$group == "Lithodidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$group == "Lithodidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$group == "Lithodidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$group == "Lithodidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$group == "Lithodidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$group == "Lithodidae" & anobra.scores$dev == "zoea"])


# Nephropidae

max(anobra.scores$PC1[anobra.scores$fine_group == "Nephropidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Nephropidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Nephropidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Nephropidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Nephropidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Nephropidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Nephropidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Nephropidae" & anobra.scores$dev == "zoea"])


# Paguridae

max(anobra.scores$PC1[anobra.scores$fine_group == "Paguridae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Paguridae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Paguridae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Paguridae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Paguridae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Paguridae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Paguridae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Paguridae" & anobra.scores$dev == "zoea"])


# Porcellanidae

max(anobra.scores$PC1[anobra.scores$group == "Porcellanidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$group == "Porcellanidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$group == "Porcellanidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$group == "Porcellanidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$group == "Porcellanidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$group == "Porcellanidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$group == "Porcellanidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$group == "Porcellanidae" & anobra.scores$dev == "zoea"])


# Portunidae

max(anobra.scores$PC1[anobra.scores$group == "Portunidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$group == "Portunidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$group == "Portunidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$group == "Portunidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$group == "Portunidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$group == "Portunidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$group == "Portunidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$group == "Portunidae" & anobra.scores$dev == "zoea"])


# Raninidae

max(anobra.scores$PC1[anobra.scores$group == "Raninidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$group == "Raninidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$group == "Raninidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$group == "Raninidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$group == "Raninidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$group == "Raninidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$group == "Raninidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$group == "Raninidae" & anobra.scores$dev == "zoea"])

# Xanthidae

max(anobra.scores$PC1[anobra.scores$fine_group == "Xanthidae" & anobra.scores$dev == "adult"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Xanthidae" & anobra.scores$dev == "adult"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Xanthidae" & anobra.scores$dev == "juvenile"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Xanthidae" & anobra.scores$dev == "juvenile"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Xanthidae" & anobra.scores$dev == "megalopa"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Xanthidae" & anobra.scores$dev == "megalopa"])
max(anobra.scores$PC1[anobra.scores$fine_group == "Xanthidae" & anobra.scores$dev == "zoea"])
min(anobra.scores$PC1[anobra.scores$fine_group == "Xanthidae" & anobra.scores$dev == "zoea"])


## BOXPLOTS OF SUM VARIANCE VALUES ---------------------------------------------
# Gives ranges for average displacements for figure 2

# Carcinidae
data.car <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$group == "Carcinidae"))
data.car <- droplevels(data.car)
data <- as.data.frame(data.car[, 7:23])
group.v <- data.frame(data.car$dev)
rownames(data) <- 1:nrow(data)
car_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(car_subsets))
car_boot <- boot.matrix(data = car_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
car_pos <- dispRity(data = car_boot, metric = c(mean, displacements))
summary(car_pos)



# Dromiidae
data.dro <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$group == "Dromiidae"))
data.dro <- droplevels(data.dro)
data <- as.data.frame(data.dro[, 7:23])
group.v <- data.frame(data.dro$dev)
rownames(data) <- 1:nrow(data)
dro_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(dro_subsets))
dro_boot <- boot.matrix(data = dro_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
dro_pos <- dispRity(data = dro_boot, metric = c(mean, displacements))
summary(dro_pos)

# Earliest NO SUBSETS
data.ear <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$group == "Earliest"))
data.ear <- droplevels(data.ear)
data <- as.data.frame(data.ear[, 7:23])
rownames(data) <- 1:nrow(data)
ear_boot <- boot.matrix(data = data, bootstraps = 10000)
ear_pos <- dispRity(data = ear_boot, metric = c(mean, displacements))
summary(ear_pos)

# Diogenidae
data.dio <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$fine_group == "Diogenidae"
                              & anobra.scores$dev != "juvenile"))
# not enough juveniles for data analysis
data.dio <- droplevels(data.dio)
data <- as.data.frame(data.dio[, 7:23])
group.v <- data.frame(data.dio$dev)
rownames(data) <- 1:nrow(data)
dio_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(dio_subsets))
dio_boot <- boot.matrix(data = dio_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
dio_pos <- dispRity(data = dio_boot, metric = c(mean, displacements))
summary(dio_pos)

# Galatheidae
data.gal <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$group == "Galatheidae"))
data.gal <- droplevels(data.gal)
data <- as.data.frame(data.gal[, 7:23])
group.v <- data.frame(data.gal$dev)
rownames(data) <- 1:nrow(data)
gal_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(gal_subsets))
gal_boot <- boot.matrix(data = gal_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
gal_pos <- dispRity(data = gal_boot, metric = c(mean, displacements))
summary(gal_pos)

# Grapsidae
data.gra <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$group == "Grapsidae"))
data.gra <- droplevels(data.gra)
data <- as.data.frame(data.gra[, 7:23])
group.v <- data.frame(data.gra$dev)
rownames(data) <- 1:nrow(data)
gra_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(gra_subsets))
gra_boot <- boot.matrix(data = gra_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
gra_pos <- dispRity(data = gra_boot, metric = c(mean, displacements))
summary(gra_pos)

# Hippidae
data.hip <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$fine_group == "Hippidae"))
data.hip <- droplevels(data.hip)
data <- as.data.frame(data.hip[, 7:23])
group.v <- data.frame(data.hip$dev)
rownames(data) <- 1:nrow(data)
hip_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(hip_subsets))
hip_boot <- boot.matrix(data = hip_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
hip_pos <- dispRity(data = hip_boot, metric = c(mean, displacements))
summary(hip_pos)


# Lithodidae
data.lit <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$group == "Lithodidae"))
data.lit <- droplevels(data.lit)
data <- as.data.frame(data.lit[, 7:23])
group.v <- data.frame(data.lit$dev)
rownames(data) <- 1:nrow(data)
lit_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(lit_subsets))
lit_boot <- boot.matrix(data = lit_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
lit_pos <- dispRity(data = lit_boot, metric = c(mean, displacements))
summary(lit_pos)


# Nephropidae
data.nep <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$fine_group == "Nephropidae" 
                              & anobra.scores$dev != "juvenile"))
data.nep <- droplevels(data.nep)
data <- as.data.frame(data.nep[, 7:23])
group.v <- data.frame(data.nep$dev)
rownames(data) <- 1:nrow(data)
nep_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(nep_subsets))
nep_boot <- boot.matrix(data = nep_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
nep_pos <- dispRity(data = nep_boot, metric = c(mean, displacements))
summary(nep_pos)


# Paguridae
data.pag <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$fine_group == "Paguridae" 
))
data.pag <- droplevels(data.pag)
data <- as.data.frame(data.pag[, 7:23])
group.v <- data.frame(data.pag$dev)
rownames(data) <- 1:nrow(data)
pag_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(pag_subsets))
pag_boot <- boot.matrix(data = pag_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
pag_pos <- dispRity(data = pag_boot, metric = c(mean, displacements))
summary(pag_pos)


# Porcellanidae
data.por <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$group == "Porcellanidae" 
))
data.por <- droplevels(data.por)
data <- as.data.frame(data.por[, 7:23])
group.v <- data.frame(data.por$dev)
rownames(data) <- 1:nrow(data)
por_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(por_subsets))
por_boot <- boot.matrix(data = por_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
por_pos <- dispRity(data = por_boot, metric = c(mean, displacements))
summary(por_pos)


# Portunidae
data.port <- data.frame(subset(anobra.scores, 
                               subset = anobra.scores$group == "Portunidae" 
))
data.port <- droplevels(data.port)
data <- as.data.frame(data.port[, 7:23])
group.v <- data.frame(data.port$dev)
rownames(data) <- 1:nrow(data)
port_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(port_subsets))
port_boot <- boot.matrix(data = port_subsets, bootstraps = 10000, 
                         rarefaction = minimum_size)
port_pos <- dispRity(data = port_boot, metric = c(mean, displacements))
summary(port_pos)


# Raninidae
data.ran <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$group == "Raninidae" 
                              & anobra.scores$dev != "juvenile"))
data.ran <- droplevels(data.ran)
data <- as.data.frame(data.ran[, 7:23])
group.v <- data.frame(data.ran$dev)
rownames(data) <- 1:nrow(data)
ran_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(ran_subsets))
ran_boot <- boot.matrix(data = ran_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
ran_pos <- dispRity(data = ran_boot, metric = c(mean, displacements))
summary(ran_pos)


# Xanthidae
data.xan <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$fine_group == "Xanthidae" 
                              & anobra.scores$dev != "juvenile"))
data.xan <- droplevels(data.xan)
data <- as.data.frame(data.xan[, 7:23])
group.v <- data.frame(data.xan$dev)
rownames(data) <- 1:nrow(data)
xan_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(xan_subsets))
xan_boot <- boot.matrix(data = xan_subsets, bootstraps = 10000, 
                        rarefaction = minimum_size)
xan_pos <- dispRity(data = xan_boot, metric = c(mean, displacements))
summary(xan_pos)



## DISPARITY ANALYSIS FOR ALL DEV STAGES FOR BRACHYURA AND ANOMALA--------------

data <- data.frame(subset(anobra.scores, 
                          subset = anobra.scores$major != "Astacidea"))
data <- droplevels(data)
data <- as.data.frame(data[, 7:23])
group.v <- data.frame(indi = subset(df$indi, 
                                    subset = df$major != "Astacidea"))
rownames(data) <- 1:nrow(data)
anobra_subsets <- custom.subsets(data, group = group.v)
minimum_size <- min(size.subsets(anobra_subsets))
minimum_size
anobra_bootstrapped <- boot.matrix(data = anobra_subsets, 
                                   bootstraps = 10000)
anobra_bootrare <- boot.matrix(data = anobra_subsets, 
                               bootstraps = 10000,
                               rarefaction = minimum_size)

anobra_size <- dispRity(data = anobra_bootrare, 
                        metric = c(sum, variances))
summary(anobra_size)
plot(anobra_size, rarefaction = minimum_size)

disp.metrics <- test.metric(anobra_bootrare, metric = c(mean, dispRity::displacements), 
                            shifts = c("random", "size"))
summary(disp.metrics)
plot(disp.metrics)
#Random and outer should be low R2 and inner should high R2
anobra_position <- dispRity(data = anobra_bootrare, 
                            metric = c(mean, dispRity::displacements))
summary(anobra_position)
plot(anobra_position)

# CHECK IF DATA NORMAL DISTRIBUTED
qqnorm(anobra_size$disparity$indi.adultAnomala$elements)
qqline(anobra_size$disparity$indi.adultAnomala$elements)

#Parametric tests, becasue data is bootstrapped and sample size correceted
test.dispRity(anobra_position, test = adonis.dispRity, rarefaction = minimum_size)
test.dispRity(anobra_position, test = t.test, correction = "bonferroni", 
              rarefaction = minimum_size)
test.dispRity(anobra_size, test = t.test, correction = "bonferroni", 
              rarefaction = minimum_size)


## BETWEEN GROUPS POSITION TEST ------------------------------------------------

data.bet <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$dev == "megalopa" & 
                                anobra.scores$major != "Astacidea"))
data.bet <- data.frame(subset(anobra.scores, 
                              subset = anobra.scores$major != "Astacidea"))
data.bet <- droplevels(data.bet)
data <- as.data.frame(data.bet[, 7:23])
group.v <- data.frame(data.bet$indi)
rownames(data) <- 1:nrow(data)
bet_subsets <- custom.subsets(data, group = group.v)
bet_boot <- boot.matrix(data = bet_subsets, 
                        bootstraps = 10000)
bet_pos <- dispRity(data = bet_boot, metric = group.dist, 
                    between.groups = TRUE, probs = c(0.5))
summary(bet_pos)
plot(bet_pos)


#Plot for figure 2 in text

a <- as.data.frame(unlist(bet_pos$disparity$`data.bet.indi.megalopaAnomala:data.bet.indi.megalopaBrachyura`))
a <- cbind(a, as.data.frame(unlist(bet_pos$disparity$`data.bet.indi.adultAnomala:data.bet.indi.juvenileAnomala`)))
a <- cbind(a, as.data.frame(unlist(bet_pos$disparity$`data.bet.indi.adultAnomala:data.bet.indi.megalopaAnomala`)))
a <- cbind(a, as.data.frame(unlist(bet_pos$disparity$`data.bet.indi.adultBrachyura:data.bet.indi.zoeaBrachyura`)))
a <- cbind(a, as.data.frame(unlist(bet_pos$disparity$`data.bet.indi.adultBrachyura:data.bet.indi.zoeaAnomala`)))
a <- cbind(a, as.data.frame(unlist(bet_pos$disparity$`data.bet.indi.juvenileBrachyura:data.bet.indi.zoeaAnomala`)))
names(a)[1] <- "megano:megbra"
names(a)[2] <- "aduano:juvano"
names(a)[3] <- "aduano:megano"
names(a)[4] <- "adubra:zoebra"
names(a)[5] <- "adubra:zoeano"
names(a)[6] <- "juvbra:zoeano"
boxplot(a, col = "white")


