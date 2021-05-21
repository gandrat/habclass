
#Load packages---------
library(raster)
library(factoextra)
library(RStoolbox)
library(reshape2)
library(corrplot)
library(fpc)
library(RColorBrewer)
library(dplyr)

#read and stack rasters-------------
files<-list.files('input_data/SA',full.names = T)

var<-stack(files)
plot(var)

names(var)<-c('bbpi','bioregion','curspeed','depth','oxigen','ligth','nitrate','phosphate','poc','salinity','silicate','slope','temperature')

save(var,file='south_atlantic_dataset.Rda')
#Set the theme for plots----------
theme_set(
  theme_bw(base_size = 14)+
    theme(text=element_text(family="Arial"),
          plot.title = element_text(hjust = 0.5, face='bold',size=14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


#PCA------------------------
gc()
var.df<-as.data.frame(var)
var.df<-var.df[complete.cases(var.df),]

res.pca<-prcomp(var.df,scale = TRUE)


summary(res.pca)
plot(res.pca)


#Figuras PCA----
fviz_screeplot(res.pca, choice='eigenvalue', geom='line')+
  ggtitle(NULL)+
  geom_hline(yintercept=1, linetype='dashed')

ggsave('figures/scree_plot_SA.jpg', width=13, heigh=6, units='in',dpi=150)

#PCA Wheel
fviz_pca_var(res.pca, axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title='')
  
ggsave("figures/pca_wheel_12.jpg", dpi=150, units='in', width=7, height=7)

fviz_pca_var(res.pca, axes = c(3,4),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title='')

ggsave("figures/pca_wheel_34.jpg", dpi=150, units='in', width=7, height=7)

#Raster PCA---------------
pca<-rasterPCA(var,spca=TRUE)

#Comparing results from 2 PCA
summary(pca$model)
summary(res.pca)
plot(pca$map)

pc<-stack(subset(pca$map,c(1:4)))
# pc<-pca$map
plot(pc)
writeRaster(pca$map, "output_data/pca.tif")

#Mapas PCA-----
mapWorld <- borders("world", colour="gray50", fill="gray50", alpha=.7) # create a layer of borders

pc_scale<-(pc-cellStats(pc,"min"))/(cellStats(pc,"max")-cellStats(pc,"min"))
# pc_wgs<-projectRaster(pc_scale,crs="+init=epsg:4326")
# plot(pc_wgs)
pc.p <- rasterToPoints(pc_scale)

#Make the points a dataframe for ggplot
pc.df <- data.frame(pc.p)

minlon<-min(pc.df$x-0.2)
maxlon<-max(pc.df$x+0.2)
minlat<-min(pc.df$y-0.2)
maxlat<-max(pc.df$y+0.2)

d<-melt(pc.df,id=c('x','y'))

ggplot()+geom_raster(data=d,aes(y=y,x=x, fill=value))+theme_bw()+
  # geom_path(data = isobath_df, 
  #           aes(x = long, y = lat, group = group),
  #           linetype = "dashed", size = .2)+
  facet_wrap(~variable, ncol=2)+
  coord_equal(xlim = c(minlon, maxlon),ylim = c(minlat, maxlat))+
  xlab(NULL)+ylab(NULL)+
  scale_fill_distiller(palette = "Spectral")+
  theme(legend.position="none",
        strip.text.x = element_text(size = 8),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        line=element_blank())+
  mapWorld

ggsave('figures/pca_maps.jpg', dpi=150, 
       units='in', width=13, height=7.5)


#Correlation among descriptors and PCs--------
var_pca<-stack(var, pc)
plot(var_pca)

pca.df<-as.data.frame(var_pca)
pca.df<-pca.df[complete.cases(pca.df),]
c<-cor(pca.df, method="pearson")
write.table(c, file='figures/correlation.csv', dec = '.', sep=';')

#Correlogram
png(height=500, width=700, file="figures/corplot.png", type = "cairo")
corrplot(c,method='square',type='lower')
dev.off()

#CH-index------

df <- as.data.frame(scale(pc))

df<-as.data.frame(df[complete.cases(df),])

ci=data.frame(1,NA)
names(ci)=c('i','cii')

for(i in c(2:20)){
  c<-kmeans(df, i, iter.max = 100000)
  cii<-calinhara(df,c$cluster)
  ci<-rbind(ci,data.frame(i,cii))
}

#Grafico C-H
ggplot(ci,aes(x=i,y=cii/10000))+geom_point(size=.5)+geom_line(size=.3)+
  theme_bw()+xlab('Clusters')+ylab('Calinski-Harabasz Index')+
  geom_vline(xintercept = 8, linetype=2)+
  geom_vline(xintercept = 14, linetype=2)+
  scale_x_continuous(breaks = c(2,4,6,8,10,12,14,16,18,20))


ggsave("figures/ch_index.jpg", 
       width = 13 , height = 6, units = 'in', dpi = 200)

#Cluster analysis
full_pc<-pc

full_pc[is.na(full_pc[])] <- -20


df<-as.data.frame(full_pc)
cbig<-kmeans(df, 14, iter.max = 100000)
csmall<-kmeans(df, 8, iter.max = 100000)

max(var$depth)
mask<-reclassify(var$depth,c(-8000,0,1))
plot(mask)
cbig.r <- setValues(mask, cbig$cluster)*mask
plot(cbig.r)
csmall.r <- setValues(mask, csmall$cluster)*mask
plot(csmall.r)


table(cbig$cluster)
table(csmall$cluster)

#Filtro espacial para eliminar ruidos----------

f <- function(x){
  tab <- table(x)
  # I am using the first value here, maybe you want to use the mean, 
  # if 2 or more values occur equally often.
  names(tab)[which.max(tab)][1]
}

cbig.focal <- focal(cbig.r, w=matrix(1,5,5), fun=f)*mask
csmall.focal<-focal(csmall.r, w=matrix(1,5,5), fun=f)*mask

clust<-stack(cbig.focal, csmall.focal)

plot(clust)
names(clust)<-c('c14','c8')
save(clust,file='output_data/habitatsv1.Rda')

#Plot clusters-----------
c.p <- rasterToPoints(clust)

#Make the points a dataframe for ggplot
c.df <- data.frame(c.p)

names(c.df)<-c('x','y','14_clusters','8_clusters')

d<-melt(c.df,id=c('x','y'))

d<-d[complete.cases(d),]
unique(d$value)


# Define the number of colors you want
nb.cols <- 14
mycolors <- colorRampPalette(brewer.pal(10, "Set3"))(nb.cols)

ggplot()+geom_raster(data=d%>%filter(variable=='8_clusters'),aes(y=y,x=x, fill=factor(value)))+theme_bw()+
  # facet_wrap(~variable,ncol=1)+
  coord_equal(xlim = c(minlon, maxlon),ylim = c(minlat, maxlat))+
  xlab(NULL)+ylab(NULL)+
  scale_fill_manual(values=mycolors, name='Zones')+
  mapWorld

ggsave("figures/habitat_clusters_8c.jpg", 
       width = 14 , height = 7, units = 'in', dpi = 200)

ggplot()+geom_raster(data=d%>%filter(variable=='14_clusters'),aes(y=y,x=x, fill=factor(value)))+theme_bw()+
  # facet_wrap(~variable,ncol=1)+
  coord_equal(xlim = c(minlon, maxlon),ylim = c(minlat, maxlat))+
  xlab(NULL)+ylab(NULL)+
  scale_fill_manual(values=mycolors, name='Zones')+
  mapWorld
  

ggsave("figures/habitat_clusters_14c.jpg", 
       width = 14 , height = 7, units = 'in', dpi = 200)
table(d$value,d$variable)


#Boxplots to describe clusters------------------
s <- stack(var, clust$c14)

s.df<-as.data.frame(s)
s.df<-s.df[complete.cases(s.df),]

unique(s.df$c14)
s.m<-melt(s.df,id='c14')

ggplot(s.m,aes(x=factor(c14),y=value,fill=factor(c14)))+
  geom_boxplot(outlier.size=.1,lwd=.1)+
  facet_wrap(~variable,scales='free',ncol=5)+
  theme_bw()+
  scale_fill_manual(values=mycolors)+
  theme(legend.position="none")+
  xlab(NULL)+ylab(NULL)


ggsave('figures/boxplots.jpg',
       width=16, height=6.5, unit='in', dpi=200)

#Decision Tree--------------
load('output_data/habitatsv1.Rda')
load('south_atlantic_dataset.Rda')
s<-stack(var,clust)
names(s)<-c('bbpi','bioregion','curspeed','depth','oxigen','ligth','nitrate','phosphate','poc','salinity','silicate','slope','temperature','clust14','clust8')
plot(s)
s.df<-as.data.frame(s)

s.df<-s.df[complete.cases(s.df),]
s.df$zone<-factor(s.df$clust14)

s.df<-s.df%>%select(-clust8,-clust14)
bfit <- rpart(zone~., method="class", data=s.df)
printcp(bfit) # display the results 
plotcp(bfit) # visualize cross-validation results 
summary(bfit) # detailed summary of splits

jpeg('figures/decision_tree.jpg', width=13, height = 7, units = 'in', res = 300)
rpart.plot(bfit, type=0, nn=FALSE, family = 'Arial',box.palette = 'Blues')
dev.off()
