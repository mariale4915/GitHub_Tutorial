#Authors: Liana T. Burghardt, Regina B. Bledsoe, Maria Gil Polo, Jenn Harris, Rebecca Fudge

library(tidyverse)
library(tidyr)
library(ggplot2)

#test this file new comments
setwd()

meso_all<- read.csv("S&R3_Meso_freq.txt",sep= "\t", header=TRUE)

# This is a big dataset comprised of 106 rows (102 mesocosms + 4 initial communities) and 69 columns.

dim(meso_all)

# We need to do some data filtering to focus in on the temperature samples and wrangling to get our reponse and explanatory data into a form we can use it! 
# Let's turn meso_all into a tibble so it easier to work with. This is similar to a data.frame but has a number of improvements. It is also the format that tidyR works on.
#TidyR is amazing in all sorts of ways See here for info on the functions I am about to use (https://r4ds.had.co.nz/transform.html)

meso_all<-as_tibble(meso_all)
meso_all

# The first column is "pool" which has all the sample and treatment information coded into a single string separated by underscores.
# Each additional column is an Ensifer strain. R puts an X in front of each name because a column name can't be a number  
# Next, we will extract all the information from this column into three new columns to make it easy to filter to the data of interest
# We use the dplyr function 'separate' and save the result to a overwrite the original dataset (meso_all). 
# The first argument is the data, the second the column we want to separate, the third the names for those three new columns, fourth the character seperating the information (in this case an underscore), and the fifth is a logical arguement telling R to keep the original column.
#

meso_all <- separate(data = meso_all,col = pool,into = c("Trt","Time","Rep"), sep="_", remove=FALSE)
meso_all

# You will see that we still have the pool column but now we have three additional columns
# We are interested in analyzing the Temperature treatments at 2.5 months here. So lets subset down to just those using the filter argument

### The "pool" column contains all the information about each sample separated by underscores ## 
# 3 temperature treatments: F= ~22C, F32 = 32C, F4 = 4C
# All mesocosms were sampled at 1 timepoint: 2.5 months
# Replicate numbers (between 11 and 15)

#### First we want only 2.5 month timepoint and save it to a new tibble
meso_sub <- filter(meso_all,Time =="2.5m")
meso_sub

# Now we want to subset to the treatments (Trt) of interest (F, F32, and F4). We use the %in% notation because we want more than one Trt 
meso_sub<-filter(meso_sub,Trt %in% c("F","F32","F4")) 
meso_sub

# We can condense these two lines of code into a single one using the "&" (and) notation. What do you think the | (or) notation will give you?
meso_sub <- filter(meso_all,Trt %in% c("F","F32","F4") & Time =="2.5m")
meso_sub

# We can also "pipe" in our data frame into the filter function so we can string commands together. 
# This is also identical and can come in handy such as in making the graph below.
meso_all %>% filter(Trt %in% c("F","F32","F4") & Time =="2.5m") -> meso_sub

#### Quick plot of these frequencies just so you can see it...
freqs<-meso_sub %>% group_by(Trt) %>% summarise_if(is.numeric,mean) %>% pivot_longer(-Trt,names_to= "strain",values_to = "freq")
mycols<-c("gold","orangered","darkblue")
ggplot(data = freqs, aes(x=freq,color=Trt))+
  geom_histogram(fill="white", position="identity", alpha=.25,bins=50,lwd=1.2)+
  scale_color_manual(values = mycols)+
  theme_bw()

# Okay, now we have strain frequencies for our 14 Treatments and timepoints of interest. Next step is to convert those frequencies into log2 fold changes from the initial frequencies (or fitness)
# To do that we need to grab our initial strain frequencies from the main data frame. All of these have a "Time" of "initial" 
# Let's use the piping method I showed above

meso_all %>% filter(Time =="initial") -> initial
initial
# Let's look at the initial data. There are four sequencing replicates of the initial community used to inoculate each mesocosm
# You'll notice that the estimates are very similar across the replicates. 
# Next, we need to calculate the mean frequency of each strain.We will use the apply function to do this. 

initial<-apply(data.frame(initial[,c(-1:-4)]), 2, mean)

##### Quick plot of the initial frequencies ########
ggplot(data = data.frame(Freq=initial),aes(x=Freq))+geom_histogram(fill="white",color="black")+theme_bw()

#### Now we want to divide every strain frequency by the initial mean frequency and take the log 2 () to calculate "fitness": log2(selected_freq/initial_means). 
# log2(selected_freq/initial_means). 
#I can't for the life of me at the moment figure out how to do this function in tidyR 
#so converting to a dataframes that work with vectorized function and then reassembling the tibble.
# Can one of you figure out a solution in TidyR?

meso_fit<-as_tibble(cbind(meso_sub[,c(1:4)],log2(as.matrix(data.frame(meso_sub[,c(-1:-4)]))/initial)))

#Put temperatures in order and rename for graphing

meso_fit %>% mutate(Trt = factor(Trt, levels= c("F4","F","F32"),labels= c("4C","22C","32C"))) -> meso_fit

# Define a color scheme
mycols<-c("darkblue","goldenrod","orangered")

##### Quick Histogram ####
fit<-meso_fit %>% group_by(Trt) %>% summarise_if(is.numeric,mean) %>% pivot_longer(-Trt,names_to= "strain",values_to = "fitness")
ggplot(data = fit, aes(x=fitness,color=Trt))+
  geom_histogram(fill="white", position="identity", alpha=.25,bins=50,lwd=1.2)+
  scale_color_manual(values = mycols)+
  theme_bw()

ggplot(data = fit%>%filter(Trt=="4C"), aes(x=fitness,color=Trt))+
  geom_histogram(fill="white", position="identity", alpha=.25,bins=50,lwd=1.2)+
  scale_color_manual(values = mycols)+
  theme_bw()



### NOW lets run FINALLY run a PCA #####
library(FactoMineR)
library(factoextra)

#Quick Graphs
PCA(meso_fit[,c(-1:-4)],graph = TRUE)

# Save results to an object
res.PCA<-PCA(meso_fit[,c(-1:-4)],graph = FALSE)
# Access the Eigenvalues
get_eigenvalue(res.PCA)
# Graph the Eigenvalues. In this case, the vast majority (~80%) of the variation in the data can be explained by the first two dimensions. 
fviz_eig(res.PCA)
#Quick combined graph
fviz_pca_biplot(res.PCA, axes = c(3,4), habillage = meso_fit$Trt, palette = mycols)

#### You can also plot samples (individuals) and strains (variables) seperately
fviz_pca_ind(res.PCA,habillage = meso_fit$Trt, palette = mycols)
fviz_pca_var(res.PCA)

### We can get an idea of which axis each variable is loaded on with the cos2 values in the PCA variable results. 
# Higher values mean a greater contribution to that dimension.
var <- get_pca_var(res.PCA)
var$cos2

### We can quickly visulize this data with corrplot
library(corrplot)
corrplot(var$cos2, is.corr=FALSE,method = "shade",tl.cex = .5, tl.col = "black",cl.cex=.5,cl.align.text="l",cl.ratio= .75)
# There are only a few strains highly loaded on the second dimension (e.g. strain 1660). This means that this strain increases/decreases along Dimension 2. 
# What treatments are differentiated along Dim2? See if this is also the case in the RDA analysis below! 

#For a much more comprensive explaination of the the outcomes from these function see: 
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

##### RDA! ######
library(vegan)

# Run the rda model rda(y~x,data,scale) 
rda1<-rda(meso_fit[,c(5:72)]~Trt,meso_fit, scale=TRUE)

# Quick plot of rda results
ordiplot(rda1,type = "text")
# This is called a triplot. 

# Access all of the results
summary(rda1)

# Rsquared gives you an overall idea of well the overall model (your explanatory variables) did at explaining the variation in the data
RsquareAdj(rda1)
# In this example,treatments explain ~ 2/3rds of the variation in strain fitness. That is really high!

# Run a permutational significance test for the explanatory variables in the model
anova(rda1, step=1000, perm.max=1000, by= "terms")
# Temperature treatment has a significant effect on strain fitness (P<.001)

# Run a permutational significance test for each axis
anova(rda1, step=1000, perm.max=1000, by= "axis")
# In this case both of the RDA axis are highly significant.

####### A much prettier visualization ########

#Make the window view wider
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))

# Initiate the plot
rdaplot <- ordiplot(rda1, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",
                    xlim=c(-2,2),ylim=c(-3,6),scaling=1, xlab=paste("RDA 1* (",round(summary(rda1)$cont$importance[2,1],2)*100,"% var.)",sep=""), 
                    ylab=paste("RDA 2* (",round(summary(rda1)$cont$importance[2,2],2)*100,"% var.)",sep=""))
# Add the points to it
points(rda1,"wa", cex=0.8,pch=16,col=mycols[meso_fit$Trt])
# Add a circle around the points
ordiellipse(ord = rda1, mycols[meso_fit$Trt], kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols), border=paste(mycols),lty=c(1) ,alpha=63)
# Add lines from the Centroid 
ordispider(rda1, meso_fit$Trt, lwd=1,label =FALSE,col=paste(mycols),cex=.5)
# Add a legend
legend("topleft", legend = levels(meso_fit$Trt), bty = "n",
       col = mycols, pch = 21, pt.bg = mycols,)

#####--------Question 6-------
#Which strains are generalists? well in all environments? 
#Jenn
#these are the species that increase under 4c
increase4c<-meso_fit%>%filter(Trt=="4C")%>%
  filter_at(vars(5:72), any_vars(. > 0))

#species that increase under 32C
increase32c<-meso_fit%>%filter(Trt=="32C")%>%
  filter_at(vars(5:72), any_vars(. > 0))

#which species increase?
meso_fit[colSums(meso_fit[,c(-1:-4)] > 0)]
#which are generalists?
generalists<-fit%>%filter(fitness>0)%>%ungroup() %>% #species that increase
  group_by(strain) %>% 
  filter(n() == 3)
unique(generalists$strain)
# these strains are all generalists because they increase under all conditions

specialists<-fit%>%filter(fitness>2)%>%ungroup() %>%
  group_by(strain) %>%
  filter(n() == 1)
specialists
#these strains are all specialists because they increased more than (double?) but only in one condition

#9. How do these results differ from the inferences we would make from ordination methods that do not require a normally distributed outcome. What would happen if we used NMDS or PCoA to analyze the raw frequencies?

#Principal coordinates analysis (metric multidimensional scaling) or PCoA is a way to visually represent compositional dissimilarity between sites represented by multiple count variables (species) at each site. PCoA is a way to represent multivariate data in a low-dimension Euclidean (metric, based on eigenvalues) space. PCoA produces a set of uncorrelated (orthogonal) axes to summarize the variability in the data set. This is similar to PCA, but PCoA does not assume that data is of a certain distribution and uses a dissimilarity matrix such as Bray-Curtis, Sorensen, or Jaccard as input. 

#Using the meso_sub dataset, first convert the dataset to a dissimilarity matrix
#We use the vegan function for this... and will use the Bray-Curtis method however there are more options for distance matrix. https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/vegdist

#First we need and site by species matrix so using meso-sub, first we need to set treatments (the pool column in meso_sub dataset) as the rownames
#Assign pool column to rownames
rownames(meso_sub)=meso_sub$pool
#Calculate bray-curtis matrix
meso_sub_dist<-vegdist(meso_sub[,-c(1:4)], method="bray")

#Bray-Curtis ranks site similarity from 0 to 1, with 0 meaning the two sites are similar and 1 meaning the two sites are dissimilar. This value is used as a distance in PCoA analysis. 

#Calculate PCoA components
meso_sub_dist_pcoa <- cmdscale(meso_sub_dist, k=3, eig=TRUE, add=FALSE)
#Calculate variation for each axis
axis1 <- round(meso_sub_dist_pcoa$eig[1] / sum(meso_sub_dist_pcoa$eig), 3) * 100
axis1

axis2 <- round(meso_sub_dist_pcoa$eig[2] / sum(meso_sub_dist_pcoa$eig), 3) * 100
axis2 

#This is just the sum... how much variation in composition is explained by both axes.
sum.eiga <- sum(axis1, axis2)
sum.eiga

#make dataframe with points and add row names, this is what we will use to make plots with
pcoa_points <- as.data.frame(meso_sub_dist_pcoa$points)
rownames(pcoa_points) <- meso_sub$pool
#Add meta data to points and make factors
pcoa_points <- cbind(meso_sub[,c(2:4)], pcoa_points)
pcoa_points$Trt <- as.factor(pcoa_points$Trt)

#Use ggplot to plot pcoa_points
pcoa_temp <- ggplot(pcoa_points, aes(x=V1, y=V2, color=Trt)) + 
  geom_point(aes(color=Trt)) + 
  scale_color_manual(values = c("blue", "green", "red"), labels = c("4C", "22C", "32C"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  xlab("PCoA 1 (74.6%)") + ylab("PCoA 2 (19.7%)") + 
  labs(colour="Treatment") +
  ggtitle("PCoA of Temperature Treatments")+
  theme(legend.position="right")

pcoa_temp

#Save a copy of the plot. 
ggsave("figures/pcoa_temp.png", plot=last_plot(), device=NULL, path=NULL, scale=1, width=6.3, height=4.4, dpi=900, limitsize=TRUE,bg="transparent")


###Here I add time as a factor to see how over time the community composition changes within temperature treatments. 
meso_all %>% filter(Trt %in% c("F","F32","F4")) -> meso_sub_time
meso_sub_time<-as.data.frame(meso_sub_time)

#Assign pool to rownames
rownames(meso_sub_time)=meso_sub_time$pool
#Calculate bray-curtis matrix
meso_sub_time_dist<-vegdist(meso_sub_time[,-c(1:4)], method="bray")

#Calculate PCoA components
meso_sub_time_dist_pcoa <- cmdscale(meso_sub_time_dist, k=3, eig=TRUE, add=FALSE)
#Calculate variation for each axis 
axis1 <- round(meso_sub_time_dist_pcoa$eig[1] / sum(meso_sub_time_dist_pcoa$eig), 3) * 100
axis1

axis2 <- round(meso_sub_time_dist_pcoa$eig[2] / sum(meso_sub_time_dist_pcoa$eig), 3) * 100
axis2 

sum.eiga <- sum(axis1, axis2)
sum.eiga

#make dataframe with points and add row names 
pcoa_points <- as.data.frame(meso_sub_time_dist_pcoa$points)
rownames(pcoa_points) <- meso_sub_time$pool
#Add meta data to points and make factors
pcoa_points <- cbind(meso_sub_time[,c(2:4)], pcoa_points)
pcoa_points$Trt <- as.factor(pcoa_points$Trt)
pcoa_points$Time <- as.factor(pcoa_points$Time)

#Use ggplot to plot pcoa_points
pcoa_time <- ggplot(pcoa_points, aes(x=V1, y=V2, color=Trt, shape=Time)) + 
  geom_point(aes(color=Trt, shape=Time), size=4, stroke=1) + 
  scale_color_manual(values = c("blue", "green", "red"), labels = c("4C", "22C", "32C"))+
  scale_shape_manual(values = c(1,12,19))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  xlab("PCoA 1 (62.2%)") + ylab("PCoA 2 (23.1%)") + 
  labs(colour="Treatment", shape="Time") +
  ggtitle("PCoA with Temp and Time")+
  theme(legend.position="right")

pcoa_time

ggsave("figures/pcoa_time.png", plot=last_plot(), device=NULL, path=NULL, scale=1, width=6.3, height=4.4, dpi=900, limitsize=TRUE,bg="transparent")

#To validate to what we see, I usually pair this type of plot with a PERMANOVA analysis
meso_sub_time_ad <- adonis(meso_sub_time[,-c(1:4)] ~Trt+Time, method = "bray", data = meso_sub_time, perm=1000, set.seed=42)

meso_sub_time_ad 

# Question 5 ####
install.packages ("moments")
library (moments)


# Here, I filter to just include the 22C treatment
med <- filter(meso_all,Trt =="22C")

#Then, I find the range.
range(med[,c(5:50)])

# Then I repeat the process with the 4C treatment
low <- filter(meso_fit,Trt =="4C")
range(low[,c(5:50)])

# Then the 32C treatment
high <- filter(meso_fit,Trt =="32C")
range(high[,c(5:50)])

# The widest range occurred at 4C, indicating cold was most selective!
# Good to know for my research :)

# I couldn't figure out how to get the standard deviation or skewness to work.
# Once I figure those out, I'll run the ANOVA.

# I accidentally answered part of this question too -Jenn
# I'll just include what I did here

# 5. which environment is most selective?

ggplot(data = fit%>%filter(Trt=="4C"), aes(x=fitness,color=Trt))+
  geom_histogram(fill="white", position="identity", alpha=.25,bins=50,lwd=1.2)+
  scale_color_manual(values = mycols)+
  theme_bw()
# right skewed distribution shows that most species decrease in abundance at 4C, and few are increases greatly

fit%>%group_by(Trt)%>%summarise("min"=min(fitness), "max"=max(fitness), "sd"=sd(fitness), "mean"=mean(fitness), "median"=median(fitness))%>%
  mutate("range"= max-min, "skew"= mean-median)
#4c is right skewed because the median is smaller than the mean.
#22c and 32c is also  right skewed because the median is smaller than the mean.
#-Jenn

# Question 8 ####

# Similar to the code above, I filter by time, then find the range.
t1 <- filter(meso_all,Time =="2wk")
range(t1[,c(5:50)])

t2 <- filter(meso_all,Time =="2.5m")
range(t2[,c(5:50)])

t3 <- filter(meso_all,Time =="6m")
range(t3[,c(5:50)])


# Here is a summary of the ranges of the temperature treatments
# and the time treatments.
my_summary <- list("Temp Fitness Range" =
                     list("22C" = range(med[,c(5:50)]),
                          "4C"  = range(low[,c(5:50)]),
                          "32C" = range(high[,c(5:50)])),
                   "Time Fitness Range" = 
                     list("2wk"  = range(t1[,c(5:50)]),
                          "2.5m" = range(t2[,c(5:50)]),
                          "6m"   = range(t3[,c(5:50)]))
)

# Temperature had a bigger effect on strain fitness than time.

# Should I run an ANOVA? I couldn't figure out how to do it, but will work on it.

#------Question 10-----

# First we want only 2.5 month timepiont and our treatments of interest.
# I'm picking clay and salt added
meso_sub <- filter(meso_all,Trt %in% c("F","FNa","CY") & Time =="2.5m")
meso_sub

# Okay, now we have strain frequencies for our 14 Treatments and timepoints of interest.
# Next step is to convert those frequencies into log2 fold changes from the initial frequencies (or fitness)
# To do that we need to grab our initial strain frequencies from the main data frame. All of these have a "Time" of "initial" 
meso_all %>% filter(Time =="initial") -> initial
initial
# Next, we need to calculate the mean frequency of each strain.We will use the apply function to do this. 
initial<-apply(data.frame(initial[,c(-1:-4)]), 2, mean)

# Quick plot of the initial frequencies
ggplot(data = data.frame(Freq=initial),aes(x=Freq))+geom_histogram(fill="white",color="black")+theme_bw()

#### Now we want to divide every strain frequency by the initial mean frequency and take the log 2 () to calculate "fitness": log2(selected_freq/initial_means). 
# log2(selected_freq/initial_means). 
# I can't for the life of me at the moment figure out how to do this function in tidyR 
# so converting to a dataframes that work with vectorized function and then reassembling the tibble.
# Can one of you figure out a solution in TidyR?

meso_fit<-as_tibble(cbind(meso_sub[,c(1:4)],log2(as.matrix(data.frame(meso_sub[,c(-1:-4)]))/initial)))

#Put temperatures in order and rename for graphing

meso_fit <-meso_fit %>% mutate(Trt = factor(Trt, levels= c("F","FNa","CY"), labels= c("F","FNa","CY")))

# Define a color scheme
mycols<-c("darkblue","goldenrod","orangered")

##### Quick Histogram (again with new treatments) ####
fit<-meso_fit %>% group_by(Trt) %>% summarise_if(is.numeric,mean) %>% pivot_longer(-Trt,names_to= "strain",values_to = "fitness")
ggplot(data = fit, aes(x=fitness,color=Trt))+
  geom_histogram(fill="white", position="identity", alpha=.25,bins=50,lwd=1.2)+
  scale_color_manual(values = mycols)+
  theme_bw()

ggplot(data = fit%>%filter(Trt=="F"), aes(x=fitness,color=Trt))+
  geom_histogram(fill="white", position="identity", alpha=.25,bins=50,lwd=1.2)+
  scale_color_manual(values = mycols)+
  theme_bw()

### NOW lets run FINALLY run a PCA (again) #####
library(FactoMineR)
library(factoextra)

#Quick Graphs
PCA(meso_fit[,c(-1:-4)],graph = TRUE)

# Save results to an object
res.PCA<-PCA(meso_fit[,c(-1:-4)],graph = FALSE)
# Access the Eigenvalues
get_eigenvalue(res.PCA)
# Graph the Eigenvalues. In this case, the vast majority (~80%) of the variation in the data can be explained by the first two dimensions. 
fviz_eig(res.PCA)
#Quick combined graph
fviz_pca_biplot(res.PCA, axes = c(3,4), habillage = meso_fit$Trt, palette = mycols)

#### You can also plot samples (individuals) and strains (variables) seperately
fviz_pca_ind(res.PCA,habillage = meso_fit$Trt, palette = mycols)
fviz_pca_var(res.PCA)

### We can get an idea of which axis each variable is loaded on with the cos2 values in the PCA variable results. 
# Higher values mean a greater contribution to that dimension.
var <- get_pca_var(res.PCA)
var$cos2

### We can quickly visulize this data with corrplot
library(corrplot)
corrplot(var$cos2, is.corr=FALSE,method = "shade",tl.cex = .5, tl.col = "black",cl.cex=.5,cl.align.text="l",cl.ratio= .75)
# There are only a few strains highly loaded on the second dimension (e.g. strain 1660). This means that this strain increases/decreases along Dimension 2. 
# What treatments are differentiated along Dim2? See if this is also the case in the RDA analysis below! 

##### RDA! (again)######
library(vegan)

# Run the rda model rda(y~x,data,scale) 
rda1<-rda(meso_fit[,c(5:72)]~Trt,meso_fit, scale=TRUE)

# Quick plot of rda results
ordiplot(rda1,type = "text")
# This is called a triplot. 

# Access all of the results
summary(rda1)

# Rsquared gives you an overall idea of well the overall model (your explanatory variables) did at explaining the variation in the data
RsquareAdj(rda1)
# In this example,treatments explain ~ 2/3rds of the variation in strain fitness. That is really high!

# Run a permutational significance test for the explanatory variables in the model
anova(rda1, step=1000, perm.max=1000, by= "terms")
# Temperature treatment has a significant effect on strain fitness (P<.001)

# Run a permutational significance test for each axis
anova(rda1, step=1000, perm.max=1000, by= "axis")
# In this case both of the RDA axis are highly significant.

####### A much prettier visualization (again) ########

#Make the window view wider
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))

# Initiate the plot
rdaplot <- ordiplot(rda1, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",
                    xlim=c(-2,2),ylim=c(-3,6),scaling=1, xlab=paste("RDA 1* (",round(summary(rda1)$cont$importance[2,1],2)*100,"% var.)",sep=""), 
                    ylab=paste("RDA 2* (",round(summary(rda1)$cont$importance[2,2],2)*100,"% var.)",sep=""))
# Add the points to it
points(rda1,"wa", cex=0.8,pch=16,col=mycols[meso_fit$Trt])
# Add a circle around the points
ordiellipse(ord = rda1, mycols[meso_fit$Trt], kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols), border=paste(mycols),lty=c(1) ,alpha=63)
# Add lines from the Centroid 
ordispider(rda1, meso_fit$Trt, lwd=1,label =FALSE,col=paste(mycols),cex=.5)
# Add a legend
legend("topleft", legend = levels(meso_fit$Trt), bty = "n",
       col = mycols, pch = 21, pt.bg = mycols,)

# Question 7:

# Grab the temperature treatments
meso_sub<-filter(meso_all,Trt %in% c("F","F32","F4"))

# Calculate fitness
meso_all %>% filter(Time =="initial") -> initial
initial<-apply(data.frame(initial[,c(-1:-4)]), 2, mean)
meso_fit<-as_tibble(cbind(meso_sub[,c(1:4)],log2(as.matrix(data.frame(meso_sub[,c(-1:-4)]))/initial)))

#Put temperatures and times in order and rename for graphing and create a column the combined treatment

meso_fit <-meso_fit %>% mutate(Trt = factor(Trt, levels= c("F4","F","F32"), labels= c("F4","F22","F32")),
                               Time = factor(Time, levels= c("2wk","2.5m","6m"), labels= c("2wk","10wk","24wk")),
                               Trt.Time = factor(paste(Trt,Time,sep="_"),levels=c("F4_10wk","F4_24wk","F22_2wk","F22_10wk","F22_24wk","F32_2wk","F32_10wk","F32_24wk")))
meso_fit <-meso_fit %>% select(Trt.Time,everything())

# Define a color and shape scheme
mycols<-c("lightblue","darkblue","gold","goldenrod","darkgoldenrod","orange","orangered","red")
myshapes<-c(15,16,17)

# Save results to an object
res.PCA<-PCA(meso_fit[,c(-1:-5)],graph = FALSE)
# Access the Eigenvalues
get_eigenvalue(res.PCA)

#### You can also plot samples (individuals) and strains (variables) seperately
fviz_pca_ind(res.PCA,axes = c(1,2),habillage = meso_fit$Trt.Time, palette = mycols)
fviz_pca_ind(res.PCA,axes = c(3,4),habillage = meso_fit$Trt.Time, palette = mycols)

#### Write a few different RDA Models
rda1<-rda(meso_fit[,c(6:73)]~Trt*Time,meso_fit, scale=TRUE) #Full Model w/ interactions
rda2<-rda(meso_fit[,c(6:73)]~Trt+Time,meso_fit, scale=TRUE) # Additive Model
rda3<-rda(meso_fit[,c(6:73)]~Trt,meso_fit, scale=TRUE) # Treatment only
rda4<-rda(meso_fit[,c(6:73)]~Time,meso_fit, scale=TRUE) # Time only

# Rsquared gives you an overall idea of well the overall model adjusted for the number of predictors you have
RsquareAdj(rda1)$adj.r.squared
RsquareAdj(rda2)$adj.r.squared
RsquareAdj(rda3)$adj.r.squared
RsquareAdj(rda4)$adj.r.squared
# Here the adjusted r-squared for the model allowing for an interaction between Time and Treatment has the highest rsquared even when penalized for the additional predictors

# Run a permutational significance test for the explanatory variables in the model
anova(rda1, step=1000, perm.max=1000, by= "terms") # The interaction is significant, but temperature has the largest effect
anova(rda2, step=1000, perm.max=1000, by= "terms")

#### Plot the FULL model #######

#Make the window view wider
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))

# Initiate the plot
rdaplot <- ordiplot(rda1, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",
                    xlim=c(-2,2),ylim=c(-3,6),scaling=1, xlab=paste("RDA 1 (",round(summary(rda1)$cont$importance[2,1],2)*100,"% var.)",sep=""), 
                    ylab=paste("RDA 2 (",round(summary(rda1)$cont$importance[2,2],2)*100,"% var.)",sep=""))
# Add the points to it
points(rda1,"wa", cex=0.8,col=mycols[meso_fit$Trt.Time],pch=myshapes[meso_fit$Time])
# Add a circle around the points
ordiellipse(ord = rda1, groups = meso_fit$Trt.Time, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols), border=paste(mycols),lty=c(1) ,alpha=63)
# Add lines from the Centroid 
ordispider(rda1, meso_fit$Trt.Time, lwd=1,label =FALSE,col=paste(mycols),cex=.5)
# Add a legend
legend("topright", legend = levels(meso_fit$Trt.Time), bty = "n",
       col = mycols, pch = c(16,17,15,16,17,15,16,17), pt.bg = mycols)

##### Plot the Additive Model #######

rdaplot <- ordiplot(rda2, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",
                    xlim=c(-2,2),ylim=c(-3,6),scaling=1, xlab=paste("RDA 1 (",round(summary(rda2)$cont$importance[2,1],2)*100,"% var.)",sep=""), 
                    ylab=paste("RDA 2 (",round(summary(rda2)$cont$importance[2,2],2)*100,"% var.)",sep=""))
# Add the points to it
points(rda2,"wa", cex=0.8,col=mycols[meso_fit$Trt.Time],pch=myshapes[meso_fit$Time])
# Add a circle around the points
ordiellipse(ord = rda2, groups = meso_fit$Trt.Time, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols), border=paste(mycols),lty=c(1) ,alpha=63)
# Add lines from the Centroid 
ordispider(rda2, meso_fit$Trt.Time, lwd=1,label =FALSE,col=paste(mycols),cex=.5)
# Add a legend
legend("topright", legend = levels(meso_fit$Trt.Time), bty = "n",
       col = mycols, pch = c(16,17,15,16,17,15,16,17), pt.bg = mycols)
