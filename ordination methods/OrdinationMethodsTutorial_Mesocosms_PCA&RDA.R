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

# thinking about questions
increase4c<-meso_fit%>%filter(Trt=="4C")%>%
  filter_at(vars(5:72), any_vars(. > 0))

increase32c<-meso_fit%>%filter(Trt=="32C")%>%
  filter_at(vars(5:72), any_vars(. > 0))
#this way also works 
meso_fit[colSums(meso_fit[,c(-1:-4)] > 0)]

# which environment is most selective?

ggplot(data = fit%>%filter(Trt=="4C"), aes(x=fitness,color=Trt))+
  geom_histogram(fill="white", position="identity", alpha=.25,bins=50,lwd=1.2)+
  scale_color_manual(values = mycols)+
  theme_bw()

fit%>%group_by(Trt)%>%summarise("min"=min(fitness), "max"=max(fitness), "sd"=sd(fitness), "mean"=mean(fitness), "median"=median(fitness))%>%
  mutate("range"= max-min, "skew"= mean-median)
  
# which strains are generalists? well in all environments?

generalists<-fit%>%filter(fitness>0)%>%ungroup() %>%
  group_by(strain) %>%
  filter(n() == 3)

generalists[order(winners$strain),]
generalists[1:21,]

specialists<-fit%>%filter(fitness>2)%>%ungroup() %>%
  group_by(strain) %>%
  filter(n() == 1)
specialists
