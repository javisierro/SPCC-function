rm(list = ls())
library(tuneR);library(data.table);library(ggplot2);library(seewave)

maxtoff <- 0.01; t.resolution <- 0.001;wlcorr <- 512;ovlp=80; printspectro=F
## The default maxtoff might be too small for certain cases. This should be adjusted depending on the length of notes being compared. ## i.e. maxtoff <- 0.2

note1 <- readWave('/Users/javiersierro/Biology/2017-21 Lancaster/bioacoustic r protocols/bioacoustic protocol blue tit/SPCC-function/p4_c7_bird_01_bg_060_tag_07_file_3_(2021_07_15-05_54_59)_ASWMUX221163@221163@2021_07_15-05_54_59@45@12851422-12856177.wav') # Directory to a wave file containing a single, isolated note
note2 <- readWave('/Users/javiersierro/Biology/2017-21 Lancaster/bioacoustic r protocols/bioacoustic protocol blue tit/SPCC-function/p4_c7_bird_01_bg_060_tag_07_file_3_(2021_07_15-05_54_59)_ASWMUX221163@221163@2021_07_15-05_54_59@45@12957077-12964082.wav') # Directory to a wave file containing a single, isolated note. Then one to be compared to note1


samp.rate=note1@samp.rate
n1 <- normalize(note1, unit='1', rescale=T, pcm=note1@pcm)
n1 <- as.matrix(n1@left)
n1 <- as.matrix(c(rep(0, length=round(samp.rate*maxtoff)),n1))
n2<- normalize(note2, unit='1', rescale=T, pcm=note2@pcm)
n2 <- as.matrix(n2@left)

max.duration <- as.integer((max(length(n1)/samp.rate, length(n2)/samp.rate)+maxtoff*2)*samp.rate)
if(max.duration > length(n1)){n1 <- as.matrix(c(n1, rep(0, length=samp.rate*5)))}
n1.adj <- as.matrix(n1[1:max.duration,])
if(max.duration > length(n2)){n2 <- as.matrix(c(n2, rep(0, length=samp.rate*5)))}
n2 <- as.matrix(n2[1:max.duration,])

t.displace <- seq(0, maxtoff*2, by=t.resolution);t.loc <- t.displace-maxtoff
n1_mat <- seewave::sspectro(n1.adj, samp.rate,wl=wlcorr, ovlp=ovlp) # Transform sound wave into spectrogram matrix

## This is the SPCC function strictly speaking, previous codes are necessary steps to prepare the sound files
crosscor <- numeric()
for (k in t.displace){
  n2_off <- as.matrix(c(rep(0,length = k*samp.rate), n2, rep(0,length= samp.rate)));n2_off <- as.matrix(n2_off[1:nrow(n1.adj),])
  n2_mat <- seewave::sspectro(n2_off, samp.rate,wl=wlcorr, ovlp=ovlp)
  cor.k <- cor(c(n1_mat), c(n2_mat), method='pearson', use ='all.obs')				
  crosscor <- c(crosscor, cor.k)
  if(printspectro==T){
    longData<-melt(n1_mat + n2_mat);longData$Var1 <- longData$Var1*samp.rate/wlcorr/1000; longData$Var2 <- longData$Var2/wlcorr
    print(ggplot(longData, aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) + scale_fill_gradient(low="white", high="black") + theme_classic()+ theme(legend.position='') + coord_cartesian(ylim=c(0,11)))
  }
}		

## Plotting of SPCC curve with red note demarking the highest, maximum correlation which is taking as the SPCC index for similarity between note1 and note2
plot(x=c(-100,-100), xlim=c(min(t.loc),max(t.loc)),ylim=c(0,1), ylab='Spectrogram correlation', xlab='t offset')
lines(y= crosscor, x= t.loc)
points(y= crosscor, x= t.loc, pch=21, bg='khaki')
points(y= max(crosscor), x= t.loc[which(crosscor==max(crosscor))], pch=21, bg='firebrick', cex=2)
