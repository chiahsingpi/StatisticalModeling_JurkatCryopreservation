setwd("D:/07_Statistical Analysis/SGI")
data = read.table("SGI.csv",header=F,sep=',')
colnames(data)<-c("S","G","I","R") #assign column name

colors = c("purple", "blue", "green", "yellow", "orange", "red")

#function to make transparent colors
transp = function(color, pctTransp=0.25) {
  rgb.num = grDevices::col2rgb(color)
  rgb(rgb.num[1], rgb.num[2], rgb.num[3], pctTransp*256, max=256)
}

#function to make the plot
#main is the x-axis osmolyte
#color is the osmolyte to which the coloring will correspond
#fixed is the other osmolyte
#Recovery is the post-thaw recovery vector
#colorName is the name of the "color" osmolyte
#mainName is the name of the "main" osmolyte
plotSpaghetti = function(main, color, fixed, Recovery, colorName, mainName, legendLocation="topleft") {
  plot(main,Recovery,type="n",
       ylab="Post-Thaw Recovery",ylim=c(0,100), xlab=paste0("Level of ",mainName),
       main=paste0("Coloring by Level of ",colorName),
       xlim=c(-0.5,5.5), cex.axis=4.0,cex.lab=4.0,cex.main=3.0)
  #
  for (i in 0:5) {
    for (j in 0:5) {
      points(main[which(fixed==i & color==j)], Recovery[which(fixed==i & color==j)], 
             type="l", lty=1, lwd=3, col=transp(colors[j+1],0.5))
    }
  }
  points(main[which(fixed==0 & color==0)], Recovery[which(fixed==0 & color==0)], 
         type="l", lty=2, lwd =3, col="black")
  legend(legendLocation, paste("Level", 0:5), col=colors, lty=1, bty="n",cex=3.0,ncol=2)
}


#glycerol
jpeg('G_I.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plotSpaghetti(data$G, data$I, data$S, data$R, "Isoleucine", "Glycerol")
dev.off()

jpeg('G_S.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plotSpaghetti(data$G, data$S, data$I, data$R, "Sucrose", "Glycerol")
dev.off()

#sucrose
jpeg('S_I.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plotSpaghetti(data$S, data$I, data$G, data$R, "Isoleucine", "Sucrose", "topright")
dev.off()

jpeg('S_G.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plotSpaghetti(data$S, data$G, data$I, data$R, "Glycerol", "Sucrose", "topright")
dev.off()

#isoleucine
jpeg('I_S.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plotSpaghetti(data$I, data$S, data$G, data$R, "Sucrose", "Isoleucine")
dev.off()

jpeg('I_G.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plotSpaghetti(data$I, data$G, data$S, data$R, "Glycerol", "Isoleucine")
dev.off()