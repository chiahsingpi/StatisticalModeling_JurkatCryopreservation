# This code was written by Dr. Ashley Petersen and adjusted by Chia-Hsing Pi
# This code used (1)quasi-binomial, (2) sucrose as factor, (3) one simple model and one interaction models
####################################################################################################################################

# import data
setwd("D:/07_Statistical Analysis/SGI2")
data_SGI<-read.table('SGI_log.csv',header=F,sep=',') #input data
colnames(data_SGI) = c("S","G","I","R") 

#collapse data into counts
collapsed = aggregate(data_SGI$R, data_SGI[,1:3], FUN=sum)
temp = aggregate(data_SGI$R, data_SGI[,1:3], FUN=length)
collapsed = cbind(collapsed, temp[,4]-collapsed[,4])
names(collapsed)[4:5] = c("success", "failure")

#coloring for plots
colors = c("purple", "blue", "green", "yellow", "orange", "red")

####################################################################################################################################
## "main effects" model (no interactions)
####################################################################################################################################

mainEffects = glm(cbind(success, failure) ~ as.factor(S) + G + I, data=collapsed, family="quasibinomial")
summary(mainEffects)

#p-value for effect of S
mainEffects_noS = glm(cbind(success, failure) ~ G + I, data=collapsed, family="quasibinomial")
anova(mainEffects, mainEffects_noS, test = "F")

#predicted values for log odds
fit = predict(mainEffects, type = "link")

#plot summarizing impact of additive effects
#glycerol
jpeg('MainEffects_G.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plot(collapsed$S, fit, type="n", xlab="Sucrose Level", ylab="Log Odds of Post-Thaw Recovery", 
     cex=2.25,cex.axis=4.0,cex.lab=4.0, ylim=c(-2,1))
for (i in 0:5) points(collapsed$S[collapsed$G==i & collapsed$I==0], fit[collapsed$G==i & collapsed$I==0], 
                      type="l", col=colors[i+1])
legend("topright", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0,ncol = 3)
dev.off()

#isoleucine
jpeg('MainEffects_I.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plot(collapsed$S, fit, type="n", xlab="Sucrose Level", ylab="Log Odds of Post-Thaw Recovery", 
     cex=2.25,cex.axis=4.0,cex.lab=4.0)
for (i in 0:5) points(collapsed$S[collapsed$I==i & collapsed$G==0], fit[collapsed$I==i & collapsed$G==0], 
                      type="l", col=colors[i+1])
legend("topright", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0)
dev.off()

####################################################################################################################################
## interaction model
####################################################################################################################################

intxn = glm(cbind(success, failure) ~ as.factor(S)*G + as.factor(S)*I + G*I, data=collapsed, family="quasibinomial")
summary(intxn)

#p-value for interaction between S & G
intxn_noSG = glm(cbind(success, failure) ~ as.factor(S) + as.factor(S)*I + G*I, data=collapsed, family="quasibinomial")
anova(intxn, intxn_noSG, test = "F")

#p-value for interaction between S & I
intxn_noSI = glm(cbind(success, failure) ~ as.factor(S) + as.factor(S)*G + G*I, data=collapsed, family="quasibinomial")
anova(intxn, intxn_noSI, test = "F")

#p-value for interaction between G & I
intxn_noGI = glm(cbind(success, failure) ~ as.factor(S)*G + as.factor(S)*I, data=collapsed, family="quasibinomial")
anova(intxn, intxn_noGI, test = "F")

#predicted values for log odds
fitIntxn = predict(intxn, type = "link")


#effect of S (interaction)
#isoleucine
for (g in 0:5) {
  file_name=paste("I_G_",g,".jpeg",sep = "")
  jpeg(file_name,width=900,height=900)
  par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
  plot(collapsed$S, fitIntxn, type="n", xlab="Sucrose Level", ylab="Log Odds of Post-Thaw Recovery", 
       cex=2.25,cex.axis=4.0,cex.lab=4.0)
  #main=paste0("Color by isoleucine, glycerol=",g),
  if (g==0) legend("topright", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0)
  for (i in 0:5) points(collapsed$S[collapsed$I==i & collapsed$G==g], fitIntxn[collapsed$I==i & collapsed$G==g], 
                        type="l", col=colors[i+1])
  dev.off()
}

#glycerol
for (i in 0:5) {
  file_name=paste("G_I_",i,".jpeg",sep = "")
  jpeg(file_name,width=900,height=900)
  par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
  plot(collapsed$S, fitIntxn, type="n", xlab="Sucrose Level", ylab="Log Odds of Post-Thaw Recovery", 
       cex=2.25,cex.axis=4.0,cex.lab=4.0)
  #main=paste0("Color by glycerol, isoleucine=",i),
  if (i==0) legend("topright", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0)
  for (g in 0:5) points(collapsed$S[collapsed$G==g & collapsed$I==i], fitIntxn[collapsed$G==g & collapsed$I==i], 
                        type="l", col=colors[g+1])
  dev.off()
}

#effect of G (interaction)
#sucrose
for (i in 0:5) {
  file_name=paste("S_I_",i,".jpeg",sep = "")
  jpeg(file_name,width=900,height=900)
  par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
  plot(collapsed$G, fitIntxn, type="n", xlab="Glycerol Level", ylab="Log Odds of Post-Thaw Recovery", 
       cex=2.25,cex.axis=4.0,cex.lab=4.0)
  #main=paste0("Color by sucrose, isoleucine=",i),
  if (i==0) legend("topright", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0)
  for (s in 0:5) points(collapsed$G[collapsed$S==s & collapsed$I==i], fitIntxn[collapsed$S==s & collapsed$I==i], 
                        type="l", col=colors[s+1])
  dev.off()
}

#isoleucine
for (s in 0:5) {
  file_name=paste("I_S_",s,".jpeg",sep = "")
  jpeg(file_name,width=900,height=900)
  par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
  plot(collapsed$G, fitIntxn, type="n", xlab="Glycerol Level", ylab="Log Odds of Post-Thaw Recovery", 
       cex=2.25,cex.axis=4.0,cex.lab=4.0)
  #main=paste0("Color by isoleucine, sucrose=",s),
  if (s==0) legend("topleft", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0)
  for (i in 0:5) points(collapsed$G[collapsed$I==i & collapsed$S==s], fitIntxn[collapsed$I==i & collapsed$S==s], 
                        type="l", col=colors[i+1])
  dev.off()
}

#effect of I (interaction)
#glycerol
for (s in 0:5) {
  file_name=paste("G_S_",s,".jpeg",sep = "")
  jpeg(file_name,width=900,height=900)
  par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
  plot(collapsed$I, fitIntxn, type="n", xlab="Isoleucine Level", ylab="Log Odds of Post-Thaw Recovery", 
       cex=2.25,cex.axis=4.0,cex.lab=4.0)
  #main=paste0("Color by glycerol, sucrose=",s),
  if (s==0) legend("topleft", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0)
  for (g in 0:5) points(collapsed$I[collapsed$G==g & collapsed$S==s], fitIntxn[collapsed$G==g & collapsed$S==s], 
                        type="l", col=colors[g+1])
  dev.off()
}

#sucrose
for (g in 0:5) {
  file_name=paste("S_G_",g,".jpeg",sep = "")
  jpeg(file_name,width=900,height=900)
  par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
  plot(collapsed$I, fitIntxn, type="n", xlab="Isoleucine Level", ylab="Log Odds of Post-Thaw Recovery", 
       cex=2.25,cex.axis=4.0,cex.lab=4.0)
  #main=paste0("Color by sucrose, glycerol=",g),
  if (g==0) legend("topright", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0)
  for (s in 0:5) points(collapsed$I[collapsed$S==s & collapsed$G==g], fitIntxn[collapsed$S==s & collapsed$G==g], 
                        type="l", col=colors[s+1])
  dev.off()
}

####################################################################################################################################
## predicted vs. actual plots
####################################################################################################################################

fitp = predict(mainEffects, type = "response")
fitIntxnp = predict(intxn, type = "response")
obs = collapsed$success/(collapsed$success+collapsed$failure)

jpeg('pva_principal.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plot(obs*100, fitp*100, xlim=c(0,100), ylim=c(0,100), ylab="Predicted Post-Thaw Recovery (%)", xlab="Actual Post-Thaw Recovery (%)", 
     pch=19, cex=2, cex.axis=4.0,cex.lab=4.0)
abline(0,1)
#plot(fitIntxnp, fitIntxnp-obs)
dev.off()

jpeg('pva_interaction.jpeg',width=900,height=900)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
plot(obs*100, fitIntxnp*100, xlim=c(0,100), ylim=c(0,100), ylab="Predicted Post-Thaw Recovery (%)", xlab="Actual Post-Thaw Recovery (%)", 
     pch=19, cex=2, cex.axis=4.0,cex.lab=4.0)
abline(0,1)
#plot(fitIntxnp, fitIntxnp-obs)
dev.off()

##########################################################################################################

###########################################################################################################
#effect of S (interaction)
#glycerol
file_name=paste("S_G.jpeg",sep = "")
jpeg(file_name, width = 3600, height = 3600, units = 'px', res = 300)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
for (i in 0:5) {
  if(i==0){
    plot(collapsed$S, fitIntxn, type="n", xlab="Sucrose Level", ylab="Log Odds of Post-Thaw Recovery", 
         main=paste0("Coloring by Level of Glycerol"), cex.main=3.0, cex=2.25,cex.axis=4.0,cex.lab=4.0,ylim=c(-2,1.5))
    legend("topright", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0)
  }
  else{
    points(collapsed$S, fitIntxn, type="n")
  }
  for (g in 0:5) points(collapsed$S[collapsed$G==g & collapsed$I==i], fitIntxn[collapsed$G==g & collapsed$I==i], 
                        type="l", col=colors[g+1])
}
dev.off()

#isoleucine
file_name=paste("S_I.jpeg",sep = "")
jpeg(file_name, width = 3600, height = 3600, units = 'px', res = 300)
par(mar=c(10, 10, 4, 2),mgp=c(6,2,0),lwd=3)
for (g in 0:5) {
  if(g==0){
    plot(collapsed$S, fitIntxn, type="n", xlab="Sucrose Level", ylab="Log Odds of Post-Thaw Recovery", 
         main=paste0("Coloring by Level of Isoleucine"), cex.main=3.0, cex=2.25,cex.axis=4.0,cex.lab=4.0,ylim=c(-2,1.5))
    legend("topright", col=colors, paste0("Level ",0:5), lty=1, bty="n",cex=3.0)
  }
  else{
    points(collapsed$S, fitIntxn, type="n") 
  }
  for (i in 0:5) points(collapsed$S[collapsed$I==i & collapsed$G==g], fitIntxn[collapsed$I==i & collapsed$G==g], 
                        type="l", col=colors[i+1])
}
dev.off()