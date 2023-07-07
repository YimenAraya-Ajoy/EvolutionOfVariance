require(colorRamps)
require(lattice)
require(ggpubr)
library(gridExtra)
library(grid)
library(gridGraphics)

setwd("/home/yi/Dropbox/EvolutionOfVariance/CodeForSubmission/")

theme.novpadding1 <-list(axis.line = list(col = "transparent"), 
                         layout.heights = list(top.padding = 2,
                                               main.key.padding = 0,
                                               key.axis.padding = 0,
                                               axis.xlab.padding = 0,
                                               xlab.key.padding = 0,
                                               key.sub.padding = 0,
                                               bottom.padding = 0),
                         clip =list(panel="off"),
                         layout.widths =
                           list(left.padding = 0,
                                key.ylab.padding = 0,
                                ylab.axis.padding = 0,
                                axis.key.padding = 0,
                                right.padding = 0))

col<-matlab.like(100)

t<-theme_bw() + theme(text = element_text(size = 10), panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")   



I<-seq(-0.5,0.5,0.1)
resI<-seq(0, 1.16, by=0.04)

dI<-expand.grid(I,resI)
colnames(dI)<-c("I", "resI")

beta<-0
alpha<-0/var(I^2)
eta<--0.3/var(resI)
beta/-(2*alpha)

dI$w<-(beta*(dI$I-mean(dI$I)) + alpha*(dI$I-mean(dI$I))^2 + eta*dI$resI)
dI$w2<-(dI$w+(min(dI$w)*-1))/mean(dI$w+(min(dI$w)*-1))

p1a<-wireframe(dI$w2~dI$resI+dI$I, zlab="w", ylab=expression("i"), xlab=expression(bar("r"^2)),  colorkey = FALSE,
               screen = list(z = -50, x = -60), 
               scales = list( arrows = TRUE, col="black", cex.axis=0.5, axis=list(text=list(cex=0.7))), 
               aspect = c(1.1, 1),
               drape = TRUE, 
               col.regions = col,
               par.settings = theme.novpadding1, main="(A)")


var(I)/(var(I) + mean(resI))
vIp<-var(I) + cov(dI$w,(dI$I-mean(dI$I))^2)
vRp<-mean(resI) + cov(dI$w, dI$resI)
vIp/(vRp+vIp)

n=4000
dfa1<-data.frame(z=rnorm(n, 0, sqrt(var(I)+mean(resI))), Selection="A", V="Total")
dfa2<-data.frame(z=rnorm(n, 0, sqrt(var(I))), Selection="A", V="Among")
dfa3<-data.frame(z=rnorm(n, 0, sqrt(mean(resI))), Selection="A", V="Within")

dfb1<-data.frame(z=rnorm(n, 0, sqrt(vIp+vRp)), Selection="B", V="Total")
dfb2<-data.frame(z=rnorm(n, 0, sqrt(vIp)), Selection="B", V="Among")
dfb3<-data.frame(z=rnorm(n, 0, sqrt(vRp)), Selection="B", V="Within")


df1<-rbind(dfa2,dfb2,dfa3,dfb3)


p1b <- ggplot(df1, aes(x=z, linetype=Selection, color=V)) + geom_density(aes(y = after_stat(count)))+ t +   scale_color_manual(values=c("orange", "darkgreen")) + ylim(0,6000) + xlab("Phenotype (z)") + ylab("Frequency")



#####
beta<-0
alpha<- 0.2/var(I^2)
eta<- 0/var(resI)
beta/-(2*alpha)

dI$w<-(beta*(dI$I-mean(dI$I)) + alpha*(dI$I-mean(dI$I))^2 + eta*dI$resI)
dI$w2<-(dI$w+(min(dI$w)*-1))/mean(dI$w+(min(dI$w)*-1))

p2a<-wireframe(dI$w2~dI$resI+dI$I, zlab="w", ylab="i", xlab=expression(bar("r"^2)),  colorkey = FALSE,
               screen = list(z = -50, x = -60), 
               scales = list( arrows = TRUE, col="black", cex.axis=0.5, axis=list(text=list(cex=0.7))), 
               aspect = c(1.1, 1),
               drape = TRUE, 
               col.regions = col,
               par.settings = theme.novpadding1, main="(B)")


var(I)/(var(I) + mean(resI))
vIp<-var(I) + cov(dI$w,(dI$I-mean(dI$I))^2)
vRp<-mean(resI) + cov(dI$w, dI$resI)
vIp/(vRp+vIp)


dfa1<-data.frame(z=rnorm(n, 0, sqrt(var(I)+mean(resI))), Selection="A", V="Total")
dfa2<-data.frame(z=rnorm(n, 0, sqrt(var(I))), Selection="A", V="Among")
dfa3<-data.frame(z=rnorm(n, 0, sqrt(mean(resI))), Selection="A", V="Within")

dfb1<-data.frame(z=rnorm(n, 0, sqrt(vIp+vRp)), Selection="B", V="Total")
dfb2<-data.frame(z=rnorm(n, 0, sqrt(vIp)), Selection="B", V="Among")
dfb3<-data.frame(z=rnorm(n, 0, sqrt(vRp)), Selection="B", V="Within")


df1<-rbind(dfa2,dfb2,dfa3,dfb3)


p2b<-ggplot(df1, aes(x=z, linetype=Selection, color=V)) +
  geom_density(aes(y = after_stat(count)))+ t +   scale_color_manual(values=c("orange", "darkgreen")) + ylim(0,6000) + xlab("Phenotype (z)") + ylab("Frequency")
#####
beta<-0
alpha<-0.2/var(I^2)
eta<--0.2/var(resI)
beta/-(2*alpha)

dI$w<-(beta*(dI$I-mean(dI$I)) + alpha*(dI$I-mean(dI$I))^2 + eta*dI$resI)
dI$w2<-(dI$w+(min(dI$w)*-1))/mean(dI$w+(min(dI$w)*-1))

p3a<-wireframe(dI$w2~dI$resI+dI$I, zlab="w", ylab="i", xlab=expression(bar("r"^2)),  colorkey = FALSE,
               screen = list(z = -50, x = -60), 
               scales = list( arrows = TRUE, col="black", cex.axis=0.5, axis=list(text=list(cex=0.7))), 
               aspect = c(1.1, 1),
               drape = TRUE, 
               col.regions = col,
               par.settings = theme.novpadding1, main="(B)")


var(I)/(var(I) + mean(resI))
vIp<-var(I) + cov(dI$w,(dI$I-mean(dI$I))^2)
vRp<-mean(resI) + cov(dI$w, dI$resI)
vIp/(vRp+vIp)


dfa1<-data.frame(z=rnorm(n, 0, sqrt(var(I)+mean(resI))), Selection="A", V="Total")
dfa2<-data.frame(z=rnorm(n, 0, sqrt(var(I))), Selection="A", V="Among")
dfa3<-data.frame(z=rnorm(n, 0, sqrt(mean(resI))), Selection="A", V="Within")

dfb1<-data.frame(z=rnorm(n, 0, sqrt(vIp+vRp)), Selection="B", V="Total")
dfb2<-data.frame(z=rnorm(n, 0, sqrt(vIp)), Selection="B", V="Among")
dfb3<-data.frame(z=rnorm(n, 0, sqrt(vRp)), Selection="B", V="Within")


df1<-rbind(dfa2,dfb2,dfa3,dfb3)


p3b<-ggplot(df1, aes(x=z, linetype=Selection, color=V)) +
  geom_density(aes(y = after_stat(count)))+ t +   scale_color_manual(values=c("orange", "darkgreen")) + ylim(0,6000) + xlab("Phenotype (z)") + ylab("Frequency")


#####
beta<-0
alpha<--0.05/var(I^2)
eta<--0.2/var(resI)
beta/-(2*alpha)

dI$w<-(beta*(dI$I-mean(dI$I)) + alpha*(dI$I-mean(dI$I))^2 + eta*dI$resI)
dI$w2<-(dI$w+(min(dI$w)*-1))/mean(dI$w+(min(dI$w)*-1))

p4a<-wireframe(dI$w2~dI$resI+dI$I, zlab="w", ylab="i", xlab=expression(bar("r"^2)),  colorkey = FALSE,
               screen = list(z = -50, x = -60), 
               scales = list( arrows = TRUE, col="black", cex.axis=0.5, axis=list(text=list(cex=0.7))), 
               aspect = c(1.1, 1),
               drape = TRUE, 
               col.regions = col,
               par.settings = theme.novpadding1, main="(C)")


var(I)/(var(I) + mean(resI))
vIp<-var(I) + cov(dI$w,(dI$I-mean(dI$I))^2)
vRp<-mean(resI) + cov(dI$w, dI$resI)
vIp/(vRp+vIp)


dfa1<-data.frame(z=rnorm(n, 0, sqrt(var(I)+mean(resI))), Selection="A", V="Total")
dfa2<-data.frame(z=rnorm(n, 0, sqrt(var(I))), Selection="A", V="Among")
dfa3<-data.frame(z=rnorm(n, 0, sqrt(mean(resI))), Selection="A", V="Within")

dfb1<-data.frame(z=rnorm(n, 0, sqrt(vIp+vRp)), Selection="B", V="Total")
dfb2<-data.frame(z=rnorm(n, 0, sqrt(vIp)), Selection="B", V="Among")
dfb3<-data.frame(z=rnorm(n, 0, sqrt(vRp)), Selection="B", V="Within")


df1<-rbind(dfa2,dfb2,dfa3,dfb3)


p4b<-ggplot(df1, aes(x=z, linetype=Selection, color=V)) +
  geom_density(aes(y = after_stat(count)))+ t +   scale_color_manual(values=c("orange", "darkgreen")) + ylim(0,6000) + xlab("Phenotype (z)") + ylab("Frequency")

pdf("fig2.pdf", height= 4, width=6)
grid.arrange(p1a, p3a, p4a, p1b, p3b, p4b, ncol=3, nrow=2, heights=c(1.25,1))
dev.off()
