require(colorRamps)
require(lattice)
require(ggpubr)
library(gridExtra)
library(grid)
library(gridGraphics)

setwd("/home/yi/Dropbox/EvolutionOfVariance/CodeForSubmission/")

s=2
n.ind<-50
n.reps=4000
I<-rnorm(n.ind,20,sqrt(0.5))
resI<-rnorm(n.ind,1,sqrt(0.01))
d<-data.frame(ID=rep(1:n.ind,n.reps), I=rep(I,n.reps), resI=rep(resI,n.reps))
d$z<-d$I + rnorm(nrow(d),0,sqrt(d$resI))
d$type<-"A"
d$selection<-"A"

d2<-data.frame(ID=10000, I=20, resI=NA ,z=rnorm(n.reps*4,20,sqrt(1.5)), type="B", selection="A")
d3<-data.frame(ID=20000, I=20, resI=NA ,z=rnorm(n.reps,  20, sqrt(1)), type="C",  selection="A")
d4<-data.frame(ID=30000, I=20, resI=NA ,z=rnorm(n.reps,  20, sqrt(0.5)), type="D", selection="A")
d5<-rbind(d,d2,d3,d4)


I<-rnorm(n.ind,25,sqrt(0.1))
resI<-rnorm(n.ind,0.4,sqrt(0.01))
db<-data.frame(ID=rep(n.ind:(2*n.ind-1), n.reps), I=rep(I,n.reps), resI=rep(resI,n.reps))
db$z<-db$I + rnorm(nrow(db),0,sqrt(db$resI))
db$type<-"A"
db$selection<-"B"

db2<-data.frame(ID=40000, I=25, resI=NA , z=rnorm(n.reps*4, 25, sqrt(0.5)), type="B", selection="B")
db3<-data.frame(ID=60000, I=25, resI=NA , z=rnorm(n.reps, 25, sqrt(0.4)), type="C", selection="B")
db4<-data.frame(ID=80000, I=25, resI=NA , z=rnorm(n.reps, 25, sqrt(0.1)), type="D", selection="B")

db5<-rbind(db,db2,db3,db4)
d6<-rbind(db5,d5)

an1<- annotate("text", x = 20, y=6000, label = "Before selection", size=s)
an2<- annotate("text", x = 25, y=10000, label = "After selection", size=s)

p1 <-ggplot(d6, aes(x=z, group=as.factor(ID),linetype=selection)) +
  geom_density(aes(y = after_stat(count), colour=type)) + scale_color_manual(values=c("grey80", "black", "darkgreen",  "orange")) + t + ylim(0,10000) + xlim(15,30) + ggtitle("(A)") + xlab("Phenotype (z)") + ylab("Frequency") + an1 + an2

###
I<-rnorm(n.ind,25,sqrt(1.5))
resI<-rnorm(n.ind,1.5,sqrt(0.01))
dc<-data.frame(ID=rep(n.ind:(2*n.ind-1), n.reps), I=rep(I,n.reps), resI=rep(resI,n.reps))
dc$z<-dc$I + rnorm(nrow(dc),0,sqrt(dc$resI))
dc$type<-"A"
dc$selection<-"B"

dc2<-data.frame(ID=40000, I=25, resI=NA , z=rnorm(n.reps*4, 25, sqrt(3)), type="B", selection="B")
dc3<-data.frame(ID=60000, I=25, resI=NA , z=rnorm(n.reps, 25, sqrt(1.5)), type="C", selection="B")
dc4<-data.frame(ID=80000, I=25, resI=NA , z=rnorm(n.reps, 25, sqrt(1.5)), type="D", selection="B")

dc5<-rbind(dc,dc2,dc3,dc4)
d7<-rbind(dc5,d5)

an1<- annotate("text", x = 20, y=6000, label = "Before selection", size=s)
an2<- annotate("text", x = 25, y=4500, label = "After selection", size=s)

p2 <-ggplot(d7, aes(x=z, group=as.factor(ID),linetype=selection)) +
  geom_density(aes(y = after_stat(count), colour=type)) + scale_color_manual(values=c("grey80", "black", "darkgreen",  "orange")) + t + ylim(0,10000) + xlim(15,30) + ggtitle("(B)") + xlab("Phenotype (z)") + ylab("Frequency") + an1 + an2

p2

###
I<-rnorm(n.ind,25,sqrt(1))
resI<-rnorm(n.ind,0.5,sqrt(0.02))
dd<-data.frame(ID=rep(n.ind:(2*n.ind-1), n.reps), I=rep(I,n.reps), resI=rep(resI,n.reps))
dd$z<-dd$I + rnorm(nrow(dd),0,sqrt(dd$resI))
dd$type<-"A"
dd$selection<-"B"

dd2<-data.frame(ID=40000, I=25, resI=NA , z=rnorm(n.reps*4, 25, sqrt(1.5)), type="B", selection="B")
dd3<-data.frame(ID=60000, I=25, resI=NA , z=rnorm(n.reps, 25, sqrt(0.5)), type="C", selection="B")
dd4<-data.frame(ID=80000, I=25, resI=NA , z=rnorm(n.reps, 25, sqrt(1)), type="D", selection="B")

dd5<-rbind(dd,dd2,dd3,dd4)
d8<-rbind(dd5,d5)

an1<- annotate("text", x = 20, y=6000, label = "Before selection", size=s)
an2<- annotate("text", x = 25, y=6000, label = "After selection", size=s)
p3 <-ggplot(d8, aes(x=z, group=as.factor(ID),linetype=selection)) +
  geom_density(aes(y = after_stat(count), colour=type)) + scale_color_manual(values=c("grey80", "black", "darkgreen",  "orange")) + t +  ylim(0,10000) + xlim(15,30) + ggtitle("(C)") + xlab("Phenotype (z)") + ylab("Frequency") + an1 + an2
p3

pdf("fig1.pdf", height= 6, width=3)
grid.arrange(p1,p2,p3, ncol=1, nrow=3)
dev.off()


