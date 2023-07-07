
setwd("/home/yi/Dropbox/EvolutionOfVariance/CodeForSubmission/")

load("EvoVariance.RData")
pdf("fig3.pdf", height= 6, width=12)
par(mfrow=c(1,2), mar=c(5,5,1,1))
beeswarm(res_data_frame$delta~as.factor(res_data_frame$reps), ylim=c(-0.2,0.3),xaxt="n", ylab=expression(delta), xlab="Repeated measures per individual",
         corral = "wrap", 
         method = "hex")
points(c(delta,delta,delta)~c(1:3), pch=19, cex=1.4, col="red")
points(c(mean_delta)~c(1:3), pch=19, cex=1.4)
abline(h=0, lty=2)
axis(1, c(1:3), c(5,10, 20))

beeswarm(res_data_frame$eta~as.factor(res_data_frame$reps), ylim=c(-1,0.2), xaxt="n", ylab=expression(eta), xlab="Repeated measures per individual",
         corral = "wrap",
         method = "hex")
points(c(eta,eta,eta)~c(1:3), pch=19, cex=1.4, col="red")
points(c(mean_eta)~c(1:3), pch=19, cex=1.4)
abline(h=0, lty=2)
axis(1, c(1:3), c(5,10,20))
dev.off()


