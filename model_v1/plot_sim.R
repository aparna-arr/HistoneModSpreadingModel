library(ggplot2)
library(reshape2)

sim<-read.delim("results.txt",header=T)
sim.m<-melt(sim,id.vars=c('Event'))

#divs<-c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000)
divs<-c(0, 200000, 400000, 600000, 800000)

pdf("sim.pdf", width=20)
ggplot(sim.m, aes(Event,value,col=variable)) + geom_line() + geom_vline(data=melt(divs),aes(xintercept=value))
dev.off()
