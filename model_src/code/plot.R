library(ggplot2)
library(reshape2)

files<-commandArgs(trailingOnly=TRUE)
dat5k<-read.delim(files[1], header=T)
dat3k<-read.delim(files[2], header=T)
dat1k<-read.delim(files[3], header=T)
outfilebase<-files[4]
#head(dat)
## Gap score plot over time
## 1000 RT
recruit <- c(30000, 31000)
recruit.df <- data.frame(recruit)
pdf(paste0(outfilebase, "rt_1000_gapscore.pdf"))

ggplot(dat1k[dat1k$RecruitTime == 1000,], aes(Timestep, AvgGapScore, col=factor(FValue))) + geom_line() + geom_vline(aes(xintercept=recruit), data=recruit.df, col="red", linetype="dashed") + geom_hline(aes(yintercept=0), linetype="dashed") + ylab("AvgGapScore (M - A) / (M + A)") + ggtitle("Gap Scores Over Time: RT = 1000 ( 1 cell cycle )")

dev.off()

## 3000 RT
recruit <- c(30000, 33000)
recruit.df <- data.frame(recruit)
pdf(paste0(outfilebase, "rt_3000_gapscore.pdf"))

ggplot(dat3k[dat3k$RecruitTime == 3000,], aes(Timestep, AvgGapScore, col=factor(FValue))) + geom_line() + geom_vline(aes(xintercept=recruit), data=recruit.df, col="red", linetype="dashed") + geom_hline(aes(yintercept=0), linetype="dashed") + ylab("AvgGapScore (M - A) / (M + A)") + ggtitle("Gap Scores Over Time: RT = 3000 ( 3 cell cycles )")

dev.off()

## 5000 RT
recruit <- c(30000, 35000)
recruit.df <- data.frame(recruit)

pdf(paste0(outfilebase, "rt_5000_gapscore.pdf"))

ggplot(dat5k[dat5k$RecruitTime == 5000,], aes(Timestep, AvgGapScore, col=factor(FValue))) + geom_line() + geom_vline(aes(xintercept=recruit), data=recruit.df, col="red", linetype="dashed") + geom_hline(aes(yintercept=0), linetype="dashed") + ylab("AvgGapScore (M - A) / (M + A)") + ggtitle("Gap Scores Over Time: RT = 5000 ( 5 cell cycles )")

dev.off()


## Transcription score plot over time
## 1000 RT
recruit <- c(30000, 31000)
recruit.df <- data.frame(recruit)
pdf(paste0(outfilebase, "rt_1000_transcriptionscore.pdf"))

ggplot(dat1k[dat1k$RecruitTime == 1000,], aes(Timestep, PercOff, col=factor(FValue))) + geom_line() + geom_vline(aes(xintercept=recruit), data=recruit.df, col="red", linetype="dashed") + geom_hline(aes(yintercept=0), linetype="dashed") + ylab("% cells off (>=60% M)") + ggtitle("% Off Over Time: RT = 1000 ( 1 cell cycle )")

dev.off()

## 3000 RT
recruit <- c(30000, 33000)
recruit.df <- data.frame(recruit)
pdf(paste0(outfilebase, "rt_3000_transcriptionscore.pdf"))

ggplot(dat3k[dat3k$RecruitTime == 3000,], aes(Timestep, PercOff, col=factor(FValue))) + geom_line() + geom_vline(aes(xintercept=recruit), data=recruit.df, col="red", linetype="dashed") + geom_hline(aes(yintercept=0), linetype="dashed") + ylab("% cells off (>=60% M)") + ggtitle("% Off Over Time: RT = 3000 ( 3 cell cycles )")

dev.off()

## 5000 RT
recruit <- c(30000, 35000)
recruit.df <- data.frame(recruit)

pdf(paste0(outfilebase, "rt_5000_transcriptionscore.pdf"))

ggplot(dat5k[dat5k$RecruitTime == 5000,], aes(Timestep, PercOff, col=factor(FValue))) + geom_line() + geom_vline(aes(xintercept=recruit), data=recruit.df, col="red", linetype="dashed") + geom_hline(aes(yintercept=0), linetype="dashed") + ylab("% cells off (>= 60% M)") + ggtitle("% Off Over Time: RT = 5000 ( 5 cell cycles )")

dev.off()


