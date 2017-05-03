library(ggplot2)

## run 1 ##
inhibv1<-read.delim("trail_inhib_data_v1.tsv", header=T)

pdf("inhibitor_assay_v1.pdf")

ggplot(inhibv1, aes(factor(conc), tot_count, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv1, aes(factor(conc), live_perc, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv1, aes(factor(conc), tot_count*live_perc, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv1, aes(factor(conc), tot_count/control_count, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv1, aes(factor(conc), live_perc/control_perc, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv1, aes(factor(conc), (tot_count*live_perc)/(control_count*control_perc), col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

dev.off()

## run 2 ##
inhibv2<-read.delim("trail_inhib_data_v2.tsv", header=T)

pdf("inhibitor_assay_v2.pdf")

ggplot(inhibv2, aes(factor(conc), tot_count, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv2, aes(factor(conc), live_perc, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv2, aes(factor(conc), tot_count*live_perc, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv2, aes(factor(conc), tot_count/control_count, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv2, aes(factor(conc), live_perc/control_perc, col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

ggplot(inhibv2, aes(factor(conc), (tot_count*live_perc)/(control_count*control_perc), col=factor(rep), group=factor(rep))) + geom_point() + geom_line() + facet_grid(.~drug, scale="free_x")

dev.off()
