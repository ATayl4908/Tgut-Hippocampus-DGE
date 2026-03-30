#Mineralocorticoid receptor
MR.counts <- zbnorm$counts["NR3C2",]
names(MR.counts) <- groups
MR.counts <- stack(MR.counts)[2:1]
MR.counts <- setNames(MR.counts, c("Group","CPM"))
MR.counts <- MR.counts[sort.list(as.character(MR.counts$Group)),]
MR.counts$Group <- factor(MR.counts$Group, levels = c(as.character(unique(MR.counts$Group))))

plot(MR.counts, main = "NR3C2 (Mineralocorticoid Receptor)")
write.csv(MR.counts, file = "MR counts.csv", row.names = FALSE)

#Glucocorticoid receptor
GR.counts <- zbnorm$counts["NR3C1",]
names(GR.counts) <- groups
GR.counts <- stack(GR.counts)[2:1]
GR.counts <- setNames(GR.counts, c("Group","CPM"))
GR.counts <- GR.counts[sort.list(as.character(GR.counts$Group)),]
GR.counts$Group <- factor(GR.counts$Group, levels = c(as.character(unique(GR.counts$Group))))

plot(GR.counts, main = "NR3C1 (Glucocorticoid Receptor)")
ggplot(data = GR.counts, aes(x=Group, y = CPM)) + 
  geom_boxplot()+
  scale_x_discrete(guide=guide_axis(n.dodge = 2)) + 
  labs(title = "NR3C1 (Glucocorticoid Receptor)") +
  theme( plot.title=element_text( face = 'bold', hjust=0.5))
write.csv(MR.counts, file = "GR counts.csv", row.names = FALSE)

#Follistatin
FST.counts <- zbnorm$counts["FST",]
names(FST.counts) <- groups
FST.counts <- stack(FST.counts)[2:1]
FST.counts <- setNames(FST.counts, c("Group","CPM"))
FST.counts <- FST.counts[sort.list(as.character(FST.counts$Group)),]
FST.counts$Group <- factor(FST.counts$Group, levels = c(as.character(unique(FST.counts$Group))))

plot(FST.counts, main = "FST")
ggplot(data = FST.counts, aes(x=Group, y = CPM)) + 
  geom_boxplot()+
  scale_x_discrete(guide=guide_axis(n.dodge = 2)) + 
  labs(title = "Follistatin") +
  theme( plot.title=element_text( face = 'bold', hjust=0.5))
write.csv(FST.counts, file = "FST counts.csv", row.names = FALSE)

#Activin receptor 2A
ACVR2A.counts <- zbnorm$counts["ACVR2A",]
names(ACVR2A.counts) <- groups
ACVR2A.counts <- stack(ACVR2A.counts)[2:1]
ACVR2A.counts <- setNames(ACVR2A.counts, c("Group","CPM"))
ACVR2A.counts <- ACVR2A.counts[sort.list(as.character(ACVR2A.counts$Group)),]
ACVR2A.counts$Group <- factor(ACVR2A.counts$Group, levels = c(as.character(unique(ACVR2A.counts$Group))))

plot(ACVR2A.counts, main = "ACVR2A (Activin receptor 2A)")
ggplot(data = ACVR2A.counts, aes(x=Group, y = CPM)) + 
  geom_boxplot()+
  scale_x_discrete(guide=guide_axis(n.dodge = 2)) + 
  labs(title = "ACVR2A (Activin Receptor 2A)") +
  theme( plot.title=element_text( face = 'bold', hjust=0.5))
write.csv(FST.counts, file = "ACVR2A counts.csv", row.names = FALSE)
