# draws figures
library(ggplot2)
library(ggpubr)
library(tidyverse)
setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/HPC_Analyses/')

#Read in the data and the tree
df = readRDS('RScripts/R_Objects/metadata.RDS')


df<- subset(df, !is.na(GermTimeSimp))

cell_type_fig<- ggplot(df, aes(x=GermTimeSimp, y=log(Types))) + geom_boxplot() + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_minimal() + xlab(NULL) + ylab("Number of Cell Types (log)") + 
  scale_x_discrete(labels=c(
    "adult"="Late",
    "early"="Early",
    "no_germline"="No Germline"
  ))
cell_number_fig<- ggplot(df, aes(x=GermTimeSimp, y=log(Number))) + geom_boxplot() + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_minimal() + xlab(NULL) + ylab("Number of Cells (log)") + 
  scale_x_discrete(labels=c(
    "adult"="Late",
    "early"="Early",
    "no_germline"="No Germline"
  ))
FigureGerm<- ggarrange(cell_number_fig, cell_type_fig, labels = "AUTO", common.legend = T, legend='right')
FigureGerm<- annotate_figure(FigureGerm, bottom=text_grob("Germline Timing",hjust=0.8, vjust=0.5))
FigureGerm

pdf("figures/FigureGerm.pdf")
FigureGerm
dev.off()

cell_type_fig<- ggplot(df, aes(x=Fission, y=log(Types))) + geom_boxplot() + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_minimal() + xlab(NULL) + ylab("Number of Cell Types (log)") + 
  scale_x_discrete(labels=c(
    "adult"="Late",
    "early"="Early",
    "no_germline"="No Germline"
  ))
cell_number_fig<- ggplot(df, aes(x=Fission, y=log(Number))) + geom_boxplot() + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_minimal() + xlab(NULL) + ylab("Number of Cells (log)") + 
  scale_x_discrete(labels=c(
    "adult"="Late",
    "early"="Early",
    "no_germline"="No Germline"
  ))

FigureFission<- ggarrange(cell_number_fig, cell_type_fig, labels = "AUTO", common.legend = T, legend='right')
FigureFission<- annotate_figure(FigureFission, bottom=text_grob("Presence of Fissiparous Reproduction",
                                                                hjust=0.8, vjust=0.5))
FigureFission

pdf("figures/FigureFission.pdf")
FigureFission
dev.off()

# trying to build the figure Ashleigh is after
df <- df %>%
  mutate(TypesControlledForNumber = Types/log(Number), species_ordered = fct_reorder(species_name, TypesControlledForNumber)) 


bar_fission<- ggplot(df, aes(x=species_ordered, y=1, fill = Fission)) +
  geom_bar(stat='identity') + theme_minimal() + 
  theme(
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  coord_flip()

bar_germ<- ggplot(df, aes(x=species_ordered, y=1, fill = GermTimeSimp)) +
  geom_bar(stat='identity') + theme_minimal() + 
  theme(
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) + 
  coord_flip()

bar_number <- ggplot(df, aes(x=species_ordered, y=log(Number))) +
  geom_bar(stat='identity', fill = 'blue') + theme_minimal() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) + 
  coord_flip()

bar_type<- ggplot(df, aes(x=species_ordered, y=log(Types))) +
  geom_bar(stat='identity', fill = 'blue') + theme_minimal() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  coord_flip()

bar_type_number<- ggplot(df, aes(x=species_ordered, y=TypesControlledForNumber)) +
  geom_bar(stat='identity', fill = 'blue') + theme_minimal() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  coord_flip()

FigureBars<- ggarrange(bar_fission, bar_germ, bar_type, bar_number, bar_type_number, nrow=1, ncol = 5, align = 'h')

pdf("figures/FigureBars.pdf", width = 15, height = 10)
FigureBars
dev.off()
