# draws figures
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/HPC_Analyses/')

#Read in the data and the tree
df = readRDS('RScripts/R_Objects/metadata.RDS')

## FIGURE 1 is from anvio, so not generated here. 

## organise data
df<- subset(df, !is.na(GermTimeSimp))
df <- df %>%
  mutate(TypesControlledForNumber = Types/log(Number), species_ordered = fct_reorder(species_name, TypesControlledForNumber)) %>%
  mutate(Bottleneck = abs(as.numeric(Fission)-2))



# Figure 1b
df_fission<- df %>%
  mutate(Bottleneck = (Fission!=1))

bar_type_number<- ggplot(df_fission, aes(x=species_ordered, y=TypesControlledForNumber)) +
  geom_point(aes(shape = as.factor(Bottleneck), colour = GermNumeric), size = 2) + theme_few(base_family = "Helvetica") + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab('Number of Cell Types / log(Number of Cells)') + 
  scale_shape_discrete(
    labels = c('Absent', 'Present'),
    name = 'Strict Bottleneck') + 
  scale_colour_discrete(
    labels = c('None', 'Early', 'Late'),
    name = 'Germline Timing') + 
  coord_flip() + xlab('Species') 


pdf("figures/Figure1b.pdf", width = 8, height = 8)
bar_type_number
dev.off()

bar_type_number
ggsave('figures/figure1b.jpeg', width = 8, height = 8)


## Germline Figure (Figure 2)
cell_type_fig<- ggplot(df, aes(x=GermTimeSimp, y=log(Types))) + geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_few() + xlab(NULL) + ylab("Number of Cell Types (log)") + 
  scale_x_discrete(labels=c(
    "adult"="Late",
    "early"="Early",
    "no_germline"="No Germline"
  ))
cell_number_fig<- ggplot(df, aes(x=GermTimeSimp, y=log(Number))) + geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_few() + xlab(NULL) + ylab("Number of Cells (log)") + 
  scale_x_discrete(labels=c(
    "adult"="Late",
    "early"="Early",
    "no_germline"="No Germline"
  ))

cell_type_number_fig<- ggplot(df, aes(x=GermTimeSimp, y=TypesControlledForNumber)) + geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_few() + xlab(NULL) + ylab("Number of Cell Types / log(Number of Cells)") + 
  scale_x_discrete(labels=c(
    "adult"="Late",
    "early"="Early",
    "no_germline"="No Germline"
  ))

FigureGerm<- ggarrange(cell_number_fig, cell_type_fig,cell_type_number_fig, labels = "AUTO", common.legend = T, legend='right', nrow = 1)
FigureGerm<- annotate_figure(FigureGerm, bottom=text_grob("Germline Timing",hjust=0.8, vjust=0.5))
FigureGerm

pdf("figures/Figure2.pdf", width = 10, height = 8, useDingbats=FALSE )
FigureGerm
dev.off()

FigureGerm
ggsave('figures/figure2.jpeg', width = 10, height = 8)


## Bottleneck Figure (Figure 3)
cell_type_fig<- ggplot(df, aes(x=as.factor(Bottleneck), y=log(Types))) + geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_few() + xlab(NULL) + ylab("Number of Cell Types (log)") 

cell_number_fig<- ggplot(df, aes(x=as.factor(Bottleneck), y=log(Number))) + geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_few() + xlab(NULL) + ylab("Number of Cells (log)") 


cell_type_number_fig<- ggplot(df, aes(x=as.factor(Bottleneck), y=TypesControlledForNumber)) + geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(colour = Kingdom), position=position_jitter(width=0.175), alpha= 0.7) + 
  theme_few() + xlab(NULL) + ylab("Number of Cell Types / log(Number of Cells)")

FigureFission<- ggarrange(cell_number_fig, cell_type_fig,cell_type_number_fig, labels = "AUTO", common.legend = T, legend='right', nrow = 1)
FigureFission<- annotate_figure(FigureFission, bottom=text_grob("Presence of Generational Single-Cell Bottleneck",
                                                                hjust=0.8, vjust=0.5))
FigureFission

pdf("figures/Figure3.pdf", width = 10, height = 8)
FigureFission
dev.off()

FigureFission
ggsave('figures/figure3.jpeg', width = 10, height = 8)



# Correlation figure (Figure 4)
library(binom)

df_prob_fig<- subset(df, GermNumeric != 0) %>%
  mutate(GermNumeric = as.numeric(!(as.numeric(GermNumeric) -1))) %>%
  mutate(Bottleneck = abs(as.numeric(Fission)-2))

df_summary<- df_prob_fig %>%
#  mutate(GermNumeric = as.numeric(!(as.numeric(GermNumeric) -1))) %>%
#  mutate(Bottleneck = abs(as.numeric(Fission)-2)) %>%
  group_by(Bottleneck) %>%
  summarise(conf = binom.confint(sum(GermNumeric), length(GermNumeric), methods = 'ac')) 

df_summary<- (as.data.frame(df_summary))


Correlation_fig<- ggplot() + 
  geom_point(data = df_prob_fig, aes(x = as.factor(Bottleneck), y = GermNumeric),
             position = position_jitter(width = 0.13, height = 0.1),
             alpha = 0.8, size = 3, colour = 'grey40') + 
  geom_pointrange(data = df_summary, aes(x = as.factor(Bottleneck), y = conf$mean, ymin = conf$lower, ymax = conf$upper), size = 1) + 
  theme_few() + xlab('Strict Generational Single-Cell Bottleneck') + ylab('Probability(Early Germline Segregation)')

pdf("figures/Figure4.pdf")
Correlation_fig
dev.off()


Correlation_fig
ggsave('figures/figure4.jpeg', width = 8, height = 8)
