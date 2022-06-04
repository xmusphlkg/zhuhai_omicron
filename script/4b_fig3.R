
# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(ggsci)
library(gridExtra)
library(showtext)

font_add('Helvetica', 
         "./script/fonts/Helvetica.ttf",
         bold = "./script/fonts/Helvetica-Bold.ttf")

remove(list = ls())

load('./data/sars_2_cov.Rdata')

fill_color <- rev(pal_nejm()(4))[-3]

set.seed(202202)

# function ----------------------------------------------------------------

datafile_info_BA1$lineage <- 'BA.1'
datafile_info_BA2$lineage <- 'BA.2'
datafile_info_BA2$location1 <- '厦门'
datafile_info_BA2$location2 <- '厦门'
datafile_info_Delta$lineage <- 'Delta'

datafile_info <- rbind(datafile_info_Delta,
                       datafile_info_BA1,
                       mutate(datafile_info_BA2, gender = NA))|> 
        mutate(vaccine_mix = case_when(
                is.na(vaccine_type) ~ 'U',
                vaccine_type == "Combind" ~ 'M',
                vaccine_type != "Combind" ~ 'U',
                TRUE ~ as.character(vaccine_type)
        ),
        vaccine = factor(vaccine, 
                         levels = c(3, 2, 1, 0),
                         labels = c('A', 'B', 'C', 'C')),
        lineage = factor(lineage,
                         levels = c('Delta', 'BA.1', 'BA.2')),
        onset_positive = round(as.numeric(datepositive - dateonset), 0)) |> 
        filter(type != 'Asymptomatic')

# plot --------------------------------------------------------------------

fig_a <- ggplot(data = datafile_info,
                mapping = aes(x = onset_positive,
                              color = lineage))+
        stat_density(adjust = 1.8,
                     geom = 'line',
                     position="identity")+
        scale_color_manual(values = fill_color,
                           na.translate = F)+
        scale_y_continuous(expand = c(0, 0),
                           limits = c(0, 0.5),
                           labels = label_number(accuracy = 0.01))+
        scale_x_continuous(breaks = seq(0, 8, 2),
                           expand = c(0, 0))+
        theme_classic(base_family = 'Helvetica')+
        theme(axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 0.5, face = 'plain', color = 'black'),
              axis.text.y = element_text(size = 10, hjust = 1, vjust = .5, face = 'plain', color = 'black'),
              axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
              axis.title.y = element_text(size = 12, hjust = .5, vjust = 1, face = 'bold'),
              strip.text.x = element_text(size = 12, face = 'bold'),
              strip.text.y = element_text(size = 12, face = 'bold'),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.position = c(0.99, 0.99),
              legend.justification = c(1, 1),
              plot.margin = margin(5, 5, 5, 5),
              plot.background = element_blank(),
              plot.title = element_text(size = 16, hjust = 0, vjust = 0, face = 'bold'))+
        labs(x = 'Time from onset to testing (days)',
             y = 'Probability density functions',
             color = 'VOC',
             title = '')

# ggplot(data = datafile_info,
#        mapping = aes(x = onset_positive,
#                      y = Throat_N,
#                      color = lineage))+
#      geom_point()+
#      geom_smooth(method = 'glm')

fig_a

ggsave(filename = './outcome/publish/extend/Figure S3.pdf', 
       height = 3.5, width = 3.5)

ggsave(filename = './outcome/publish/extend/Figure S3.tiff', 
       fig_a,
       height = 3.5, width = 3.5)
