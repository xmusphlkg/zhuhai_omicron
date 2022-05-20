
# packages ----------------------------------------------------------------

library(tidyverse)
library(ggraph)
library(cowplot)
library(patchwork)
library(ggsci)
library(showtext)

font_add('Helvetica', 
         "./script/fonts/Helvetica.ttf",
         bold = "./script/fonts/Helvetica-Bold.ttf")

remove(list = ls())

fill_color <- c("#20854EFF", alpha("#20854EFF", 0.3),
                "#E18727FF", alpha("#E18727FF", 0.3), 
                "#BC3C29FF", alpha("#BC3C29FF", 0.3))

load('./data/sars_2_cov.Rdata')

datafile_demo_delta <- data.frame(
     age_g = c('[0,15)', '[15,60)', '[60,Inf)'),
     prop = c(0.1952, 0.6060, 0.1988),
     lineage = 'Delta',
     type = 'demo'
     )

datafile_demo_ba1 <- data.frame(
     age_g = c('[0,15)', '[15,60)', '[60,Inf)'),
     prop = c(0.1731, 0.7866, 0.0403),
     lineage = 'BA.1',
     type = 'demo'
     )

datafile_demo_ba2 <- data.frame(
     age_g = c('[0,15)', '[15,60)', '[60,Inf)'),
     prop = c(0.1716, 0.7328, 0.0956),
     lineage = 'BA.2',
     type = 'demo'
     )

# contact age distribution ------------------------------------------------


datafile_info_BA1$lineage <- 'BA1'
datafile_info_BA2$lineage <- 'BA2'
datafile_info_BA2$location1 <- '厦门'
datafile_info_BA2$location2 <- '厦门'
datafile_info_Delta$lineage <- 'Delta'

datafile_age_delta <- datafile_cont_Delta |> 
     mutate(age_g = cut(age,
                breaks = c(0, 15, 60, Inf),
                right = F)) |> 
     filter(!is.na(age_g)) |> 
     group_by(age_g) |> 
     count() |> 
     ungroup() |> 
     mutate(prop = n/sum(n),
            lineage = 'Delta',
            type = 'contact') |> 
     select(-n)

datafile_age_ba1 <- datafile_cont_BA1 |> 
     mutate(age_g = cut(age,
                        breaks = c(0, 15, 60, Inf),
                        right = F)) |> 
     filter(!is.na(age_g)) |> 
     group_by(age_g) |> 
     count() |> 
     ungroup() |> 
     mutate(prop = n/sum(n),
            lineage = 'BA.1',
            type = 'contact') |> 
     select(-n)

datafile_age_ba2 <- datafile_cont_BA2 |> 
     mutate(age_g = cut(age,
                        breaks = c(0, 15, 60, Inf),
                        right = F)) |> 
     filter(!is.na(age_g)) |> 
     group_by(age_g) |> 
     count() |> 
     ungroup() |> 
     mutate(prop = n/sum(n),
            lineage = 'BA.2',
            type = 'contact') |> 
     select(-n)

datafile_age <- rbind(datafile_age_ba1,
                      datafile_age_ba2,
                      datafile_age_delta,
                      datafile_demo_delta,
                      datafile_demo_ba1,
                      datafile_demo_ba2)

# plot --------------------------------------------------------------------

ggplot(data = datafile_age)+
     geom_bar(mapping = aes(x = age_g,
                            y = prop,
                            fill = type),
              color = 'black',
              stat="identity",
              position = "dodge")+
        facet_wrap(vars(factor(lineage, levels = c('BA.2', 'BA.1', 'Delta'))),
                   nrow = 1,
                   scales = 'free',
                   labeller = as_labeller(c(
                           Delta = 'c',
                           BA.1 = 'b',
                           BA.2 = 'a'
                   )))+
        scale_y_continuous(n.breaks = 6,
                           expand =  expansion(add = c(0, 0)),
                           limits = c(0, 1))+
        scale_x_discrete(labels = c('0-14', '15-59', '60-'))+
        scale_fill_nejm(labels = c('Contacts', 'Population'))+
        theme_classic(base_family = 'Helvetica')+
        theme(plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
              plot.margin = margin(0, 1, 0, 1, "cm"),
              axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
              axis.text.y = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
              axis.title.y = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
              strip.text = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
              strip.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = 'bottom')+
        labs(y = "Propotions",
             x = '',
             fill = 'Age distribution')

ggsave('./outcome/publish/extend/Figure S6.pdf',
       width = 7,
       height = 3)
