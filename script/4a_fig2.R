
# packages ----------------------------------------------------------------

library(tidyverse)
library(ggrepel)
library(ggraph)
library(cowplot)
library(patchwork)
library(ggsci)
library(prismatic, include.only = 'best_contrast')
library(showtext)

font_add('Helvetica', 
         "./script/fonts/Helvetica.ttf",
         bold = "./script/fonts/Helvetica-Bold.ttf")

remove(list = ls())

fill_color <- c("#20854EFF", alpha("#20854EFF", 0.3),
                "#E18727FF", alpha("#E18727FF", 0.3), 
               "#BC3C29FF", alpha("#BC3C29FF", 0.3))

load('./data/sars_2_cov.Rdata')


# infections vaccine and age ----------------------------------------------

datafile_info_BA1$lineage <- 'BA1'
datafile_info_BA2$lineage <- 'BA2'
datafile_info_BA2$location1 <- '厦门'
datafile_info_BA2$location2 <- '厦门'
datafile_info_Delta$lineage <- 'Delta'

datafile_info <- rbind(datafile_info_Delta,
                       datafile_info_BA1,
                       mutate(datafile_info_BA2, gender = NA))|> 
     select(age, vaccine, vaccine_type, lineage) |> 
     mutate(age_g = if_else(age <=18, 'c', 'a'),
            age_g = if_else(age >=65, 'o', age_g),
            age_g = factor(age_g, 
                           levels = c('c', 'a', 'o')),
            vaccine_mix = case_when(
                 is.na(vaccine_type) ~ 'U',
                 vaccine_type == "Combind" ~ 'M',
                 vaccine_type != "Combind" ~ 'U',
                 TRUE ~ as.character(vaccine_type)
            ),
            vaccine = factor(vaccine, 
                             levels = c(3, 2, 1, 0),
                             labels = c('A', 'B', 'C', 'C')))

datafile_percent <- datafile_info |> 
     group_by(age_g, lineage, vaccine, vaccine_mix) |> 
     count() |> 
     mutate(group = paste(vaccine, vaccine_mix, sep = "_"))

datafile_total <- datafile_info |> 
     group_by(age_g, lineage) |> 
     count()

datafile_label <- data.frame(x = rep(1, 3),
                             y = rep('a', 3),
                             label = rep('Number of infections', 3),
                             lineage = c('BA2', 'BA1', 'Delta'))


# infections plot ---------------------------------------------------------


fig_inf <- ggplot(filter(datafile_percent))+
     geom_bar(mapping = aes(y = age_g,
                            x = n, 
                            fill = group),
              color = 'black',
              stat="identity",
              position = "fill",
              show.legend = T)+
     geom_text(data = filter(datafile_total),
               mapping = aes(x = 1,
                             y = age_g,
                             label = n),
               vjust = -0.2,
               angle = -90,
               family = 'Helvetica')+
     geom_text(data = datafile_label,
               mapping = aes(x = x,
                             y = y,
                             label = label),
               vjust = -2,
               angle = -90,
               family = 'Helvetica')+
     facet_wrap(vars(factor(lineage, levels = c('BA2', 'BA1', 'Delta'))), 
                ncol = 1,
                scales = 'free_x',
                labeller = as_labeller(c(
                      Delta = 'c',
                      BA1 = 'b',
                      BA2 = 'a'
                )))+
     scale_y_discrete(labels = c('0-18', '19-64', '65-'),
                      expand = c(0, 0.6),
                      breaks = c('c', 'a', 'o'))+
     scale_x_continuous(n.breaks = 6,
                        expand =  expansion(add = c(0, 0)))+
     scale_fill_manual(values = fill_color,
                       labels = c('Booster Dose & Mixed', 
                                  'Booster Dose & Unmixed', 
                                  'Fully Vaccinated & Mixed', 
                                  'Fully Vaccinated & Unmixed',
                                  'Unfully Vaccinated & Mixed',
                                  'Unfully Vaccinated & Unmixed'),
                       na.translate = F)+
     coord_cartesian(clip = "off")+
     theme_classic(base_family = 'Helvetica')+
     labs(y = "",
          x = 'Proportions',
          fill = 'Vaccination status')

# contacts vaccine and age ------------------------------------------------

datafile_cont_BA1$lineage <- 'BA1'
datafile_cont_BA2$lineage <- 'BA2'
datafile_cont_Delta$lineage <- 'Delta'

datafile_cont <- rbind(datafile_cont_Delta,
                       datafile_cont_BA1,
                       mutate(datafile_cont_BA2, gender = NA))|> 
     select(age, vaccine, vaccine_mix, lineage) |> 
     mutate(age_g = if_else(age <=18, 'c', 'a'),
            age_g = if_else(age >=65, 'o', age_g),
            age_g = factor(age_g, levels = c('c', 'a', 'o')),
            vaccine_mix = if_else(is.na(vaccine_mix), 
                                  'N', 
                                  vaccine_mix)) |> 
     mutate(vaccine = factor(vaccine, 
                             levels = c(3, 2, 1, 0),
                             labels = c('A', 'B', 'C', 'C')),
            vaccine_mix = if_else(vaccine_mix == 'Y', 'M', 'U'))

datafile_percent <- datafile_cont |> 
     group_by(age_g, lineage, vaccine, vaccine_mix) |> 
     count() |> 
     filter(!is.na(age_g)) |> 
     mutate(group = paste(vaccine, vaccine_mix, sep = "_"))

datafile_total <- datafile_cont |> 
     group_by(age_g, lineage) |> 
     count() |> 
     filter(!is.na(age_g))

datafile_label <- data.frame(x = rep(1, 3),
                             y = rep('a', 3),
                             label = rep('Number of contacts', 3),
                             lineage = c('BA2', 'BA1', 'Delta'))

# contacts ----------------------------------------------------------------

fig_cont <- ggplot(filter(datafile_percent))+
     geom_bar(mapping = aes(y = age_g,
                            x = n, 
                            fill = group),
              color = 'black',
              stat="identity",
              position = "fill",
              show.legend = T)+
     geom_text(data = filter(datafile_total),
               mapping = aes(x = 1,
                             y = age_g,
                             label = n),
               vjust = -0.2,
               angle = -90,
               family = 'Helvetica')+
     geom_text(data = datafile_label,
               mapping = aes(x = x,
                             y = y,
                             label = label),
               vjust = -2,
               angle = -90,
               family = 'Helvetica')+
     facet_wrap(vars(factor(lineage, levels = c('BA2', 'BA1', 'Delta'))), 
                ncol = 1,
                scales = 'free_x',
                labeller = as_labeller(c(
                     Delta = 'f',
                     BA1 = 'e',
                     BA2 = 'd'
                )))+
     scale_y_discrete(labels = c('0-18', '19-64', '65-'),
                      expand = c(0, 0.6),
                      breaks = c('c', 'a', 'o'))+
     scale_x_continuous(n.breaks = 6,
                        expand =  expansion(add = c(0, 0)))+
     scale_fill_manual(values = fill_color,
                       labels = c('Booster Dose & Mixed', 
                                  'Booster Dose & Unmixed', 
                                  'Fully Vaccinated & Mixed', 
                                  'Fully Vaccinated & Unmixed',
                                  'Unfully Vaccinated & Mixed',
                                  'Unfully Vaccinated & Unmixed'),
                       na.translate = F)+
     coord_cartesian(clip = "off")+
     theme_classic(base_family = 'Helvetica')+
     labs(y = "",
          x = 'Proportions',
          fill = 'Vaccination status')

# combind plot ------------------------------------------------------------

fig_1 <- fig_inf + fig_cont +
     plot_layout(guides = 'collect')&
     theme(plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
           plot.margin = margin(0, 1, 0, 0, "cm"),
           axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.text.y = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
           strip.text = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
           strip.background = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.position = 'bottom')

fig_1

ggsave(filename = './outcome/publish/Figure 2.pdf',
       fig_1,
       width = 14, height = 7)
