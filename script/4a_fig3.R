
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

fill_color <- c(alpha("#BC3C29FF", 0.3), "#BC3C29FF", 
                alpha("#E18727FF", 0.3), "#E18727FF",
                alpha("#20854EFF", 0.3), "#20854EFF" 
)
fill_color_1 <- pal_nejm()(8)[c(2, 5:8)]

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
     ungroup() |> 
     mutate(group = paste(vaccine, vaccine_mix, sep = "_"),
            group = factor(group,
                           levels = c("C_U", "C_M", "B_U", "B_M", "A_U", "A_M")))

datafile_total <- datafile_info |> 
     group_by(age_g, lineage) |> 
     count()

datafile_label <- data.frame(x = rep(1, 3),
                             y = rep('a', 3),
                             label = rep('Number of infections', 3),
                             title = c('Omicron BA.2 infections',
                                       'Omicron BA.1 infections',
                                       'Delta infections'),
                             lineage = c('BA2', 'BA1', 'Delta'))


# infections plot ---------------------------------------------------------

fig_inf <- ggplot(datafile_percent)+
     geom_bar(mapping = aes(y = age_g,
                            x = n, 
                            fill = group),
              color = 'black',
              stat="identity",
              position = position_fill(reverse = T),
              show.legend = T)+
     geom_text(data = filter(datafile_total),
               mapping = aes(x = 1,
                             y = age_g,
                             label = n),
               vjust = -0.6,
               angle = -90,
               family = 'Helvetica')+
     geom_text(data = datafile_label,
               mapping = aes(x = x,
                             y = y,
                             label = label),
               vjust = -2,
               angle = -90,
               family = 'Helvetica')+
     geom_text(data = datafile_label,
               mapping = aes(x = 0.5,
                             y = Inf,
                             label = title),
               vjust = -0.6,
               family = 'Helvetica')+
     facet_wrap(vars(factor(lineage, levels = c('BA2', 'BA1', 'Delta'))), 
                ncol = 1,
                scales = 'free_x',
                labeller = as_labeller(c(
                     Delta = 'e',
                     BA1 = 'c',
                     BA2 = 'a'
                )))+
     scale_y_discrete(labels = c('0-18', '19-64', '65-'),
                      expand = c(0, 0.6),
                      breaks = c('c', 'a', 'o'))+
     scale_x_continuous(n.breaks = 6,
                        expand =  expansion(add = c(0, 0)))+
     scale_fill_manual(values = fill_color,
                       labels = c('Unfully Vaccinated & Unmixed', 
                                  'Unfully Vaccinated & Mixed',
                                  'Fully Vaccinated & Unmixed',
                                  'Fully Vaccinated & Mixed', 
                                  'Booster Dose & Unmixed', 
                                  'Booster Dose & Mixed'
                       ),
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
                       datafile_cont_BA2)|> 
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
     ungroup() |> 
     mutate(group = paste(vaccine, vaccine_mix, sep = "_"),
            group = factor(group,
                           levels = c("C_U", "C_M", "B_U", "B_M", "A_U", "A_M")))

datafile_total <- datafile_cont |> 
     group_by(age_g, lineage) |> 
     count() |> 
     filter(!is.na(age_g))

datafile_label <- data.frame(x = rep(1, 3),
                             y = rep('a', 3),
                             label = rep('Number of contacts', 3),
                             title = c('Contacts of Omicron BA.2 infections',
                                       'Contacts of Omicron BA.1 infections',
                                       'Contacts of Delta infections'),
                             lineage = c('BA2', 'BA1', 'Delta'))

# contacts ----------------------------------------------------------------

fig_cont <- ggplot(datafile_percent)+
     geom_bar(mapping = aes(y = age_g,
                            x = n, 
                            fill = group),
              color = 'black',
              stat="identity",
              position = position_fill(reverse = T),
              show.legend = T)+
     geom_text(data = filter(datafile_total),
               mapping = aes(x = 1,
                             y = age_g,
                             label = n),
               vjust = -0.6,
               angle = -90,
               family = 'Helvetica')+
     geom_text(data = datafile_label,
               mapping = aes(x = x,
                             y = y,
                             label = label),
               vjust = -2,
               angle = -90,
               family = 'Helvetica')+
     geom_text(data = datafile_label,
               mapping = aes(x = 0.5,
                             y = Inf,
                             label = title),
               vjust = -0.6,
               family = 'Helvetica')+
     facet_wrap(vars(factor(lineage, levels = c('BA2', 'BA1', 'Delta'))), 
                ncol = 1,
                scales = 'free_x',
                labeller = as_labeller(c(
                     Delta = 'f',
                     BA1 = 'd',
                     BA2 = 'b'
                )))+
     scale_y_discrete(labels = c('0-18', '19-64', '65-'),
                      expand = c(0, 0.6),
                      breaks = c('c', 'a', 'o'))+
     scale_x_continuous(n.breaks = 6,
                        expand =  expansion(add = c(0, 0)))+
     scale_fill_manual(values = fill_color,
                       labels = c('Unfully Vaccinated & Unmixed', 
                                  'Unfully Vaccinated & Mixed',
                                  'Fully Vaccinated & Unmixed',
                                  'Fully Vaccinated & Mixed', 
                                  'Booster Dose & Unmixed', 
                                  'Booster Dose & Mixed'
                       ),
                       na.translate = F)+
     coord_cartesian(clip = "off")+
     theme_classic(base_family = 'Helvetica')+
     labs(y = "",
          x = 'Proportions',
          fill = 'Vaccination status')

# combined plot ------------------------------------------------------------

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

# vaccine manufacturer ----------------------------------------------------

datafile_cont <- rbind(datafile_cont_Delta,
                       datafile_cont_BA1,
                       mutate(datafile_cont_BA2, gender = NA))|> 
     select(age, vaccine, vaccine_type, vaccine_mix, lineage) |> 
     mutate(age_g = if_else(age <=18, 'c', 'a'),
            age_g = if_else(age >=65, 'o', age_g),
            age_g = factor(age_g, levels = c('c', 'a', 'o')),
            vaccine_mix = if_else(is.na(vaccine_mix), 
                                  'N', 
                                  vaccine_mix)) |> 
     mutate(vaccine = factor(vaccine, 
                             levels = 0:3,
                             labels = c('C', 'C', 'B', 'A')),
            vaccine_mix = if_else(vaccine_mix == 'Y', 'M', 'U'),
            vaccine_g = str_remove_all(vaccine_type, '[_]'),
            vaccine_g = sapply(vaccine_g, FUN = function(x){
                 paste(sort(unique(strsplit(x, "")[[1]])), collapse = '')
            })) |> 
     group_by(vaccine, vaccine_g, lineage) |> 
     count() |> 
     ungroup() |> 
     filter(vaccine_g != "")

datafile_cont_back <- datafile_cont |> 
     group_by(lineage, vaccine) |> 
     summarise(n = sum(n),
               .groups = 'drop')

# datafile_cont_top <- datafile_cont |>
#         arrange(desc(n)) |>
#         group_by(lineage, vaccine) |>
#         slice(1:3)
datafile_cont_top <- datafile_cont |> 
     group_by(lineage, vaccine) |> 
     filter(vaccine_g %in% c('P', 'V', 'C', 'PV')) |> 
     select(vaccine, vaccine_g, lineage, n)

datafile_cout_sum <- datafile_cont_top |> 
     summarise(n = sum(n),
               .groups = 'drop') |> 
     left_join(datafile_cont_back,
               by = c("lineage", "vaccine")) |> 
     mutate(n = n.y - n.x,
            vaccine_g = 'Other') |> 
     select(lineage, vaccine, vaccine_g, n) |> 
     rbind(datafile_cont_top) |> 
     mutate(vaccine_g = if_else(vaccine_g == "",
                                'Unvaccine',
                                vaccine_g),
            vaccine_g = factor(vaccine_g,
                               levels = c('P', 'V', 'C', 'PV', 'Other', 'Unvaccine')))

# plot --------------------------------------------------------------------

datafile_label <- data.frame(x = rep('B', 3),
                             label = c('Omicron BA.2 outbreak',
                                       'Omicron BA.1 outbreak',
                                       'Delta outbreak'),
                             lineage = c('BA2', 'BA1', 'Delta'))

fig_2 <- ggplot(data = datafile_cout_sum)+
     geom_col(mapping = aes(x = vaccine,
                            y = n, 
                            fill = vaccine_g),
              color = 'black',
              position = position_dodge2(width = 0.9, preserve = "single"),
              show.legend = T)+
     geom_text(data = filter(datafile_label),
               mapping = aes(x = x,
                             y = Inf,
                             label = label),
               vjust = -0.6,
               family = 'Helvetica')+
     coord_cartesian(clip = "off")+
     facet_wrap(vars(factor(lineage, levels = c('BA2', 'BA1', 'Delta'))),
                nrow = 1,
                scales = 'free_y',
                labeller = as_labeller(c(
                     Delta = 'i',
                     BA1 = 'h',
                     BA2 = 'g'
                )))+
     scale_x_discrete(labels = c('Unfully Vaccinated', 'Fully Vaccinated', 'Booster Dose'),
                      expand = c(0, 0.6),
                      breaks = c('C', 'B', 'A'))+
     scale_y_continuous(limits = c(0, 2000),
                        expand = c(0, 0))+
     scale_fill_manual(values = fill_color_1,
                       labels = c('Sinopharm',
                                  'Sinovac',
                                  'Cansino',
                                  'Sinopharm & Sinovac',
                                  'Other',
                                  'Unvaccined'),
                       na.translate = F)+
     coord_cartesian(clip = "off")+
     theme_classic(base_family = 'Helvetica')+
     theme(plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
           plot.margin = margin(0, 1, 0, 0.25, "cm"),
           axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.text.y = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.title.y = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
           strip.text = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
           strip.background = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.position = 'bottom')+
     labs(y = "Number of contacts",
          x = '',
          fill = 'Manufacturer')

# combined plot -----------------------------------------------------------

cowplot::plot_grid(fig_1, fig_2, 
                   ncol = 1,
                   rel_heights = c(1.75, 1))

ggsave(filename = './outcome/publish/Figure 3.pdf',
       width = 14, height = 11)
