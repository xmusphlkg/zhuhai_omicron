
# packages ----------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggsci)
library(ggpubr)

remove(list = ls())

load('./data/data_case.Rdata')

extrafont::font_import('./script/fonts', prompt = F)
extrafont::loadfonts(device = "win")

# vaccine effectiveness susceptibility ------------------------------------

## data clean --------------------------------------------------------------

## inner household summary
datafile_household_sum <- datafile_cont_all %>% 
  mutate(vaccine_g = factor(vaccine, levels = 0:3,
                            labels = c(0, 0, 1, 2))) %>% 
  filter(contact == '同住') %>%
  group_by(vaccine_g, outcome) %>% 
  count() %>% 
  group_by(vaccine_g) %>% 
  mutate(freq = round(n/sum(n), digits = 4),
         type = 'Household')

## outer household summary
datafile_outer_sum <- datafile_cont_all %>% 
  mutate(vaccine_g = factor(vaccine, levels = 0:3,
                            labels = c(0, 0, 1, 2))) %>% 
  group_by(vaccine_g, outcome) %>% 
  count() %>% 
  group_by(vaccine_g) %>% 
  mutate(freq = round(n/sum(n), digits = 4),
         type = 'Entire')

datafile_plot <- rbind(datafile_household_sum, datafile_outer_sum) %>% 
  mutate(vaccine_g = factor(vaccine_g, levels = 0:2,
                            labels = c('Unfully Vaccinated',
                                       'Fully Vaccinated',
                                       'Booster Dose'))) %>% 
  filter(outcome == 1)

plot_text <- c(as.numeric(((datafile_plot[1,'freq'] - datafile_plot[2,'freq'])/datafile_plot[1,'freq'])*100),
               as.numeric(((datafile_plot[1,'freq'] - datafile_plot[3,'freq'])/datafile_plot[1,'freq'])*100))
plot_text <- round(plot_text, 2)

## plot --------------------------------------------------------------------

fig1 <- ggplot(data = datafile_plot)+
  geom_point(mapping = aes(x = vaccine_g, y = freq, 
                           group = type, color = type), 
             shape = 17)+
  geom_line(mapping = aes(x = vaccine_g, y = freq, 
                          group = type, color = type),
            show.legend = F)+
  scale_color_nejm()+
  scale_y_continuous(limits = c(0, 0.3), 
                     expand = c(0, 0),
                     labels = function(x) paste0(x*100, "%"))+
  theme_bw(base_family = 'Helvetica')+
  theme(legend.title = element_text(face = 'bold', size = 12),
        legend.text = element_text(size = 10),
        legend.box.background = element_rect(fill = "transparent", colour = 'transparent'),
        legend.background = element_rect(fill = "transparent", colour = 'transparent'),
        legend.position = c(0.3,0.8))+
  labs(x = '',
       y = 'SAR',
       color = 'Classification of contact',
       title = 'A')


# vaccine effectiveness contagious ----------------------------------------

## data clean --------------------------------------------------------------

datafile_info <- datafile_info_clean %>% 
  select(id, vaccine) %>% 
  mutate(vaccine = factor(vaccine, levels = 0:3,
                            labels = c('Unfully Vaccinated',
                                       'Unfully Vaccinated',
                                       'Fully Vaccinated',
                                       'Booster Dose')))

## inner household each case

datafile_household <- datafile_cont_all %>% 
  select(id_cases, outcome, contact) %>% 
  filter(contact == '同住') %>%
  group_by(id_cases, outcome) %>% 
  count() %>% 
  ungroup() %>% 
  complete(id_cases,
           outcome,
           fill = list(
             n = 0,
             freq = 0
           )) %>% 
  group_by(id_cases) %>% 
  mutate(freq = round(n/sum(n), digits = 4),
         type = 'Household') %>% 
  left_join(datafile_info, by = c('id_cases' = 'id')) %>% 
  filter(outcome == 1)

## inner household summary
datafile_household_sum <- datafile_cont_all %>%
  select(id_cases, outcome, contact) %>%
  filter(contact == '同住') %>%
  left_join(datafile_info, by = c('id_cases' = 'id')) %>%
  group_by(vaccine, outcome) %>%
  count() %>%
  group_by(vaccine) %>%
  mutate(freq = round(n/sum(n), digits = 4),
         type = 'Household')

## outer household
datafile_outer_sum <- datafile_cont_all %>% 
  select(id_cases, outcome, contact) %>% 
  left_join(datafile_info, by = c('id_cases' = 'id')) %>% 
  group_by(vaccine, outcome) %>% 
  count() %>% 
  group_by(vaccine) %>% 
  mutate(freq = round(n/sum(n), digits = 4),
         type = 'Entire')

datafile_plot <- rbind(datafile_household_sum, datafile_outer_sum) %>% 
  filter(outcome == 1)

## plot --------------------------------------------------------------------

fig2 <- ggplot(data = datafile_plot)+
  geom_point(mapping = aes(x = vaccine, y = freq, 
                           group = type, color = type), 
             shape = 17)+
  geom_line(mapping = aes(x = vaccine, y = freq, 
                          group = type, color = type),
            show.legend = F)+
  scale_color_nejm()+
  scale_y_continuous(limits = c(0, 0.3), 
                     expand = c(0, 0),
                     labels = function(x) paste0(x*100, "%"))+
  theme_bw(base_family = 'Helvetica')+
  theme(legend.title = element_text(face = 'bold', size = 12),
        legend.text = element_text(size = 10),
        legend.box.background = element_rect(fill = "transparent", colour = 'transparent'),
        legend.background = element_rect(fill = "transparent", colour = 'transparent'),
        legend.position = c(0.3,0.8))+
  labs(x = '',
       y = 'SAR',
       color = 'Classification of contact',
       title = 'B')
fig1+fig2

# Ct value ----------------------------------------------------------------

swab_x <- c('Nasopharyngeal_O', 'Nasopharyngeal_N', 'Throat_O', 'Throat_N')

for (i in 3:6) {
  datafile_ct_value <- datafile_info_clean %>%
    select(id, vaccine, type, Nasopharyngeal_O, Nasopharyngeal_N, Throat_O, Throat_N,
           Anal_O, Anal_N) %>% 
    pivot_longer(cols = -c(id, type, vaccine), names_to = 'swab', values_to = 'value') %>% 
    filter(swab %in% swab_x[i-2]) %>% 
    filter(!is.na(value)) %>% 
    mutate(vaccine_g = factor(vaccine,
                              levels = 0:3,
                              labels = c('Unfully Vaccinated',
                                         'Unfully Vaccinated',
                                         'Fully Vaccinated',
                                         'Booster Dose')))
  
  fig_ct <- ggboxplot(data = datafile_ct_value,
                      x = 'vaccine_g',
                      y = 'value',
                      add = 'jitter',
                      add.params = list(size = 8*5/14),
                      color = 'vaccine_g',
                      short.panel.labs = T,
                      legend = 'none')+
    stat_compare_means(comparison = list(c("Fully Vaccinated", "Booster Dose")), 
                       label.y = 40,
                       family = 'Helvetica',
                       label.x.npc = "left",
                       bracket.size = 1,
                       size = 8*5/14)+
    stat_compare_means(label.y = 15,
                       label.x = 0.8,
                       family = 'Helvetica',
                       label.x.npc = "left",
                       size = 8*5/14)+
    geom_rect(mapping = aes(xmin = -Inf, xmax = Inf, ymax = Inf, ymin = 40),
              fill = 'grey', alpha = 0.01)+
    geom_text(mapping = aes(x = 0.7, y = 46), label = 'Negative', 
              color = 'black', 
              family = 'Helvetica', size = 10*5/14, 
              nudge_x = 0)+
    scale_y_continuous(limits = c(0, 49),
                       expand = c(0, 0))+
    scale_color_nejm()+
    theme_bw(base_family = 'Helvetica')+
    theme(legend.position = 'none')+
    labs(x = '',
         y = 'Ct value',
         title = paste0(LETTERS[i]))
  assign(paste0('fig', i), fig_ct)
}

(fig1 + fig2) /
  (fig3 + fig4) / 
  (fig5 + fig6)&
  theme(axis.text.x = element_text(size = 10, color = 'black'),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.x = element_text(face = 'bold', size = 12),
        axis.title.y = element_text(face = 'bold', size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
        plot.margin = margin(1, 1, 1, 1))

ggsave(filename = './outcome/science/Figure 3.pdf', 
       device = cairo_pdf, height = 10, width = 9)

ggsave(filename = './outcome/science/Figure 3.tiff',
       width = 9, height = 10, dpi = 300)

ggsave(filename = './outcome/science/Figure 3.png',
       width = 9, height = 10, dpi = 300)


