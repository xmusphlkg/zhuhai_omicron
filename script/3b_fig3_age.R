
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

# data clean --------------------------------------------------------------

datafile_info <- datafile_info_clean %>% 
  select(id, vaccine, age) %>% 
  mutate(vaccine = factor(vaccine, levels = 0:3,
                          labels = c('Unfull Vaccine',
                                     'Unfull Vaccine',
                                     'Full Vaccine',
                                     'Booster Dose')))

datafile_sar <- datafile_cont_all %>% 
  select(id_cases, outcome, contact) %>% 
  group_by(id_cases, outcome) %>% 
  count() %>% 
  ungroup() %>% 
  complete(id_cases = 1:38,
           outcome,
           fill = list(
             n = 0
           )) %>% 
  group_by(id_cases) %>% 
  mutate(freq = round(n/sum(n), digits = 4)) %>% 
  left_join(datafile_info, by = c('id_cases' = 'id')) %>% 
  filter(outcome == 1) %>% 
  mutate(age_g = ifelse(age < 18, 1, 2))

datafile_sar$freq[is.nan(datafile_sar$freq)] <- 0


# age sar plot --------------------------------------------------------------------

ggboxplot(data = datafile_sar,
          x = 'age_g',
          y = 'freq',
          add = 'jitter',
          color = 'age_g',
          short.panel.labs = T,
          legend = 'none')+
  stat_compare_means(comparison = list(c("1", "2")),
                     label.y = 0.8,
                     label.x = 0.8,
                     family = 'Helvetica',
                     label.x.npc = "left",
                     size = 10*5/14)+
  scale_x_discrete(labels = c('<18', '≥18'))+
  scale_color_nejm()+
  theme_bw(base_family = 'Helvetica')+
  theme(axis.title.x = element_text(face = 'bold', size = 14),
        axis.title.y = element_text(face = 'bold', size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold', family = 'Times New Roman'),
        plot.margin = margin(5, 5, 5, 5),
        legend.position = 'none')+
  labs(x = 'Age of infector',
       y = 'SAR')

ggsave(filename = './outcome/publish/extend/Figure 3.pdf',
       device = cairo_pdf, height = 5.5, width = 5.5)

ggsave(filename = './outcome/publish/extend/Figure 3.tiff',
       width = 5.5, height = 5.5, dpi = 300)

ggsave(filename = './outcome/publish/extend/Figure 3.png',
       width = 5.5, height = 5.5, dpi = 300)
