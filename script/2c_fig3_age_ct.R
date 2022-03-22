
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

# plot --------------------------------------------------------------------

datafile_ct_value <- datafile_info_clean %>%
     select(vaccine, age, Nasopharyngeal_O, Nasopharyngeal_N, Throat_O, Throat_N) %>% 
     pivot_longer(cols = -c(vaccine, age), names_to = 'swab', values_to = 'value') %>% 
     filter(!is.na(value)) %>% 
     mutate(age_g = as.character(ifelse(age < 18, 1, 2)))

fig <- ggboxplot(data = datafile_ct_value,
                 x = 'age_g',
                 y = 'value',
                 add = 'jitter',
                 color = 'age_g',
                 short.panel.labs = T,
                 ggtheme = theme_bw(),
                 legend = 'none') +
     stat_compare_means(comparison = list(c("1", "2")),
                        label.y = 40,
                        label.x = 0.8,
                        family = 'Helvetica',
                        label.x.npc = "left",
                        size = 10*5/14)+
     geom_rect(mapping = aes(xmin = -Inf, xmax = Inf, ymax = Inf, ymin = 40),
               fill = 'grey', alpha = 0.01)+
     geom_text(mapping = aes(x = 0.7, y = 46), label = 'Negative', 
               color = 'black', 
               family = 'Helvetica', size = 10*5/14, 
               nudge_x = 0)+     
     scale_x_discrete(labels = c('<18', '≥18'))+
     scale_y_continuous(limits = c(0, 49),
                        expand = c(0, 0))+
     scale_color_nejm()+
     theme_bw(base_family = 'Helvetica')+
     theme(axis.title.x = element_text(face = 'bold', size = 14),
           axis.title.y = element_text(face = 'bold', size = 14),
           axis.text.x = element_text(size = 10),
           axis.text.y = element_text(size = 10),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),
           strip.text = element_text(hjust = 0, face = 'bold', size = 12),
           plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold', family = 'Times New Roman'),
           plot.margin = margin(5, 5, 5, 5),
           legend.position = 'none')+
     labs(x = 'Age of infector',
          y = 'Ct value')

facet(fig, facet.by = 'swab', ncol = 4, panel.labs = list(swab = letters[1:4]))


ggsave(filename = './outcome/publish/extend/Figure 4.pdf',
       device = cairo_pdf, height = 3, width = 12)

ggsave(filename = './outcome/publish/extend/Figure 4.png',
       width = 12, height = 3, dpi = 300)
