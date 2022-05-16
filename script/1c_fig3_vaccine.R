
# packages ----------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggsci)
library(ggpubr)
library(lme4)

remove(list = ls())

load('./data/data_case.Rdata')

library(showtext)
font_add('Helvetica', 
         "./script/fonts/Helvetica.ttf",
         bold = "./script/fonts/Helvetica-Bold.ttf")

# vaccine effectiveness ---------------------------------------------------

datafile_log_BA1 <- datafile_cont_all |> 
  select(age, vaccine, vaccine_type, vaccine_mix, outcome, seq) |> 
  mutate(lineage = 'BA.1')

datafile_log_BA1[which(datafile_log_BA1$vaccine_type == 'P_V'),'vaccine_mix'] <- 'Y'

datafile_log_BA2 <- datafile_cont_BA2 |> 
  select(age, vaccine, vaccine_type, vaccine_mix, outcome, seq) |> 
  mutate(lineage = 'BA.2',
         vaccine = if_else(is.na(vaccine),0,vaccine))

datafile_log <- rbind(datafile_log_BA1, datafile_log_BA2) |> 
  mutate(seq_g = if_else(seq >=14 & seq <=120,
                         "inner",
                         "outer"),
         vaccine_g = str_remove_all(vaccine_type, '[_]'),
         vaccine_g = sapply(vaccine_g, FUN = function(x){
           paste(sort(unique(strsplit(x, "")[[1]])), collapse = '')
         })) |> 
  filter(!(is.na(seq_g)&vaccine>0)) |> 
  filter(!is.na(age)) |> 
  mutate(vaccine = factor(vaccine, levels = 0:3),
         vaccine_g = if_else(vaccine == 0,
                             '0',
                             vaccine_g),
         vaccine_g = if_else(vaccine_mix == 'Y' & vaccine != 0,
                             "M",
                             vaccine_g),
         vaccine_g = paste(vaccine, vaccine_g, sep = "_"))

datafile_count <- datafile_log |> 
  group_by(vaccine_g) |> 
  count() |> 
  filter(n>10)

datafile_log_BA1 <- datafile_log |> 
  filter(vaccine_g %in% datafile_count$vaccine_g &
           lineage == 'BA.1')

model <- glm(outcome ~ vaccine_g + age,
             family = "binomial",
             data = datafile_log_BA1)

confint(model)


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
  print(swab_x[i-2])
  test <- datafile_ct_value %>% 
    group_by(vaccine_g) %>% 
    summarise(value = mean(value))
  print(test)
  
  fig_ct <- ggboxplot(data = datafile_ct_value,
                      x = 'vaccine_g',
                      y = 'value',
                      add = 'jitter',
                      add.params = list(size = 1),
                      color = 'vaccine_g',
                      short.panel.labs = T,
                      legend = 'none')+
    stat_compare_means(comparison = list(c("Fully Vaccinated", "Booster Dose")), 
                       label.y = 40,
                       family = 'Helvetica',
                       label.x.npc = "left",
                       bracket.size = 1,
                       size = 8*5/14)+
    stat_compare_means(label.y = 10,
                       label.x = 0.8,
                       family = 'Helvetica',
                       label.x.npc = "left",
                       size = 8*5/14)+
    geom_rect(mapping = aes(xmin = -Inf, xmax = Inf, ymax = Inf, ymin = 40),
              fill = 'grey', alpha = 0.01)+
    geom_text(mapping = aes(x = 1, y = 46), label = 'Negative', 
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
         title = paste0(letters[i]))
  assign(paste0('fig', i), fig_ct)
}

(fig1 + fig2) /
  (fig3 + fig4)/
     (fig5 + fig6)&
  theme(axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', angle = 15, color = 'black'),
        axis.text.y = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', color = 'black'),
        axis.title.y = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 12, hjust = 0, vjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
        plot.margin = margin(1, 1, 1, 1))

ggsave(filename = './outcome/publish/Figure 3.pdf', 
       device = cairo_pdf, height = 9, width = 6)

ggsave(filename = './outcome/publish/Figure 3.tiff',
       width = 6, height = 9, dpi = 300)


