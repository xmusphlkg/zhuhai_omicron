
# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(ggsci)
library(MASS)

extrafont::font_import('./script/fonts', prompt = F)
extrafont::loadfonts(device = "win")

remove(list = ls())

load('data/data_case.Rdata')

set.seed(202202)

data_expose <- datafile_info_clean %>% 
  dplyr::select(id, dateexpose1, dateexpose2)

expose_seq <- function(x){
  infector_expose_dates <- seq.Date(from = datafile_cont_2[x, 'infector_exposedate1'],
                                    to = datafile_cont_2[x, 'infector_exposedate2'],
                                    by = 'day')
  infectee_expose_dates <- seq.Date(from = datafile_cont_2[x, 'infectee_exposedate1'],
                                    to = datafile_cont_2[x, 'infectee_exposedate2'],
                                    by = 'day')
  datafile <- as.data.frame(expand.grid(infectee_expose_dates, infector_expose_dates)) %>% 
    mutate(median = as.numeric(Var2 - Var1),
           median = ifelse(median<0, 0, median))
  gt_median <- median(datafile$median)
  gt_max <- max(infectee_expose_dates) - min(infector_expose_dates)
  gt_min <- min(infectee_expose_dates) - max(infector_expose_dates)
  gt_min <- ifelse(gt_min<0, 0, gt_min)
  return(c(gt_min, gt_median, gt_max))
}

gamma_min_add <- function(x){
  gt_min <- fitdistr(as.numeric(datafile_cont_2$gt_min + x), "gamma")
  return(c(x, shape = gt_min$estimate[1], rate = gt_min$estimate[2]))
}

datafile_cont_2 <- datafile_cont %>% 
  dplyr::select(to, from) %>% 
  rename(c(
    'infectee' = 'to',
    'infector' = 'from'
  )) %>% 
  left_join(data_expose, by = c('infectee' = 'id')) %>% 
  rename(c('infectee_exposedate1' = dateexpose1,
           'infectee_exposedate2' = dateexpose2))%>% 
  left_join(data_expose, by = c('infector' = 'id')) %>% 
  rename(c('infector_exposedate1' = dateexpose1,
           'infector_exposedate2' = dateexpose2))

datafile_expose_seq <- t(sapply(rownames(datafile_cont_2), expose_seq))
datafile_cont_2 <- datafile_cont_2 %>% 
  mutate(gt_min = datafile_expose_seq[,1],
         gt_max = datafile_expose_seq[,3])

add_time <- seq(0.001, 0.020, 0.001)

df_gt_min <- lapply(add_time, gamma_min_add)
df_gt_min <- data.frame(do.call('rbind', df_gt_min))
names(df_gt_min) <- c('add', 'shape', 'rate')

# plot --------------------------------------------------------------------

gamma_lines <- plyr::alply(as.matrix(df_gt_min), 1, function(gt_select){
  stat_function(fun = dgamma, 
                n = 100,
                args = list(shape = gt_select[2], rate = gt_select[3]),
                mapping = aes(colour = '0.001-0.020')
                # mapping = aes(colour = as.character(gt_select[1]))
                )
})


fig <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
  gamma_lines+
  stat_function(fun = dgamma, n = 100,
                args = list(shape = df_gt_min$shape[3], rate = df_gt_min$rate[3]),
                mapping = aes(colour = '0.002')) +
  coord_cartesian(xlim = c(0, 5))+
  scale_colour_manual(values = c('grey', 'red'))+
  scale_x_continuous(breaks = seq(0, 5, 1),
                     expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     labels = label_number(accuracy = 0.01))+
  theme_bw(base_family = 'Helvetica', base_line_size = 0.2)+
  theme(axis.title.x = element_text(face = 'bold', size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 5, color = 'black'),
        axis.text.y = element_text(size = 5, color = 'black'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(face = 'bold', size = 6),
        plot.title = element_text(size = 9, hjust = 0, vjust = .5, face = 'bold'),
        plot.margin = margin(1, 5, 1, 1),
        legend.spacing.y = unit(0, 'pt'),
        legend.text = element_text(margin = margin(r = 10, unit = "pt"),
                                   size = 5),
        legend.background = element_rect(fill = "transparent", colour = 'transparent'),
        legend.position = c(0.8,0.7))+
  labs(x = 'Time (days)',
       y = 'Relative frequency',
       colour = 'Addition (day)')

fig

ggsave(filename = './outcome/science/Figure S1.pdf', 
       device = cairo_pdf, height = 5.5, width = 5.5, units = 'cm')

ggsave(filename = './outcome/science/Figure S1.tiff', 
       width = 5.5, height = 5.5, dpi = 300, units = 'cm')

ggsave(filename = './outcome/science/Figure S1.png', 
       width = 5.5, height = 5.5, dpi = 300, units = 'cm')
