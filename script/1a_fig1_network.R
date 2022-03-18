
# packages ----------------------------------------------------------------

library(tidyverse)
library(igraph)
library(scales)
library(RColorBrewer)
library(patchwork)
library(ggforce)
library(ggraph)
library(graphlayouts)
library(sf)
library(ggrepel)
library(ghibli)

extrafont::font_import('./script/fonts', prompt = F)
extrafont::loadfonts(device = "win")

remove(list = ls())

border_color <- brewer_pal(palette = "PuRd")(3)
fill_color <- ggsci::pal_nejm()(3)
vaccine_color <- c('#A6611A', '#018571', '#80CDC1', '#EDF8E9')
fill_color_gannt <- ghibli::ghibli_palette('PonyoLight', 4)

load('./data/data_case.Rdata')

datafile_info <- datafile_info_clean %>% 
  select(id, type, vaccine, locationcluster) %>% 
  mutate(type = factor(type, labels = LETTERS[1:3]),
         width = factor(vaccine, levels = c(3, 2, 1, 0),
                        labels = LETTERS[1:4]),
         label = id)

datafile_cont <- datafile_cont[,c('from', 'to', 'freq')]


# Epicurve ----------------------------------------------------------------

set.seed(202202)

datafile_label <- data.frame(
  date = as.Date(c('2022-1-13', '2022-1-14', '2022-1-15', '2022-1-16', '2022-1-18', '2022-1-19')),
  value = c(4, 4, 7, 3, 4, 4),
  label = c('A case reported in neighbouring city',
            'First local case reported\nFirst citywide mass COVID testing',
            'Entirely locked Nanping town',
            'Second citywide mass COVID testing',
            'Third citywide mass COVID testing',
            'Partly locked of Xiangwan subdistrict')
)

datafile_info <- datafile_info_clean %>% 
  mutate(type = factor(type, labels = LETTERS[1:3]),
         width = factor(vaccine, levels = c(3, 2, 1, 0),
                        labels = LETTERS[1:4]),
         label = id)

fig_curve <- ggplot(data = datafile_info)+
  geom_col(mapping = aes(x = dateonset, y = 1, fill = type),
           color = "white",
           width = 1,
           size = 0.1)+
  geom_text_repel(mapping = aes(x = date, y = value, label = label),
                  data = datafile_label,
                  family = 'Helvetica',
                  size = 8*5/14,
                  force_pull   = 0,
                  nudge_y      = 0.8,
                  direction    = "x",
                  angle        = 90,
                  hjust        = 0,
                  segment.size = 0.2,
                  max.iter = 1e4, max.time = 1)+
  coord_equal()+
  scale_fill_manual(labels = c('Moderate', 'Mild', 'Asymptomatic'),
                    values = fill_color)+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 12),
                     labels = scales::number_format(accuracy = 1)
  )+
  scale_x_date(expand = c(0,0),
               limits = c(as.Date('2022/1/2'), as.Date('2022/1/22')),
               breaks = '2 days',
               date_labels = '%Y-%m-%d')+
  theme_bw(base_family = 'Helvetica')+
  theme(
    axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', 
                               angle = 45, 
                               color = 'black'),
    axis.text.y = element_text(size = 10, hjust = .5, vjust = .5, face = 'plain', 
                               color = 'black'),
    axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
    axis.title.y = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
    plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin=grid::unit(c(0,0,0,0), "cm"),
    strip.text = element_text(size = 12, hjust = .5, vjust = 0, face = 'bold')
  )+
  labs(x = 'Onset Date',
       y = 'Number of infections',
       title = 'a',
       fill = 'Infections\nclassification')

fig_curve


# Network plot ------------------------------------------------------------

fig_net <- graph_from_data_frame(d = datafile_cont, 
                                 vertices = datafile_info, 
                                 directed = T)

fig_b <- ggraph(fig_net,layout = "kk")+
  geom_edge_link(mapping = aes(colour = freq), 
                 arrow = arrow(length = unit(3, 'mm')), 
                 end_cap = circle(4, 'mm'),
                 check_overlap = F,
                 width = 0.7,
                 show.legend = T)+
  geom_node_point(aes(fill = type,
                      colour = width),
                  size = 7,
                  stroke = 2,
                  shape = 21,
                  show.legend = T) +
  scale_fill_manual(name = "Infections\nclassification",
                    labels = c('Moderate', 'Mild', 'Asymptomatic'),
                    values = fill_color)+
  scale_edge_colour_manual(name = 'Contact frequency',
                           labels = c('Low', 'Mild', 'High'),
                           values = border_color)+
  scale_colour_manual(name = 'Vaccination status',
                      labels = c('Booster Dose', 'Fully Vaccinated', 'One Dose', 'Non-vaccinated'),
                      values = vaccine_color)+
  geom_node_text(aes(label = label), 
                 family = 'Helvetica') +
  coord_equal() +
  theme_graph(base_family = 'Helvetica')+
  theme(plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        legend.key.height = unit(1.5,"line"))+
  guides(fill = guide_legend(override.aes = list(color = 'white', size = 9,
                                                 stroke = 0, 
                                                 edge_color = 'transparent')),
         edge_colour = guide_legend(override.aes = list(color = 'transparent')),
         colour = guide_legend(override.aes = list(fill = '#FFDC91FF', size = 7, 
                                                   stroke = 2, edge_color = 'transparent'))
         )+
  labs(title = 'b')
fig_b


# contact tracing ---------------------------------------------------------

datafile_contact_uninfect <- datafile_cont_all %>% 
  filter(outcome == 0) %>% 
  select(id_cases, ids) %>% 
  rename(c('from' = 'id_cases',
           'to' = 'ids')) %>% 
  mutate(to = paste0('un-', to),
         type = 'Uninfect')

datafile_contact_infect <- datafile_cont %>% 
  select(from, to) %>% 
  mutate(type = 'Infect')

datafile_contact_all <- rbind(datafile_contact_infect, datafile_contact_uninfect)

datafile_info_uninfect <- datafile_cont_all %>% 
  filter(outcome == 0) %>% 
  select(ids, vaccine, gender) %>% 
  mutate(type = 'Uninfect') %>% 
  mutate(ids = paste0('un-', ids)) %>% 
  rename(c('id' = 'ids')) %>% 
  distinct(id, .keep_all= TRUE)

datafile_info_infect <- datafile_info %>% 
  select(id, vaccine, gender, type)

datafile_info_all <- rbind(datafile_info_infect, datafile_info_uninfect) %>% 
  mutate(vaccine = factor(vaccine))

fig_net <- graph_from_data_frame(d = datafile_contact_all, 
                                 vertices = datafile_info_all, 
                                 directed = T)

fig_d <- ggraph(fig_net, "kk")+
  geom_edge_link(mapping = aes(colour = type,
                               width = type,
                               alpha = type), 
                 check_overlap = T,
                 show.legend = F)+
  geom_node_point(aes(colour = vaccine),
                  stroke = 0,
                  size = 0.15,
                  shape = 19,
                  show.legend = T)+
  # scale_fill_manual(name = "Cases classification",
  #                   labels = c('Moderate', 'Mild', 'Asymptomatic', 'Uninfected'),
  #                   values = c(fill_color, 'grey'))+
  scale_colour_manual(name = 'Vaccination status',
                    labels = c('Booster Dose', 'Fully Vaccinated', 'One Dose', 'Non-vaccinated'),
                      values = vaccine_color)+
  scale_edge_alpha_manual(values = c(1, 0.1))+
  scale_edge_colour_manual(values = c('black', 'gray'))+
  scale_edge_width_manual(values = c(0.1, 0.001))+
  # scale_edge_width_manual(values = c(1, 1, 1, 0.001))+
  coord_equal() +
  theme_graph(base_family = 'Helvetica')+
  theme(plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))+
  guides(colour = guide_legend(override.aes = list(size = 0.5, shape = 20)),
         edge_colour = 'none',
         edge_width = 'none',
         edge_alpha = 'none',
         size = 'none')+
  labs(title = 'd')

# gannt -------------------------------------------------------------------

datafile_gannt <- datafile_info_clean %>% 
  select(id, dateexpose1, dateexpose2, dateonset, datepositive, datereported, locationcluster) %>% 
  mutate(cluster = factor(locationcluster,
                          exclude = NULL,
                          levels = c("Kindergarten", "Primary school",
                                     "Dental clinic", "Undistinct"),
                          labels = c('Kindergarten', 'Primary\nschool', 
                                     'Dental\nclinic', 'Indistinct')
                          )) %>% 
  arrange(cluster) %>% 
  group_by(cluster) %>% 
  mutate(index = paste(cluster, row_number(), sep = '.'))

datafile_gannt_plot <- datafile_gannt[rep(datafile_gannt$id, 4), c('id', 'index', 'cluster')]
datafile_gannt_plot$wp <- datafile_gannt_plot$cluster
datafile_gannt_plot$start_date <- as.Date(unlist(datafile_gannt[,2:5]), origin = as.Date("1970-01-01"))
datafile_gannt_plot$end_date <- as.Date(unlist(datafile_gannt[,3:6]), origin = as.Date("1970-01-01"))
datafile_gannt_plot$activity <- datafile_gannt_plot$id
datafile_gannt_plot$type <- rep(c("A", "B", "C", "D"), each = nrow(datafile_gannt))

y_labels <- paste('Case', sprintf('%02d', datafile_gannt$id))
names(y_labels) <- as.character(datafile_gannt$index)

fig_c <- ggplot(data = datafile_gannt_plot)+
  geom_linerange(mapping = aes(y = index,
                               xmin = start_date,
                               xmax = end_date,
                               colour = type),
                 size = 4)+
  scale_y_discrete(labels = y_labels)+
  scale_x_date(expand = c(0,0),
               limits = c(as.Date('2022/1/2'), as.Date('2022/1/22')),
               date_breaks = '2 days',
               date_labels = '%Y-%m-%d')+
  scale_colour_manual(values = fill_color_gannt,
                      labels = c("First Exposed\nto Last Exposed",
                                 "Last Exposed\nto Onset",
                                 "Onset\nto positive test",
                                 "Positive test\nto Report"))+
  theme_bw(base_family = 'Helvetica')+
  theme(plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', angle = 45, color = 'black'),
        axis.text.y = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', color = 'black'),
        axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
        strip.text.x = element_text(size = 12, face = 'bold'),
        strip.text.y = element_text(size = 12, face = 'bold'),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key.height = unit(2,"line"),
        panel.spacing.y = unit(0, 'mm'),
        text = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = 'bold'))+
  labs(title = 'c',
       y = "",
       x = 'Date',
       colour = 'Periods')+
  facet_grid(cluster~., scales = "free_y", switch = "y", space = "free_y")

fig_c

fig <- fig_curve + fig_b + fig_c + fig_d +
  plot_layout(ncol = 2, byrow = T)&
  theme(legend.text = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', 
                                   color = 'black'),
        legend.title = element_text(size = 12, hjust = 0, vjust = 0.5, face = 'bold'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,-7),
        legend.position = 'right',
        legend.justification = 'left')

ggsave(filename = './outcome/publish/Figure 1.pdf',
       fig,
       width = 20, height = 13, device = cairo_pdf)

ggsave(filename = './outcome/publish/Figure 1.tiff',
       fig,
       width = 18, height = 16, dpi = 300)
