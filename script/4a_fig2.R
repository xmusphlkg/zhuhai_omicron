
# packages ----------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggsci)
library(prismatic, include.only = 'best_contrast')
library(showtext)

font_add('Helvetica', 
         "./script/fonts/Helvetica.ttf",
         bold = "./script/fonts/Helvetica-Bold.ttf")

remove(list = ls())

fill_color <- pal_nejm()(8)[c(1, 3:4)]
fill_color_1 <- pal_nejm()(8)[c(2, 5:8)]

load('./data/sars_2_cov.Rdata')

# BA2 ---------------------------------------------------------------------

datafile_gannt <- datafile_info_BA2 |> 
     select(id, dateexpose1, dateexpose2, dateonset, datepositive, datereported) |> 
     mutate(index = sprintf('%02d', id))

datafile_gannt_plot <- datafile_gannt[rep(datafile_gannt$id, 4), c('id', 'index')] |> 
     mutate(start_date = as.Date(unlist(datafile_gannt[,2:5]), origin = as.Date("1970-01-01")),
            end_date = as.Date(unlist(datafile_gannt[,3:6]), origin = as.Date("1970-01-01")),
            type = rep(c("A", "B", "C", "D"), each = nrow(datafile_gannt)),
            just = if_else(start_date == end_date,
                            0.1,
                            0)) |> 
     arrange(index) |> 
     group_by(index) |> 
     mutate(index = fct_inorder(index))

y_labels <- paste('Infection', datafile_gannt$index)
names(y_labels) <- as.character(datafile_gannt$index)
limits <- c(as.Date('2022/3/7'), as.Date('2022/4/10'))

fig_a <- ggplot(data = datafile_gannt_plot)+
     geom_linerange(mapping = aes(y = index,
                                  xmin = start_date,
                                  xmax = end_date),
                    colour = '#C8C8C8FF',
                    size = 4,
                    show.legend = F)+
     geom_point(mapping = aes(y = index,
                              x = start_date-just,
                              colour = type),
                size = 6)+
     geom_point(data = filter(datafile_gannt_plot, type == 'D'),
                mapping = aes(y = index,
                              x = end_date,
                              colour = 'E'),
                size = 6)+
     annotate(geom = 'text',
              x = mean.Date(limits),
              y = Inf,
              vjust = -0.5,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Infections of Omicron BA.2')+
     coord_cartesian(clip = "off")+
     scale_colour_manual(values = fill_color_1,
                         labels = c("First Exposed",
                                    "Last Exposed",
                                    "Onset",
                                    "Positive testing",
                                    "Report"))+
     scale_y_discrete(labels = y_labels)+
     scale_x_date(expand = c(0,1),
                  limits = limits,
                  date_breaks = '2 days',
                  date_labels = '%Y-%m-%d')+
     labs(title = 'a',
          y = "",
          x = 'Date',
          colour = 'Date')

# BA1 ---------------------------------------------------------------------

datafile_gannt <- datafile_info_BA1 |> 
     select(id, dateexpose1, dateexpose2, dateonset, datepositive, datereported) |> 
     mutate(index = sprintf('%02d', id))

datafile_gannt_plot <- datafile_gannt[rep(datafile_gannt$id, 4), c('id', 'index')] |> 
     mutate(start_date = as.Date(unlist(datafile_gannt[,2:5]), origin = as.Date("1970-01-01")),
            end_date = as.Date(unlist(datafile_gannt[,3:6]), origin = as.Date("1970-01-01")),
            type = rep(c("A", "B", "C", "D"), each = nrow(datafile_gannt)),
            just = if_else(start_date == end_date,
                           0.1,
                           0)) |> 
     arrange(index) |> 
     group_by(index) |> 
     mutate(index = fct_inorder(index))

y_labels <- paste('Infection', datafile_gannt$index)
names(y_labels) <- as.character(datafile_gannt$index)
limits <- c(as.Date('2022/1/4'), as.Date('2022/2/7'))

fig_b <- ggplot(data = datafile_gannt_plot)+
     geom_linerange(mapping = aes(y = index,
                                  xmin = start_date,
                                  xmax = end_date),
                    colour = '#C8C8C8FF',
                    size = 4,
                    show.legend = F)+
     geom_point(mapping = aes(y = index,
                              x = start_date-just,
                              colour = type),
                size = 6)+
     geom_point(data = filter(datafile_gannt_plot, type == 'D'),
                mapping = aes(y = index,
                              x = end_date,
                              colour = 'E'),
                size = 6)+
     annotate(geom = 'text',
              x = mean.Date(limits),
              y = Inf,
              vjust = -0.5,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Infections of Omicron BA.1')+
     coord_cartesian(clip = "off")+
     scale_colour_manual(values = fill_color_1,
                         labels = c("First Exposed",
                                    "Last Exposed",
                                    "Onset",
                                    "Positive testing",
                                    "Report"))+
     scale_y_discrete(labels = y_labels)+
     scale_x_date(expand = c(0,1),
                  limits = limits,
                  date_breaks = '2 days',
                  date_labels = '%Y-%m-%d')+
     labs(title = 'b',
          y = "",
          x = 'Date',
          colour = 'Date')+
     guides(colour = 'none')
# delta ---------------------------------------------------------------------

datafile_gannt <- datafile_info_Delta |> 
     select(id, dateexpose1, dateexpose2, dateonset, datepositive, datereported) |> 
     mutate(index = sprintf('%03d', id))

datafile_gannt_plot <- datafile_gannt[rep(datafile_gannt$id, 4), c('id', 'index')] |> 
     mutate(start_date = as.Date(unlist(datafile_gannt[,2:5]), origin = as.Date("1970-01-01")),
            end_date = as.Date(unlist(datafile_gannt[,3:6]), origin = as.Date("1970-01-01")),
            type = rep(c("A", "B", "C", "D"), each = nrow(datafile_gannt)),
            just = if_else(start_date == end_date,
                           0.1,
                           0)) |> 
     arrange(index) |> 
     group_by(index) |> 
     mutate(index = fct_inorder(index))

y_labels <- paste('Infection', datafile_gannt$index)
names(y_labels) <- as.character(datafile_gannt$index)
limits <- c(as.Date('2021/7/14'), as.Date('2021/8/17'))

fig_c <- ggplot(data = datafile_gannt_plot)+
     geom_linerange(mapping = aes(y = index,
                                  xmin = start_date,
                                  xmax = end_date),
                    colour = '#C8C8C8FF',
                    size = 4,
                    show.legend = F)+
     geom_point(mapping = aes(y = index,
                              x = start_date-just,
                              colour = type),
                size = 6)+
     geom_point(data = filter(datafile_gannt_plot, type == 'D'),
                mapping = aes(y = index,
                              x = end_date,
                              colour = 'E'),
                size = 6)+
     annotate(geom = 'text',
              x = mean.Date(limits),
              y = Inf,
              vjust = -0.5,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Infections of Delta')+
     coord_cartesian(clip = "off")+
     scale_colour_manual(values = fill_color_1,
                         labels = c("First Exposed",
                                    "Last Exposed",
                                    "Onset",
                                    "Positive testing",
                                    "Report"))+
     scale_y_discrete(labels = y_labels)+
     scale_x_date(expand = c(0,1),
                  limits = limits,
                  date_breaks = '2 days',
                  date_labels = '%Y-%m-%d')+
     labs(title = 'c',
          y = "",
          x = 'Date',
          colour = 'Date')+
     guides(colour = 'none')

# combined plot -----------------------------------------------------------

design <- "
AACC
BBCC
"

fig_a + fig_b + fig_c +
     plot_layout(design = design)&
     theme_bw(base_family = 'Helvetica')+
     theme(plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
           plot.margin = margin(0, 0.5, 0, 0, "cm"),
           axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', angle = 45, color = 'black'),
           axis.text.y = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', color = 'black'),
           axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
           legend.position = 'bottom',
           strip.text = element_blank(),
           strip.placement = "outside",
           strip.background = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.spacing.y = unit(0, 'mm'),
           text = element_text(size = 12),
           axis.title.y = element_text(size = 12, face = 'bold'))

ggsave('./outcome/publish/Figure 2.pdf',
       height = 18,
       width = 14)
