
# packages ----------------------------------------------------------------

library(tidyverse)
library(sf)
library(ggrepel)
library(cowplot)
library(patchwork)
library(ggraph)
library(igraph, include.only = 'graph_from_data_frame')
library(ggsci)
library(prismatic, include.only = 'best_contrast')
library(showtext)

font_add('Helvetica', 
         "./script/fonts/Helvetica.ttf",
         bold = "./script/fonts/Helvetica-Bold.ttf")

remove(list = ls())

# border_color <- rgb_material('deep-orange', 4)[-1]
fill_color <- pal_nejm()(4)
type_levels <- c("Server", "Moderate", "Mild", "Asymptomatic")
text_color <- prismatic::best_contrast(fill_color)
vaccine_color <- c('#F1EB58FF', '#AA812DFF', '#283F99FF')

load('./data/sars_2_cov.Rdata')


# Map ---------------------------------------------------------------------

map_data_border <- st_read("./data/geo/border.shp")[,c('NAME', 'geometry')]
map_data_province <- st_read("./data/geo/province.shp")[,c('sheng_mc', 'geometry')] |> 
        st_make_valid()
map_data_city <- st_read("./data/geo/city.shp")[,c('sheng_mc', 'shi_mc', 'geometry')] |> 
        st_make_valid() |> 
        filter(sheng_mc %in% c('湖南省', '广东省', '福建省'))

datafile_map <- data.frame(
        lineage = c(rep('Delta', 7), 'BA1', 'BA2'),
        city = c('常德市', '湘潭市', '湘西土家族苗族自治州', 
                 '益阳市', '张家界市', '长沙市', '株洲市',
                 '珠海市', '厦门市'),
        n = c(3, 4, 2, 8, 76, 5, 31, 38, 35))

datafile_map <- map_data_city %>% 
        left_join(datafile_map, by = c('shi_mc' = 'city'))

datafile_point <- data.frame(
        lineage = c('Delta', 'Omicron BA.1', 'Omicron BA.2'),
        label_lng = rep(124, 3),
        label_lat = c(30, 25, 20),
        point_lng = c(110.5, 113.5, 118.2),
        point_lat = c(29.5, 22.2, 24.8)) |> 
        mutate(x_nudge = label_lng - point_lng,
               y_nudge = label_lat - point_lat)

fig_map <- ggplot(data = datafile_map)+
        geom_sf(mapping = aes(fill = n),
                color = 'black',
                size = 0.3,
                alpha = 1, 
                show.legend = T)+
        geom_sf(data = map_data_province,
                color = 'black',
                fill = NA,
                size = 0.5,
                show.legend = F)+
        geom_sf(data = map_data_border,
                color = 'black',
                fill = NA,
                size = 0.5,
                show.legend = F)+
        geom_text_repel(data = datafile_point,
                        mapping = aes(x = point_lng,
                                      y = point_lat,
                                      label = lineage),
                        nudge_x = datafile_point$x_nudge,
                        nudge_y = datafile_point$y_nudge,
                        hjust = 0,
                        color = fill_color[1],
                        family = 'Helvetica')+
        scale_fill_distiller(palette = "Oranges", 
                             na.value = 'white', 
                             limit = c(0, 80),
                             direction = 1)+
        coord_sf(xlim = c(108, 135), ylim = c(5, 32))+
        theme_bw(base_family = 'Helvetica')+
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
              plot.margin=grid::unit(c(0,0.2,0,0), "cm"),
              legend.position = c(0.99, 0.01),
              legend.justification = c(1, 0)
        )+
        guides(fill = guide_colorbar(frame.colour = "black", 
                                     ticks.colour = "black"))+
        labs(fill = 'Cumulative\nCases',
             title = 'a')

fig_map

# Epicurve ----------------------------------------------------------------

set.seed(202202)

theme_curve <- function(){
        theme_bw(base_family = 'Helvetica')+
                theme(legend.position = c(0.01, 0.99),
                      legend.justification = c(0, 1),
                      axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', 
                                                 color = 'black'),
                      axis.text.y = element_text(size = 10, hjust = .5, vjust = .5, face = 'plain', 
                                                 color = 'black'),
                      axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
                      axis.title.y = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
                      plot.title = element_text(size = 16, hjust = 0, vjust = 0, face = 'bold'),
                      plot.margin = margin(0, 0.2, 0, 0, "cm"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.key.size = unit(0.28, 'cm'),
                      strip.text = element_text(size = 12, hjust = .5, vjust = 0, face = 'bold'))
}

limits <- c(as.Date('2022/3/9'), as.Date('2022/4/7'))
fig_curve_ba2 <- datafile_info_BA2 |>  
     mutate(type = factor(type, 
                          levels = type_levels,
                          labels = LETTERS[1:4])
     ) |> 
     ggplot()+
     geom_col(mapping = aes(x = dateonset, y = 1, fill = type),
              color = "white",
              width = 1,
              size = 0.1,
              show.legend = T)+
     annotate(geom = 'text',
              x = mean.Date(limits),
              y = 16,
              vjust = -0.5,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Epicurve of Omicron BA.2')+
     coord_equal(clip = 'off')+
     scale_fill_manual(labels = type_levels[-1],
                       values = fill_color[-1])+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 16),
                        labels = scales::number_format(accuracy = 1)
     )+
     scale_x_date(expand = c(0,0),
                  limits = limits,
                  breaks = '4 days',
                  date_labels = '%m-%d\n%Y')+
     theme_curve()+
     labs(x = 'Onset Date',
          y = 'Number of infections',
          title = 'b',
          fill = 'Infections\nclassification')

limits <- c(as.Date('2022/1/1'), as.Date('2022/1/30'))
fig_curve_ba1 <- datafile_info_BA1 |>  
     mutate(type = factor(type, 
                          levels = type_levels,
                          labels = LETTERS[1:4])
            ) |> 
     ggplot()+
     geom_col(mapping = aes(x = dateonset, y = 1, fill = type),
              color = "white",
              width = 1,
              size = 0.1,
              show.legend = T)+
     annotate(geom = 'text',
              x = mean.Date(limits),
              y = 16,
              vjust = -0.5,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Epicurve of Omicron BA.1')+
     coord_equal(clip = 'off')+
     scale_fill_manual(labels = type_levels[-1],
                       values = fill_color[-1])+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 16),
                        labels = scales::number_format(accuracy = 1)
     )+
     scale_x_date(expand = c(0,0),
                  limits = limits,
                  breaks = seq.Date(as.Date('2022/1/1'), 
                                    as.Date('2022/1/30'), 
                                    '4 days'),
                  date_labels = '%m-%d\n%Y')+
     theme_curve()+
     labs(x = 'Onset Date',
          y = 'Number of infections',
          title = 'c',
          fill = 'Infections\nclassification')

limits <- c(as.Date('2021/7/20'), as.Date('2021/8/18'))
fig_curve_delta <- datafile_info_Delta |>  
     mutate(type = factor(type, 
                          levels = type_levels,
                          labels = LETTERS[1:4])
     ) |> 
     ggplot()+
     geom_col(mapping = aes(x = dateonset, y = 1, fill = type),
              color = "white",
              width = 1,
              size = 0.1,
              show.legend = T)+
     annotate(geom = 'text',
              x = mean.Date(limits),
              y = 16,
              vjust = -0.5,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Epicurve of Delta')+
     coord_equal(clip = 'off')+
     scale_fill_manual(values = fill_color,
                       labels = type_levels)+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 16),
                        labels = scales::number_format(accuracy = 1)
     )+
     scale_x_date(expand = c(0,0),
                  limits = limits,
                  breaks = seq.Date(as.Date('2021/7/20'), 
                                    as.Date('2021/8/18'),
                                    '4 days'),
                  date_labels = '%m-%d\n%Y')+
     theme_curve()+
     labs(x = 'Onset Date',
          y = 'Number of infections',
          title = 'd',
          fill = 'Infections\nclassification')

fig_1 <- fig_map + fig_curve_ba2 + fig_curve_ba1 + fig_curve_delta+
     plot_layout(ncol = 1)

# fig_1
# 
# ggsave(filename = './outcome/publish/test.pdf',
#        width = 4, height = 12)

# network plot ------------------------------------------------------------

fig <- graph_from_data_frame(d = datafile_chains_BA2[,c('from', 'to')],
                             vertices = mutate(datafile_info_BA2,
                                               label = id,
                                               type = factor(type, 
                                                             levels = type_levels,
                                                             labels = LETTERS[1:4]),
                                               width = factor(vaccine, levels = c(3, 2, 1, 0),
                                                              labels = c('A', 'B', 'C', 'C'))),
                             directed = T)

fig_net_ba2 <- ggraph(fig,layout = "kk")+
     geom_edge_link(arrow = arrow(length = unit(1.5, 'mm')), 
                    end_cap = circle(4, 'mm'),
                    check_overlap = F,
                    width = 0.7,
                    show.legend = F)+
     geom_node_circle(aes(fill = type,
                          linetype = width,
                          r = 0.35),
                      size = 1,
                      show.legend = F)+
     geom_node_text(aes(label = label),
                    colour = 'white',
                    show.legend = F,
                    family = 'Helvetica') +
     coord_equal(clip = 'off')+
     scale_fill_manual(name = "Infections\nclassification",
                       labels = type_levels[-1],
                       values = fill_color[-1])+
     scale_linetype_manual(name = 'Vaccination status',
                           labels = c('Booster Dose', 'Fully Vaccinated', 'Unfully Vaccinated'),
                           values = c('solid', 'longdash', 'dotted'))+
     theme_graph(base_family = 'Helvetica')+
     labs(title = 'e')

fig_net_ba2 <- fig_net_ba2+
     annotate(geom = 'text',
              x = (max(fig_net_ba2$data$x) + min(fig_net_ba2$data$x))/2,
              y = max(fig_net_ba2$data$y),
              vjust = -3.5,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Transmission chain of Omicron BA.2')

fig <- graph_from_data_frame(d = datafile_chains_BA1[,c('from', 'to')],
                             vertices = mutate(datafile_info_BA1,
                                               label = id,
                                               type = factor(type, 
                                                             levels = type_levels,
                                                             labels = LETTERS[1:4]),
                                               width = factor(vaccine, levels = c(3, 2, 1, 0),
                                                              labels = c('A', 'B', 'C', 'C'))),
                             directed = T)

fig_net_ba1 <- ggraph(fig,layout = "kk")+
     geom_edge_link(arrow = arrow(length = unit(1.5, 'mm')), 
                    end_cap = circle(4, 'mm'),
                    check_overlap = F,
                    width = 0.7,
                    show.legend = F)+
     geom_node_circle(aes(fill = type,
                          linetype = width,
                          r = 0.35),
                      size = 1,
                      show.legend = F)+
     geom_node_text(aes(label = label),
                    colour = 'white',
                    show.legend = F,
                    family = 'Helvetica') +
     coord_equal(clip = 'off')+
     scale_fill_manual(name = "Infections\nclassification",
                       labels = type_levels[-1],
                       values = fill_color[-1])+
     scale_linetype_manual(name = 'Vaccination status',
                           labels = c('Booster Dose', 'Fully Vaccinated', 'Unfully Vaccinated'),
                           values = c('solid', 'longdash', 'dotted'))+
     theme_graph(base_family = 'Helvetica')+
     labs(title = 'f')

fig_net_ba1 <- fig_net_ba1 +
     annotate(geom = 'text',
              x = (max(fig_net_ba1$data$x) + min(fig_net_ba1$data$x))/2,
              y = max(fig_net_ba1$data$y),
              vjust = -3.5,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Transmission chain of Omicron BA.1')

fig <- graph_from_data_frame(d = datafile_chains_Delta[,c('from', 'to')],
                             vertices = mutate(datafile_info_Delta,
                                               label = id,
                                               type = factor(type, 
                                                             levels = type_levels,
                                                             labels = LETTERS[1:4]),
                                               width = factor(vaccine, levels = c(2, 1, 0),
                                                              labels = c('B', 'C', 'C'))),
                             directed = T)

fig_net_delta <- ggraph(fig,layout = "kk")+
     geom_edge_link(arrow = arrow(length = unit(1.5, 'mm')), 
                    end_cap = circle(4, 'mm'),
                    check_overlap = F,
                    width = 0.7,
                    show.legend = F)+
     geom_node_circle(aes(fill = type,
                          linetype = width,
                          r = 0.35),
                      size = 1,
                      show.legend = F)+
     geom_node_text(aes(label = label),
                    colour = 'white',
                    show.legend = F,
                    family = 'Helvetica') +
     scale_fill_manual(name = "Infections classification",
                       labels = type_levels,
                       values = fill_color)+
     scale_linetype_manual(name = 'Vaccination status',
                           labels = c('Fully Vaccinated', 'Unfully Vaccinated'),
                           values = c('longdash', 'dotted'))+
     theme_graph(base_family = 'Helvetica')+
     coord_cartesian(clip = "off")+
     labs(title = 'g')

fig_net_delta <- fig_net_delta +
     annotate(geom = 'text',
              x = (max(fig_net_delta$data$x) + min(fig_net_delta$data$x))/2,
              y = max(fig_net_delta$data$y),
              vjust = -5.5,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Transmission chain of Delta')

# combined plot ------------------------------------------------------------

layout <- "
AEF
BGG
CGG
DGG
"

fig_2 <- fig_map + fig_curve_ba2 + fig_curve_ba1 + fig_curve_delta+
     fig_net_ba2 + fig_net_ba1 + fig_net_delta + 
     plot_layout(design = layout)&
     theme(
          # legend.position = 'bottom',
          # plot.title = element_text(size = 16, hjust = 0, vjust = 0, face = 'bold'),
          plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 12, hjust = .5, vjust = 0, face = 'bold')
     )

ggsave(filename = './outcome/publish/Figure 1r.pdf',
       fig_2,
       width = 14, height = 15)
