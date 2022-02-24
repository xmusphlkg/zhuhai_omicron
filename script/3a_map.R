
# packages ----------------------------------------------------------------

library(tidyverse)
library(scales)
library(RColorBrewer)
library(patchwork)
library(sf)


extrafont::font_import('./script/fonts', prompt = F)
extrafont::loadfonts(device = "win")

load('./data/data_case.Rdata')

datafile_info <- datafile_info_clean

# map plot ----------------------------------------------------------------

map_data_city1 <- st_read("./data/geo/zhuhai1.shp")[,c('Name', 'geometry')]
map_data_county <- st_read("./data/geo/zhuhai2.shp")[,c('Name', 'geometry')]
map_data_town <- st_read("./data/geo/zhuhai3.shp")[,c('Name', 'geometry')]

datafile_map <- data.frame(table(datafile_info$location2))
datafile_map <- map_data_town %>% 
  left_join(datafile_map, by = c('Name' = 'Var1'))
datafile_map$Freq[is.na(datafile_map$Freq)] <- 0


map_data_city <- map_data_city1 %>% 
  filter(Name != '珠海市')

datafile_text <- data.frame(
  long = c(113.12, 113.385, 113.992),
  lat = c(22.42, 22.5, 22.6),
  label = c('Jiangmen City', 'Zhongshan City', 'Shenzhen City')
)

datafile_text_1 <- data.frame(
  long = c(113.317, 113.471, 113.503, 113.65),
  lat = c(22.111, 22.188, 22.260, 22.305),
  label = c('Hongqi Town', 'Nanping Town', 'Qianshan Subdistrict', 'Xiangwan Subdistrict')
)

fig_min <- ggplot()+
  geom_sf(data = map_data_city,
          fill = 'grey', 
          size = 0.5, color = 'black',
          show.legend = F)+
  geom_sf(data = datafile_map, aes(fill = Freq), alpha = 1, show.legend = T)+
  geom_text(data = datafile_text, aes(x = long, y = lat, label = label),
            size = 3.5,
            fontface = "bold",
            family = 'Helvetica')+
  geom_text(data = datafile_text_1, aes(x = long, y = lat, label = label),
            size = 3.5,
            family = 'Helvetica')+
  scale_fill_distiller(palette = "Oranges", 
                       na.value = 'white', 
                       limit = c(0, 40),
                       direction = 1)+
  coord_sf(xlim = c(113.050, 114.450), ylim = c(21.255, 22.655))+
  theme_bw(base_family = 'Helvetica')+
  theme(
    legend.position = 'right',
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
    plot.margin=grid::unit(c(0,0,0,0), "cm"),
    legend.text = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', 
                               color = 'black'),
    legend.justification="left",
    legend.title = element_text(size = 12, hjust = 0, vjust = 0, face = 'bold')
  )+
  labs(fill = 'Cumulative Cases',
       title = '')
fig_min

fig_panel <- ggplot()+
  geom_sf(data = map_data_city1)+
  geom_rect(mapping = aes(xmin = 113.050, xmax = 114.450, ymin = 21.255, ymax = 22.655), 
            color = 'black', fill = 'transparent', size = 1)+
  geom_text(mapping = aes(x = 111.16, y = 25.536, label = 'Guangdong Province'),
            size = 3,
            fontface = "bold",
            family = 'Helvetica')+
  # theme_void()
  theme_bw()+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin=grid::unit(c(0,0,0,0), "cm")
  )

fig_map <- fig_min + inset_element(fig_panel,
                                   left = 0.42, bottom = 0.01, right = 1, top = 0.4, align_to = 'plot')

ggsave(filename = './outcome/science/Figure S1.pdf',
       fig_map,
       width = 7, height = 6, device = cairo_pdf)

ggsave(filename = './outcome/science/Figure S1.png',
       fig_map,
       width = 7, height = 6)

ggsave(filename = './outcome/science/Figure S1.tiff',
       fig_map,
       width = 7, height = 6)
