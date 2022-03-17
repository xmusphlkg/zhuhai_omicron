
# packages ----------------------------------------------------------------

library(tidyverse)


extrafont::font_import('./script/fonts', prompt = F)
extrafont::loadfonts(device = "win")

remove(list = ls())

load('./data/data_case.Rdata')

# data analysis -----------------------------------------------------------

datafile_cont_matrix <- datafile_cont_all %>% 
  select(age, durationexpose, id_cases) %>% 
  left_join(select(datafile_info_clean, id, age), by = c('id_cases' = 'id')) %>% 
  # select(-id_cases) %>% 
  filter(!is.na(age.x)) %>% 
  mutate(age.x = cut(age.x, 
                     breaks = 10*(0:9),
                     include.lowest = T,
                     right = T),
         # age.y = cut(age.y, 
         #             breaks = 10*(0:9),
         #             include.lowest = T,
         #             right = T),
         
         )

labels_matrix <- levels(datafile_cont_matrix$age.x)

df_matrix_x <- summary(datafile_cont_matrix$age.x)
df_matrix_y <- summary(datafile_cont_matrix$age.y)

datafile_cont_matrix <- datafile_cont_matrix %>% 
  rbind(data.frame(age.x = datafile_cont_matrix$age.y,
                   contactfreq = datafile_cont_matrix$contactfreq,
                   age.y = datafile_cont_matrix$age.x)) %>% 
  mutate(age.x = as.numeric(age.x),
         age.y = as.numeric(age.y)) %>% 
  # filter(age.x >= age.y) %>% 
  group_by(age.x, age.y) %>% 
  count() %>%
  # summarise(n = sum(contactfreq, na.rm = T)) %>%
  # mutate(n = ifelse(age.x < age.y,
  #                    n * df_matrix_x[age.x] / df_matrix_y[age.y],
  #                    n)) %>%
  ungroup() %>% 
  pivot_wider(names_from = age.x,
              values_from = n) %>% 
  as.data.frame()

rownames(datafile_cont_matrix) <- datafile_cont_matrix$age.y
datafile_cont_matrix <- as.matrix(round(datafile_cont_matrix[,-1]))

# for (c in 1:5) {
#   for (r in (c+1):6) {
#     datafile_cont_matrix[r, c] <- round(datafile_cont_matrix[c, r] * df_matrix_x[r] / df_matrix_y[c])
#   }
# }

datafile_cont_plot <- datafile_cont_matrix %>% 
  as.data.frame() %>% 
  mutate(age.y = as.numeric(rownames(.))) %>% 
  pivot_longer(!age.y, names_to = "age.x", values_to = "value") %>% 
  filter(!is.na(value)) %>% 
  mutate(age.x = age.x)

# plot --------------------------------------------------------------------

ggplot(datafile_cont_plot, aes(x = age.x, y = age.y)) + 
  geom_tile(aes(fill = value), colour = "white") + 
  scale_x_discrete(breaks = as.character(1:9), 
                   labels = labels_matrix,
                   expand = c(0, 0),
                   position = 'top') + 
  scale_y_reverse(breaks = 1:9,
                  labels = labels_matrix,
                  expand = c(0, 0))+
  scale_fill_viridis_c(direction = -1)+
  theme_bw()+
  theme(legend.position = "right", 
        legend.key.height = unit(1, "cm"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1)
        ) + 
  coord_equal()+
  labs(x = '',
       y = '',
       fill = '')
