
# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(ggsci)
library(gridExtra)
library(showtext)

font_add('Helvetica', 
         "./script/fonts/Helvetica.ttf",
         bold = "./script/fonts/Helvetica-Bold.ttf")

remove(list = ls())

load('./data/sars_2_cov.Rdata')

fill_color <- rev(pal_nejm()(4))[-3]

set.seed(202202)

# logistics regression ----------------------------------------------------

datafile_cont_BA1$lineage <- 'BA1'
datafile_cont_BA2$lineage <- 'BA2'
datafile_cont_Delta$lineage <- 'Delta'

datafile_cont_inactive <- rbind(datafile_cont_Delta,
                                datafile_cont_BA1,
                                mutate(datafile_cont_BA2, gender = NA))|> 
        select(age, vaccine, vaccine_mix, vaccine_type, lineage, outcome, 
               vaccinelastdate, lastexposedate) |> 
        mutate(vaccine_mix = if_else(is.na(vaccine_mix), 
                                     'N', 
                                     vaccine_mix),
               vaccine = if_else(vaccine == 0,
                                 0,
                                 vaccine - 1),
               vaccine_mix = if_else(vaccine_mix == 'Y', 1, 0),
               age_g = if_else(age <=18, 'c', 'a'),
               age_g = if_else(age >=65, 'o', age_g),
               age_g = factor(age_g, levels = c('c', 'a', 'o')),
               vaccine_g = str_remove_all(vaccine_type, '[_]'),
               vaccine_g = sapply(vaccine_g, FUN = function(x){
                       paste(sort(unique(strsplit(x, "")[[1]])), collapse = '')
               }),
               vaccine = factor(vaccine),
               outcome = as.factor(outcome),
               lineage = factor(lineage, levels = c('Delta', 'BA1', 'BA2'))) |> 
        filter(vaccine_g %in% c('PV', 'P', 'V', 'K', 'KP', 'KV')) |> 
        select(vaccine, age_g, outcome, lineage)

datafile_resample_delta <- filter(datafile_cont_inactive, lineage == 'Delta')

datafile_resample_ba1 <- SMOTE(vaccine~.,
                               data = filter(datafile_cont_inactive, lineage == 'BA1'),
                               perc.over = 1000,perc.under=200)

datafile_resample_ba2 <- SMOTE(vaccine~.,
                               data = filter(datafile_cont_inactive, lineage == 'BA2'),
                               perc.over = 1000,perc.under=200)

datafile_cont <- rbind(datafile_resample_delta,
                       datafile_resample_ba1,
                       datafile_resample_ba2)

datafile_cont$outcome <- as.numeric(datafile_cont$outcome) - 1
datafile_cont_inactive$outcome <- as.numeric(datafile_cont_inactive$outcome) - 1

# adjust ------------------------------------------------------------------

res_clog_delta <- clogit(formula = outcome ~ vaccine  + strata(age_g), 
                         data = filter(datafile_cont_inactive, lineage == 'Delta')) |> 
        summary() %>%
        .[["conf.int"]] |> 
        as.data.frame() |> 
        mutate(lineage = 'Delta') |> 
        select(-`exp(-coef)`) |> 
        rownames_to_column('var')
res_clog_ba1 <- clogit(formula = outcome ~ vaccine  + strata(age_g), 
                       data = filter(datafile_cont_inactive, lineage == 'BA1')) |> 
        summary() %>%
        .[["conf.int"]] |> 
        as.data.frame() |> 
        mutate(lineage = 'BA1') |> 
        select(-`exp(-coef)`) |> 
        rownames_to_column('var')
res_clog_ba2 <- clogit(formula = outcome ~ vaccine  + strata(age_g), 
                       data = filter(datafile_cont_inactive, lineage == 'BA2')) |> 
        summary() %>%
        .[["conf.int"]] |> 
        as.data.frame() |> 
        mutate(lineage = 'BA2') |> 
        select(-`exp(-coef)`) |> 
        rownames_to_column('var')

datafile_res_adjust_raw <- rbind(res_clog_delta,
                                 res_clog_ba1,
                                 res_clog_ba2)
names(datafile_res_adjust_raw)[2:4] <- c('OR', 'OR_1', 'OR_2')

res_clog_delta <- clogit(formula = outcome ~ vaccine  + strata(age_g), 
                         data = filter(datafile_cont, lineage == 'Delta')) |> 
        summary() %>%
        .[["conf.int"]] |> 
        as.data.frame() |> 
        mutate(lineage = 'Delta') |> 
        select(-`exp(-coef)`) |> 
        rownames_to_column('var')
res_clog_ba1 <- clogit(formula = outcome ~ vaccine  + strata(age_g), 
                       data = filter(datafile_cont, lineage == 'BA1')) |> 
        summary() %>%
        .[["conf.int"]] |> 
        as.data.frame() |> 
        mutate(lineage = 'BA1') |> 
        select(-`exp(-coef)`) |> 
        rownames_to_column('var')
res_clog_ba2 <- clogit(formula = outcome ~ vaccine  + strata(age_g), 
                       data = filter(datafile_cont, lineage == 'BA2')) |> 
        summary() %>%
        .[["conf.int"]] |> 
        as.data.frame() |> 
        mutate(lineage = 'BA2') |> 
        select(-`exp(-coef)`) |> 
        rownames_to_column('var')

datafile_res_adjust <- rbind(res_clog_delta,
                             res_clog_ba1,
                             res_clog_ba2)
names(datafile_res_adjust)[2:4] <- c('OR', 'OR_1', 'OR_2')

# unadjust ----------------------------------------------------------------

res_log_delta <- glm(formula = outcome ~ vaccine,
                     data = filter(datafile_cont_inactive, lineage == 'Delta'),
                     family = binomial(link = "logit")) |>
        summary() %>%
        .[["coefficients"]] |>
        as.data.frame() |>
        mutate(lineage = 'Delta') |>
        rownames_to_column('var')
res_log_ba1 <- glm(formula = outcome ~ vaccine,
                   data = filter(datafile_cont_inactive, lineage == 'BA1'),
                   family = binomial(link = "logit")) |>
        summary() %>%
        .[["coefficients"]] |>
        as.data.frame() |>
        mutate(lineage = 'BA1') |>
        rownames_to_column('var')
res_log_ba2 <- glm(formula = outcome ~ vaccine,
                   data = filter(datafile_cont_inactive, lineage == 'BA2'),
                   family = binomial(link = "logit")) |>
        summary() %>%
        .[["coefficients"]] |>
        as.data.frame() |>
        mutate(lineage = 'BA2') |>
        rownames_to_column('var')

datafile_res_unadjust <- rbind(res_log_delta,
                               res_log_ba1,
                               res_log_ba2)
datafile_res_unadjust_raw <- datafile_res_unadjust[-grep('age_|Intercept', datafile_res_unadjust$var),] |> 
        select(var, lineage, Estimate, `Std. Error`) |> 
        mutate(OR = exp(Estimate),
               OR_1 = exp(Estimate - 1.96*`Std. Error`),
               OR_2 = exp(Estimate + 1.96*`Std. Error`)) |> 
        select(var, OR, OR_1, OR_2, lineage)

res_log_delta <- glm(formula = outcome ~ vaccine,
                     data = filter(datafile_cont, lineage == 'Delta'),
                     family = binomial(link = "logit")) |>
        summary() %>%
        .[["coefficients"]] |>
        as.data.frame() |>
        mutate(lineage = 'Delta') |>
        rownames_to_column('var')
res_log_ba1 <- glm(formula = outcome ~ vaccine,
                   data = filter(datafile_cont, lineage == 'BA1'),
                   family = binomial(link = "logit")) |>
        summary() %>%
        .[["coefficients"]] |>
        as.data.frame() |>
        mutate(lineage = 'BA1') |>
        rownames_to_column('var')
res_log_ba2 <- glm(formula = outcome ~ vaccine,
                   data = filter(datafile_cont, lineage == 'BA2'),
                   family = binomial(link = "logit")) |>
        summary() %>%
        .[["coefficients"]] |>
        as.data.frame() |>
        mutate(lineage = 'BA2') |>
        rownames_to_column('var')

datafile_res_unadjust <- rbind(res_log_delta,
                               res_log_ba1,
                               res_log_ba2)
datafile_res_unadjust <- datafile_res_unadjust[-grep('age_|Intercept', datafile_res_unadjust$var),] |> 
        select(var, lineage, Estimate, `Std. Error`) |> 
        mutate(OR = exp(Estimate),
               OR_1 = exp(Estimate - 1.96*`Std. Error`),
               OR_2 = exp(Estimate + 1.96*`Std. Error`)) |> 
        select(var, OR, OR_1, OR_2, lineage)

# plot --------------------------------------------------------------------

datafile_res_adjust_raw$just <- 'Age'
datafile_res_unadjust_raw$just <- 'No'
datafile_res_raw <- rbind(datafile_res_adjust_raw, datafile_res_unadjust_raw) |> 
        mutate(var = factor(var,
                            levels = c(paste0('vaccine', 1:2))),
               lineage = factor(lineage,
                                levels = c('Delta', 'BA1', 'BA2'))) |> 
        filter(!is.na(OR_1) & OR_1>0 & OR_2<100) |> 
        mutate_at(vars(OR, OR_1, OR_2), round, digits = 4)

datafile_res_adjust$just <- 'Age'
datafile_res_unadjust$just <- 'No'
datafile_res <- rbind(datafile_res_adjust, datafile_res_unadjust) |> 
        mutate(var = factor(var,
                            levels = c(paste0('vaccine', 1:2))),
               lineage = factor(lineage,
                                levels = c('Delta', 'BA1', 'BA2'))) |> 
        filter(!is.na(OR_1) & OR_1>0 & OR_2<100) |> 
        mutate_at(vars(OR, OR_1, OR_2), round, digits = 4)

fig_log_unjust <- ggplot(data = filter(datafile_res, just == 'No'),
                         mapping = aes(x = lineage,
                                       y = OR,
                                       color = lineage))+
        geom_point()+
        geom_pointrange(mapping = aes(ymin = OR_1, 
                                      ymax = OR_2))+
        geom_hline(yintercept = 1,
                   linetype = 'dashed',
                   color = 'black')+
        facet_grid(cols = vars(var),
                   switch = 'both',
                   labeller = as_labeller(c('vaccine1' = 'Fully vaccinated vs.\nUnfully vaccinated',
                                            'vaccine2' = 'Booster dose vs.\nUnfully vaccinated')))+
        scale_x_discrete(expand = c(0, 0.6))+
        scale_y_continuous(limits = c(0, 3),
                           breaks = seq(0, 3, 1),
                           expand = c(0, 0))+
        scale_color_manual(values = fill_color,
                           labels = c('Delta', 'BA.1', 'BA.2'),
                           na.translate = F)+
        theme_classic(base_family = 'Helvetica')+
        guides(color=guide_legend(direction='horizontal',
                                  title.position = 'top',
                                  title.hjust = 0.5))+
        labs(x = '',
             y = 'Odds ratio',
             title = 'a',
             color = 'COVID-19 Variants of Concern')

fig_log_just <- ggplot(data = filter(datafile_res, just == 'Age'),
                       mapping = aes(x = lineage,
                                     y = OR,
                                     color = lineage))+
        geom_point()+
        geom_pointrange(mapping = aes(ymin = OR_1, 
                                      ymax = OR_2))+
        geom_hline(yintercept = 1,
                   linetype = 'dashed',
                   color = 'black')+
        facet_grid(cols = vars(var),
                   switch = 'both',
                   labeller = as_labeller(c('vaccine1' = 'Fully vaccinated vs.\nUnfully vaccinated',
                                            'vaccine2' = 'Booster dose vs.\nUnfully vaccinated')))+
        scale_x_discrete(expand = c(0, 0.6))+
        scale_y_continuous(limits = c(0, 3),
                           breaks = seq(0, 3, 1),
                           expand = c(0, 0))+
        scale_color_manual(values = fill_color,
                           labels = c('Delta', 'BA.1', 'BA.2'),
                           na.translate = F)+
        theme_classic(base_family = 'Helvetica')+
        guides(color=guide_legend(direction='horizontal',
                                  title.position = 'top',
                                  title.hjust = 0.5))+
        labs(x = '',
             y = 'Adjusted odds ratio',
             title = 'b',
             color = 'COVID-19 Variants of Concern')

fig_log_unjust_raw <- ggplot(data = filter(datafile_res_raw, just == 'No'),
                         mapping = aes(x = lineage,
                                       y = OR,
                                       color = lineage))+
        geom_point()+
        geom_pointrange(mapping = aes(ymin = OR_1, 
                                      ymax = OR_2))+
        geom_hline(yintercept = 1,
                   linetype = 'dashed',
                   color = 'black')+
        facet_grid(cols = vars(var),
                   switch = 'both',
                   labeller = as_labeller(c('vaccine1' = 'Fully vaccinated vs.\nUnfully vaccinated',
                                            'vaccine2' = 'Booster dose vs.\nUnfully vaccinated')))+
        scale_x_discrete(expand = c(0, 0.6))+
        scale_y_continuous(limits = c(0, 3),
                           breaks = seq(0, 3, 1),
                           expand = c(0, 0))+
        scale_color_manual(values = fill_color,
                           labels = c('Delta', 'BA.1', 'BA.2'),
                           na.translate = F)+
        theme_classic(base_family = 'Helvetica')+
        guides(color=guide_legend(direction='horizontal',
                                  title.position = 'top',
                                  title.hjust = 0.5))+
        labs(x = '',
             y = 'Odds ratio',
             title = 'c',
             color = 'COVID-19 Variants of Concern')

fig_log_just_raw <- ggplot(data = filter(datafile_res_raw, just == 'Age'),
                       mapping = aes(x = lineage,
                                     y = OR,
                                     color = lineage))+
        geom_point()+
        geom_pointrange(mapping = aes(ymin = OR_1, 
                                      ymax = OR_2))+
        geom_hline(yintercept = 1,
                   linetype = 'dashed',
                   color = 'black')+
        facet_grid(cols = vars(var),
                   switch = 'both',
                   labeller = as_labeller(c('vaccine1' = 'Fully vaccinated vs.\nUnfully vaccinated',
                                            'vaccine2' = 'Booster dose vs.\nUnfully vaccinated')))+
        scale_x_discrete(expand = c(0, 0.6))+
        scale_y_continuous(limits = c(0, 3),
                           breaks = seq(0, 3, 1),
                           expand = c(0, 0))+
        scale_color_manual(values = fill_color,
                           labels = c('Delta', 'BA.1', 'BA.2'),
                           na.translate = F)+
        theme_classic(base_family = 'Helvetica')+
        guides(color=guide_legend(direction='horizontal',
                                  title.position = 'top',
                                  title.hjust = 0.5))+
        labs(x = '',
             y = 'Adjusted odds ratio',
             title = 'd',
             color = 'COVID-19 Variants of Concern')

fig_log <- fig_log_unjust + fig_log_just + fig_log_unjust_raw + fig_log_just_raw+
        plot_layout(guides = 'collect')&
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              plot.title = element_text(size = 16, hjust = 0, vjust = 0, face = 'bold'),
              plot.margin = margin(0, 0.1, 0, 0.1, "cm"),
              axis.text.y = element_text(size = 10, hjust = 1, vjust = .5, face = 'plain', color = 'black'),
              axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold', color = 'black'),
              axis.title.y = element_text(size = 12, hjust = .5, vjust = 0, face = 'bold', color = 'black'),
              strip.text.x = element_text(size = 10, face = 'plain'),
              strip.text.y = element_text(size = 10, face = 'plain'),
              strip.placement = "outside",
              strip.background = element_rect(color = NA),
              panel.spacing.x = unit(0, 'mm'),
              legend.position = 'bottom',
              legend.box.margin = margin(0, 0, 0 ,0 , "cm"))

# ggsave(filename = './outcome/publish/extend/Figure S3.pdf', 
#        fig_log,
#        height = 3.5, width = 7)
