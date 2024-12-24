
# packages ----------------------------------------------------------------
# devtools::install_github('cran/DMwR')

library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggsci)
library(ggpubr)
library(survival)
library(DMwR)

remove(list = ls())

load('./data/sars_2_cov.Rdata')

library(showtext)
font_add('Helvetica', 
	    "./script/fonts/Helvetica.ttf",
	    bold = "./script/fonts/Helvetica-Bold.ttf")

fill_color <- rev(pal_nejm()(4))[-3]
fill_color_1 <- c('#F1EB58FF', '#AA812DFF', '#283F99FF', "#E64B35FF")

set.seed(202205)


# datafile modified -------------------------------------------------------

datafile_cont_BA1 <- datafile_cont_BA1 %>% 
	mutate(vaccine_dose = as.numeric(str_length(vaccine_type)),
		  vaccine_dose = if_else(is.na(vaccine_dose), 0, vaccine_dose),
		  vaccine_dose = if_else(vaccine > 0 & 
		  				   	vaccinelastdate >= lastexposedate - 14 &
		  				   	!is.na(vaccinelastdate) &
		  				   	!is.na(lastexposedate),
		  				   vaccine_dose - 1,
		  				   vaccine_dose))

datafile_cont_BA2 <- datafile_cont_BA2 %>% 
	mutate(vaccine_dose = as.numeric(str_length(vaccine_type)),
		  vaccine_dose = if_else(is.na(vaccine_dose), 0, vaccine_dose),
		  vaccine_dose = if_else(vaccine > 0 & 
		  				   	vaccinelastdate >= lastexposedate - 14 &
		  				   	!is.na(vaccinelastdate) &
		  				   	!is.na(lastexposedate),
		  				   vaccine_dose - 1,
		  				   vaccine_dose))

datafile_cont_Delta <- datafile_cont_Delta %>% 
	mutate(vaccine_dose = as.numeric(str_length(vaccine_type)),
		  vaccine_dose = if_else(is.na(vaccine_dose), 0, vaccine_dose),
		  vaccine_dose = if_else(vaccine > 0 & 
		  				   	vaccinelastdate >= lastexposedate - 14 &
		  				   	!is.na(vaccinelastdate) &
		  				   	!is.na(lastexposedate),
		  				   vaccine_dose - 1,
		  				   vaccine_dose),
		  vaccine_dose = if_else(is.na(vaccine_type),
		  				   vaccine,
		  				   vaccine_dose))

# Ct value ----------------------------------------------------------------

datafile_info <- rbind(mutate(datafile_info_Delta, 
						lineage = 'Delta'),
						mutate(datafile_info_BA2, 
							  lineage = 'BA.2', 
							  location1 = '厦门',
							  location2 = '厦门'),
						mutate(datafile_info_BA1, 
							  lineage = 'BA.1')) |> 
	mutate(vaccine = factor(vaccine, 
					    levels = c(3, 2, 1, 0),
					    labels = c('A', 'B', 'C', 'D')),
		  vaccine = factor(vaccine, levels = c('D', 'C', 'B', 'A')),
		  lineage = factor(lineage,
		  			  levels = c('Delta', 'BA.1', 'BA.2')),
		  tag = '')


datafile_info_ci <- datafile_info |> 
	group_by(vaccine, lineage) |> 
	summarise(Throat_N_median = median(Throat_N, na.rm = T),
			Throat_N_q1 = quantile(Throat_N, 0.25, na.rm = T),
			Throat_N_q3 = quantile(Throat_N, 0.75, na.rm = T),
			Throat_O_median = median(Throat_O, na.rm = T),
			Throat_O_q1 = quantile(Throat_O, 0.25, na.rm = T),
			Throat_O_q3 = quantile(Throat_O, 0.75, na.rm = T),
			.groups = 'drop')


# plot --------------------------------------------------------------------

fig_ct_voc_n <- ggplot(datafile_info,
				   mapping = aes(x = lineage,
				   		    y = Throat_N,
				   		    fill = lineage)) +
	geom_boxplot(alpha = 0.3,
			   show.legend = F)+
	geom_boxplot(aes(color = lineage),
			   fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0.3,
			   show.legend = F)+
	stat_compare_means(label.y = 42,
				    label.x = 1,
				    family = 'Helvetica',
				    label.x.npc = "center",
				    size = 10*5/14)+
	geom_jitter(aes(color = lineage),
			  alpha = 0.3,
			  na.rm = T,
			  width = 0.2,
			  height = 0,
			  show.legend = F) +
	scale_x_discrete(expand = c(0, 0.6),
				  breaks = c('Delta', 'BA.1', 'BA.2'))+
	scale_y_continuous(limits = c(10, 45),
				    expand = c(0, 0))+
	scale_fill_manual(values = fill_color,
				   na.translate = F)+
	scale_color_manual(values = fill_color,
				    na.translate = F)+
	theme_classic(base_family = 'Helvetica')+
	labs(x = '',
		y = '',
		title = 'b')

fig_ct_vaccine_n <- ggplot(datafile_info,
					  mapping = aes(x = vaccine,
					  		    y = Throat_N,
					  		    fill = vaccine)) +
	geom_boxplot(alpha = 0.3,
			   show.legend = F)+
	geom_boxplot(aes(color = vaccine),
			   fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0.3,
			   show.legend = F)+
	stat_compare_means(label.y = 42,
				    label.x = 1,
				    family = 'Helvetica',
				    label.x.npc = "center",
				    size = 10*5/14)+
	geom_jitter(aes(color = vaccine),
			  alpha = 0.3,
			  na.rm = T,
			  width = 0.2,
			  height = 0,
			  show.legend = F) +
	scale_x_discrete(expand = c(0, 0.6),
				  breaks = c('D', 'C', 'B', 'A'),
				  labels = 0:3)+
	scale_y_continuous(limits = c(10, 45),
				    expand = c(0, 0))+
	scale_fill_manual(values = fill_color_1,
				   na.translate = F)+
	scale_color_manual(values = fill_color_1,
				    na.translate = F)+
	theme_classic(base_family = 'Helvetica')+
	theme(axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 0.5, angle = 0, face = 'plain', color = 'black'))+
	labs(x = 'Vaccine doses',
		y = 'N gene Ct value',
		title = 'a')

fig_ct_voc_o <- ggplot(datafile_info,
				   mapping = aes(x = lineage,
				   		    y = Throat_O,
				   		    fill = lineage)) +
	geom_boxplot(alpha = 0.3,
			   show.legend = F)+
	geom_boxplot(aes(color = lineage),
			   fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0.3,
			   show.legend = F)+
	stat_compare_means(label.y = 42,
				    label.x = 1,
				    family = 'Helvetica',
				    label.x.npc = "center",
				    size = 10*5/14)+
	geom_jitter(aes(color = lineage),
			  alpha = 0.3,
			  na.rm = T,
			  width = 0.2,
			  height = 0,
			  show.legend = F) +
	scale_x_discrete(expand = c(0, 0.6),
				  breaks = c('Delta', 'BA.1', 'BA.2'))+
	scale_y_continuous(limits = c(10, 45),
				    expand = c(0, 0))+
	scale_fill_manual(values = fill_color,
				   na.translate = F)+
	scale_color_manual(values = fill_color,
				    na.translate = F)+
	theme_classic(base_family = 'Helvetica')+
	labs(x = '',
		y = '',
		title = 'd')

fig_ct_vaccine_o <- ggplot(datafile_info,
					  mapping = aes(x = vaccine,
					  		    y = Throat_O,
					  		    fill = vaccine)) +
	geom_boxplot(alpha = 0.3,
			   show.legend = F)+
	geom_boxplot(aes(color = vaccine),
			   fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0.3,
			   show.legend = F)+
	stat_compare_means(label.y = 42,
				    label.x = 1,
				    family = 'Helvetica',
				    label.x.npc = "center",
				    size = 10*5/14)+
	geom_jitter(aes(color = vaccine),
			  alpha = 0.3,
			  na.rm = T,
			  width = 0.2,
			  height = 0,
			  show.legend = F) +
	scale_x_discrete(expand = c(0, 0.6),
				  breaks = c('D', 'C', 'B', 'A'),
				  labels = 0:3)+
	scale_y_continuous(limits = c(10, 45),
				    expand = c(0, 0))+
	scale_fill_manual(values = fill_color_1,
				   na.translate = F)+
	scale_color_manual(values = fill_color_1,
				    na.translate = F)+
	theme_classic(base_family = 'Helvetica')+
	theme(axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 0.5, angle = 0, face = 'plain', color = 'black'))+
	labs(x = 'Vaccine doses',
		y = 'ORF gene Ct value',
		title = 'c')

# Nested ANOVA ------------------------------------------------------------

nest <- aov(datafile_info$Throat_N ~ datafile_info$vaccine/datafile_info$lineage)
nest <- summary(nest)
nest_vaccine <- paste0('Vaccination status: Nested, p = ', sprintf('%.2f', nest[[1]]$`Pr(>F)`[1]))
nest_lineage <- paste0('Vaccination status~Lineage: Nested, p = ', sprintf('%.2f', nest[[1]]$`Pr(>F)`[2]))
label = paste0(nest_vaccine, '\n', nest_lineage)

df_text <- data.frame(x = 2, y = 42,
				  label = label,
				  vaccine = factor('A', levels = c('C', 'B', 'A')))

fig_ct_n <- ggplot(datafile_info,
			    mapping = aes(x = lineage,
			    		    y = Throat_N,
			    		    fill = lineage)) +
	geom_boxplot(alpha = 0.3,
			   show.legend = F)+
	geom_boxplot(aes(color = lineage),
			   fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0.3,
			   show.legend = F)+
	geom_jitter(aes(color = lineage),
			  alpha = 0.3,
			  na.rm = T,
			  width = 0.2,
			  height = 0,
			  show.legend = F) +
	geom_text(data = df_text,
			mapping = aes(x = x,
					    y = y,
					    label = label,
					    fill = NULL),
			hjust = 1)+
	scale_x_discrete(expand = c(0, 0.6),
				  breaks = c('Delta', 'BA.1', 'BA.2'))+
	scale_y_continuous(limits = c(10, 45),
				    expand = c(0, 0))+
	scale_fill_manual(values = fill_color,
				   na.translate = F)+
	scale_color_manual(values = fill_color,
				    na.translate = F)+
	coord_cartesian(clip="off")+
	theme_classic(base_family = 'Helvetica')+
	theme(plot.title = element_text(size = 16, hjust = 0, vjust = 0, face = 'bold'),
		 plot.margin = margin(0, 0.2, 0, 0, "cm"),
		 axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 0.5, face = 'plain', color = 'black'),
		 axis.text.y = element_text(size = 10, hjust = 1, vjust = .5, face = 'plain', color = 'black'),
		 axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold', color = 'black'),
		 axis.title.y = element_text(size = 12, hjust = .5, vjust = 0, face = 'bold', color = 'black'),
		 strip.text.x = element_text(size = 10, face = 'plain'),
		 strip.text.y = element_text(size = 10, face = 'plain'),
		 strip.placement = "outside",
		 strip.background = element_rect(color = NA),
		 panel.spacing.x = unit(0, 'mm'))+
	facet_grid(. ~ factor(vaccine, levels = c('D', 'C', 'B', 'A')),
			 switch = 'both',
			 labeller = as_labeller(c('D' = 'Unvaccine',
			 					'C' = '1 dose',
			 					'B' = '2 doses',
			 					'A' = '3 doses')))+
	labs(x = '',
		y = 'N gene Ct value',
		title = 'e')

nest <- aov(datafile_info$Throat_O ~ datafile_info$vaccine/datafile_info$lineage)
nest <- summary(nest)
nest_vaccine <- paste0('Vaccination status: Nested, p = ', sprintf('%.2f', nest[[1]]$`Pr(>F)`[1]))
nest_lineage <- paste0('Vaccination status~Lineage: Nested, p = ', sprintf('%.2f', nest[[1]]$`Pr(>F)`[2]))
label = paste0(nest_vaccine, '\n', nest_lineage)

df_text <- data.frame(x = 2, y = 42,
				  label = label,
				  vaccine = factor('A', levels = c('C', 'B', 'A')))

fig_ct_o <- ggplot(datafile_info,
			    mapping = aes(x = lineage,
			    		    y = Throat_O,
			    		    fill = lineage)) +
	geom_boxplot(alpha = 0.3,
			   show.legend = F)+
	geom_boxplot(aes(color = lineage),
			   fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0.3,
			   show.legend = F)+
	geom_jitter(aes(color = lineage),
			  alpha = 0.3,
			  na.rm = T,
			  width = 0.2,
			  height = 0,
			  show.legend = F) +
	geom_text(data = df_text,
			mapping = aes(x = x,
					    y = y,
					    label = label,
					    fill = NULL),
			hjust = 1)+
	scale_x_discrete(expand = c(0, 0.6),
				  breaks = c('Delta', 'BA.1', 'BA.2'))+
	scale_y_continuous(limits = c(10, 45),
				    expand = c(0, 0))+
	scale_fill_manual(values = fill_color,
				   na.translate = F)+
	scale_color_manual(values = fill_color,
				    na.translate = F)+
	coord_cartesian(clip="off")+
	theme_classic(base_family = 'Helvetica')+
	theme(plot.title = element_text(size = 16, hjust = 0, vjust = 0, face = 'bold'),
		 plot.margin = margin(0, 0.2, 0, 0, "cm"),
		 axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 0.5, face = 'plain', color = 'black'),
		 axis.text.y = element_text(size = 10, hjust = 1, vjust = .5, face = 'plain', color = 'black'),
		 axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold', color = 'black'),
		 axis.title.y = element_text(size = 12, hjust = .5, vjust = 0, face = 'bold', color = 'black'),
		 strip.text.x = element_text(size = 10, face = 'plain'),
		 strip.text.y = element_text(size = 10, face = 'plain'),
		 strip.placement = "outside",
		 strip.background = element_rect(color = NA),
		 panel.spacing.x = unit(0, 'mm'))+
	facet_grid(. ~ factor(vaccine, levels = c('D', 'C', 'B', 'A')),
			 switch = 'both',
			 labeller = as_labeller(c('D' = 'Unvaccine',
			 					'C' = '1 dose',
			 					'B' = '2 doses',
			 					'A' = '3 doses')))+
	labs(x = '',
		y = 'ORF gene Ct value',
		title = 'f')

# logistics regression ----------------------------------------------------

datafile_cont_BA1$lineage <- 'BA1'
datafile_cont_BA2$lineage <- 'BA2'
datafile_cont_Delta$lineage <- 'Delta'

datafile_cont <- rbind(datafile_cont_Delta,
				   datafile_cont_BA1,
				   mutate(datafile_cont_BA2, gender = NA))|> 
	select(age, vaccine_dose, vaccine_mix, lineage, outcome, 
		  vaccinelastdate, lastexposedate) |> 
	mutate(vaccine_mix = if_else(is.na(vaccine_mix), 
						    'N', 
						    vaccine_mix),
		  vaccine_mix = if_else(vaccine_mix == 'Y', 1, 0),
		  age_g = if_else(age <=18, 'c', 'a'),
		  age_g = if_else(age >=65, 'o', age_g),
		  age_g = factor(age_g, levels = c('c', 'a', 'o')),
		  vaccine_g = vaccine_dose,
		  vaccine_g = as.factor(vaccine_g),
		  outcome = as.factor(outcome),
		  lineage = factor(lineage, levels = c('Delta', 'BA1', 'BA2'))) |> 
	select(vaccine_g, age_g, outcome, lineage)
# 
# datafile_resample_delta <- SMOTE(vaccine_g~.,
#                                  data = filter(datafile_cont, lineage == 'Delta'))
# datafile_resample_delta <- filter(datafile_cont, lineage == 'Delta')
# 
# datafile_resample_ba1 <- SMOTE(vaccine_g~.,
# 						 data = filter(datafile_cont, lineage == 'BA1'))
# 
# datafile_resample_ba2 <- SMOTE(vaccine_g~.,
# 						 data = filter(datafile_cont, lineage == 'BA2'))
# 
# datafile_resample <- rbind(datafile_resample_delta,
# 					  datafile_resample_ba1,
# 					  datafile_resample_ba2)
# 
# # datafile_resample <- SMOTE(vaccine_g ~ .,
# #                            data = datafile_cont)
# 
# table(datafile_resample$lineage)
# 
# datafile_cont <- datafile_resample
datafile_cont$outcome <- as.numeric(datafile_cont$outcome) - 1

# adjust ----------------------------------------------------------------

res_clog_delta <- clogit(formula = outcome ~ vaccine_g  + strata(age_g), 
					data = filter(datafile_cont, lineage == 'Delta')) |> 
	summary() %>%
	.[["conf.int"]] |> 
	as.data.frame() |> 
	mutate(lineage = 'Delta') |> 
	select(-`exp(-coef)`) |> 
	rownames_to_column('var')
res_clog_ba1 <- clogit(formula = outcome ~ vaccine_g  + strata(age_g), 
				   data = filter(datafile_cont, lineage == 'BA1')) |> 
	summary() %>%
	.[["conf.int"]] |> 
	as.data.frame() |> 
	mutate(lineage = 'BA1') |> 
	select(-`exp(-coef)`) |> 
	rownames_to_column('var')
res_clog_ba2 <- clogit(formula = outcome ~ vaccine_g  + strata(age_g), 
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

res_log_delta <- glm(formula = outcome ~ vaccine_g,
				 data = filter(datafile_cont, lineage == 'Delta'),
				 family = binomial(link = "logit")) |>
	summary() %>%
	.[["coefficients"]] |>
	as.data.frame() |>
	mutate(lineage = 'Delta') |>
	rownames_to_column('var')
res_log_ba1 <- glm(formula = outcome ~ vaccine_g,
			    data = filter(datafile_cont, lineage == 'BA1'),
			    family = binomial(link = "logit")) |>
	summary() %>%
	.[["coefficients"]] |>
	as.data.frame() |>
	mutate(lineage = 'BA1') |>
	rownames_to_column('var')
res_log_ba2 <- glm(formula = outcome ~ vaccine_g,
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

datafile_res_adjust$just <- 'Age'
datafile_res_unadjust$just <- 'No'
datafile_res <- rbind(datafile_res_adjust, datafile_res_unadjust) |> 
	mutate(var = factor(var,
					levels = c(paste0('vaccine_g', 1:3))),
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
			 labeller = as_labeller(c('vaccine_g1' = '1 dose vs.\nUnvaccinated',
			 					'vaccine_g2' = '2 doses vs.\nUnvaccinated',
			 					'vaccine_g3' = '3 doses vs.\nUnvaccinated')))+
	scale_x_discrete(expand = c(0, 0.6))+
	scale_y_continuous(limits = c(0, 14),
				    breaks = seq(0, 14, 2),
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
		title = 'g',
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
			 labeller = as_labeller(c('vaccine_g1' = '1 dose vs.\nUnvaccinated',
			 					'vaccine_g2' = '2 doses vs.\nUnvaccinated',
			 					'vaccine_g3' = '3 doses vs.\nUnvaccinated')))+
	scale_x_discrete(expand = c(0, 0.6))+
	scale_y_continuous(limits = c(0, 14),
				    breaks = seq(0, 14, 2),
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
		title = 'h',
		color = 'COVID-19 Variants of Concern')


# combined -----------------------------------------------------------------

gg_ct_o <- ggplotGrob(fig_ct_o)
gg_ct_n <- ggplotGrob(fig_ct_n)

for (i in 1:4) {
	grob.i <- grep("strip-b", gg_ct_n$layout$name)[i]
	gg_ct_n$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- alpha(fill_color_1[i], 0.3)
	grob.i <- grep("strip-b", gg_ct_o$layout$name)[i]
	gg_ct_o$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- alpha(fill_color_1[i], 0.3)
}

fig_log <- fig_log_unjust + fig_log_just&
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
		 legend.position = c(0.5, 1.1),
		 legend.justification = c(0.5, 1),
		 legend.box.margin = margin(0, 0, 0 ,0 , "cm"))

layout <- "
ABCD
EEFF
"

fig_ct <- fig_ct_vaccine_n + fig_ct_voc_n + fig_ct_vaccine_o + fig_ct_voc_o + 
	gg_ct_n + gg_ct_o+
	plot_layout(nrow = 1, design = layout, heights = c(1, 1.4))&
	theme(plot.title = element_text(size = 16, hjust = 0, vjust = 0, face = 'bold'),
		 plot.margin = margin(0, 0.1, 0, 0.1, "cm"),
		 axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
		 axis.text.y = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
		 axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold', color = 'black'),
		 axis.title.y = element_text(size = 12, hjust = .5, vjust = 0, face = 'bold', color = 'black'),
		 strip.text = element_text(size = 10, hjust = .5, vjust = 1, face = 'plain', color = 'black'),
		 strip.placement = "outside",
		 strip.background = element_rect(color = NA),
		 panel.spacing.x = unit(0, 'mm'),
		 legend.position = 'none')

cowplot::plot_grid(fig_ct, fig_log, 
			    ncol = 1,
			    rel_heights = c(2, 1))

ggsave(filename = './outcome/publish/extend/Figure S4.tiff',
	  width = 9, height = 10)
