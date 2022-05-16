
# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(SimInf)
library(ggsci)
library(scales)

remove(list = ls())
load('data/data_case.Rdata')
load('outcome/base/para.RData')

extrafont::font_import('./script/fonts', prompt = F)
extrafont::loadfonts(device = "win")

# parameter ---------------------------------------------------------------

set.seed(202202)

## Simulate Times
n <- 1000

## Zhuhai City Population
para_pop <- 2439585

## Zhuhai City Population Accepted Vaccine
para_vac <- round(para_pop * 0.901)
# para_vac <- 141.7 * 10e3

## Vaccine Effectiveness
para_sus <- 0.5267
para_trans <- 0

## Serial Interval
time_si <- si
time_si_shape <- (time_si$mean/time_si$sd)^2
time_si_rate <- time_si_shape / time_si$mean

## Incubation Period
time_ib <- ib
time_ib_shape <- time_ib$distribution$parameters$shape
time_ib_scale <- time_ib$distribution$parameters$scale

## Infectious Period Asymptomatic
time_ip_asy <- rgamma(n, shape = time_ib_shape, scale = time_ib_scale)
gamma_2 <- 1 / time_ip_asy

## Infectious Period Pre-symptomatic
time_ip_pre <- rgamma(n, shape = time_ib_shape - time_si_shape, rate = time_si_rate)
gamma_1 <- 1 / time_ip_pre

## Infectious Period Symptomatic
time_ip_sym <- rgamma(n, shape = time_si_shape, rate = time_si_rate)
gamma_3 <- 1 / time_ip_sym

## Symptomatic cases communicability period

# github_url <- 'https://raw.githubusercontent.com/blythejane/covid_safety/main/discordant_data.csv'
# github_data <- read_csv(github_url, show_col_types = FALSE) %>% 
#   select(case, day, infectious) %>% 
#   filter(infectious & 
#            case %in% c('B', 'C', 'F', 'H*', 'L', 'P', 'Q', 'S', 'CC')) %>% ## filter case contain the end of communicability period
#   group_by(case) %>% 
#   summarise(d = max(day))

time_cp_sym <- rep(4.5, n)
gamma_5 <- 1 / time_cp_sym

## Asymptomatic Cases Communicability Period
time_cp_asy <- rep(4.5, n)
gamma_4 <- 1 / time_cp_asy

## Percent of Asymptomatic Case
p <- rep(0.31, n)

## Beta (the transmission rate)
R <- 5
# beta <- rep(0.5, n)
# data clean --------------------------------------------------------------

data_clean <- datafile_info_clean %>% 
  group_by(dateonset) %>% 
  count() %>% 
  rename(c(
    'cases' = 'n',
    'date' = 'dateonset'
  )) %>% 
  as.data.frame() %>% 
  complete(
    date = seq.Date(from = min(date), to = max(date), by = 'day'),
    fill = list(cases = 0)
  ) %>% 
  mutate(city = 'Zhuhai',
         intervention = rep(c(0,1), c(8, 6)),
         pop = para_pop)

# simulate ----------------------------------------------------------------

ldata <- data.frame(gamma_1 = gamma_1,
                    gamma_2 = gamma_2,
                    gamma_3 = gamma_3,
                    gamma_4 = gamma_4,
                    gamma_5 = gamma_5,
                    p = p,
                    vsi = rep(para_sus, n),
                    vei = rep(para_trans, n),
                    N1 = rep(para_pop - para_vac, n),
                    N2 = rep(para_vac, n))

ldata <- ldata %>% 
  mutate(beta = R * ((1-p)*gamma_1 + p*gamma_2)/
           ((1-p)*gamma_1/gamma_3 + p*gamma_2/gamma_4 + (1-p)*gamma_1/gamma_5))

model  <- mparse(transitions = c("S1 -> beta*S1*((I2_p + I2_s + I2_a)*(1 - vei) + I1_p + I1_s + I1_a)/N1 -> E1",
                                 "E1 -> gamma_1*(1-p)*E1 -> I1_p",
                                 "E1 -> gamma_2*p*E1 -> I1_a",
                                 "I1_p -> gamma_3*I1_p -> I1_s",
                                 "I1_a -> gamma_4*I1_a -> R1",
                                 "I1_s -> gamma_5*I1_s -> R1",
                                 ## accept vaccine
                                 "S2 -> beta*S2*((I2_p + I2_s + I2_a)*(1 - vei) + I1_p + I1_s + I1_a)*(1 - vsi)/N2 -> E2",
                                 "E2 -> gamma_1*(1-p)*E2 -> I2_p",
                                 "E2 -> gamma_2*p*E2 -> I2_a",
                                 "I2_p -> gamma_3*I2_p -> I2_s",
                                 "I2_a -> gamma_4*I2_a -> R2",
                                 "I2_s -> gamma_5*I2_s -> R2"),
                 compartments = c("S1", "E1", "I1_a", "I1_p", "I1_s", "R1",
                                  "S2", "E2", "I2_a", "I2_p", "I2_s", "R2"),
                 ldata = ldata,
                 u0 = data.frame(S1 = rep(para_pop - para_vac, n),
                                 E1 = rep(0, n),
                                 I1_a = rep(0, n),
                                 I1_p = rep(0, n),
                                 I1_s = rep(0, n),
                                 R1 = rep(0, n),
                                 S2 = rep(para_vac, n),
                                 E2 = rep(0, n),
                                 I2_a = rep(0, n),
                                 I2_p = rep(1, n),
                                 I2_s = rep(0, n),
                                 R2 = rep(0, n)),
                 tspan = 1:200)
result <- run(model)

df_nto1 <- function(x){
  datafile_select <- outcome %>% 
    .[,seq((x-1)*12+1, x*12)]
  
  datafile_outcome <- data.frame(
    t         = 1:nrow(datafile_select),
    S         = datafile_select[,1] + datafile_select[,7],
    E         = datafile_select[,2] + datafile_select[,8],
    I_a_unvac = datafile_select[,3],
    I_a_vac   = datafile_select[,9],
    I_a       = datafile_select[,3] + datafile_select[,9],
    I_unvac_p = datafile_select[,4],
    I_vac_p   = datafile_select[,10],
    I_p       = datafile_select[,4] + datafile_select[,10],
    I_unvac_s = datafile_select[,5],
    I_vac_s   = datafile_select[,11],
    I_s       = datafile_select[,5] + datafile_select[,11],
    R         = datafile_select[,6] + datafile_select[,12]
  ) %>% 
    mutate(I_total       = I_a + I_p + I_s,
           I_total_unvac = I_a_unvac + I_unvac_p + I_unvac_s,
           I_total_vac   = I_a_vac + I_vac_p + I_vac_s,
           I_c           = I_p + I_s,
           I_c_unvac     = I_unvac_p + I_unvac_s,
           I_c_vac       = I_vac_p + I_vac_s,
           n             = x)
  return(datafile_outcome)
}

outcome <- as.data.frame(t(result@U))
outcome <- lapply(seq(1, ncol(outcome)/12), df_nto1)

save(result, outcome, file = './outcome/base/simu_vac.RData')

load('./outcome/base/simu_vac.RData')

# data modify -------------------------------------------------------------

simulate_summary <- function(x){
  # x <- 2
  datafile_select <- outcome[[x]] %>% 
    select(t, I_a, I_p, I_s, R, I_total)
  
  if(max(datafile_select$R) == 1){
    return(c(0, 0, 0, 0))
  }else{
    total_attack_rate <- max(datafile_select$R)/para_pop
    I_total <- datafile_select[datafile_select$R>0, 'I_total']
    duration_time <- min(which(I_total == 0))
    duration_time <- ifelse(is.infinite(duration_time), 150, duration_time)
    peak_time <- which(I_total == max(I_total))[1]
    peak_attack_rate <- datafile_select[datafile_select$t == peak_time, "I_total"]/para_pop
    return(c(total_attack_rate, duration_time, peak_time, peak_attack_rate))
  }
}

datafile_summary <- lapply(1:length(outcome), simulate_summary)
datafile_summary <- do.call('rbind', datafile_summary)

names(datafile_summary) <- c('TAR', 'Duration', 'Peak', 'PAR')
datafile_summary_outcome <- data.frame(
  name = c('TAR', 'Duration', 'Peak', 'PAR'),
  median = apply(datafile_summary, 2, median),
  Q1 = apply(datafile_summary, 2, quantile, probs = 0.25),
  Q3 = apply(datafile_summary, 2, quantile, probs = 0.75)
)

write.csv(datafile_summary_outcome, file = './outcome/base/simulate_vaccine.csv')

# plot --------------------------------------------------------------------

outcome <- do.call('rbind', outcome)

fill_color <- pal_nejm()(3)
case_list <- list(c('I_a', 'I_a_unvac', 'I_a_vac'), 
                  c('I_c', 'I_c_unvac', 'I_c_vac'),
                  c('I_total', 'I_total_unvac', 'I_total_vac'))
y_label <- c("Number of infections",
             "Number of infections",
             "Number of infections")
limit_value <- rep(c(2e5, 10e5, 10e5), 2)

scientific_10 <- function(x) {
  # parse(text=gsub("e\\+", " %*% 10^", scales::scientific_format()(x)))
  ifelse(x == 0, 0, sprintf('%.1f', x/1e5))
}

datafile_ci <- outcome %>% 
  group_by(t) %>% 
  summarize(I_a_min    = quantile(I_a, 0.025),
            I_a_q1     = quantile(I_a, 0.25),
            I_a_median = quantile(I_a, 0.5),
            I_a_mean   = mean(I_a),
            I_a_q3     = quantile(I_a, 0.75),
            I_a_max    = quantile(I_a, 0.975),
            
            I_a_unvac_min    = quantile(I_a_unvac, 0.025),
            I_a_unvac_q1     = quantile(I_a_unvac, 0.25),
            I_a_unvac_median = quantile(I_a_unvac, 0.5),
            I_a_unvac_mean   = mean(I_a_unvac),
            I_a_unvac_q3     = quantile(I_a_unvac, 0.75),
            I_a_unvac_max    = quantile(I_a_unvac, 0.975),
            
            I_a_vac_min    = quantile(I_a_vac, 0.025),
            I_a_vac_q1     = quantile(I_a_vac, 0.25),
            I_a_vac_median = quantile(I_a_vac, 0.5),
            I_a_vac_mean   = mean(I_a_vac),
            I_a_vac_q3     = quantile(I_a_vac, 0.75),
            I_a_vac_max    = quantile(I_a_vac, 0.975),
            
            I_p_min    = quantile(I_p, 0.025),
            I_p_q1     = quantile(I_p, 0.25),
            I_p_median = quantile(I_p, 0.5),
            I_p_mean   = mean(I_p),
            I_p_q3     = quantile(I_p, 0.75),
            I_p_max    = quantile(I_p, 0.975),
            
            I_unvac_p_min    = quantile(I_unvac_p, 0.025),
            I_unvac_p_q1     = quantile(I_unvac_p, 0.25),
            I_unvac_p_median = quantile(I_unvac_p, 0.5),
            I_unvac_p_mean   = mean(I_unvac_p),
            I_unvac_p_q3     = quantile(I_unvac_p, 0.75),
            I_unvac_p_max    = quantile(I_unvac_p, 0.975),
            
            I_vac_p_min    = quantile(I_vac_p, 0.025),
            I_vac_p_q1     = quantile(I_vac_p, 0.25),
            I_vac_p_median = quantile(I_vac_p, 0.5),
            I_vac_p_mean   = mean(I_vac_p),
            I_vac_p_q3     = quantile(I_vac_p, 0.75),
            I_vac_p_max    = quantile(I_vac_p, 0.975),
            
            I_c_min    = quantile(I_c, 0.025),
            I_c_q1     = quantile(I_c, 0.25),
            I_c_median = quantile(I_c, 0.5),
            I_c_mean   = mean(I_c),
            I_c_q3     = quantile(I_c, 0.75),
            I_c_max    = quantile(I_c, 0.975),
            
            I_c_unvac_min    = quantile(I_c_unvac, 0.025),
            I_c_unvac_q1     = quantile(I_c_unvac, 0.25),
            I_c_unvac_median = quantile(I_c_unvac, 0.5),
            I_c_unvac_mean   = mean(I_c_unvac),
            I_c_unvac_q3     = quantile(I_c_unvac, 0.75),
            I_c_unvac_max    = quantile(I_c_unvac, 0.975),
            
            I_c_vac_min    = quantile(I_c_vac, 0.025),
            I_c_vac_q1     = quantile(I_c_vac, 0.25),
            I_c_vac_median = quantile(I_c_vac, 0.5),
            I_c_vac_mean   = mean(I_c_vac),
            I_c_vac_q3     = quantile(I_c_vac, 0.75),
            I_c_vac_max    = quantile(I_c_vac, 0.975),
            
            I_total_min      = quantile(I_total, 0.025),
            I_total_q1       = quantile(I_total, 0.25),
            I_total_median   = median(I_total),
            I_total_mean     = quantile(I_total, 0.5),
            I_total_q3       = quantile(I_total, 0.75),
            I_total_max      = quantile(I_total, 0.975),
            
            I_total_unvac_min      = quantile(I_total_unvac, 0.025),
            I_total_unvac_q1       = quantile(I_total_unvac, 0.25),
            I_total_unvac_median   = median(I_total_unvac),
            I_total_unvac_mean     = quantile(I_total_unvac, 0.5),
            I_total_unvac_q3       = quantile(I_total_unvac, 0.75),
            I_total_unvac_max      = quantile(I_total_unvac, 0.975),
            
            I_total_vac_min      = quantile(I_total_vac, 0.025),
            I_total_vac_q1       = quantile(I_total_vac, 0.25),
            I_total_vac_median   = median(I_total_vac),
            I_total_vac_mean     = quantile(I_total_vac, 0.5),
            I_total_vac_q3       = quantile(I_total_vac, 0.75),
            I_total_vac_max      = quantile(I_total_vac, 0.975)
            )

for (x in 1:3) {
  # x <- 1
  datafile_raw_plot <- outcome[,c('t', 'n', case_list[[x]])]
  names(datafile_raw_plot)[3:5] <- c('var', 'var_unvac', 'var_vac')
  
  datafile_ci_plot <- datafile_ci %>% 
    select(c(t, contains(all_of(case_list[[x]]))))
  names(datafile_ci_plot) <- str_remove_all(names(datafile_ci_plot), paste0(case_list[[x]][1], '_'))
  
  fig <- ggplot()+
    geom_line(data = datafile_raw_plot,
              mapping = aes(x = t, y = var, group = n),
              alpha = 0.03, color = fill_color[x],
              show.legend = F)+
    geom_line(data = datafile_ci_plot,
              mapping = aes(x = t, y =median,
                            group = 'Total',
                            linetype = 'Total'),
              size = 0.7, color = fill_color[x],
              show.legend = F)+
    geom_line(data = datafile_ci_plot,
              mapping = aes(x = t, y =unvac_median, 
                            group = 'Non-vaccinated',
                            linetype = 'Non-vaccinated'),
              size = 0.7, color = fill_color[x],
              show.legend = F)+
    geom_line(data = datafile_ci_plot,
              mapping = aes(x = t, y =vac_median, 
                            group = 'Fully vaccinated',
                            linetype = 'Fully vaccinated'),
              size = 0.7, color = fill_color[x],
              show.legend = F)+
    coord_cartesian(ylim = c(0, limit_value[x]),
                    xlim = c(0, 120))+
    scale_x_continuous(expand = c(0,0),
                       breaks = c(0, 30, 60, 90, 120))+
    scale_y_continuous(expand = c(0, NA),
                       labels = scientific_10)+
    scale_linetype_manual(values = c('dashed', 'dotted', 'solid'))+
    theme_bw(base_family = 'Helvetica')+
    theme(legend.text = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', 
                                     color = 'black'),
          legend.title = element_text(size = 12, hjust = 0, vjust = 0.5, face = 'bold'),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,-7),
          legend.box.background = element_rect(fill = "transparent", colour = 'transparent'),
          legend.background = element_rect(fill = "transparent", colour = 'transparent'),
          legend.position = "none",
          axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain',
                                     color = 'black'),
          axis.text.y = element_text(size = 10, hjust = 1, vjust = .5, face = 'plain', 
                                     color = 'black'),
          axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
          axis.title.y = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
          plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(5, 5, 5, 5))+
    labs(x = 'Days',
         # y = paste0(y_label[x], parse(text = "(10^5)")),
         y = expr(paste(!!y_label[x], " (", 10^5, ")", sep = '')),
         title = letters[x])
  
  assign(paste0('fig_temp_', x), fig)
}

fig_temp_1 + fig_temp_2 + fig_temp_3

# non vaccine -------------------------------------------------------------

# parameter ---------------------------------------------------------------

para_sus <- 0
para_trans <- 0

# simulate ----------------------------------------------------------------

ldata <- data.frame(gamma_1 = gamma_1,
                    gamma_2 = gamma_2,
                    gamma_3 = gamma_3,
                    gamma_4 = gamma_4,
                    gamma_5 = gamma_5,
                    p = p,
                    vsi = rep(para_sus, n),
                    vei = rep(para_trans, n),
                    N1 = rep(para_pop - para_vac, n),
                    N2 = rep(para_vac, n))

ldata <- ldata %>% 
  mutate(beta = R * ((1-p)*gamma_1 + p*gamma_2)/
           ((1-p)*gamma_1/gamma_3 + p*gamma_2/gamma_4 + (1-p)*gamma_1/gamma_5))

model  <- mparse(transitions = c("S1 -> beta*S1*((I2_p + I2_s + I2_a)*(1 - vei) + I1_p + I1_s + I1_a)/N1 -> E1",
                                 "E1 -> gamma_1*(1-p)*E1 -> I1_p",
                                 "E1 -> gamma_2*p*E1 -> I1_a",
                                 "I1_p -> gamma_3*I1_p -> I1_s",
                                 "I1_a -> gamma_4*I1_a -> R1",
                                 "I1_s -> gamma_5*I1_s -> R1",
                                 ## accept vaccine
                                 "S2 -> beta*S2*((I2_p + I2_s + I2_a)*(1 - vei) + I1_p + I1_s + I1_a)*(1 - vsi)/N2 -> E2",
                                 "E2 -> gamma_1*(1-p)*E2 -> I2_p",
                                 "E2 -> gamma_2*p*E2 -> I2_a",
                                 "I2_p -> gamma_3*I2_p -> I2_s",
                                 "I2_a -> gamma_4*I2_a -> R2",
                                 "I2_s -> gamma_5*I2_s -> R2"),
                 compartments = c("S1", "E1", "I1_a", "I1_p", "I1_s", "R1",
                                  "S2", "E2", "I2_a", "I2_p", "I2_s", "R2"),
                 ldata = ldata,
                 u0 = data.frame(S1 = rep(para_pop - para_vac, n),
                                 E1 = rep(0, n),
                                 I1_a = rep(0, n),
                                 I1_p = rep(0, n),
                                 I1_s = rep(0, n),
                                 R1 = rep(0, n),
                                 S2 = rep(para_vac, n),
                                 E2 = rep(0, n),
                                 I2_a = rep(0, n),
                                 I2_p = rep(1, n),
                                 I2_s = rep(0, n),
                                 R2 = rep(0, n)),
                 tspan = 1:200)
result <- run(model)

outcome <- as.data.frame(t(result@U))
outcome <- lapply(seq(1, ncol(outcome)/12), df_nto1)

save(result, outcome, file = './outcome/base/simu_unvac.RData')
load('./outcome/base/simu_unvac.RData')

# data modify -------------------------------------------------------------

datafile_summary <- lapply(1:length(outcome), simulate_summary)
datafile_summary <- do.call('rbind', datafile_summary)

names(datafile_summary) <- c('TAR', 'Duration', 'Peak', 'PAR')
datafile_summary_outcome <- data.frame(
  name = c('TAR', 'Duration', 'Peak', 'PAR'),
  median = apply(datafile_summary, 2, median),
  Q1 = apply(datafile_summary, 2, quantile, probs = 0.25),
  Q3 = apply(datafile_summary, 2, quantile, probs = 0.75)
)

write.csv(datafile_summary_outcome, file = './outcome/base/simulate_unvaccine.csv')

# plot --------------------------------------------------------------------

outcome <- do.call('rbind', outcome)

datafile_ci <- outcome %>% 
  group_by(t) %>% 
  summarize(I_a_min    = quantile(I_a, 0.025),
            I_a_q1     = quantile(I_a, 0.25),
            I_a_median = quantile(I_a, 0.5),
            I_a_mean   = mean(I_a),
            I_a_q3     = quantile(I_a, 0.75),
            I_a_max    = quantile(I_a, 0.975),
            
            I_a_unvac_min    = quantile(I_a_unvac, 0.025),
            I_a_unvac_q1     = quantile(I_a_unvac, 0.25),
            I_a_unvac_median = quantile(I_a_unvac, 0.5),
            I_a_unvac_mean   = mean(I_a_unvac),
            I_a_unvac_q3     = quantile(I_a_unvac, 0.75),
            I_a_unvac_max    = quantile(I_a_unvac, 0.975),
            
            I_a_vac_min    = quantile(I_a_vac, 0.025),
            I_a_vac_q1     = quantile(I_a_vac, 0.25),
            I_a_vac_median = quantile(I_a_vac, 0.5),
            I_a_vac_mean   = mean(I_a_vac),
            I_a_vac_q3     = quantile(I_a_vac, 0.75),
            I_a_vac_max    = quantile(I_a_vac, 0.975),
            
            I_p_min    = quantile(I_p, 0.025),
            I_p_q1     = quantile(I_p, 0.25),
            I_p_median = quantile(I_p, 0.5),
            I_p_mean   = mean(I_p),
            I_p_q3     = quantile(I_p, 0.75),
            I_p_max    = quantile(I_p, 0.975),
            
            I_unvac_p_min    = quantile(I_unvac_p, 0.025),
            I_unvac_p_q1     = quantile(I_unvac_p, 0.25),
            I_unvac_p_median = quantile(I_unvac_p, 0.5),
            I_unvac_p_mean   = mean(I_unvac_p),
            I_unvac_p_q3     = quantile(I_unvac_p, 0.75),
            I_unvac_p_max    = quantile(I_unvac_p, 0.975),
            
            I_vac_p_min    = quantile(I_vac_p, 0.025),
            I_vac_p_q1     = quantile(I_vac_p, 0.25),
            I_vac_p_median = quantile(I_vac_p, 0.5),
            I_vac_p_mean   = mean(I_vac_p),
            I_vac_p_q3     = quantile(I_vac_p, 0.75),
            I_vac_p_max    = quantile(I_vac_p, 0.975),
            
            I_c_min    = quantile(I_c, 0.025),
            I_c_q1     = quantile(I_c, 0.25),
            I_c_median = quantile(I_c, 0.5),
            I_c_mean   = mean(I_c),
            I_c_q3     = quantile(I_c, 0.75),
            I_c_max    = quantile(I_c, 0.975),
            
            I_c_unvac_min    = quantile(I_c_unvac, 0.025),
            I_c_unvac_q1     = quantile(I_c_unvac, 0.25),
            I_c_unvac_median = quantile(I_c_unvac, 0.5),
            I_c_unvac_mean   = mean(I_c_unvac),
            I_c_unvac_q3     = quantile(I_c_unvac, 0.75),
            I_c_unvac_max    = quantile(I_c_unvac, 0.975),
            
            I_c_vac_min    = quantile(I_c_vac, 0.025),
            I_c_vac_q1     = quantile(I_c_vac, 0.25),
            I_c_vac_median = quantile(I_c_vac, 0.5),
            I_c_vac_mean   = mean(I_c_vac),
            I_c_vac_q3     = quantile(I_c_vac, 0.75),
            I_c_vac_max    = quantile(I_c_vac, 0.975),
            
            I_total_min      = quantile(I_total, 0.025),
            I_total_q1       = quantile(I_total, 0.25),
            I_total_median   = median(I_total),
            I_total_mean     = quantile(I_total, 0.5),
            I_total_q3       = quantile(I_total, 0.75),
            I_total_max      = quantile(I_total, 0.975),
            
            I_total_unvac_min      = quantile(I_total_unvac, 0.025),
            I_total_unvac_q1       = quantile(I_total_unvac, 0.25),
            I_total_unvac_median   = median(I_total_unvac),
            I_total_unvac_mean     = quantile(I_total_unvac, 0.5),
            I_total_unvac_q3       = quantile(I_total_unvac, 0.75),
            I_total_unvac_max      = quantile(I_total_unvac, 0.975),
            
            I_total_vac_min      = quantile(I_total_vac, 0.025),
            I_total_vac_q1       = quantile(I_total_vac, 0.25),
            I_total_vac_median   = median(I_total_vac),
            I_total_vac_mean     = quantile(I_total_vac, 0.5),
            I_total_vac_q3       = quantile(I_total_vac, 0.75),
            I_total_vac_max      = quantile(I_total_vac, 0.975)
  )

# x <- 3

# write.csv(datafile_ci, file = './outcome/science/simulate_unvaccine.csv')

for (x in 1:3) {
  # x <- 1
  datafile_raw_plot <- outcome[,c('t', 'n', case_list[[x]])]
  names(datafile_raw_plot)[3:5] <- c('var', 'var_unvac', 'var_vac')
  
  datafile_ci_plot <- datafile_ci %>% 
    select(c(t, contains(all_of(case_list[[x]]))))
  names(datafile_ci_plot) <- str_remove_all(names(datafile_ci_plot), paste0(case_list[[x]][1], '_'))
  
  fig <- ggplot()+
    geom_line(data = datafile_raw_plot,
              mapping = aes(x = t, y = var, group = n),
              alpha = 0.03, color = fill_color[x],
              show.legend = F)+
    geom_line(data = datafile_ci_plot,
              mapping = aes(x = t, y =median,
                            group = 'Total',
                            linetype = 'Total'),
              size = 0.7, color = fill_color[x],
              show.legend = T)+
    geom_line(data = datafile_ci_plot,
              mapping = aes(x = t, y =unvac_median, 
                            group = 'Non-vaccinated',
                            linetype = 'Non-vaccinated'),
              size = 0.7, color = fill_color[x],
              show.legend = T)+
    geom_line(data = datafile_ci_plot,
              mapping = aes(x = t, y =vac_median, 
                            group = 'Fully vaccinated',
                            linetype = 'Fully vaccinated'),
              size = 0.7, color = fill_color[x],
              show.legend = T)+
    coord_cartesian(ylim = c(0, limit_value[x]),
                    xlim = c(0, 120))+
    scale_x_continuous(expand = c(0,0),
                       breaks = c(0, 30, 60, 90, 120))+
    scale_y_continuous(expand = c(0, NA),
                       labels = scientific_10)+
    scale_linetype_manual(values = c('dashed', 'dotted', 'solid'))+
    theme_bw(base_family = 'Helvetica')+
    theme(legend.text = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', 
                                     color = 'black'),
          legend.title = element_text(size = 12, hjust = 0, vjust = 0.5, face = 'bold'),
          legend.margin=margin(0,0,0,1),
          legend.box.margin=margin(0,0,0,-10),
          legend.box.background = element_rect(fill = "transparent", colour = 'transparent'),
          legend.background = element_rect(fill = "transparent", colour = 'transparent'),
          legend.position = "right",
          axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain',
                                     color = 'black'),
          axis.text.y = element_text(size = 10, hjust = 1, vjust = .5, face = 'plain', 
                                     color = 'black'),
          axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
          axis.title.y = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
          plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(5, 5, 5, 5))+
    guides(linetype = guide_legend(override.aes = list(color = 'black')))+
    labs(x = 'Days',
         # y = paste0(y_label[x], parse(text = "(10^5)")),
         y = expr(paste(!!y_label[x], " (", 10^5, ")", sep = '')),
         title = letters[x+3],
         linetype = 'Vaccinated status')
  
  assign(paste0('fig_temp_', x+3), fig)
}

# plot combind ------------------------------------------------------------

fig_temp_1$labels$x <- NULL
fig_temp_2$labels$x <- NULL
fig_temp_3$labels$x <- NULL

fig_temp_2$labels$y <- NULL
fig_temp_3$labels$y <- NULL
fig_temp_5$labels$y <- NULL
fig_temp_6$labels$y <- NULL

(fig_temp_1 + fig_temp_2 + fig_temp_3)/
  (fig_temp_4 + fig_temp_5 + fig_temp_6)+
  plot_layout(guides = 'collect')&
  theme(legend.position = 'bottom',
        plot.margin = margin(5, 10, 5, 5))

ggsave(filename = './outcome/publish/Figure 4.pdf', 
       device = cairo_pdf, height = 6, width = 10)

ggsave(filename = './outcome/publish/Figure 4.tiff', 
       height = 6, width = 10,
       dpi = 300)
