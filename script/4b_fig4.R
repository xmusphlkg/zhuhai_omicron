res_clog <- clogit(formula = outcome ~ lineage + strata(age),
                   data = datafile_cont) |> 
     summary() %>%
     .[["conf.int"]] |> 
     as.data.frame() |> 
     mutate(lineage = 'All') |> 
     select(-`exp(-coef)`) |> 
     rownames_to_column('var')


datafile_cont_smote <- SMOTE(lineage~., 
                             data = datafile_cont[,c('lineage', 'vaccine_g', 'outcome')],
                             perc.over = 600,perc.under=100)

res_log <- glm(formula = outcome ~ lineage + vaccine_g,
               data = datafile_cont_smote,
               family = binomial(link = "logit")) |>
     summary() %>%
     .[["coefficients"]] |>
     as.data.frame() |>
     mutate(lineage = 'All') |>
     rownames_to_column('var')