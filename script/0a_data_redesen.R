

datafile_BA1 <- datafile_info_clean |> 
     left_join(datafile_cont[,2:3],
               by = c('id' = 'to'))

datafile_BA1 <- datafile_BA1 |> 
     select(id, from, gender, age, type,
            location1, location2,
            dateexpose1, dateexpose2,
            dateonset, datepositive,
            datereported, vaccine,
            vaccine_type, vaccinelastdate,
            Nasopharyngeal_O, Nasopharyngeal_N,
            Throat_O, Throat_N)

datafile_cont_BA1 <- datafile_cont_all |> 
     select(ids, age, gender, vaccine,
            vaccine_type, vaccine_mix, 
            outcome, vaccinelastdate, dateexpose,
            datecontact, id_cases)
