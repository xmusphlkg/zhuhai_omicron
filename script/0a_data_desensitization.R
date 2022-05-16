# packages ----------------------------------------------------------------

library(tidyverse)
library(openxlsx)


# data --------------------------------------------------------------------

## read data
datafile_info_raw <- read.xlsx('./data/raw/case info.xlsx')
datafile_cont_raw <- read.xlsx('./data/raw/case contact.xlsx')

## identify
datafile_ids <- datafile_cont_raw %>% 
  select(姓名, 序列号, 年龄, 性别) %>% 
  mutate(ids = paste(姓名, 性别,年龄, sep = '-')) %>% 
  filter(!duplicated(ids)) %>% 
  arrange(序列号) %>% 
  select(-c(序列号, ids))

x <- datafile_info_raw[!datafile_info_raw$姓名 %in% datafile_ids$姓名,c('姓名', '性别', '年龄')]

datafile_ids <- rbind(datafile_ids, x) %>% 
  mutate(ids = 1:nrow(.),
         年龄 = as.numeric(年龄)) %>% 
  rename(c(
    'name' = '姓名',
    'gender' = '性别',
    'age' = '年龄'
  )) %>% 
  mutate(gender = factor(gender, levels = c("男", "女"),
                         labels = c('Male', 'Female')),
         match = paste(name, gender, age, sep = '-')) %>% 
  select(name, match, ids)
  


# information datafile ----------------------------------------------------

swap_fun <- function(x){
  x <- str_replace_all(x, '/', '')
  x <- str_replace_all(x, '阴性', '45')
  return(as.numeric(x))
}

datafile_info_clean <- datafile_info_raw %>% 
  select(-接种时间) %>% 
  rename(c(
    'id' = '序号',
    'name' = '姓名',
    'gender' = '性别',
    'age' = '年龄',
    'type' = '病例类型',
    'type1' = '病例分型',
    'location1' = '家庭地址区',
    'location2' = '家庭地址乡镇/街道',
    'location' = '家庭所在村',
    'work' = '职业',
    'class' = '班级',
    'vaccine' = '疫苗剂次',
    'vaccinetype' = '接种类型',
    'vaccinelastdate' = '末次接种时间',
    'datepositive' = '初筛阳性采样日期',
    'datereported' = '第一次复采报告日期',
    'Nasopharyngeal_O' = '鼻咽拭子O', 
    'Nasopharyngeal_N' = '鼻咽拭子N', 
    'Throat_O' = '咽拭子O', 
    'Throat_N' = '咽拭子N', 
    'Anal_O' = '肛拭子O', 
    'Anal_N' = '肛拭子N',
    'dateonset' = '发病时间',
    'locationcluster' = '聚集地',
    'reporttype' = '发现途径',
    'dateexpose1' = '可能暴露时段（起）',
    'dateexpose2' = '可能暴露时段（止）',
    'singleexpose' = '是否单次暴露',
    'exposetype' = '暴露方式',
    'exposefreq' = '暴露频次',
    'exposenum' = '密接人数',
  )) %>% 
  arrange(id) %>% 
  mutate_at(vars(contains('date')), convertToDate) %>% 
  mutate_at(c("Nasopharyngeal_O", "Nasopharyngeal_N", "Throat_O", "Throat_N", "Anal_O", "Anal_N"), swap_fun) %>% 
  mutate(gender = factor(gender, levels = c("男", "女"),
                         labels = c('Male', 'Female')),
         id = 1:nrow(.),
         type = factor(type, levels = c("普通型", "轻型", "无症状感染者"),
                       labels = c('Moderate', 'Mild', 'Asymptomatic')),
         singleexpose = factor(singleexpose, levels = c('是', '否'),
                               labels = c(T, F)),
         reporttype = factor(reporttype, levels = c("区域筛查", "密接筛查", "发热门诊筛查（已判密接待管控）"),
                             labels = c('Mass Testing', 'Contact Tracing', 'Fever Clinic')),
         locationcluster = factor(locationcluster,
                                  exclude = NULL,
                                  levels = c("广生榕园幼儿园", "广生小学", "拜瑞口腔", NA),
                                  labels = c('Kindergarten', 'Primary school', 
                                             'Dental clinic', 'Undistinct')
         ),
         vaccine = ifelse(vaccinelastdate > dateexpose1 - 14,
                          vaccine - 1,
                          vaccine),
         vaccine = ifelse(is.na(vaccine), 0, vaccine)) %>% 
  # left_join(datafile_ids, by = c('name')) %>% 
  select(-c(type1, work, class, vaccinetype,
            exposetype, exposefreq, location, IgG, IgM, exposenum))

# write.csv(datafile_info_clean, file = '../data/raw/information_case.csv',
#           fileEncoding = 'UTF-8', row.names = F)

# contact file ------------------------------------------------------------

datafile_case_id <- datafile_info_clean %>% select(name, id)

datafile_cont <- datafile_cont_raw %>% 
  filter(姓名 %in% datafile_case_id$name) %>% 
  select(关联病例, 姓名, 接触频率) %>% 
  rename(c(
    'to' = 姓名,
    'from' = 关联病例,
    'freq' = 接触频率
  )) %>% 
  mutate(freq = factor(freq,
                       levels = c('一般', '偶尔', '经常'),
                       labels = LETTERS[1:3])) %>% 
  left_join(datafile_case_id, by = c('to' = 'name')) %>% 
  select(-to) %>% 
  rename(c('to' = 'id')) %>% 
  left_join(datafile_case_id, by = c('from' = 'name')) %>% 
  select(-from) %>% 
  rename(c('from' = 'id'))

# write.csv(datafile_cont, file = '../data/raw/contact_case.csv',
#           fileEncoding = 'UTF-8', row.names = F)
# write.xlsx(datafile_cont, file = '../data/raw/contact_case.xlsx', overwrite = T)

datafile_cont_all <- datafile_cont_raw %>% 
  rename(c(
    'id' = '序列号',
    'location1' = '现所在省',
    'location2' = '现所在市',
    'location3' = '珠海市内分区',
    'location4' = '现所在镇街',
    'casename' = '关联病例',
    'name' = '姓名',
    'gender' = '性别',
    'age' = '年龄',
    'work' = '职业',
    'contact' = '接触',
    'dateexpose' = '暴露日期',
    'datecontact' = '最后接触日期',
    'durationexpose' = 'days',
    'contactfreq' = '接触频率',
    'vaccine' = '新冠疫苗接种剂次',
    'vaccine1date' = '第一针接种日期',
    'vaccine1make' = '第一针厂家',
    'vaccine2date' = '第二针接种日期',
    'vaccine2make' = '第二针厂家',
    'vaccine3date' = '第三针接种日期',
    'vaccine3make' = '第三针厂家',
    'virustest' = '核酸检测次数（自2022年1月13日以来）',
    'virustesttime' = '末次采样日期'
  )) %>% 
  mutate_at(vars(contains('date')), convertToDate) %>% 
  mutate(age = as.numeric(age)) %>% 
  mutate(vaccinelastdate = apply(.[,c('vaccine1date', 'vaccine2date', 'vaccine3date')], 1, max, na.rm = T),
         vaccine = ifelse(vaccinelastdate > as.Date('2022/01/08')-14,
                          vaccine - 1,
                          vaccine),
         vaccine = ifelse(is.na(vaccine), 0, vaccine),
         vaccine = case_when(
           vaccine_type %in% c('CCC', 'PCC') ~ 2,
           vaccine_type %in% c('C', 'CC') ~ 1,
           vaccine_type == 'A' ~ 2,
           TRUE ~ vaccine
         ))

## add id

datafile_cont_all <- datafile_cont_all %>% 
  left_join(datafile_case_id, by = c('casename' = 'name')) %>% 
  rename('id_cases' = 'id.y') %>% 
  mutate(outcome = ifelse(name %in% datafile_case_id$name, 1, 0),
         gender = factor(gender, levels = c("男", "女"),
                         labels = c('Male', 'Female')),
         match = paste(name, gender, age, sep = '-')) %>% 
  left_join(datafile_ids, by = c('name' = 'name', 'match' = 'match')) %>% 
  select(-c(casename, name, match, vaccine1date, vaccine1make, vaccine2date,
            vaccine2make, vaccine3date, vaccine3make, work, id.x,
            location1, location2, location3, location4))

datafile_info_clean <- datafile_info_clean %>% select(-name)

save(datafile_info_clean, datafile_cont, datafile_cont_all,file = 'data/data_case.Rdata')

# summary -----------------------------------------------------------------

datafile_s1 <- datafile_info_clean %>% 
  mutate_at(vars(contains('date')), format, '%Y/%m/%d') %>% 
  mutate(exposedate = paste(dateexpose1, 'to', dateexpose2)) %>% 
  select(c("id", "gender", "age", "type", "locationcluster", 'vaccine',
           "vaccine_type", "datepositive", "datereported", "dateonset", 
           "exposedate"))
names(datafile_s1) <- c('ID', 'Gender', 'Age', 'Infections classification', 
                        'Cluster of infections', 'Vaccine dose','Vaccine type', 
                        'Date of positive test', 'Date of infections report',
                        'Date of onset', 'Date of possibly exposed')
write.xlsx(datafile_s1, file = 'outcome/science/Table S1.xlsx')

df_conts <- datafile_cont_all %>% 
     group_by(id_cases) %>% 
     count()

# BA2 ---------------------------------------------------------------------

datafile_cont_BA2 <- read.xlsx('./data/Omicron_BA2_contact.xlsx')
datafile_cont_BA2$vac1date <- convertToDate(datafile_cont_BA2$vac1date)
datafile_cont_BA2$vac2date <- convertToDate(datafile_cont_BA2$vac2date)
datafile_cont_BA2$vac3date <- convertToDate(datafile_cont_BA2$vac3date)
datafile_cont_BA2$vaccinelastdate <- convertToDate(datafile_cont_BA2$vaccinelastdate)
datafile_cont_BA2$lastexposedate <- convertToDate(datafile_cont_BA2$lastexposedate)

datafile_cont_BA2 <- datafile_cont_BA2 |> 
  mutate(seq = lastexposedate - vaccinelastdate,
         vaccine = if_else(seq <14, 
                           vaccine-1,
                           vaccine))

datafile_cont_all <- datafile_cont_all |> 
  mutate(vaccinelastdate = as.Date(vaccinelastdate),
         seq = datecontact - vaccinelastdate,
         seq = ifelse(is.na(seq),
                      as.Date(virustesttime) - vaccinelastdate,
                      seq))

save.image('./data/data_case.Rdata')
