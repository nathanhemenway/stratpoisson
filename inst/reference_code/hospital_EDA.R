#hospital EDA
# hospital <- read_csv('covid19_hospitalizations.csv')
#
# hospital %>%
#   rename('hospital_days' = 'Hospitalization time in days',
#          'age_group' = 'Age group',
#          'num_comorb' = 'Total number of comorbidities',
#          'risk_class' = 'Risk Classification Protocol',
#          'pulm_impair' = 'Pulmonary impairment') %>%
#   drop_na() -> hospital
#
#
# # Find average number of days in hospital for each age group
# hospital %>%
#   mutate(age_group = as.factor(age_group)) %>%
#   group_by(age_group) %>%
#   summarize(avg_days = mean(hospital_days, na.rm = TRUE),
#             num_subjects = n()) -> avg_days_age_group
#
# #Find average nummber of days in hospital for each comorbidity number, and number of subjects with each number of comorbidities
# hospital %>%
#   group_by(num_comorb) %>%
#   summarize(avg_days_comorb = mean(hospital_days, na.rm = TRUE),
#             num_subjects = n()) -> avg_days_comorb
#
#
# #find average number of days in hospital for each risk class
# hospital %>%
#   group_by(risk_class) %>%
#   summarize(avg_days_risk = mean(hospital_days, na.rm = TRUE),
#             num_subjects = n()) -> avg_days_risk
#
# #Find average number of days in hospital for each pulmonary impairment status
# hospital %>%
#   group_by(pulm_impair) %>%
#   summarize(avg_days_pulm = mean(hospital_days, na.rm = TRUE),
#             num_subjects = n()) -> avg_days_pulm
#
#
# avg_days_pulm %>%
#   kable(digits = 2) %>%
#   kable_styling(full_width = F)
#
#
# avg_days_risk %>%
#   kable(digits = 2) %>%
#   kable_styling(full_width = F)
#
# avg_days_comorb %>%
#   kable(digits = 2) %>%
#   kable_styling(full_width = F)
#
# avg_days_age_group %>%
#   kable(digits = 2) %>%
#   kable_styling(full_width = F)
