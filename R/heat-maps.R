library(ggplot2)
library(tidyverse)

load('data-saves/lassoversions.RData')
load('data-saves/Results.RData')


set.seed(101010)
doParallel::registerDoParallel(12)
future::plan(future::multisession, workers = 8)

strain.data <- haven::read_dta('~/Dropbox (Partners HealthCare)/Teams/drshah_shared/ARIC V5 LA strain/Master_ARIC_LA_Analysis.dta')
strain.vars <- names(strain.data)[-1]
strain.data <- strain.data %>% mutate_at(strain.vars, ~ scale(.x, T, T)[,1])

echo.vars <- c('lvedvi', 'lvedd', 'mlvwt',
               'lvrwt', 'lvmi', 'lvef',
               'lvesvi', 'lvesd', 'peakcs',
               'peakls', 'awave', 'eeplatratio',
               'eepsepratio', 'eprimelat',
               'eprimesep', 'ewave', 'lavi',
               'maxlaapd', 'pasp', 'trvel',
               'rvfac', 'tapsmv', strain.vars)

data.st <- strain.data %>%
  right_join(fifth.visit.echo.scaled, by = "id") %>%
  select(all_of(c(echo.vars, adjust)), starts_with("Seq"))

echo.data <- data.st %>% select(all_of(c(strain.vars, echo.vars)))
eigen.data <- bind_cols(eigenData %>% select(-fuptime, -hfdiag), echo.data)
mods <- grep("^ME", names(eigen.data), value = TRUE)

eigen.all.scaled <- lin.arry.aggr(mods,  eigen.data, echo.vars, adjust, labels)
echod.all.scaled <- lin.arry.aggr(proteins, data.st, echo.vars, adjust, labels)

echo.order <- echo.vars

col.order <- c('darkturquoise', 'salmon', 'darkgrey',
               'steelblue', 'pink', 'blue', 'darkorange',
               'magenta', 'cyan', 'lightgreen', 'royalblue',
               'darkgreen', 'darkred', 'brown', 'yellow',
               'orange', 'lightyellow', 'midnightblue',
               'grey60', 'black', 'red', 'turquoise',
               'greenyellow', 'tan', 'purple', 'skyblue',
               'saddlebrown', 'green', 'lightcyan', 'white', 'grey')

heatmapper <- eigen.all.scaled$all.res %>%
  select(term, outcome, estm) %>%
  mutate(outcome = factor(outcome, levels = echo.order),
         term = factor(substr(term, 3, 1000), levels = col.order),
         category = recode(outcome,
                           lvedvi = "Structure",
                           lvedd = "Structure",
                           mlvwt = "Structure",
                           lvrwt = "Structure",
                           lvmi = "Structure",
                           lvef = "Systolic Function",
                           lvesvi = "Systolic Function",
                           lvesd = "Systolic Function",
                           peakcs = "Systolic Function",
                           peakls = "Systolic Function",
                           awave = "Diastolic Function",
                           eeplatratio = "Diastolic Function",
                           eepsepratio = "Diastolic Function",
                           eprimelat = "Diastolic Function",
                           eprimesep = "Diastolic Function",
                           ewave = "Diastolic Function",
                           lavi = "Diastolic Function",
                           maxlaapd = "Diastolic Function",
                           pasp = "Diastolic Function",
                           trvel = "Diastolic Function",
                           rvfac = "RV Fxn",
                           tapsmv = "RV Fxn",
                           LAReservoir = "Strain",
                           LAConduit = "Strain",
                           LAContraction = "Strain",
                           LAIndex = "Strain",
                           maxvolume = "Strain",
                           minvolume = "Strain",
                           preavolume = "Strain"))


ggplot(data = heatmapper, aes(term, outcome, fill = estm))+
  geom_tile(color = "gray")+
  scale_fill_gradient2(low = "blue", high = "red",
                       midpoint = 0, space = "Lab",
                       name="Linear Regr. Coefficient") +
  theme_minimal()+
  labs(title = "Heatmap Echo Outcome vs WGCNA Module",
       x = "WGCNA Module",
       y = "Echo Outcome",
       fill = "Linear Regr. Estimate") +
  facet_grid(category ~ ., scales = "free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 8, hjust = 1))

ggsave(filename = "~/Desktop/modules.png", width = 12, height = 6, device='png', dpi=700)


echo.all.sign <- echod.all.scaled$all.res %>%
  filter(pval <= 0.05 / 4877)

freq <- echo.all.sign$term %>% table() %>% sort(decreasing = TRUE)
keeps <- names(freq)[freq > 10]

final.echo.sign <- echo.all.sign %>% filter(term %in% keeps) %>%
  select(term, outcome, estm) %>%
  mutate(outcome = factor(outcome, levels = echo.order),
         term1 = substr(labels[term, 'name'],1,25),
         term = factor(term, levels = keeps),
         category = recode(outcome,
                           lvedvi = "Structure",
                           lvedd = "Structure",
                           mlvwt = "Structure",
                           lvrwt = "Structure",
                           lvmi = "Structure",
                           lvef = "Systolic Function",
                           lvesvi = "Systolic Function",
                           lvesd = "Systolic Function",
                           peakcs = "Systolic Function",
                           peakls = "Systolic Function",
                           awave = "Diastolic Function",
                           eeplatratio = "Diastolic Function",
                           eepsepratio = "Diastolic Function",
                           eprimelat = "Diastolic Function",
                           eprimesep = "Diastolic Function",
                           ewave = "Diastolic Function",
                           lavi = "Diastolic Function",
                           maxlaapd = "Diastolic Function",
                           pasp = "Diastolic Function",
                           trvel = "Diastolic Function",
                           rvfac = "RV Fxn",
                           tapsmv = "RV Fxn",
                           LAReservoir = "Strain",
                           LAConduit = "Strain",
                           LAContraction = "Strain",
                           LAIndex = "Strain",
                           maxvolume = "Strain",
                           minvolume = "Strain",
                           preavolume = "Strain"))

ggplot(data = final.echo.sign, aes(term, outcome, fill = estm))+
  geom_tile(color = "gray")+
  scale_fill_gradient2(low = "blue", high = "red",
                       midpoint = 0, space = "Lab",
                       name="Linear Regr. Coefficient") +
  theme_minimal()+
  labs(title = "Heatmap Echo Outcome vs Significant Proteins in Linear Outcomes (>10x)",
       x = "Protein",
       y = "Echo Outcome",
       fill = "Linear Regr. Estimate") +
  facet_grid(category ~ ., scales = "free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 8, hjust = 1))+
  scale_x_discrete(labels = substr(labels[keeps, 'name'],1,25))

ggsave(filename = "~/Desktop/proteins.png", width = 15, height = 10, device='png', dpi=700)


cand <- c('SeqId_7655_11', 'SeqId_11178_21', 'SeqId_11109_56', 'SeqId_2602_2',
          'SeqId_16751_15', 'SeqId_4297_62', 'SeqId_2900_53', 'SeqId_4721_54',
          'SeqId_8368_102', 'SeqId_9526_3', 'SeqId_3152_57', 'SeqId_16851_50',
          'SeqId_5462_62', 'SeqId_15640_54', 'SeqId_9197_4', 'SeqId_13126_52',
          'SeqId_2201_17', 'SeqId_17337_1', 'SeqId_7211_2', 'SeqId_16614_27',
          'SeqId_13392_13', 'SeqId_12014_19', 'SeqId_6281_51', 'SeqId_2597_8',
          'SeqId_16308_14', 'SeqId_18925_24', 'SeqId_18930_28', 'SeqId_8464_31',
          'SeqId_8089_173', 'SeqId_13463_1', 'SeqId_3216_2', 'SeqId_4374_45',
          'SeqId_8556_5', 'SeqId_15565_102', 'SeqId_17761_2', 'SeqId_5661_15',
          'SeqId_9266_1', 'SeqId_3339_33', 'SeqId_5124_69', 'SeqId_4469_78',
          'SeqId_9468_8', 'SeqId_7994_41', 'SeqId_8978_30', 'SeqId_4413_3',
          'SeqId_2643_57', 'SeqId_6517_14', 'SeqId_18382_109', 'SeqId_18429_10',
          'SeqId_15492_1', 'SeqId_6462_12', 'SeqId_10748_216', 'SeqId_5599_88',
          'SeqId_13941_82', 'SeqId_8841_65', 'SeqId_15533_97', 'SeqId_9188_119',
          'SeqId_13118_5', 'SeqId_9715_15', 'SeqId_6464_40', 'SeqId_5701_81',
          'SeqId_12632_14', 'SeqId_7863_50', 'SeqId_16317_20', 'SeqId_13095_51',
          'SeqId_3348_49', 'SeqId_13597_20', 'SeqId_2765_4', 'SeqId_18821_9',
          'SeqId_12671_35', 'SeqId_3057_55', 'SeqId_9296_15', 'SeqId_3362_61',
          'SeqId_18882_7', 'SeqId_6103_70', 'SeqId_15422_12', 'SeqId_9126_171',
          'SeqId_16770_3', 'SeqId_5755_29', 'SeqId_8597_1', 'SeqId_12759_47',
          'SeqId_8252_2', 'SeqId_15315_64', 'SeqId_4959_2', 'SeqId_13133_73',
          'SeqId_6390_18', 'SeqId_16612_28', 'SeqId_6520_87')


echo.all.cand <- echod.all.scaled$all.res %>%
  filter(term %in% cand) %>%
  arrange(pval)


final.echo.cand <- echo.all.cand %>%
  select(term, outcome, estm) %>%
  mutate(outcome = factor(outcome, levels = echo.vars),
         term = factor(term, levels = echo.all.cand$term %>% unique()),
         category = recode(outcome,
                           lvedvi = "Structure",
                           lvedd = "Structure",
                           mlvwt = "Structure",
                           lvrwt = "Structure",
                           lvmi = "Structure",
                           lvef = "Systolic Function",
                           lvesvi = "Systolic Function",
                           lvesd = "Systolic Function",
                           peakcs = "Systolic Function",
                           peakls = "Systolic Function",
                           awave = "Diastolic Function",
                           eeplatratio = "Diastolic Function",
                           eepsepratio = "Diastolic Function",
                           eprimelat = "Diastolic Function",
                           eprimesep = "Diastolic Function",
                           ewave = "Diastolic Function",
                           lavi = "Diastolic Function",
                           maxlaapd = "Diastolic Function",
                           pasp = "Diastolic Function",
                           trvel = "Diastolic Function",
                           rvfac = "RV Fxn",
                           tapsmv = "RV Fxn",
                           LAReservoir = "Strain",
                           LAConduit = "Strain",
                           LAContraction = "Strain",
                           LAIndex = "Strain",
                           maxvolume = "Strain",
                           minvolume = "Strain",
                           preavolume = "Strain"))

ggplot(data = final.echo.cand, aes(term, outcome, fill = estm))+
  geom_tile(color = "gray")+
  scale_fill_gradient2(low = "blue", high = "red",
                       midpoint = 0, space = "Lab",
                       name="Linear Regr. Coefficient") +
  theme_minimal()+
  labs(title = "Heatmap Echo Outcome vs Candidate Proteins",
       x = "Protein",
       y = "Echo Outcome",
       fill = "Linear Regr. Estimate") +
  facet_grid(category ~ ., scales = "free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 8, hjust = 1))+
  scale_x_discrete(labels = substr(labels[cand, 'name'],1,25))

ggsave(filename = "~/Desktop/cand-proteins.png", width = 20, height = 15, device='png', dpi=700)
