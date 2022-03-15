
##########################################################################
##########################################################################
### Correcting for reporting delay: an example using dengue fever data ###
##########################################################################
##########################################################################

### Developed by Dr. Lauren J. Beesley, PhD
### Contact: lvandervort@lanl.gov
### Last Updated: 3/10/2022

### Note: This script provides example code for handling delayed case reporting for an example dataset.
### This code uses data provided as part of the R package NobBS (Version 0.1.0) on CRAN.
### This analysis models the log(cases+0.1) using an ARMA(2,2) model.
### These methods are not specific to ARMA(2,2) and can be applied in many modeling settings. 

#######################################
### Read in Libraries and Functions ###
#######################################
library(ggplot2)
library(NobBS)
library(dplyr)
library(MCMCglmm)
detach(package:plyr) #ignore warnings, just ensuring no package clashes

### NEED TO SOURCE FUNCTIONS ###
current_directory = dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(current_directory,'/Example_Reporting_Delay_Functions.R'))


#############################
### Get Dengue Fever Data ###
#############################
denguedat$onset_week = as.character(denguedat$onset_week)
denguedat$report_week = as.character(denguedat$report_week)
DENGUE_DATA = NULL
DAT = data.frame( Weeks = names(table(denguedat$onset_week)), True_Counts = as.numeric(table(denguedat$onset_week)))
DAT$Years = format(as.Date(DAT$Weeks, format="%Y-%m-%d"),"%Y")
### Calculate the counts reported initially and up to 6 weeks after initial report for each week
for(d in 0:6){
  denguedat$prev_reported = as.numeric(as.Date(denguedat$onset_week, format="%Y-%m-%d")  >= (as.Date(denguedat$report_week, format="%Y-%m-%d") - 7*d))
  reported_week = aggregate(prev_reported~onset_week,  data = denguedat, FUN = sum)
  DENGUE_DATA = rbind(DENGUE_DATA, data.frame(DAT, Reported_Counts = reported_week$prev_reported, lag = d))
}
DENGUE_DATA$transformed_weeks = as.numeric(format(as.Date(DENGUE_DATA$Weeks, format="%Y-%m-%d"),"%W"))
unique_weeks = unique(sort(DENGUE_DATA$Weeks))
DENGUE_DATA = merge(DENGUE_DATA, data.frame(Weeks = unique_weeks, week_index = c(1:length(unique_weeks))), by = 'Weeks', all.x = T)
DENGUE_DATA = DENGUE_DATA[order(DENGUE_DATA$Years, DENGUE_DATA$Weeks, DENGUE_DATA$lag),]


######################################################
### Estimate pi(d) and f(N_ts(d); X) for each week ###
######################################################
YEARS = c(1992:2009) 
PREEST = replicate(18, list(NULL))
for(y in c(1:length(YEARS))){
  PREEST[[y]] = Get_PreEstimates(YEARS[y])
  print(y)
}


###########################
### Visualize Reporting ###
###########################

### Plot Observed Data ###
TO_PLOT = DENGUE_DATA[DENGUE_DATA$lag %in% c(0,1) & DENGUE_DATA$Years %in% c(2001:2009),]
ggplot(data = TO_PLOT)+
  geom_line(aes(x=transformed_weeks, y = Reported_Counts, color = 'red'), lwd = 1, data = TO_PLOT[TO_PLOT$lag == 0,])+
  geom_line(aes(x=transformed_weeks, y = Reported_Counts, color = 'darkgreen'), lwd = 1, data = TO_PLOT[TO_PLOT$lag == 1,])+
  geom_line(aes(x=transformed_weeks, y = True_Counts, color = 'black'), lwd = 1, data = TO_PLOT[TO_PLOT$lag == 0,])+
  theme_classic()+
  labs(title = 'Puerto Rico Dengue Fever Case Data')+
  facet_wrap(. ~ Years,shrink = TRUE, scales = "free_y")+
  ylab('Dengue fever cases')+xlab('Weeks')+
  scale_color_identity(name = '', breaks = c('red','darkgreen', 'black'), labels = c('Cases Initially Reported (Lag = 0)','Cases Reported after 1 Week (Lag = 1)', 'Validation Cases'), guide = 'legend')+
  theme(legend.position = 'top')+
  theme(legend.text=element_text(size=12))

### Plot estimated pi(d) [LAG ESTIMATION] ###
PILONG = bind_rows(simplify2array(PREEST)[1,], .id = 'Version' )
ggplot(PILONG[PILONG$Years %in% c(2001:2009),])+
  geom_line(aes(x=CURRENT_WEEK, y =cumul_pi, color = factor(lag)), size = 1.2)+
  scale_colour_viridis_d(name = 'Lag')+
  theme_classic()+
  xlab('Weeks')+ylab('Cumulative Reporting Probability')+
  labs(title = 'Estimated Cumulative Reporting Probabilities')+
  facet_wrap(Years~.)+
  theme(legend.position = 'top')+
  guides(color = guide_legend( nrow = 1, byrow = TRUE))+
  theme(legend.text=element_text(size=12))

### Compare Lag and Local pi(0) ###
PILONG = do.call("rbind", simplify2array(PREEST)['CUMUL_PI',])
PILONG$cumul_pi_local = do.call("rbind", simplify2array(PREEST)['CUMUL_PI_LOCAL',])$cumul_pi
PILONG = merge(PILONG, data.frame(CURRENT_WEEK = DENGUE_DATA$transformed_weeks, lag = DENGUE_DATA$lag, Years = DENGUE_DATA$Years, cumul_pi_obs = DENGUE_DATA$Reported_Counts/DENGUE_DATA$True_Counts), by = c('CURRENT_WEEK','lag', 'Years'), all.x = T, all.y = F)
PILONG = PILONG[PILONG$lag ==0  ,]
year_limits = aggregate(week_index~Years, data = PILONG, FUN = function(x){cbind(min(x), max(x))})
year_limits = data.frame(year_limits[,1], year_limits[,2])
names(year_limits) = c('Years', 'min_week', 'max_week')
ggplot(PILONG)+
  geom_vline(xintercept =unique(c(year_limits$min_week, year_limits$max_week)), color = 'gray90', linetype = 2)+
  geom_line(aes(x=week_index, y =cumul_pi_obs, color = 'gray'),size = 1.2)+
  geom_line(aes(x=week_index, y =cumul_pi_local, color = 'blue'),size = 1.2)+
  geom_line(aes(x=week_index, y =cumul_pi, color = 'red'),size = 1.2)+
  scale_x_continuous(breaks =year_limits$min_week + 25, labels = year_limits$Years)+
  scale_color_identity(name = '', breaks=c('gray', 'blue', 'red'),labels = c('Observed Data', 'Local (prev. 6 weeks)', 'Lag (prev. 2 yrs)'), guide = 'legend')+
  theme_classic()+
  guides(size = 'none')+
  xlab('Season')+ylab(expression("Estimated inverse reporting factor, " ~ pi[ts](0)))+
  labs(title = 'Estimated Lag-0 Reporting Probabilities for Dengue Fever')+
  theme(legend.position = 'top', legend.text= element_text(size = 12))


#############################
### Specify Forecast Year ###
#############################
### To obtain forecasts for different years, change y:
### Corresponds to year in 1992:2009 to make predictions for (y=18 -> 2009)
y = 18 

###################
###################
### Forecasting ### 
###################
###################
DENGUE_SUB = DENGUE_DATA[(DENGUE_DATA$Years == YEARS[y] & DENGUE_DATA$transformed_weeks <=53) | DENGUE_DATA$Years < YEARS[y] , ]
DENGUE_SUB = DENGUE_SUB[order(DENGUE_SUB$Years, DENGUE_SUB$transformed_weeks, DENGUE_SUB$lag),]

### Analyze Validation Data ###
FORECASTS_Validated = ARIMA_Vanilla_Dengue(DENGUE_SUB, OUTCOME = 'True_Counts', PREEST, y)

### Analyze Observed Data without Correction ###
FORECASTS_Observed = ARIMA_Vanilla_Dengue(DENGUE_SUB, OUTCOME = 'Reported_Counts', PREEST, y)

### Analyze Observed Data after Excluding Recent-Reported Weeks ###
FORECASTS_Exclude1 = ARIMA_Exclude_Dengue(DENGUE_SUB, OUTCOME = 'Reported_Counts', PREEST, y, EXCLUDE_THRESHOLD = 1)
FORECASTS_Exclude2 = ARIMA_Exclude_Dengue(DENGUE_SUB, OUTCOME = 'Reported_Counts', PREEST, y, EXCLUDE_THRESHOLD = 2)
FORECASTS_Exclude3 = ARIMA_Exclude_Dengue(DENGUE_SUB, OUTCOME = 'Reported_Counts', PREEST, y, EXCLUDE_THRESHOLD = 3)

### Rescaling: Lag ###
OUTCOME_FUNCTION = function(DENGUE_CURRENT, PREEST, CURRENT_WEEK, y){
  CUMUL_PI = data.frame(PREEST[[y]][['CUMUL_PI']])
  CUMUL_PI$cumul_pi[CUMUL_PI$lag >= 6] = 1
  CUMUL_PI = CUMUL_PI[CUMUL_PI$CURRENT_WEEK == CURRENT_WEEK,]
  CUMUL_PI = CUMUL_PI[order(CUMUL_PI$lag),]
  results = round(DENGUE_CURRENT$Reported_Counts/unlist(CUMUL_PI[DENGUE_CURRENT$lag+1,'cumul_pi']),0)
  results[DENGUE_CURRENT$Years < max(DENGUE_CURRENT$Years, na.rm=T)] = DENGUE_CURRENT$True_Counts[DENGUE_CURRENT$Years < max(DENGUE_CURRENT$Years, na.rm=T)]
  return(results)
}
FORECASTS_RescalingLag = ARIMA_Vanilla_Dengue(DENGUE_SUB, OUTCOME = 'Counts_corr', PREEST, y, OUTCOME_FUNC = OUTCOME_FUNCTION)

### Rescaling: Model ###
OUTCOME_FUNCTION = function(DENGUE_CURRENT, PREEST, CURRENT_WEEK, y){
  fit_valid = PREEST[[y]][['fit_valid']]
  results = round(predict(fit_valid[[CURRENT_WEEK]], newdata =DENGUE_CURRENT, type = 'response' ),0)
  results[DENGUE_CURRENT$lag >= 6] = DENGUE_CURRENT$Reported_Counts[DENGUE_CURRENT$lag >= 6]
  results[DENGUE_CURRENT$Years < max(DENGUE_CURRENT$Years, na.rm=T)] = DENGUE_CURRENT$True_Counts[DENGUE_CURRENT$Years < max(DENGUE_CURRENT$Years, na.rm=T)]
  return(results)
}
FORECASTS_RescalingModel = ARIMA_Vanilla_Dengue(DENGUE_SUB, OUTCOME = 'Counts_corr', PREEST, y, OUTCOME_FUNC = OUTCOME_FUNCTION)

### Rescaling: Local ###
OUTCOME_FUNCTION = function(DENGUE_CURRENT, PREEST, CURRENT_WEEK, y){
  CUMUL_PI_LOCAL = data.frame(PREEST[[y]][['CUMUL_PI_LOCAL']])
  CUMUL_PI_LOCAL$cumul_pi[CUMUL_PI_LOCAL$lag >= 6] = 1
  CUMUL_PI_LOCAL = CUMUL_PI_LOCAL[CUMUL_PI_LOCAL$CURRENT_WEEK == CURRENT_WEEK,]
  CUMUL_PI_LOCAL = CUMUL_PI_LOCAL[order(CUMUL_PI_LOCAL$lag),]
  results = round(DENGUE_CURRENT$Reported_Counts/unlist(CUMUL_PI_LOCAL[DENGUE_CURRENT$lag+1,'cumul_pi']),0)
  results[DENGUE_CURRENT$Years < max(DENGUE_CURRENT$Years, na.rm=T)] = DENGUE_CURRENT$True_Counts[DENGUE_CURRENT$Years < max(DENGUE_CURRENT$Years, na.rm=T)]
  return(results)
}
FORECASTS_RescalingLocal = ARIMA_Vanilla_Dengue(DENGUE_SUB, OUTCOME = 'Counts_corr', PREEST, y, OUTCOME_FUNC = OUTCOME_FUNCTION)

### Offset: Lag ###
FORECASTS_OffsetLag = ARIMA_Offset_Dengue(DENGUE_SUB, OUTCOME = 'Reported_Counts', PREEST, y, pi_type = 'Lag')

### Offset: Model ###
FORECASTS_OffsetModel = ARIMA_Offset_Dengue(DENGUE_SUB, OUTCOME = 'Reported_Counts', PREEST, y, pi_type = 'Model')

### Offset: Local ###
FORECASTS_OffsetLocal = ARIMA_Offset_Dengue(DENGUE_SUB, OUTCOME = 'Reported_Counts', PREEST, y, pi_type = 'Local')

### Imputation: Lag ### (MAY TAKE A SEC)
FORECASTS_ImputationLag = ARIMA_Imp_Dengue(DENGUE_SUB, PREEST, y, pi_type = 'Lag')

### Imputation: Model ### (MAY TAKE A SEC)
FORECASTS_ImputationModel = ARIMA_Imp_Dengue(DENGUE_SUB, PREEST, y, pi_type = 'Model')

### Imputation: Model ### (MAY TAKE A SEC)
FORECASTS_ImputationLocal = ARIMA_Imp_Dengue(DENGUE_SUB, PREEST, y, pi_type = 'Local')

####################################
### Calculate Summary Statistics ### 
####################################

### Get Validation Values ###
validation_values = replicate(18, list(NULL))
for(year in 1:length(YEARS)){
  DENGUE_SUB_YEAR = DENGUE_DATA[DENGUE_DATA$Years == YEARS[year], ]
  DENGUE_SUB_YEAR = DENGUE_SUB_YEAR[order(DENGUE_SUB_YEAR$transformed_weeks, DENGUE_SUB_YEAR$lag),]
  DENGUE_SUB_YEAR$current_week = DENGUE_SUB_YEAR$transformed_weeks + DENGUE_SUB_YEAR$lag
  validation_values[[year]] = DENGUE_SUB_YEAR$True_Counts[!duplicated(DENGUE_SUB_YEAR$transformed_weeks)]
}

### Apply Rubin's Combining Rules for Multiple Imputation Methods
FORECASTS_ImputationLagRUBIN = ARIMA_Rubin_Est(FORECASTS_ImputationLag)
FORECASTS_ImputationModelRUBIN = ARIMA_Rubin_Est(FORECASTS_ImputationModel)
FORECASTS_ImputationLocalRUBIN = ARIMA_Rubin_Est(FORECASTS_ImputationLocal)

### Calculate Performance Summaries
RESULTS_Validated = Summarize_ARIMA(FORECASTS_Validated, validation_values[[y]])
RESULTS_Observed = Summarize_ARIMA(FORECASTS_Observed, validation_values[[y]])
RESULTS_Exclude1 = Summarize_ARIMA(FORECASTS_Exclude1, validation_values[[y]])
RESULTS_Exclude2 = Summarize_ARIMA(FORECASTS_Exclude2, validation_values[[y]])
RESULTS_Exclude3 = Summarize_ARIMA(FORECASTS_Exclude3, validation_values[[y]])
RESULTS_RescalingLag = Summarize_ARIMA(FORECASTS_RescalingLag, validation_values[[y]])
RESULTS_RescalingModel = Summarize_ARIMA(FORECASTS_RescalingModel, validation_values[[y]])
RESULTS_RescalingLocal = Summarize_ARIMA(FORECASTS_RescalingLocal, validation_values[[y]])
RESULTS_OffsetLag = Summarize_ARIMA(FORECASTS_OffsetLag, validation_values[[y]])
RESULTS_OffsetModel = Summarize_ARIMA(FORECASTS_OffsetModel, validation_values[[y]])
RESULTS_OffsetLocal = Summarize_ARIMA(FORECASTS_OffsetLocal, validation_values[[y]])
RESULTS_ImputationLag = Summarize_ARIMA(FORECASTS_ImputationLagRUBIN[,names(FORECASTS_Validated)], validation_values[[y]])
RESULTS_ImputationModel = Summarize_ARIMA(FORECASTS_ImputationModelRUBIN[,names(FORECASTS_Validated)], validation_values[[y]])
RESULTS_ImputationLocal = Summarize_ARIMA(FORECASTS_ImputationLocalRUBIN[,names(FORECASTS_Validated)], validation_values[[y]])

### Performance for Equal-Weight Ensemble (excluding validated data analysis from ensemble)
### - Prediction intervals for each method were calculated on the log(cases+0.1) scale using normal distribution approximations
### - We calculate ensemble predictions using the empirical distribution of stacked draws from each method's normal prediction distribution
TO_AVG = rbind(FORECASTS_Observed, FORECASTS_Exclude1, FORECASTS_Exclude2,FORECASTS_Exclude3,FORECASTS_RescalingLag, 
               FORECASTS_RescalingModel, FORECASTS_RescalingLocal, FORECASTS_OffsetLag, FORECASTS_OffsetModel, 
               FORECASTS_OffsetLocal, FORECASTS_ImputationLagRUBIN[,names(FORECASTS_Observed)], 
               FORECASTS_ImputationModelRUBIN[,names(FORECASTS_Observed)], 
               FORECASTS_ImputationLocalRUBIN[,names(FORECASTS_Observed)])
FORECASTS_Ensemble = Ensemble_Estimates(TO_AVG, ndraws = 100)
RESULTS_Ensemble= Summarize_ARIMA(FORECASTS_Ensemble, validation_values[[y]])


RESULTS_ALL = rbind(data.frame(RESULTS_Validated, transformed_weeks = c(1:50), Method = "Validation Data (Truth)", Method_Num = 1),
                    data.frame(RESULTS_Observed, transformed_weeks = c(1:50), Method = "Observed Data (Naive)", Method_Num = 2),
                    data.frame(RESULTS_Exclude1, transformed_weeks = c(1:50), Method = "Exclusion Method: 1 Week", Method_Num = 3),
                    data.frame(RESULTS_Exclude2, transformed_weeks = c(1:50), Method = "Exclusion Method: 2 Weeks", Method_Num = 4),
                    data.frame(RESULTS_Exclude3, transformed_weeks = c(1:50), Method = "Exclusion Method: 3 Weeks", Method_Num = 5),
                    data.frame(RESULTS_RescalingLag, transformed_weeks = c(1:50), Method = "Rescaling: Lag", Method_Num = 6),
                    data.frame(RESULTS_RescalingModel, transformed_weeks = c(1:50), Method = "Rescaling: Model", Method_Num = 7),
                    data.frame(RESULTS_RescalingLocal, transformed_weeks = c(1:50), Method = "Rescaling: Local", Method_Num = 8),
                    data.frame(RESULTS_OffsetLag, transformed_weeks = c(1:50), Method = "Offset: Lag", Method_Num = 9),
                    data.frame(RESULTS_OffsetModel, transformed_weeks = c(1:50), Method = "Offset: Model", Method_Num = 10),
                    data.frame(RESULTS_OffsetLocal, transformed_weeks = c(1:50), Method = "Offset: Local", Method_Num = 11),
                    data.frame(RESULTS_ImputationLag, transformed_weeks = c(1:50), Method = "Imputation: Lag", Method_Num = 12),
                    data.frame(RESULTS_ImputationModel, transformed_weeks = c(1:50), Method = "Imputation: Model", Method_Num = 13),
                    data.frame(RESULTS_ImputationLocal, transformed_weeks = c(1:50), Method = "Imputation: Local", Method_Num = 14),
                    data.frame(RESULTS_Ensemble, transformed_weeks = c(1:50), Method = "Equal Weight Ensemble", Method_Num = 15))
METHODS = c("Validation Data (Truth)", "Observed Data (Naive)", 
            "Exclusion Method: 1 Week", "Exclusion Method: 2 Weeks", "Exclusion Method: 3 Weeks",
            "Rescaling: Lag",'Rescaling: Model','Rescaling: Local',"Offset: Lag","Offset: Model","Offset: Local",
            "Imputation: Lag", "Imputation: Model","Imputation: Local","Equal Weight Ensemble")


#####################################################
### Aggregate performance results across 50 weeks ###
#####################################################
TEMP = aggregate(cbind(nowcast_bias,forecast_bias,WIS)~Method, FUN = function(x){mean(abs(x), na.rm=T)}, data = RESULTS_ALL, drop = FALSE, na.action = na.pass)
(TEMP = TEMP[order(match(TEMP$Method,METHODS)),])


#####################################################
### Plot relative forecast accuracy for each week ###
#####################################################
# Accuracy defined as abs(forecast bias)^{-1}
counts_max = max(DENGUE_DATA[DENGUE_DATA$lag == 0 & DENGUE_DATA$Years == YEARS[y],'True_Counts'])
METHODS_ACCURACY = c("Observed Data (Naive)", "Exclusion Method: 1 Week", "Exclusion Method: 2 Weeks", "Exclusion Method: 3 Weeks",
                     "Rescaling: Local","Offset: Local","Imputation: Local")
TO_PLOT = RESULTS_ALL[RESULTS_ALL$Method %in% METHODS_ACCURACY,] %>% group_by(transformed_weeks) %>% mutate(accuracy_forecast = (1/abs(forecast_bias))/sum((1/abs(forecast_bias)), na.rm=T))
TO_PLOT$Method = factor(TO_PLOT$Method, levels = METHODS_ACCURACY)
ggplot(TO_PLOT[TO_PLOT$transformed_weeks < 50,])+
  geom_bar(aes(x=transformed_weeks, y = accuracy_forecast, fill = Method ),position="stack", stat="identity" )+
  theme_classic()+theme(legend.position = 'top')+
  scale_fill_brewer(name = '', limits = METHODS_ACCURACY, palette = 'Spectral')+
  labs(title = 'Relative Accuracy of 1 Week Forecast')+
  xlab('Weeks')+ylab('Relative Accuracy')+
  geom_point(aes(x=1, y = accuracy_forecast, group = Method),position=position_fill(vjust = 0.5), stat='identity', data =TO_PLOT[TO_PLOT$transformed_weeks == 1,], pch = 21, size = 5, fill = 'white')+
  geom_point(aes(x=1, y = accuracy_forecast, shape=Method, fill = Method),position=position_fill(vjust = 0.5), stat='identity', data =TO_PLOT[TO_PLOT$transformed_weeks == 1,], size = 3)+
  scale_shape_manual(name = "", breaks = METHODS_ACCURACY, values = as.character(c(1:length(METHODS_ACCURACY))))+
  coord_cartesian(ylim=c(0,1))+
  geom_line(aes(x=transformed_weeks, y = True_Counts/max(True_Counts)), lwd = 1.5, data = DENGUE_DATA[DENGUE_DATA$lag == 0 & DENGUE_DATA$Years == YEARS[y] & DENGUE_DATA$transformed_weeks < 50,])+
  scale_y_continuous( sec.axis = sec_axis( trans=~(.*counts_max), name="Validation Cases"))


########################################################
### Plot proportion of weeks each method ranked best ###
########################################################
METHODS_RANK = c("Observed Data (Naive)", "Exclusion Method: 1 Week", "Exclusion Method: 2 Weeks", "Exclusion Method: 3 Weeks",
                     "Rescaling: Local","Offset: Local","Imputation: Local")
RESULTS_LONG = rbind(data.frame(est = RESULTS_ALL$nowcast_bias, RESULTS_ALL[,c('transformed_weeks', 'Method', 'Method_Num')], Metric = 'Nowcast Bias'),
                data.frame(est = RESULTS_ALL$forecast_bias, RESULTS_ALL[,c('transformed_weeks', 'Method', 'Method_Num')], Metric = 'Forecast Bias'),
                data.frame(est = RESULTS_ALL$WIS, RESULTS_ALL[,c('transformed_weeks', 'Method', 'Method_Num')], Metric = 'Forecast WIS'))
RESULTS_LONG = RESULTS_LONG[order(RESULTS_LONG$Metric, RESULTS_LONG$transformed_weeks, RESULTS_LONG$Method_Num),]
RANKINGS = aggregate(abs(est)~transformed_weeks + Metric, FUN =function(x){cbind(which.min(x), length(x))}, data = RESULTS_LONG[RESULTS_LONG$Method %in% METHODS_RANK,])
RANKINGS = data.frame(RANKINGS[,1:2], data.frame(RANKINGS[,3]))
RANKINGS = RANKINGS[RANKINGS$X2 == length(METHODS_RANK),] #only calculate rankings if all methods have estimate
AGG_RANK = aggregate(X1~Metric, FUN =function(x){cbind(as.vector(prop.table(table(factor(x, levels = c(1:length(METHODS_RANK)))))))}, data = RANKINGS)
AGG_RANK = as.matrix(data.frame(AGG_RANK[,2]))
TO_PLOT = rbind(data.frame(Metric = 'Forecast Bias', est =as.vector(AGG_RANK[1,]), Method = METHODS_RANK ),
                data.frame(Metric = 'Forecast WIS', est =as.vector(AGG_RANK[2,]), Method = METHODS_RANK ),
                data.frame(Metric = 'Nowcast Bias', est =as.vector(AGG_RANK[3,]), Method = METHODS_RANK ))
TO_PLOT$Method = factor(TO_PLOT$Method, levels = METHODS_RANK)
TO_PLOT$Metric = factor(TO_PLOT$Metric, levels = c('Nowcast Bias', 'Forecast Bias', 'Forecast WIS'))
ggplot(TO_PLOT)+
  geom_bar(aes(x=Metric, y = est, fill = Method ),position="stack", stat="identity" )+
  theme_classic()+theme(legend.position = 'top')+
  scale_fill_brewer(name = '', limits = METHODS_RANK, palette = 'Spectral')+
  labs(title = 'Relative Accuracy of 1 Week Forecast')+
  xlab('Metric')+ylab('Proportion of Weeks Method Performed Best')+
  geom_point(aes(x=0.8, y = est, group = Method),position=position_fill(vjust = 0.5), stat='identity', data =TO_PLOT[TO_PLOT$Metric == 'Nowcast Bias',], pch = 21, size = 5, fill = 'white')+
  geom_point(aes(x=0.8, y = est, shape=Method, fill = Method),position=position_fill(vjust = 0.5), stat='identity', data =TO_PLOT[TO_PLOT$Metric == 'Nowcast Bias',], size = 3)+
  scale_shape_manual(name = "", breaks = METHODS_RANK, values = as.character(c(1:length(METHODS_RANK))))+
  geom_text(aes(x=Metric, y = est, label=paste0(round(100*est,1),'%'), fill = Method),
            position=position_fill(vjust = 0.5), stat='identity', size = 3)+
  coord_cartesian(ylim=c(0,1))




