################################################
################################################
### Functions for Modeling Dengue Fever Data ###
################################################
################################################

### Developed by Dr. Lauren J. Beesley, PhD
### Contact: lvandervort@lanl.gov
### Last Updated: 3/10/2022

### Note: This script provides example functions for handling delayed case reporting for ARMA(2,2) models. 
### To avoid issues with 0 case counts, the ARMA models are applied to the log(cases + 0.1) rather than log(cases)


########################################
### Estimate pi(d) and f(N_ts(d); X) ###
########################################
Get_PreEstimates = function(YEAR){
  indices = sort(unique(DENGUE_DATA$week_index[DENGUE_DATA$Years == YEAR]))
  CUMUL_PI = CUMUL_PI_LOCAL = NULL
  fit_valid = replicate(length(indices), list(NULL))
  for(i in c(1:length(indices))){
    ### (1) Use previous 2 years of data to estimate pi(d), excluding recent 6 weeks [LAG ESTIMATION] ###
    temp = DENGUE_DATA[DENGUE_DATA$week_index < as.numeric(indices[i]-6) & DENGUE_DATA$week_index > as.numeric(indices[i]-52*2)& !is.na(DENGUE_DATA$Reported_Counts),]
    temp <- temp %>% group_by(lag) %>% mutate(cumul_pi = sum(Reported_Counts, na.rm=T)/sum(True_Counts, na.rm=T))
    CUMUL_PI = rbind(CUMUL_PI, data.frame(data.frame(temp[!duplicated(temp$lag)  ,c('lag', 'cumul_pi')]), week_index = indices[i], CURRENT_WEEK = i, Years = YEAR))
    ### (2) Use previous 2 years of data to estimate pi(d), excluding recent 6 weeks [MODEL ESTIMATION] ###
    fit_valid[[i]] = glm(True_Counts~splines::ns(transformed_weeks,df=3) + 
                           as.numeric(lag == 0) + as.numeric(lag == 1) + as.numeric(lag == 2) + 
                           as.numeric(lag == 3) + as.numeric(lag == 4) + as.numeric(lag == 5), 
                         offset = as.numeric(log(Reported_Counts+0.001)), data = temp, family = 'poisson')
    ### (3) Use recent K=6 weeks' data to estimate pi(d) [LOCAL ESTIMATION]
    ### A correction to zero-counts is added to improve estimator stability and avoid divide-by-zero errors
    temp = DENGUE_DATA[DENGUE_DATA$week_index >= as.numeric(indices[i]-6) & DENGUE_DATA$week_index < as.numeric(indices[i])& !is.na(DENGUE_DATA$Reported_Counts),]
    temp = temp[temp$week_index+temp$lag <= indices[i],]
    temp$Reported_Counts[temp$Reported_Counts == 0] = 0.5 
    most_recent = temp[order(temp$week_index, temp$lag, decreasing = TRUE),]
    most_recent = most_recent[!duplicated(most_recent$week_index), ]
    temp = merge(temp, data.frame(True_Counts_local = most_recent$Reported_Counts, week_index = most_recent[,c('week_index')]), by = c('week_index'), all.x = T, all.y = F)
    temp <- temp %>% group_by(lag) %>% mutate(cumul_pi = sum(Reported_Counts, na.rm=T)/sum(True_Counts_local, na.rm=T))
    CUMUL_TEMP = temp[!duplicated(temp$lag)  ,c('lag', 'cumul_pi')]
    CUMUL_PI_LOCAL = rbind(CUMUL_PI_LOCAL, data.frame(data.frame(CUMUL_TEMP), week_index = indices[i], CURRENT_WEEK = i, Years = YEAR))
  }
  CUMUL_PI$cumul_pi[CUMUL_PI$lag >= 6] = 1 ### Assume pi(d) = 1 for all d >= 6
  CUMUL_PI_LOCAL$cumul_pi[CUMUL_PI_LOCAL$lag >= 6] = 1 ### Assume pi(d) = 1 for all d >= 6
  return(list(CUMUL_PI = CUMUL_PI, fit_valid = fit_valid, CUMUL_PI_LOCAL = CUMUL_PI_LOCAL))
}




############################
### Standard ARIMA Model ###
############################
ARIMA_Vanilla_Dengue = function(DENGUE_SUB, OUTCOME, PREEST, y, OUTCOME_FUNC = NULL){
  ESTIMATES = NULL
  for(CURRENT_WEEK in c(1:50)){
    DENGUE_CURRENT = rbind(DENGUE_SUB[(DENGUE_SUB$transformed_weeks + DENGUE_SUB$lag) <= CURRENT_WEEK & as.numeric(DENGUE_SUB$Years) == max(as.numeric(DENGUE_SUB$Years)),],
                           DENGUE_SUB[ DENGUE_SUB$lag == 6 &  as.numeric(DENGUE_SUB$Years) < max(as.numeric(DENGUE_SUB$Years)),])
    DENGUE_CURRENT = DENGUE_CURRENT[order(DENGUE_CURRENT$Years, DENGUE_CURRENT$transformed_weeks, DENGUE_CURRENT$lag, decreasing = TRUE),]
    DENGUE_CURRENT$most_recent = as.numeric(!duplicated(DENGUE_CURRENT[,c('Years', 'transformed_weeks')] ))
    DENGUE_CURRENT = DENGUE_CURRENT[DENGUE_CURRENT$most_recent == 1,]
    DENGUE_CURRENT = DENGUE_CURRENT[order(DENGUE_CURRENT$Years, DENGUE_CURRENT$transformed_weeks, DENGUE_CURRENT$lag),]
    ### Scaling and other methods have trouble with 0 case counts
    DENGUE_CURRENT$Reported_Counts[DENGUE_CURRENT$Reported_Counts == 0 & DENGUE_CURRENT$Years == max(DENGUE_CURRENT$Years)] = 0.5 
    
    ### Obtain rescaling estimate estimate of validation predictions, if applicable
    if(!is.null(OUTCOME_FUNC)){
      DENGUE_CURRENT[,OUTCOME] = OUTCOME_FUNC(DENGUE_CURRENT, PREEST, CURRENT_WEEK, y)
    }
    ### Fit ARMA Model (ARIMA(p,0,q) = ARMA(p,q)) and Forecast
    ts_input = ts(log(DENGUE_CURRENT[,OUTCOME]+0.1), start = 1, end = length(DENGUE_CURRENT[,1]))
    fit = try(arima(ts_input , order = c(2,0,2), seasonal = list(order = c(0,0,0))), silent = T)
    
    if(class(fit) != 'try-error'){
      nowcasting = exp(ts_input[length(ts_input)] - fit$residuals[length(fit$residuals)])-0.1
      forecasting = forecast::forecast(fit,h=1, level = c(10,20,30,40,50,60,70,80,90,95,98), lambda = NULL)
      ESTIMATES = rbind(ESTIMATES,data.frame(forecasting = as.numeric(exp(forecasting$mean))-0.1, 
                                             forecasting_lower = as.numeric(exp(forecasting$lower))-0.1, 
                                             forecasting_upper = as.numeric(exp(forecasting$upper))-0.1,
                                             nowcasting = nowcasting,
                                             LEVELS = c(10,20,30,40,50,60,70,80,90,95,98),
                                             CURRENT_WEEK = CURRENT_WEEK))
    }else{
      ESTIMATES = rbind(ESTIMATES,data.frame(forecasting = rep(NA,11), 
                                             forecasting_lower = rep(NA,11), 
                                             forecasting_upper = rep(NA,11), 
                                             nowcasting = rep(NA,11),
                                             LEVELS = c(10,20,30,40,50,60,70,80,90,95,98),
                                             CURRENT_WEEK = CURRENT_WEEK))
    }
  }
  return(ESTIMATES)
}


#########################################
### Excluding Recently-Reported Weeks ###
#########################################
ARIMA_Exclude_Dengue = function(DENGUE_SUB, OUTCOME, PREEST, y, EXCLUDE_THRESHOLD = 1){
  ESTIMATES = NULL
  for(CURRENT_WEEK in c(1:50)){
    DENGUE_CURRENT = rbind(DENGUE_SUB[(DENGUE_SUB$transformed_weeks + DENGUE_SUB$lag) <= CURRENT_WEEK & as.numeric(DENGUE_SUB$Years) == max(as.numeric(DENGUE_SUB$Years)),],
                           DENGUE_SUB[ DENGUE_SUB$lag == 6 &  as.numeric(DENGUE_SUB$Years) < max(as.numeric(DENGUE_SUB$Years)),])
    DENGUE_CURRENT = DENGUE_CURRENT[order(DENGUE_CURRENT$Years, DENGUE_CURRENT$transformed_weeks, DENGUE_CURRENT$lag, decreasing = TRUE),]
    DENGUE_CURRENT$most_recent = as.numeric(!duplicated(DENGUE_CURRENT[,c('Years', 'transformed_weeks')] ))
    DENGUE_CURRENT = DENGUE_CURRENT[DENGUE_CURRENT$most_recent == 1,]
    DENGUE_CURRENT = DENGUE_CURRENT[order(DENGUE_CURRENT$Years, DENGUE_CURRENT$transformed_weeks, DENGUE_CURRENT$lag),]
    ### Scaling and other methods have trouble with 0 case counts
    DENGUE_CURRENT$Reported_Counts[DENGUE_CURRENT$Reported_Counts == 0 & DENGUE_CURRENT$Years == max(DENGUE_CURRENT$Years)] = 0.5 
    ### Exclude most recently reported EXCLUDE_THRESHOLD weeks
    DENGUE_CURRENT = DENGUE_CURRENT[DENGUE_CURRENT$lag >= EXCLUDE_THRESHOLD,]
    
    ### Fit ARMA Model (ARIMA(p,0,q) = ARMA(p,q)) and Forecast
    ts_input = ts(log(DENGUE_CURRENT[,OUTCOME]+0.1), start = 1, end = length(DENGUE_CURRENT[,1]))
    fit = try(arima(ts_input , order = c(2,0,2), seasonal = list(order = c(0,0,0))), silent = T)
    if(class(fit) != 'try-error'){
      forecasting = forecast::forecast(fit,h=EXCLUDE_THRESHOLD+1, level = c(10,20,30,40,50,60,70,80,90,95,98), lambda = NULL) 
      ESTIMATES = rbind(ESTIMATES,data.frame(forecasting = as.numeric(exp(forecasting$mean[EXCLUDE_THRESHOLD+1]))-0.1, 
                                             forecasting_lower = as.numeric(exp(forecasting$lower[EXCLUDE_THRESHOLD+1,]))-0.1, 
                                             forecasting_upper = as.numeric(exp(forecasting$upper[EXCLUDE_THRESHOLD+1,]))-0.1,
                                             nowcasting = as.numeric(exp(forecasting$mean[EXCLUDE_THRESHOLD]))-0.1,
                                             LEVELS = c(10,20,30,40,50,60,70,80,90,95,98),
                                             CURRENT_WEEK = CURRENT_WEEK))
    }else{
      ESTIMATES = rbind(ESTIMATES,data.frame(forecasting = rep(NA,11), 
                                             forecasting_lower = rep(NA,11), 
                                             forecasting_upper = rep(NA,11), 
                                             nowcasting = rep(NA,11),
                                             LEVELS = c(10,20,30,40,50,60,70,80,90,95,98),
                                             CURRENT_WEEK = CURRENT_WEEK))
    }
  }
  return(ESTIMATES)
}

########################################
### Fit ARIMA with Mean Model Offset ###
########################################
ARIMA_Offset_Dengue = function(DENGUE_SUB, OUTCOME, PREEST, y, pi_type = 'Lag'){
  ESTIMATES = NULL
  for(CURRENT_WEEK in c(1:50)){
    DENGUE_CURRENT = rbind(DENGUE_SUB[(DENGUE_SUB$transformed_weeks + DENGUE_SUB$lag) <= CURRENT_WEEK & DENGUE_SUB$Years == max(as.numeric(DENGUE_SUB$Years)),],
                           DENGUE_SUB[ DENGUE_SUB$lag == 6 &  as.numeric(DENGUE_SUB$Years) < max(as.numeric(DENGUE_SUB$Years)),])
    DENGUE_CURRENT = DENGUE_CURRENT[order(DENGUE_CURRENT$Years, DENGUE_CURRENT$transformed_weeks, DENGUE_CURRENT$lag, decreasing = TRUE),]
    DENGUE_CURRENT$most_recent = as.numeric(!duplicated(DENGUE_CURRENT[,c('Years', 'transformed_weeks')] ))
    DENGUE_CURRENT = DENGUE_CURRENT[DENGUE_CURRENT$most_recent == 1,]
    DENGUE_CURRENT = DENGUE_CURRENT[order(DENGUE_CURRENT$Years, DENGUE_CURRENT$transformed_weeks, DENGUE_CURRENT$lag),]
    ### Scaling and other methods have trouble with 0 case counts
    DENGUE_CURRENT$Reported_Counts[DENGUE_CURRENT$Reported_Counts == 0 & DENGUE_CURRENT$Years == max(DENGUE_CURRENT$Years)] = 0.5
    ### Obtain pi_ts(d)
    if(pi_type == 'Lag'){
      CUMUL_PI = data.frame(PREEST[[y]][['CUMUL_PI']])
      CUMUL_PI$cumul_pi[CUMUL_PI$lag >= 6] = 1
      CUMUL_PI = CUMUL_PI[CUMUL_PI$CURRENT_WEEK == CURRENT_WEEK,]
      CUMUL_PI = CUMUL_PI[order(CUMUL_PI$lag),]
      DENGUE_CURRENT$CUMUL_PI = unlist(CUMUL_PI[DENGUE_CURRENT$lag+1,'cumul_pi'])
    }else if(pi_type == 'Model'){
      fit_valid = PREEST[[y]][['fit_valid']]
      MODEL_PREDICTED = predict(fit_valid[[CURRENT_WEEK]], newdata =DENGUE_CURRENT, type = 'response' )
      DENGUE_CURRENT$CUMUL_PI = DENGUE_CURRENT$Reported_Counts/MODEL_PREDICTED
      ### Restrict inverse reporting factors to be equal to or less than 1
      DENGUE_CURRENT$CUMUL_PI[DENGUE_CURRENT$CUMUL_PI>1 | DENGUE_CURRENT$lag >= 6] = 1
    }else if(pi_type == 'Local'){
      CUMUL_PI_LOCAL = data.frame(PREEST[[y]][['CUMUL_PI_LOCAL']])
      CUMUL_PI_LOCAL$cumul_pi[CUMUL_PI_LOCAL$lag >= 6] = 1
      CUMUL_PI_LOCAL = CUMUL_PI_LOCAL[CUMUL_PI_LOCAL$CURRENT_WEEK == CURRENT_WEEK,]
      CUMUL_PI_LOCAL = CUMUL_PI_LOCAL[order(CUMUL_PI_LOCAL$lag),]
      DENGUE_CURRENT$CUMUL_PI = unlist(CUMUL_PI_LOCAL[DENGUE_CURRENT$lag+1,'cumul_pi'])
    }else{
      stop('invalid valid for pi_type argument')
    }
    ### Calculate fixed offset
    DENGUE_CURRENT$OFFSET = log(DENGUE_CURRENT$CUMUL_PI)
    ### Fit ARMA Model (ARIMA(p,0,q) = ARMA(p,q)) and Forecast
    ts_input = ts(log(DENGUE_CURRENT[,OUTCOME]+0.1), start = 1, end = length(DENGUE_CURRENT[,1]))
    fit = try(arima(ts_input , order = c(2,0,2), seasonal = list(order = c(0,0,0)), 
                    xreg = matrix(DENGUE_CURRENT$OFFSET), fixed = c(NA,NA,NA,NA,NA,1)), silent = T)
    if(class(fit)[1] != 'try-error'){
      nowcasting = exp(ts_input[length(ts_input)] - fit$residuals[length(fit$residuals)] - DENGUE_CURRENT$OFFSET[length(ts_input)])-0.1
      quantile_vals = qnorm(p =(100-c(10,20,30,40,50,60,70,80,90,95,98))*0.01/2, lower.tail = F )
      A = stats:::predict.Arima(fit, n.ahead = 4,newxreg = matrix(c(0,0,0,0)))
      B = data.frame(pred = as.numeric(A$pred), se = as.numeric(A$se))
      forecasting = data.frame(mean = rep(B$pred[1],11), lower = B$pred[1]-quantile_vals*B$se[1], upper = B$pred[1]+quantile_vals*B$se[1])
      ESTIMATES = rbind(ESTIMATES,data.frame(forecasting = as.numeric(exp(forecasting$mean))-0.1, 
                                             forecasting_lower = as.numeric(exp(forecasting$lower))-0.1, 
                                             forecasting_upper = as.numeric(exp(forecasting$upper))-0.1, 
                                             nowcasting = nowcasting,
                                             LEVELS = c(10,20,30,40,50,60,70,80,90,95,98),
                                             CURRENT_WEEK = CURRENT_WEEK))
    }else{
      ESTIMATES = rbind(ESTIMATES,data.frame(forecasting = rep(NA,11), 
                                             forecasting_lower = rep(NA,11), 
                                             forecasting_upper = rep(NA,11), 
                                             nowcasting = rep(NA,11),
                                             LEVELS = c(10,20,30,40,50,60,70,80,90,95,98),
                                             CURRENT_WEEK = CURRENT_WEEK))
    }
  }
  return(ESTIMATES)
}


##########################################
### Fit ARIMA with Multiple Imputation ###
##########################################
ARIMA_Imp_Dengue = function(DENGUE_SUB, PREEST, y, pi_type = 'Lag'){
  ESTIMATES = NULL
  for(CURRENT_WEEK in c(1:50)){
    DENGUE_CURRENT = rbind(DENGUE_SUB[(DENGUE_SUB$transformed_weeks + DENGUE_SUB$lag) <= CURRENT_WEEK & DENGUE_SUB$Years == max(as.numeric(DENGUE_SUB$Years)),],
                           DENGUE_SUB[ DENGUE_SUB$lag == 6 &  as.numeric(DENGUE_SUB$Years) < max(as.numeric(DENGUE_SUB$Years)),])
    DENGUE_CURRENT = DENGUE_CURRENT[order(DENGUE_CURRENT$Years, DENGUE_CURRENT$transformed_weeks, DENGUE_CURRENT$lag, decreasing = TRUE),]
    DENGUE_CURRENT$most_recent = as.numeric(!duplicated(DENGUE_CURRENT[,c('Years', 'transformed_weeks')] ))
    DENGUE_CURRENT = DENGUE_CURRENT[DENGUE_CURRENT$most_recent == 1,]
    DENGUE_CURRENT = DENGUE_CURRENT[order(DENGUE_CURRENT$Years, DENGUE_CURRENT$transformed_weeks, DENGUE_CURRENT$lag),]
    
    ### Scaling and other methods have trouble with 0 case counts
    DENGUE_CURRENT$Reported_Counts[DENGUE_CURRENT$Reported_Counts == 0 & DENGUE_CURRENT$Years == max(DENGUE_CURRENT$Years,na.rm=T)] = 0.5
    ### Obtain pi_ts(d)
    if(pi_type == 'Lag'){
      CUMUL_PI = data.frame(PREEST[[y]][['CUMUL_PI']])
      CUMUL_PI$cumul_pi[CUMUL_PI$lag >= 6] = 1
      CUMUL_PI = CUMUL_PI[CUMUL_PI$CURRENT_WEEK == CURRENT_WEEK,]
      CUMUL_PI = CUMUL_PI[order(CUMUL_PI$lag),]
      DENGUE_CURRENT$CUMUL_PI = unlist(CUMUL_PI[DENGUE_CURRENT$lag+1,'cumul_pi'])
    }else if(pi_type == 'Model'){
      fit_valid = PREEST[[y]][['fit_valid']]
      MODEL_PREDICTED = predict(fit_valid[[CURRENT_WEEK]], newdata =DENGUE_CURRENT, type = 'response' )
      DENGUE_CURRENT$CUMUL_PI = DENGUE_CURRENT$Reported_Counts/MODEL_PREDICTED
      ### Restrict inverse reporting factors to be equal to or less than 1
      DENGUE_CURRENT$CUMUL_PI[DENGUE_CURRENT$CUMUL_PI>1 | DENGUE_CURRENT$lag >= 6] = 1
    }else if(pi_type == 'Local'){
      CUMUL_PI_LOCAL = data.frame(PREEST[[y]][['CUMUL_PI_LOCAL']])
      CUMUL_PI_LOCAL$cumul_pi[CUMUL_PI_LOCAL$lag >= 6] = 1
      CUMUL_PI_LOCAL = CUMUL_PI_LOCAL[CUMUL_PI_LOCAL$CURRENT_WEEK == CURRENT_WEEK,]
      CUMUL_PI_LOCAL = CUMUL_PI_LOCAL[order(CUMUL_PI_LOCAL$lag),]
      DENGUE_CURRENT$CUMUL_PI = unlist(CUMUL_PI_LOCAL[DENGUE_CURRENT$lag+1,'cumul_pi'])
    }else{
      stop('invalid valid for pi_type argument')
    }
    ### Obtain imputation distribution mean and standard deviation
    DENGUE_CURRENT$mu = as.numeric(DENGUE_CURRENT$Reported_Counts/DENGUE_CURRENT$CUMUL_PI)
    DENGUE_CURRENT$sd = sqrt(DENGUE_CURRENT$mu*abs(1-DENGUE_CURRENT$CUMUL_PI)/DENGUE_CURRENT$CUMUL_PI)
    ### Truncated distributions throw errors with 0 sd. Change to very small sd. 
    DENGUE_CURRENT$sd[DENGUE_CURRENT$sd == 0] = 0.001
    
    for(impnum in c(1:10)){
      ### Obtain imputed version of validation data using normal distribution approximation
      DENGUE_CURRENT$Counts_imp = apply(cbind(data.frame(mu = DENGUE_CURRENT$mu,sd = DENGUE_CURRENT$sd, trunc_left = DENGUE_CURRENT$Reported_Counts) ),1,
                                        FUN = function(x){ rtnorm(n=1, mean = x[1],sd = x[2],lower=x[3])})
      DENGUE_CURRENT$Counts_imp[DENGUE_CURRENT$Years < max(as.numeric(DENGUE_CURRENT$Years))] = DENGUE_CURRENT$True_Counts[DENGUE_CURRENT$Years < max(as.numeric(DENGUE_CURRENT$Years))]
      ### Fit ARMA Model (ARIMA(p,0,q) = ARMA(p,q)) and Forecast
      ts_input = ts(log(DENGUE_CURRENT$Counts_imp+0.1), start = 1, end = length(DENGUE_CURRENT[,1]))
      fit = try(arima(ts_input , order = c(2,0,2), seasonal = list(order = c(0,0,0))), silent = T)
      if(class(fit) != 'try-error'){
        nowcasting = exp(ts_input[length(ts_input)] - fit$residuals[length(fit$residuals)])-0.1
        forecasting = forecast::forecast(fit,h=1, level = c(10,20,30,40,50,60,70,80,90,95,98), lambda = NULL) 
        ESTIMATES = rbind(ESTIMATES,data.frame(forecasting = as.numeric(exp(forecasting$mean))-0.1, 
                                               forecasting_lower = as.numeric(exp(forecasting$lower))-0.1, 
                                               forecasting_upper = as.numeric(exp(forecasting$upper))-0.1, 
                                               nowcasting = nowcasting,
                                               LEVELS = c(10,20,30,40,50,60,70,80,90,95,98),
                                               CURRENT_WEEK = CURRENT_WEEK, 
                                               impnum = impnum))
      }else{
        ESTIMATES = rbind(ESTIMATES,data.frame(forecasting = rep(NA,11), 
                                               forecasting_lower = rep(NA,11), 
                                               forecasting_upper = rep(NA,11), 
                                               nowcasting = rep(NA,11),
                                               LEVELS = c(10,20,30,40,50,60,70,80,90,95,98),
                                               CURRENT_WEEK = CURRENT_WEEK, 
                                               impnum = impnum))
      }
    }
  }
  return(ESTIMATES)
}



#############################
### Aggregation Functions ###
#############################
Interval_score = function(QUANTILES, alpha, y){
  interval = c(QUANTILES[QUANTILES$alpha == alpha & QUANTILES$type == 'lower','quantile' ], QUANTILES[QUANTILES$alpha == alpha & QUANTILES$type == 'upper','quantile'])
  IS = (interval[2]-interval[1])+ ((2/alpha)*(interval[1]-y)*as.numeric(y<interval[1]))  + ((2/alpha)*(y-interval[2])*as.numeric(y>interval[2]))
  return(IS)
}

### Calculate Performance Statistics
Summarize_ARIMA = function(ESTIMATES, validation_values){
  NUM_WEEKS = length(unique(ESTIMATES$CURRENT_WEEK))
  RESULTS = replicate(6,list(matrix(NA, nrow = NUM_WEEKS, ncol = 1)))
  names(RESULTS) = c('nowcast_bias', 'forecast_bias', 'WIS', 'interval_width', 'lower_forecast', 'upper_forecast')
  for(CURRENT_WEEK in c(1:NUM_WEEKS)){
    DATA_CURRENT = ESTIMATES[ESTIMATES$CURRENT_WEEK == CURRENT_WEEK,]
    if(length(DATA_CURRENT[,1])!=0){
      if(CURRENT_WEEK<NUM_WEEKS){
        DATA_CURRENT$log_SE = (log(DATA_CURRENT$forecasting_upper+0.1) - log(DATA_CURRENT$forecasting+0.1))/qnorm(p = (1-0.01*DATA_CURRENT$LEVELS)/2, mean = 0, sd = 1, lower = F)
      }
      ### BIAS
      RESULTS[['nowcast_bias']][CURRENT_WEEK,1] = DATA_CURRENT$nowcasting[1]-  validation_values[CURRENT_WEEK ]
      if(CURRENT_WEEK<NUM_WEEKS){
        RESULTS[['forecast_bias']][CURRENT_WEEK,1] = DATA_CURRENT$forecasting[1]-  validation_values[CURRENT_WEEK+1 ]
        upper = exp(log(DATA_CURRENT$forecasting[1]+0.1) + 1.96*DATA_CURRENT$log_SE[1])-0.1
        lower = exp(log(DATA_CURRENT$forecasting[1]+0.1) - 1.96*DATA_CURRENT$log_SE[1])-0.1
        RESULTS[['interval_width']][CURRENT_WEEK,1] =as.numeric( upper - lower)
        RESULTS[['lower_forecast']][CURRENT_WEEK,1] = as.numeric(lower)
        RESULTS[['upper_forecast']][CURRENT_WEEK,1] = as.numeric(upper)
      }
      ### WEIGHTED INTERVAL SCORE
      if(CURRENT_WEEK<NUM_WEEKS){
        QUANTILES = data.frame(quantile = log(DATA_CURRENT$forecasting_lower+0.1),
                               alpha = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.10, 0.05, 0.02),
                               type = rep('lower',11))
        QUANTILES = rbind(QUANTILES, data.frame(quantile = log(DATA_CURRENT$forecasting_upper+0.1),
                                                alpha = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.10, 0.05, 0.02),
                                                type = rep('upper',11)))
        QUANTILES = rbind(QUANTILES, data.frame(quantile = log(DATA_CURRENT$forecasting[1]+0.1),
                                                alpha = 1, type = 'median'))
        median_val = QUANTILES[QUANTILES$type == 'median','quantile' ]
        
        alphas = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.10, 0.05, 0.02)
        IS = apply(cbind(alphas,log(validation_values[CURRENT_WEEK+1]+0.1)), 1, FUN = function(x,QUANTILES){Interval_score(QUANTILES,alpha = x[1], y=x[2])}, QUANTILES)
        RESULTS[['WIS']][CURRENT_WEEK,1] = (1/(length(alphas)+0.5))* (0.5*abs(log(validation_values[CURRENT_WEEK+1]+0.1) - median_val) + sum( (alphas/2)*IS))
      }
    }
  }
  return(RESULTS)
}



### Apply Rubin's Multiple Imputation Combining Rules to Forecasts (on Log Scale)
ARIMA_Rubin_Est = function(ESTIMATES){
  ESTIMATES$forecasting = log(ESTIMATES$forecasting+0.1)
  ESTIMATES$forecasting_lower = log(ESTIMATES$forecasting_lower+0.1)
  ESTIMATES$forecasting_upper = log(ESTIMATES$forecasting_upper+0.1)
  ESTIMATES_SHORT = subset(ESTIMATES[ESTIMATES$impnum == 1,], select = -c(forecasting, impnum ,forecasting_lower,forecasting_upper, nowcasting ))
  ESTIMATES$log_SE = (ESTIMATES$forecasting_upper - ESTIMATES$forecasting)/qnorm(p = (1-0.01*ESTIMATES$LEVELS)/2, mean = 0, sd = 1, lower = F)

  MEANS = aggregate(cbind(forecasting, nowcasting)~CURRENT_WEEK+LEVELS, data = ESTIMATES, FUN = mean, na.rm=T, drop = FALSE)
  ESTIMATES_SHORT = merge(ESTIMATES_SHORT, MEANS, by = c('CURRENT_WEEK', 'LEVELS'), all.x = T)
  
  VAR_WITHIN = aggregate(cbind(log_SE^2)~CURRENT_WEEK+LEVELS, data = ESTIMATES, FUN = mean, na.rm=T, drop = FALSE)
  VAR_BETWEEN = aggregate(cbind(forecasting)~CURRENT_WEEK+LEVELS, data = ESTIMATES, FUN = var, na.rm=T, drop = FALSE)
  SE_TOTAL = VAR_WITHIN
  SE_TOTAL[,3] = sqrt(VAR_WITHIN[,3] + (1+0.1)*VAR_BETWEEN[,3])
  names(SE_TOTAL) = c('CURRENT_WEEK', 'LEVELS', 'log_SE')
  ESTIMATES_SHORT = merge(ESTIMATES_SHORT, SE_TOTAL, by = c('CURRENT_WEEK', 'LEVELS'), all.x = T)
  ESTIMATES_SHORT$forecasting_lower = ESTIMATES_SHORT$forecasting-qnorm(as.numeric( (1-0.01*ESTIMATES_SHORT$LEVELS)/2), mean = 0, sd = 1, lower = F)*ESTIMATES_SHORT$log_SE
  ESTIMATES_SHORT$forecasting_upper = ESTIMATES_SHORT$forecasting+qnorm(as.numeric( (1-0.01*ESTIMATES_SHORT$LEVELS)/2), mean = 0, sd = 1, lower = F)*ESTIMATES_SHORT$log_SE
  
  ESTIMATES_SHORT$forecasting = exp(ESTIMATES_SHORT$forecasting)-0.1
  ESTIMATES_SHORT$forecasting_upper = exp(ESTIMATES_SHORT$forecasting_upper)-0.1
  ESTIMATES_SHORT$forecasting_lower = exp(ESTIMATES_SHORT$forecasting_lower)-0.1
  ESTIMATES_SHORT$interval_width = ESTIMATES_SHORT$forecasting_upper - ESTIMATES_SHORT$forecasting_lower
  return(ESTIMATES_SHORT)
}


### Construct Ensemble Nowcasts/Forecasts using the empirical distribution of stacked draws from each method's normal prediction distribution
Ensemble_Estimates = function(STACK, ndraws = 100){
  ### Prediction intervals for each method were calculated on the log(cases+0.1) scale using normal distribution approximations
  STACK = STACK[STACK$LEVELS == 95,]
  STACK$log_forecast = log(STACK$forecasting+0.1)
  STACK$log_SE = (log(STACK$forecasting_upper+0.1) - log(STACK$forecasting+0.1))/qnorm(p = (1-0.01*STACK$LEVELS)/2, mean = 0, sd = 1, lower = F)
  my_func = function(x,ndraws){
    stand_dev = x[2]
    stand_dev[is.na(stand_dev)] = 0.001
    return(rnorm(n=ndraws, mean = x[1], sd = stand_dev))
  }
  DRAWS = apply(STACK[,c('log_forecast', 'log_SE')],1, FUN = my_func, ndraws=ndraws)
  LONG_DRAWS = data.frame(forecasting = exp(as.vector(DRAWS))-0.1,
             draw_num = rep(c(1:ndraws), length(STACK[,1])), 
             CURRENT_WEEK = rep(STACK$CURRENT_WEEK, each = ndraws), 
             nowcasting = rep(STACK$nowcasting, each = ndraws))
  ### Gets averages and quantiles of nowcasts and forecasts across methods for each week
  nowcast_mean = aggregate(nowcasting~CURRENT_WEEK, data = LONG_DRAWS,FUN = function(x){quantile(x,probs = c(0.5))})
  forecasting_mean = aggregate(forecasting~CURRENT_WEEK, data = LONG_DRAWS,FUN = function(x){quantile(x,probs = c(0.5))})
  forecasting_lower = aggregate(forecasting~CURRENT_WEEK, data = LONG_DRAWS,FUN = function(x){quantile(x,probs = c(0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0.025,0.01))})
  forecasting_upper = aggregate(forecasting~CURRENT_WEEK, data = LONG_DRAWS,FUN = function(x){quantile(x,probs = c(0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99))})
  ### Output dataset with same structure
  M = length(nowcast_mean[,1])
  OUTPUT = data.frame(forecasting = rep(forecasting_mean$forecasting, 11),
                      forecasting_lower = as.vector(forecasting_lower[,2]), 
                      forecasting_upper = as.vector(forecasting_upper[,2]), 
                      nowcasting = rep(nowcast_mean$nowcasting, 11),
                      LEVELS = rep(c(10,20,30,40,50,60,70,80,90,95,98), each = M),
                      CURRENT_WEEK = rep(c(1:M),11))
  return(OUTPUT)
}






