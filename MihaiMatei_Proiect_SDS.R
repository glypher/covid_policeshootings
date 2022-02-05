cat("\f")
options(scipen=999)
options(stringsAsFactors = FALSE)

load_packages <- function(package_names) {
  list.of.packages <- as.vector(package_names)
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  sapply(list.of.packages, require, character.only = TRUE)
}

load_packages(c('rstudioapi', 'tidyverse', 'plotmo',
                'astsa', 'xts', 'gridExtra', 'stats4', 'forecast', 'stringi'))


wait<-function() {invisible(readline(prompt="Press [enter] to continue"))}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


data_c19 <- as.data.frame(read.csv('./data/covid_19_clean_complete.csv'))
summary(data_c19)

data_c19 <- data_c19 %>% filter(Country.Region == 'Romania') %>%
  dplyr::select(-c(Province.State, Country.Region, Lat, Long, WHO.Region)) %>%
  filter(Confirmed > 0) %>%
  mutate(Date=as.Date(Date, "%Y-%m-%d"))

data_c19 <- data_c19 %>% mutate(Daily_New_Confirmed = Confirmed - lag(Confirmed),
                   Daily_New_Deaths = Deaths - lag(Deaths),
                   Daily_New_Recovered = Recovered - lag(Recovered),
                   Daily_New_Active = Active - lag(Active))
data_c19[is.na(data_c19)] <- 0
nrow(data_c19)
head(data_c19, 2)
tail(data_c19, 2)

data_c19_toplot <- data_c19 %>% dplyr::select(c(Date, Confirmed, Deaths, Recovered, Active)) %>% gather(case_type, number,-Date)

ggplot(data_c19_toplot, aes(x=Date)) +
  geom_line(aes(y=number, color=case_type)) + 
  geom_point(aes(y=number, color=case_type),size=1) +
  scale_color_manual(values=c("Confirmed"="darkblue","Deaths"="red", "Recovered"="green", "Active"="yellow")) +
  labs(title="Covid impact in Romania", x="Date", y="Number of Cases", colour='Cases') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

data_c19_toplot<-data_c19 %>% dplyr::select(c(Date, Daily_New_Confirmed, Daily_New_Deaths, Daily_New_Recovered)) %>%
  gather(case_type, number,-Date)

data_c19_toplot

ggplot(data_c19_toplot, aes(x=Date)) +
  geom_line(aes(y=number, color=case_type)) + 
  geom_point(aes(y=number, color=case_type),size=1) +
  scale_color_manual(values=c("Daily_New_Confirmed"="darkblue","Daily_New_Deaths"="red",
                              "Daily_New_Recovered"="green")) +
  labs(title="Covid impact in Romania", x="Date", y="Number of New Cases per day", colour='New Cases') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


ts_to_dataframe <- function(ts, log=FALSE) {
  tt <- as.vector(time(ts))
  df <- data.frame(Cases=matrix(ts), Date=tt)
  if (log==TRUE) {
    df$Cases <- exp(df$Cases)
  }
  return (df)
}

analize_ts <- function(ts, max.lag = 30) {
  par(mfrow=c(1,2))
  acf <- acf(ts, lag.max=max.lag, na.action = na.pass, plot=TRUE)
  pacf <- pacf(ts, lag.max=max.lag, na.action = na.pass, plot=TRUE)
  wait()
}

diagnostic <- function(t) {
  tsdiag(t)
  wait()
  rs<-t$residuals
  stdres <- rs/sqrt(t$sigma2)
  qqnorm(stdres,  main="Normal Q-Q Plot of Std Residuals")
  wait()
}

grid_search<- function (ts, AR, MA, D, performance_function) {
  best_model <- NULL
  best_score <- Inf
  for (ar in AR) {
    for (ma in MA) {
      for (d in D) {
        tryCatch({
          model <- Arima(ts, order=c(ar, d, ma))
          score <- performance_function(model, ar, d, ma)
          if (score < best_score) {
            best_model = model
          }
        }, error = function(e) {
          print(sprintf("ARIMA[p=%d, d=%d, q=%d]: ERROR=%s", ar, d, ma, e))
        }, warning = function(e) {
          print(sprintf("ARIMA[p=%d, d=%d, q=%d]: WARN=%s", ar, d, ma, e))
        })
      }
    }
  }
  return (best_model)
}


do_time_series <- function(time_series, data1, log=FALSE, preditions=20, confidence=95) {
  analize_ts(time_series)
  model <- grid_search(time_series, rep(20:25), c(0, 1 ,2), c(0, 1),
                       function(model, ar, d, ma) {
                         print(sprintf("ARIMA[p=%d, d=%d, q=%d]: AIC=%f", ar, d, ma, model$aic))
                         return(model$aic)
                       })
  #model <- model <- Arima(time_series, order=c(15, 1, 0))
  #wait()
  summary(model)
  diagnostic(model)

  forecast <- forecast::forecast(model, h=preditions, level=confidence)
  forecast1 <- predict(model, n.ahead=preditions)
  
  
  data_df <- time_series %>% ts_to_dataframe(log=log) %>% mutate(Date=as.Date(Date))

  first_date = head(time(time_series), n=1) - 1
  
  
  forecast_mean1 <- forecast1$pred %>% ts_to_dataframe(log=log) %>% mutate(Date=as.Date(Date, origin=first_date)) %>% rename(Pred=Cases)
  
  forecast_mean <- forecast$mean %>% ts_to_dataframe(log=log) %>% mutate(Date=as.Date(Date, origin=first_date)) %>% rename(Pred=Cases)
  forecast_lower <- forecast$lower %>% ts_to_dataframe(log=log) %>% mutate(Date=as.Date(Date, origin=first_date)) %>% rename(Lower=Cases)
  forecast_upper <- forecast$upper %>% ts_to_dataframe(log=log) %>% mutate(Date=as.Date(Date, origin=first_date)) %>% rename(Upper=Cases)
  conf_df <-merge(forecast_lower, forecast_upper, all=TRUE)
  
  y_max <- ifelse(log, exp(max(time_series)), max(time_series))
  return(ggplot(data_df, aes(Date, Cases)) +
    geom_point(color="blue") + geom_line() +
    #geom_point(data=forecast_mean, aes(Date, Pred), color="red") +
    #geom_line(data=forecast_mean, aes(Date, Pred), color="red") +
    geom_point(data=data1, aes(x=Date, y=Daily_New_Deaths), color="blue") +
      geom_point(data=forecast_mean1, aes(Date, Pred), color="yellow") +
    geom_ribbon(data=conf_df, aes(x=Date, y=NULL, ymin=Lower,ymax=Upper), color="red" ,alpha=0.3) +
    #lims(y=c(0, 1.5*y_max)) +
    labs(title="Predicted Covid impact in Romania", x="Date", y="Cases per day") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  )
}

deaths_df <- data_c19 %>% filter(Deaths > 0) %>% head(60)
deaths_df1 <- data_c19 %>% filter(Deaths > 0) %>% tail(19)

nrow(data_c19 %>% filter(Deaths > 0))
tail(deaths_df, 1)
time_series <- xts(deaths_df$Daily_New_Deaths, order.by=deaths_df$Date)
do_time_series(time_series, deaths_df1, preditions = 50, confidence = 95)





deaths_df <- data_c19 %>% filter(Deaths > 0)
time_series <- xts(deaths_df$Deaths, order.by=deaths_df$Date)
do_time_series(time_series, preditions = 20, confidence = 95)

deaths_df_log <- data_c19 %>% filter(Deaths > 0) %>% mutate(Deaths=log(Deaths)) %>% head(60)

deaths_df_log1 <- data_c19 %>% filter(Deaths > 0) %>% mutate(Deaths=log(Deaths)) %>% tail(19)

time_series <- xts(deaths_df_log$Deaths, order.by=deaths_df_log$Date)
do_time_series(time_series, log=TRUE, preditions = 20, confidence = 95)


confirmed_df <- data_c19 %>% filter(Confirmed > 0)
time_series <- xts(confirmed_df$Confirmed, order.by=confirmed_df$Date)
do_time_series(time_series, preditions = 50, confidence = 95)

confirmed_df_log <- data_c19 %>% filter(Confirmed > 0) %>% mutate(Confirmed=log(Confirmed))
time_series <- xts(confirmed_df_log$Confirmed, order.by=confirmed_df_log$Date)
do_time_series(time_series, log=TRUE, preditions = 20, confidence = 95)


deaths_df <- data_c19 %>% filter(Deaths > 0) %>% head(60)
deaths_df1 <- data_c19 %>% filter(Deaths > 0) %>% tail(19)

nrow(data_c19 %>% filter(Deaths > 0))
tail(deaths_df, 1)
time_series <- xts(deaths_df$Daily_New_Deaths, order.by=deaths_df$Date)
do_time_series(time_series, deaths_df1, preditions = 50, confidence = 95)

deaths_df_log <- data_c19 %>% filter(Daily_New_Deaths > 0) %>% mutate(Daily_New_Deaths=log(Daily_New_Deaths+1))
time_series <- xts(deaths_df_log$Daily_New_Deaths, order.by=deaths_df_log$Date)
do_time_series(time_series, log=TRUE, preditions = 20, confidence = 95)


confirmed_df <- data_c19 %>% filter(Confirmed > 0)
time_series <- xts(confirmed_df$Daily_New_Confirmed, order.by=confirmed_df$Date)
do_time_series(time_series, preditions = 20, confidence = 95)

confirmed_df_log <- data_c19 %>% filter(Confirmed > 0) %>% mutate(Daily_New_Confirmed=log(Daily_New_Confirmed+1))
time_series <- xts(confirmed_df_log$Daily_New_Confirmed, order.by=confirmed_df_log$Date)
do_time_series(time_series, log=TRUE, preditions = 20, confidence = 95)




us_police_df <- as.data.frame(read.csv('./data/fatal-police-shootings-data.csv'))
nrow(us_police_df)
head(us_police_df,1)
tail(us_police_df,1)

anova_df <- us_police_df %>% dplyr::select(c(date, age, gender, race, state, threat_level)) %>%
  mutate_all(list(~na_if(.,""))) %>% drop_na() %>%
  mutate(date=as.Date(date), gender=droplevels(as.factor(gender)), race=droplevels(as.factor(race)),
         state=as.factor(state), threat_level=as.factor(threat_level)) %>%
  mutate(race=recode(race, A="American Indian", B="Black", H="Hispanic", N="Not Hispanic", O="Unknown", W="White")) %>%
  filter(race != "Unknown") %>% mutate(race=droplevels(race))

nrow(anova_df)
min(anova_df$date)
max(anova_df$date)
head(anova_df,5)
levels(anova_df$gender)
levels(anova_df$race)
levels(anova_df$state)
levels(anova_df$threat_level)


age_race_histogram <- function(df, title="") {
  g1 <- ggplot(data=df, aes(x=age)) + 
    geom_histogram(breaks=seq(20, 50, by=2), col="red", fill="green", alpha = .2) + 
    labs(title=paste0("Histogram for Age for suspect", title), x="Age", y="Count")  +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  g2<- ggplot(data=df, aes(x=race)) + 
    geom_histogram(stat="count",col="darkblue", fill="blue", alpha = .2) + 
    labs(title=paste0("Histogram for Race for suspect", title), x="Race", y="Count")  +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))

  return (list(g1=g1, g2=g2))
}

plots <- age_race_histogram(anova_df)
grid.arrange(plots$g1, plots$g2, ncol = 2, nrow = 1)


# https://www.businessinsider.com/regions-of-united-states-2018-5
anova_df <- anova_df %>% mutate(Region=fct_collapse(state,
                                                    Northeast=c("ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA", "DC"),
                                                    Midwest=c("OH", "MI", "IN", "WI", "IL", "MN", "IA", "MO", "ND", "SD", "NE", "KS"),
                                                    South=c("DE", "MD", "VA", "WV", "KY", "NC", "SC", "TN", "GA", "FL", "AL", "MS", "AR", "LA", "TX", "OK"),
                                                    West=c("MT", "ID", "WY", "CO", "NM", "AZ", "UT", "NV", "CA", "OR", "WA", "AK", "HI")
                                                    ))

head(anova_df, 4)


group1 <- anova_df %>% filter(Region=="Northeast")
group2 <- anova_df %>% filter(Region=="Midwest")
group3 <- anova_df %>% filter(Region=="South")
group4 <- anova_df %>% filter(Region=="West")
bartlett.test(list(group1$age, group2$age, group3$age, group4$age))


plotsN <- age_race_histogram(group1, " - Northeast")
plotsM <- age_race_histogram(group2, " - Midwest")
plotsS <- age_race_histogram(group3, " - South")
plotsW <- age_race_histogram(group4, " - West")
grid.arrange(plotsN$g1, plotsN$g2, plotsM$g1, plotsM$g2, plotsS$g1, plotsS$g2, plotsW$g1, plotsW$g2, ncol = 2, nrow = 4)


anova_df_18_40 <- anova_df %>% filter(age > 18 & age <= 40)
group1 <- anova_df_18_40 %>% filter(Region=="Northeast")
group2 <- anova_df_18_40 %>% filter(Region=="Midwest")
group3 <- anova_df_18_40 %>% filter(Region=="South")
group4 <- anova_df_18_40 %>% filter(Region=="West")

bartlett.test(list(group1$age, group2$age, group3$age, group4$age))

summary(anova_df_18_40)
nrow(anova_df_18_40)

a<- aov(formula=age~race*Region, data=anova_df_18_40)
summary(a)

a<- aov(formula=age~race*threat_level, data=anova_df_18_40)
summary(a)

anova_df_40_60 <- anova_df %>% filter(age > 40 & age <= 60)
group1 <- anova_df_40_60 %>% filter(Region=="Northeast")
group2 <- anova_df_40_60 %>% filter(Region=="Midwest")
group3 <- anova_df_40_60 %>% filter(Region=="South")
group4 <- anova_df_40_60 %>% filter(Region=="West")

bartlett.test(list(group1$age, group2$age, group3$age, group4$age))

a<- aov(formula=age~race*Region, data=anova_df_40_60)
summary(a)

a<- aov(formula=age~race*threat_level, data=anova_df_40_60)
summary(a)


anova_df_no_race <- anova_df %>% group_by(Region, state, race) %>%
  summarise(no_kills_per_race = length(age),
            avg_age = sum(age) / length(age)) %>% ungroup()

head(anova_df_no_race,5)

levels(anova_df_no_race$Region)
group1 <- anova_df_no_race %>% filter(Region=="Northeast")
group2 <- anova_df_no_race %>% filter(Region=="Midwest")
group3 <- anova_df_no_race %>% filter(Region=="South")
group4 <- anova_df_no_race %>% filter(Region=="West")

bartlett.test(list(group1$no_kills_per_race, group2$no_kills_per_race, group3$no_kills_per_race, group4$no_kills_per_race))

a<- aov(formula=no_kills_per_race~state*avg_age, data=anova_df_no_race)
summary(a)

bartlett.test(list(group1$avg_age, group2$avg_age, group3$avg_age, group4$avg_age))

a<- aov(formula=avg_age~state*no_kills_per_race, data=anova_df_no_race)
summary(a)




us_police_df <- as.data.frame(read.csv('./data/fatal-police-shootings-data.csv'))

lm_police_df <- us_police_df %>% dplyr::select(c(date, age, gender, race, state, threat_level)) %>%
  mutate_all(list(~na_if(.,""))) %>% drop_na() %>%
  mutate(date=as.Date(date), gender=droplevels(as.factor(gender)), race=droplevels(as.factor(race)),
         state=as.factor(state), threat_level=as.factor(threat_level)) %>%
  mutate(race=recode(race, A="American Indian", B="Black", H="Hispanic", N="Not Hispanic", O="Unknown", W="White")) %>%
  filter(race != "Unknown") %>% mutate(race=droplevels(race))

head(lm_police_df,1)
cor.test(lm_police_df$age, as.numeric(lm_police_df$gender))
cor.test(lm_police_df$age, as.numeric(lm_police_df$race))
cor.test(lm_police_df$age, as.numeric(lm_police_df$state))


split_data <- function(data, test) {
  spec = c(train = 1 - test, test = test)
  N = nrow(data)
  
  set.seed(Sys.time())
  data_split = sample(cut( seq(N), N * cumsum(c(0, spec)), labels = names(spec) ))
  
  return(split(data, data_split))
}

rsquared <- function(true, predicted, has_intercept) {
  sse <- sum((predicted - true)^2)
  sst <- ifelse(has_intercept, sum((true - mean(true))^2), sum(true^2))
  rsq <- 1 - sse / sst
  
  return (rsq)
}

mean_squared_error <- function(true, predicted) {
  return( mean((predicted - true) ^ 2) )
}

model_info <- function(y_true, y_models, has_intercept) {
  data = as.data.frame(list(R2=c(), MSE=c()))
  for (name in names(y_models)) {
    y_hats = y_models[[name]]
    model_info <- cbind(R2=c(rsquared(y_true, y_hats, has_intercept)), MSE=c(mean_squared_error(y_true, y_hats)))
    data <- rbind(data, model_info)
  }
  rownames(data) <- names(y_models)
  return(data)
}

summaries <- function(models) {
  for (name in names(models)) {
    model = models[[name]]
    print(summary(model))
    wait()
  }
}

plot_residuals <- function(models) {
  for (name in names(models)) {
    model = models[[name]]
    print(paste('Ploiting residual information for', name))
    plotres(model)
    wait()
  }
}

goodness_of_fit <- function(models) {
  for (name in names(models)) {
    model = models[[name]]
    print(anova(model))
    wait()

    residuals <- residuals(model)
    print(shapiro.test(residuals))
    wait()
    
    print(ggplot(data.frame(Residuals=abs(residuals), Predictions=predict(model)), aes(Predictions, Residuals)) +
             geom_point(color="blue") +
             labs(title=paste0("Homoscedascity for residuals model=", name)) +
             theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    )
    wait()
  }
}



lm_police_df <- lm_police_df %>% split_data(0.1)
nrow(lm_police_df$train)
nrow(lm_police_df$test)
head(lm_police_df$train, 2)


model<- lm(formula=age~date+gender+race+state+threat_level-1, data = lm_police_df$train)
LM_models <- list(
  LM_ALL=model,
  AIC=step(model, trace=0),
  BIC=step(model, trace=0, k=log(nrow(lm_police_df$train))),
  LM_STATE=lm(formula = age~state+gender+race-1, data = lm_police_df$train)
)


summaries(LM_models)

plot_residuals(LM_models)

goodness_of_fit(LM_models)



LM_train_yhats <- list(LM_ALL=predict(LM_models$LM_ALL, lm_police_df$train),
                       AIC=predict(LM_models$AIC, lm_police_df$train),
                       BIC=predict(LM_models$BIC, lm_police_df$train),
                       LM_STATE=predict(LM_models$BIC, lm_police_df$train))

LM_test_yhats <- list(LM_ALL=predict(LM_models$LM_ALL, lm_police_df$test),
                      AIC=predict(LM_models$AIC, lm_police_df$test),
                      BIC=predict(LM_models$BIC, lm_police_df$test),
                      LM_STATE=predict(LM_models$BIC, lm_police_df$test))

print(model_info(lm_police_df$test$age, LM_test_yhats, FALSE))

toplot_df <- data.frame(data_type="TRAIN DATA", date=lm_police_df$train$date, values=lm_police_df$train$age)
toplot_df <- merge(toplot_df, data.frame(data_type="TEST DATA", date=lm_police_df$test$date, values=lm_police_df$test$age), all=TRUE)
toplot_df <- merge(toplot_df, data.frame(data_type="LM_ALL TRAIN", date=lm_police_df$train$date, values=LM_train_yhats$LM_ALL), all=TRUE)
toplot_df <- merge(toplot_df, data.frame(data_type="LM_AIC TRAIN", date=lm_police_df$train$date, values=LM_train_yhats$AIC), all=TRUE)
toplot_df <- merge(toplot_df, data.frame(data_type="LM_BIC TRAIN", date=lm_police_df$train$date, values=LM_train_yhats$BIC), all=TRUE)
toplot_df <- merge(toplot_df, data.frame(data_type="LM_STATE TRAIN", date=lm_police_df$train$date, values=LM_train_yhats$LM_STATE), all=TRUE)
toplot_df <- merge(toplot_df, data.frame(data_type="LM_ALL TEST", date=lm_police_df$test$date, values=LM_test_yhats$LM_ALL), all=TRUE)
toplot_df <- merge(toplot_df, data.frame(data_type="LM_AIC TEST", date=lm_police_df$test$date, values=LM_test_yhats$AIC), all=TRUE)
toplot_df <- merge(toplot_df, data.frame(data_type="LM_BIC TEST", date=lm_police_df$test$date, values=LM_test_yhats$BIC), all=TRUE)
toplot_df <- merge(toplot_df, data.frame(data_type="LM_STATE TEST", date=lm_police_df$test$date, values=LM_test_yhats$LM_STATE), all=TRUE)


ggplot(toplot_df, aes(x=date)) +
  geom_point(aes(y=values, shape=data_type, color=data_type),size=1) +
  scale_shape_manual(values=c("TRAIN DATA"=0, "TEST DATA"=5, "LM_ALL TRAIN"=3, "LM_ALL TEST"=4,
                              "LM_AIC TRAIN"=12, "LM_AIC TEST"=13, "LM_BIC TRAIN"=15,
                              "LM_BIC TEST"=16, "LM_STATE TRAIN"=22, "LM_STATE TEST"=21)) +
  scale_color_manual(values=c("TRAIN DATA"="darkblue","TEST DATA"="blue",
                              "LM_ALL TRAIN"="darkgreen", "LM_ALL TEST"="green",
                              "LM_AIC TRAIN"="darkred", "LM_AIC TEST"="red",
                              "LM_BIC TRAIN"="darkorange", "LM_BIC TEST"="orange",
                              "LM_STATE TRAIN"="darkgray", "LM_STATE TEST"="gray")) +
  labs(title="Age at death in Police shootings Linear Model", x="Observered", y="Age", color='Linear Models', shape='Linear Models') +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

