# randomly sample n runids from the indicated scenario
#  input: key.filename is table mapping scenario ids to concentration and duration,
#			  concentration is {0.01,0.1,1.0}, duration is {72,120,168,240,360,480},
#         n is the number of outrbreaks, between 1 and 100
#  output: a vector of runids
get.runids <- function(key.filename, concentration, duration, n) {
  key = read.csv(key.filename, colClasses=c("character", "numeric", "integer"))
  scenarioids = substr(key$file_name, 1, (nchar(key$file_name)-6))
  scenarioid = scenarioids[which(key$concentration_per_litre==concentration & key$duration_in_hours==duration)]
  
  # generate runids to load for scenario
  ids = seq(0,99)
  if (n > 0 & n <= 100) {
    runids = as.character(sample(ids, n))
    
    # if single digit, prepend a "0" to match file name
    singledigit.runids = which(nchar(runids)==1)
    if (length(singledigit.runids > 0)) {
      runids[singledigit.runids] = paste("0",runids[singledigit.runids],sep="")
    } # if - singledigit.runids
    
  } else {
    runids = NULL
  } # if - n
  
  return(paste(scenarioid, runids, sep=""))
} # get.runids


# load the runs for the indicated ids
#  input: directory is the realtive path to the datafiles, runids is a vector
#         of runids to load.
#  output: a list of data.frames contain data for inidcated runs
load.runs.old <- function(directory, runids, os="mac") {
  
  if (os == "mac") parse = 43
  else parse = 39
  
  filenames = list.files(directory, pattern= "*.csv", full.names=TRUE)
  
  filenames.toread = filenames[substr(filenames, parse, nchar(filenames)-4) %in% runids]
  
  return(lapply(filenames.toread, read.csv))
} # load.runs.old


# load the runs for the indicated ids without having to rely on directory structure
#  input: directory is the realtive path to the datafiles, runids is a vector
#         of runids to load.
#  output: a list of data.frames contain data for inidcated runs
load.runs <- function(directory, runids) {
  filenames = list.files(directory, pattern="*.csv", full.names=TRUE)
  
  if (length(filenames) == 0)
    stop("load.runs: Cannot access filenames for runs - check directory path - directory: " + directory)
  
  names_no_ext = gsub(".csv", '', basename(filenames))
  return(lapply(filenames[names_no_ext %in% runids], read.csv))
} # load.runs


# calculate summary information on outbreak
#  input df has 3 col: date, baseline, baseline + outbreak
#  output is list with elements: state vector, start date, end date, length and max height of outbreak
o.summary <- function(df) {
  
  date = as.Date(df[,1])
  o = df[,3] - df[,2]
  is.o = rep(FALSE, nrow(df))
  is.o[o > 0] = TRUE
  
  o.start = min(date[is.o == TRUE])
  o.end = max(date[is.o == TRUE])
  o.length = as.numeric(o.end - o.start)
  
  o.max = max(o)
  
  return(list(state=is.o, start=o.start, end=o.end, length=o.length, height=o.max))
  
} # o.summary


# calculate false alert rate during non-outbreak intervals
#  input df is data frame with a column "alarm" that contains binary indicator of alarm from algorithm
#  input truth is list with element: state, which is binary vector indicating presence of outbreak
#  output is a single number, the false alert rate per unit time
a.far <- function(df, truth) {
  state = truth$state
  alarm = df$alarm
  
  # Alarm is a numeric vector, with value 1 for an alarm, 0 for no alarm
  # Summing alarm vector when state is FALSE gives false alerts
  fa = sum(alarm[state==FALSE & !is.na(alarm)])
  # Length of alarm vector when state is FALSE gives number of days with an alarm value
  far = fa / length(alarm[state == FALSE & !is.na(alarm)])
  
  return(far)
  
} # a.far


# calculate timeliness for a single outbreak
#  input df is data frame with 3 col: date, baseline, baseline + outbreak
#  input truth is list with element: state, which is binary vector indicating presence of outbreak
#  output is proportion of outbreak "prevented" (from 1 for first day to 0 for last day of outbreak, or later)
o.prevented <- function(df, truth) {
  state = truth$state
  alarm = df$alarm
  length = sum(truth$state)
  prevented = 0
  
  detect = sum(alarm[state==TRUE]) > 0
  if (detect) {
    first.alarm = min(which(alarm[state==TRUE]==TRUE))
    prevented = (length - first.alarm) / length
  } # if
  
  return(prevented)
} # o.prevented


o.detected <- function(df, truth) {
  return(sum(df$alarm[truth$state==TRUE]) > 0)
}


# derive points on an ROC curve
# input: predictions, truth
# return: data.frame with threshold, fpr, tpr, prevented

# sort unique values of prediction, add "anchor" values at either end
# at each prediction value
# apply threshold to determine alarm values
# score far, detection, prevented
# add to data.frame evaluation results along with threshold
a.performance <- function(df, cutoffs) {
  state = df$state
  prediction = df$prediction
  
  performance = data.frame(cutoff=cutoffs,
                           far=rep(NA, length(cutoffs)),
                           detected=rep(NA, length(cutoffs)),
                           prevented=rep(NA, length(cutoffs)))
  
  for (cutoff in cutoffs) {
    alarms = data.frame(alarm=rep(NA, length(prediction)))
    alarms$alarm[prediction > cutoff] = TRUE
    alarms$alarm[prediction <= cutoff] = FALSE
    
    truth = data.frame(state = df$state)
    far = a.far(alarms, truth)
    detected = o.detected(alarms, truth)
    prevented = o.prevented(alarms, truth)
    
    performance[performance$cutoff == cutoff,] = c(cutoff, far, detected, prevented)
  } # for
  
  return(as.matrix(performance))
} # a.performance


# get cutoffs that cover range of predictions
get.cutoffs <- function(results, N) {
  # get range in each run
  ranges = lapply(results, function(df) range(df$prediction, na.rm=TRUE))
  # get range of range
  range = range(unlist(ranges, recursive=TRUE))
  
  return(cutoffs = c(-Inf, seq(range[1], range[2], length.out=N-2), Inf))
}


# add truth vector to algorithm results
add.truth <- function(res, outbreaks) {
  
  new.res = list(rep(NA, length(res)))
  states = sapply(outbreaks, "[", "state")
  
  for (i in 1:length(res)) {
    result = res[[i]]
    result$state = states[[i]]
    new.res[[i]] = result
  }
  
  return(new.res)
} # add.truth


a.performance.all <- function(results, outbreaks, n.cutoffs) {
  a.cutoffs = get.cutoffs(results, n.cutoffs)
  results.truth = add.truth(results, outbreaks)
  performance.all = lapply(results.truth, a.performance, cutoffs=a.cutoffs)
  
  return(performance.all)
} # a.performance.all


a.performance.avg <- function(performance.all) {
  n.runs = length(performance.all)
  performance.avg = Reduce('+', performance.all) / n.runs
  
  return(as.data.frame(performance.avg))
} # a.performance.avg


# calculate Poisson regression, limit, and alarm value
#  for last value in the time series 'x'.
poisson <- function(x, dates=NA, gap=2, window=56, interval=14) {
  observed = NA
  expected = NA
  p = NA
  alarm.value = NA
  
  x.l = length(x)
  
  if (x.l > (gap + window)) {
    
    
    # add indicator variable
    full.interval = floor(x.l / interval)
    part.interval = x.l - full.interval*interval
    # note indicator values are reveresed to ensure level exits for predicted data
    if (part.interval > 0) {
      indicator = as.factor(rev(c(rep(1:full.interval, each=interval), rep(full.interval+1, times=part.interval))))
    } else {
      indicator = as.factor(rev(rep(1:full.interval, each=interval)))
    } # if - part.interval
    
    # if dates not provided, fit model with only indicator
    if (length(dates) == 1) {
      x.df = data.frame(x=x, i=indicator)
      x.fit = x.df[(x.l-gap-1-window):(x.l-gap-1),]
      x.predict = data.frame(i=tail(indicator, 1))
      
      # fit glm to x[1:(x.l-gap)]
      x.glm = glm(x ~ i, family="quasipoisson", data=x.fit)
      
    } else {
      # if dates provided, fit model with indicator and dow
      dow = as.factor(format(dates, "%a"))
      
      x.df = data.frame(o=x, i=indicator, d=dow)
      x.fit = x.df[(x.l-gap-1-window):(x.l-gap-1),]
      x.predict = data.frame(i=tail(indicator, 1), d=x.df$d[x.l])
      
      # fit glm
      x.glm = glm(o ~ i + d, family="quasipoisson", data=x.fit)
    } # if - length dates
    
    expected = predict(x.glm, x.predict, type="response")
    observed = x[length(x)]
  } # if
  
  return(c(observed, expected))
} # poisson


# apply poisson regression to whole dataset
poisson_all <- function(df, dow=TRUE, gap=2, window=56, interval=14, threshold=0.05) {
  x = df$baseline.outbreak_count
  
  if (dow) {
    dates = as.Date(df$date)
  } else {
    dates = NA
  } # if
  
  ndays = length(x)
  
  results = data.frame(observed = rep(NA, ndays),
                       expected = rep(NA, ndays),
                       prediction = rep(NA, ndays),
                       alarm = rep(NA, ndays))
  
  # Apply to data series
  for (day in (window+gap+1):ndays) {
    
    if (dow) {
      results.day = poisson(x[1:day], dates=dates[1:day], gap, window, interval)
    } else {
      results.day = poisson(x[1:day], dates=dates, gap, window, interval)
    } # if
    
    observed = c(results$observed[(day-(window+gap)):(day-1)], results.day[1])
    expected = c(results$expected[(day-(window+gap)):(day-1)], results.day[2])
    
    # compute sd to this point
    if (length(observed > 3)) {
      sd = sqrt((sum((observed-expected)^2) / (length(observed) - 2 - 1)))
      
      # calculated z for latest observation
      z = (results.day[1] - results.day[2]) / sd
      p = pnorm(z)
      
      # determine if alarm
      alarm = p < threshold
      
    } # if - enough points for sd
    
    results[day,] = c(results.day[1], results.day[2], p=p, alarm=alarm)
    
  } # for
  
  return(results)
} # poisson_all

# calculate Poisson regression, limit, and alarm value
#  for last value in the time series 'x'.
poisson.aman <- function(y, full.mm, gap=2, window=56, interval=14) {
  y.l = length(y)
  if (y.l <= (gap + window)) return(c(NA,NA))
  
  mm <- full.mm[(y.l-gap-1-window):(y.l-gap-1),]
  x.predict <- full.mm[nrow(full.mm),]
  
  x.glm = glm.fit(x = mm,
                  y = y[(y.l-gap-1-window):(y.l-gap-1)],
                  family=quasipoisson())
  
  coefs <- ifelse(is.na(coef(x.glm)), 0, coef(x.glm))
  c(observed = y[length(y)], expected = exp(x.predict %*% coefs))
} # poisson


# apply poisson regression to whole dataset
poisson_all.aman <- function(df, dow=TRUE, gap=2, window=56, interval=14, threshold=0.05) {
  y <- df$baseline.outbreak_count
  if (dow) d <- format(as.Date(df$date), '%a')
  n.days = length(y)
  
  i <- as.factor(rep(1:ceiling(n.days / interval), each=interval)[1:n.days])
  indicator <- model.matrix(~i) # Include intercept here.
  if(dow) dates <- model.matrix(~d-1)
  
  make.mm <- function(day)
    cbind(indicator[day:1,], if(dow) dates[1:day,] else NULL)
  
  results.day <- lapply((window+gap+1):n.days, function(day)
    poisson.aman(y[1:day], make.mm(day),
                 gap, window, interval))
  
  results.day <- data.frame(do.call(rbind, results.day))
  
  # Rolling sd
  diff <- results.day$observed - results.day$expected
  sd <- sapply(1:nrow(results.day), function(x) sd(diff[1:x], na.rm=TRUE))
  
  z = diff / sd
  results.day$p = dnorm(z)
  results.day$alarm = as.numeric(results.day$p < threshold)
  
  results.day[c(rep(NA,window+gap),1:nrow(results.day)),]
} # poisson_all




# calculate C2 test statistic, x limit, and alarm value
#  for last value in the time series 'x'.
c2 <- function (x, gap=2, window=28, threshold=3.9) {
  test.statistic = NaN
  alarm.value = NaN
  x.limit = NaN
  
  x.l = length(x)
  
  if (x.l > (gap + window)) {
    # todo missing values
    x.mean = mean(x[(x.l - gap - window):(x.l - gap - 1)])
    x.sd = sd(x[(x.l - gap - window):(x.l - gap - 1)])
    
    # todo divide by zero
    test.statistic = max(0, (x[length(x)] - x.mean) / x.sd)
    x.limit = x.mean + threshold * x.sd
    alarm.value = (x[length(x)] > x.limit)
  }
  
  return(c(test.statistic, x.limit, alarm.value))
} # c2


# apply c2 to whole dataset
c2_all <- function(df, gap=2, window=28, threshold=3.9) {
  
  x = df[,3]
  ndays = length(x)
  
  results = data.frame(prediction = rep(NaN, ndays),
                       limit = rep(NaN, ndays),
                       alarm = rep(NaN, ndays))
  
  # Apply c2 to test series
  for (day in 1:ndays) {
    results[day,] = c2(x[1:day], gap=gap, window=window,
                       threshold=threshold)
  } # for
  
  return(results)
} # c2_all