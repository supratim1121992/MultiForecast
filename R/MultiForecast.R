#' Generate multiple Time Series forecasts
#'
#' Uses a time series object as input and generates multiple forecasts using different models to provide an ensembled forecast
#' @param ts The Time Series object to be used for forecasting
#' @param rdm The distinct numerical values to be used in combination as p,d,q (trend component) & P,D,Q (seasonal component) for the ARIMA models
#' @param act The actual values for the validation period based on which the models will be evaluated
#' @param fper The forecasting period
#' @param type The type of Time Series being analysed ("additive" or "multiplicative")
#' @param mod The different models to be used for forecasting ("aa" for auto.arima,"sm" for Simulated ARIMA,"hw" for HoltWinters)
#' @param weight If the ensembled forecasts should be generated based on a weighted average as per the errors from each model instead of a simple mean
#' @param xreg Data for external regressors to be used for forecasting
#' @param freg Data for external regressors for the validation period
#' @return A list containing the forecasted values from each model along with the ensembled forecast
#' @export

MultiForecast<-function(ts,rdm = c(0,1),act,fper,type = "additive",mod = c("aa","sm","hw"),weight = F,xreg = NULL,freg = NULL){
  require(forecast)
  require(data.table)
  MAPE<-function(act,pred){
    mean(abs(act-pred)/abs(act))*100
  }
  pred<-rep(0,times = 5)
  dim(pred)<-c(1,5)
  pred<-as.data.frame(pred)
  colnames(pred)<-c("Point Forecast","Lo 80","Hi 80","Lo 95","Hi 95")
  fit<-rep(0,times = length(ts))
  dim(fit)<-c(1,length(ts))
  fit<-as.data.frame(fit)

  plot(decompose(ts))

  if ("aa" %in% mod == T){
    if (type == "multiplicative"){
      ts_aa<-log(ts)
    }
    else if (type == "additive"){
      ts_aa<-ts
    }
    mod_aa<-auto.arima(ts_aa,xreg = xreg)
    cat("Completed Auto Arima model. \n")
    fore_aa<-forecast(object = mod_aa,h = fper,xreg = freg)
    if (type == "additive"){
      pred_aa<-as.integer(fore_aa$mean)
    }
    else if (type == "multiplicative"){
      pred_aa<-fore_aa$mean
      pred_aa<-as.integer(exp(pred_aa))
    }
    mape_aa<-MAPE(act,pred_aa)
    mape_aa2<-MAPE(as.numeric(ts_aa),as.numeric(mod_aa$fitted))
    cat("Training MAPE for Auto Arima model:",mape_aa2,"\n")
    cat("Validation MAPE for Auto Arima model:",mape_aa,"\n")
    rng_aa<-range(c(as.numeric(mod_aa$fitted),as.numeric(ts_aa)))
    plot(as.numeric(ts_aa),type = "l",lwd = 2,main = "Auto Arima Model Fit",xlab = "Time",ylab = "Counts",
         ylim = c(rng_aa[1]*0.8,rng_aa[2]*1.2))
    lines(as.numeric(mod_aa$fitted),col = "red")
    plot(fore_aa)
    if (type  == "multiplicative"){
      aa_dt<-as.data.frame(fore_aa)
      cnames<-colnames(aa_dt)
      aa_dt<-as.data.frame(lapply(X = aa_dt,FUN = exp))
      colnames(aa_dt)<-cnames
      fit<-as.data.frame(rbind(fit,exp(as.numeric(mod_aa$fitted))))
    }
    else if (type == "additive"){
      aa_dt<-as.data.frame(fore_aa)
      fit<-as.data.frame(rbind(fit,as.numeric(mod_aa$fitted)))
    }
    pred<-rbind(pred,aa_dt)
    cat("Model fit and forecast plots generated. \n")
    cat("\n")
  }

  if("sm" %in% mod == T){
    SimArima<-function(ts,rdm,act,fper,type,xreg,freg){
      Parameter<-function(rdm){
        par_vec<-character()
        check<-TRUE
        iter<-0
        iter_max<-1000
        len_ini<-0
        len_fin<-1
        Sampler<-function(rdm){
          par_1<-sample(x = rdm,size = 1,replace = F)
          par_2<-sample(x = rdm,size = 1,replace = F)
          par_3<-sample(x = rdm,size = 1,replace = F)
          par_4<-sample(x = rdm,size = 1,replace = F)
          par_5<-sample(x = rdm,size = 1,replace = F)
          par_6<-sample(x = rdm,size = 1,replace = F)
          par_i<-paste(par_1,par_2,par_3,par_4,par_5,par_6,sep = ",")
          return(par_i)
        }
        while (len_ini<len_fin){
          len_ini<-length(par_vec)
          while (check==T & iter<iter_max){
            par<-Sampler(rdm = rdm)
            check<-par %in% par_vec
            iter<-iter+1
          }
          if (par %in% par_vec == T){
            par_vec<-par_vec
          }
          else {
            par_vec<-c(par_vec,par)
          }
          check<-TRUE
          iter<-0
          len_fin<-length(par_vec)
        }
        return(par_vec)
      }
      prm<-Parameter(rdm)
      prm_list<-strsplit(x = prm,split = ",")
      aic_list<-numeric(length = length(prm))
      aicc_list<-numeric(length = length(prm))
      bic_list<-numeric(length = length(prm))
      MAPE_list<-numeric(length = length(prm))
      for (i in 1:length(prm_list)){
        o1<-as.numeric(prm_list[[i]][1])
        o2<-as.numeric(prm_list[[i]][2])
        o3<-as.numeric(prm_list[[i]][3])
        s1<-as.numeric(prm_list[[i]][4])
        s2<-as.numeric(prm_list[[i]][5])
        s3<-as.numeric(prm_list[[i]][6])
        tryCatch({
          mod<-Arima(y = ts,order = c(o1,o2,o3),seasonal = c(s1,s2,s3),xreg = xreg)
          if (type == "multiplicative"){
            pred<-as.integer(exp(forecast(object = mod,h = fper,xreg = freg)$mean))
          }
          else if (type == "additive"){
            pred<-as.integer(forecast(object = mod,h = fper,xreg = freg)$mean)
          }
          aic_list[i]<-mod$aic
          aicc_list[i]<-mod$aicc
          bic_list[i]<-mod$bic
          MAPE_list[i]<-MAPE(act = act,pred = pred)
        }, error=function(e){})
      }
      res<-as.data.frame(cbind("Parameters"=prm,"AIC"=aic_list,"AICC"=aicc_list,"BIC"=bic_list,
                               "MAPE"=MAPE_list))
      return(res)
    }
    if (type == "multiplicative"){
      ts_sm<-log(ts)
    }
    else if (type == "additive"){
      ts_sm<-ts
    }
    simres<-SimArima(ts = ts_sm,rdm = rdm,act = act,type = type,fper = fper,xreg = xreg,freg = freg)
    simres$AIC<-as.numeric(as.character(simres$AIC))
    simres$AICC<-as.numeric(as.character(simres$AICC))
    simres$BIC<-as.numeric(as.character(simres$BIC))
    simres$MAPE<-as.integer(as.character(simres$MAPE))
    simres<-as.data.table(simres)
    simres<-simres[MAPE != 0]
    simres<-simres[Parameters != "0,0,0,0,0,0"]
    simres<-simres[order(MAPE,AIC,BIC,AICC)]
    simpara<-as.character(simres$Parameters[1])
    simpara<-as.integer(unlist(strsplit(x = simpara,split = ",")))
    mod_sm<-Arima(y = ts_sm,order = simpara[1:3],seasonal = simpara[4:6],xreg = xreg)
    cat("Completed Simulated Arima model. \n")
    fore_sm<-forecast(object = mod_sm,h = fper,xreg = freg)
    if (type == "additive"){
      pred_sm<-as.integer(fore_sm$mean)
    }
    else if (type == "multiplicative"){
      pred_sm<-fore_sm$mean
      pred_sm<-as.integer(exp(pred_sm))
    }
    mape_sm<-MAPE(act,pred_sm)
    mape_sm2<-MAPE(as.numeric(ts_sm),as.numeric(mod_sm$fitted))
    cat("Training MAPE for Simulated Arima model:",mape_sm2,"\n")
    cat("Validation MAPE for Simulated Arima model:",mape_sm,"\n")
    rng_sm<-range(c(as.numeric(mod_sm$fitted),as.numeric(ts_sm)))
    plot(as.numeric(ts_sm),type = "l",lwd = 2,main = "Simulated Arima Model Fit",xlab = "Time",ylab = "Counts",
         ylim = c(rng_sm[1]*0.8,rng_sm[2]*1.2))
    lines(as.numeric(mod_sm$fitted),col = "red")
    plot(fore_sm)
    if (type  == "multiplicative"){
      sm_dt<-as.data.frame(fore_sm)
      cnames<-colnames(sm_dt)
      sm_dt<-as.data.frame(lapply(X = sm_dt,FUN = exp))
      colnames(sm_dt)<-cnames
      fit<-as.data.frame(rbind(fit,exp(as.numeric(mod_sm$fitted))))
    }
    else if (type == "additive"){
      sm_dt<-as.data.frame(fore_sm)
      fit<-as.data.frame(rbind(fit,as.numeric(mod_sm$fitted)))
    }
    pred<-rbind(pred,sm_dt)
    cat("Model fit and forecast plots generated. \n")
    cat("\n")

  }

  if("hw" %in% mod == T){
    if (is.null(xreg) == F){
      cat("Constructing Holt Winters model without regressors. \n")
    }
    if (type == "multiplicative"){
      ts_hw<-log(ts)
    }
    else if (type == "additive"){
      ts_hw<-ts
    }
    mod_hw<-hw(y = ts_hw,h = fper,seasonal = type)
    cat("Completed Holt Winters model. \n")
    fore_hw<-forecast(object = mod_hw,h = fper)
    if (type == "multiplicative"){
      pred_hw<-fore_hw$mean
      pred_hw<-as.integer(exp(pred_hw))
    }
    else if (type == "additive"){
      pred_hw<-as.integer(fore_hw$mean)
    }
    mape_hw<-MAPE(act,pred_hw)
    mape_hw2<-MAPE(as.integer(ts_hw),as.numeric(mod_hw$fitted))
    cat("Training MAPE for Holt Winters model:",mape_hw2,"\n")
    cat("Validation MAPE for Holt Winters model:",mape_hw,"\n")
    rng_hw<-range(c(as.numeric(mod_hw$fitted),as.numeric(ts_hw)))
    plot(as.numeric(ts_hw),type = "l",lwd = 2,main = "Holt Winters Model Fit",xlab = "Time",ylab = "Counts",
         ylim = c(rng_hw[1]*0.8,rng_hw[2]*1.2))
    lines(as.numeric(mod_hw$fitted),col = "red")
    plot(fore_hw)
    if (type  == "multiplicative"){
      hw_dt<-as.data.frame(fore_hw)
      cnames<-colnames(hw_dt)
      hw_dt<-as.data.frame(lapply(X = hw_dt,FUN = exp))
      colnames(hw_dt)<-cnames
      fit<-as.data.frame(rbind(fit,exp(as.numeric(mod_hw$fitted))))
    }
    else if (type == "additive"){
      hw_dt<-as.data.frame(fore_hw)
      fit<-as.data.frame(rbind(fit,as.numeric(mod_hw$fitted)))
    }
    pred<-rbind(pred,hw_dt)
    cat("Model fit and forecast plots generated. \n")
    cat("\n")

  }

  if (weight == T){
    wt<-numeric()
    if ("aa" %in% mod == T){
      ac_aa<-100-mape_aa
      wt<-c(wt,ac_aa)
    }
    if ("sm" %in% mod == T){
      ac_sm<-100-mape_sm
      wt<-c(wt,ac_sm)
    }
    if ("hw" %in% mod == T){
      ac_hw<-100-mape_hw
      wt<-c(wt,ac_hw)
    }
    wt_rt<-wt/sum(wt)
    fit<-fit[-1,]
    pred<-pred[-1,]
    if (length(wt_rt) == 3){
      fit<-fit[1,]*wt_rt[1] + fit[2,]*wt_rt[2] + fit[3,]*wt_rt[3]
      for (i in 1:fper){
        pred[i,]<-pred[i,]*wt_rt[1] + pred[(i+fper),]*wt_rt[2] + pred[(i+2*fper),]*wt_rt[3]
      }
    }
    else if (length(wt_rt) == 2){
      fit<-fit[1,]*wt_rt[1] + fit[2,]*wt_rt[2]
      for (i in 1:fper){
        pred[i,]<-pred[i,]*wt_rt[1] + pred[(i+fper),]*wt_rt[2]
      }
    }
    else if (length(wt_rt) == 1){
      fit<-fit[1,]*wt_rt[1]
      for (i in 1:fper){
        pred[i,]<-pred[i,]*wt_rt[1]
      }
    }
  }
  else {
    fit<-fit[-1,]
    if (length(mod) == 3){
      fit<-(fit[1,] + fit[2,] + fit[3,])/3
      pred<-pred[-1,]
      for (i in 1:fper){
        pred[i,]<-(pred[i,] + pred[(i+fper),] + pred[(i+2*fper),])/3
      }
    }
    else if (length(mod) == 2){
      fit<-(fit[1,] + fit[2,])/2
      pred<-pred[-1,]
      for (i in 1:fper){
        pred[i,]<-(pred[i,] + pred[(i+fper),])/2
      }
    }
    else if (length(mod) == 1){
      fit<-fit[1,]
      pred<-pred[-1,]
      for (i in 1:fper){
        pred[i,]<-pred[i,]
      }
    }
  }

  fore<-pred[1:fper,]
  pred_mn<-as.integer(fore[,1])
  cat("Final forecasted values:",pred_mn,"\n")
  mape_mn<-MAPE(act,pred_mn)
  cat("Validation MAPE from ensembling:",mape_mn,"\n")
  if (frequency(ts) == 12){
    mini<-min(unlist(lapply(X = fore,FUN = min)))
    maxi<-max(unlist(lapply(X = fore,FUN = max)))
    mini<-min(c(mini,as.integer(ts)))*0.8
    maxi<-max(c(maxi,as.integer(ts)))*1.2
    plot(c(as.numeric(ts),rep(NA,times = fper)),type = "l",lwd = 2,ann = F,ylim = c(mini,maxi),xaxt = "n")
    lines(as.numeric(fit),col = "red")
    lines(c(rep(NA,times = length(ts)),fore[,2]),lwd = 2,col = "grey50")
    lines(c(rep(NA,times = length(ts)),fore[,3]),lwd = 2,col = "grey50")
    lines(c(rep(NA,times = length(ts)),fore[,4]),lwd = 2,col = "grey65")
    lines(c(rep(NA,times = length(ts)),fore[,5]),lwd = 2,col = "grey65")
    polygon(x = c((length(ts)+1):(length(ts)+fper),rev((length(ts)+1):(length(ts)+fper))),
            y = c(fore[,5],rev(fore[,4])),col = "grey65",border = NA)
    polygon(x = c((length(ts)+1):(length(ts)+fper),rev((length(ts)+1):(length(ts)+fper))),
            y = c(fore[,3],rev(fore[,2])),col = "grey50",border = NA)
    lines(c(rep(NA,times = length(ts)),act),lwd = 2,col = "white",type = "b")
    lines(c(rep(NA,times = length(ts)),fore[,1]),lwd = 3,col = "grey25")
    ts2<-ts(data = c(as.numeric(ts),pred_mn),start = start(ts),frequency = frequency(ts))
    start<-paste(start(ts2),collapse = "-")
    start<-paste(start,"01",sep = "-")
    start<-as.Date(start,format = "%Y-%m-%d")
    end<-paste(paste(end(ts2),collapse = "-"),"01",sep = "-")
    end<-as.Date(end,format = "%Y-%m-%d")
    seq<-seq(from = start,to = end,by = "month")
    seq<-substr(x = seq,start = 1,stop = 7)
    id<-c(1,((1:as.integer(length(seq)/(frequency(ts)/2)))*frequency(ts)/2))
    axis(side = 1,at = c(1,seq(from = 6,to = length(ts2),by = 6)),labels = seq[id],las = 2)
    title(main = "Forecasts after ensembling",ylab = "Counts")
  }
  res<-list("Auto Arima Forecasts" = pred_aa,"Simulated Arima Forecasts" = pred_sm,"HoltWinters Forecast" = pred_hw,
            "Ensembled Forecast Values" = pred_mn,"Simulation Results" = simres,"Forecast Table" = fore)
  return(res)
}
