#' @title Dataset with manipulated (cleaned) x and y points
#'
#' @description This function manipulate x and y points according to KM curve characteristics, which are left and upward resistant. (x2 can't be less than x1, y2 can't be greater than y1)
#' @param data input data
#' @export
#' @return cleaned dataset
#' @keywords internal
#' @examples
#' data <- data.org(data)

data.org <- function(data){
  for (i in 2:nrow(data)){
    data[1,1:2] <- c(0,1)
    if (data[i,1] < data[i-1,1]) {data[i,1] <- data[(i-1),1]}
    if (data[i,2] > data[i-1,2]) {data[i,2] <- data[(i-1),2]}
  }
  return(data)
}




#' @title IPD Calculation
#'
#' @description This function calculate all the IPDs based on input data
#' @param t.points vector of time to event points; they represents the x axis values marked from the KM plot using digizeit, which is greater than 0.
#' @param s.points vector of survival rate points; they represents the y axis values marked from the KM plot using digizeit, which ranges from 0 to 1.
#' @param t.risk.T vector of time points in the number at risk table from the original KM plot, which starts from zero.
#' @param n.risk.T vector of number at risk at each time point in the number at risk table from the original KM plot, which has the same length as t.risk.T.
#' @param n.t.T number of clicks in the group
#' @param lower.T numeric vector; number of data points at start between time intervals in the number at risk table
#' @param upper.T numeric vector; number of data points at end between time intervals in the number at risk table
#' @param t.event number of events
#' @param gr.number name for the group
#' @keywords internal
#' @export
#' @return data frame with time to event, censoring, and group name information
#' 


km.cal.tab <- function(t.points, 
                       s.points, 
                       t.risk.T, n.risk.T, 
                       lower.T, upper.T, 
                       t.event="NA",gr.number="group"){
  IPD <- NULL  
  
  t.S 	<- t.points
  S 	<- s.points
  t.risk	<- t.risk.T
  n.risk	<- n.risk.T
  n.int 	<- length(n.risk)
  n.t	<- length(t.points)	    # change
  
  lower <- lower.T
  upper	<- upper.T	
  
  tot.events 	<- t.event	
  
  #Initialise vectors
  n.censor<- rep(0, n.int) 	      	       
  n.hat<-rep(n.risk[1], n.t)      
  cen<-rep(0,n.t)		              
  d<-rep(0,n.t)	
  KM.hat<-rep(1,n.t)
  last.i<-rep(1,n.int)
  sumdL<-0
  
  
  if (n.int > 1){
    #Time intervals 1,...,(n.int-1)
    for (i in 1:(n.int-1)){
      
      
      #First approximation of no. censored on interval i
      n.censor[i]<- round(n.risk[i]*S[lower[i+1]]/S[lower[i]]- n.risk[i+1])
      
      #Adjust tot. no. censored until n.hat = n.risk at start of interval (i+1)
      while((n.hat[lower[i+1]]>n.risk[i+1])||((n.hat[lower[i+1]]<n.risk[i+1])&&(n.censor[i]>0))){
        if (n.censor[i]<=0){
          cen[lower[i]:upper[i]]<-0
          n.censor[i]<-0
        }
        if (n.censor[i]>0){
          cen.t<-rep(0,n.censor[i])
          for (j in 1:n.censor[i]){
            cen.t[j]<- t.S[lower[i]] +
              j*(t.S[lower[(i+1)]]-t.S[lower[i]])/(n.censor[i]+1)
          }
          #Distribute censored observations evenly over time. Find no. censored on each time interval.
          cen[lower[i]:upper[i]]<-hist(cen.t,breaks=t.S[lower[i]:lower[(i+1)]],
                                       plot=F)$counts
        }
        #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
        n.hat[lower[i]]<-n.risk[i]
        last<-last.i[i]
        for (k in lower[i]:upper[i]){
          if (i==1 & k==lower[i]){
            d[k]<-0
            KM.hat[k]<-1
          }
          else {
            KM.hat <- ifelse(is.na(KM.hat), 0.01, KM.hat)   ## HH
            d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
            
            KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          }
          n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
          # d[is.na(d)] <- 0   ## HH
          if (d[k] != 0) last<-k
        }
        n.censor[i]<- n.censor[i]+(n.hat[lower[i+1]]-n.risk[i+1])
      } # end of while loop
      if (n.hat[lower[i+1]]<n.risk[i+1]) n.risk[i+1]<-n.hat[lower[i+1]]
      last.i[(i+1)]<-last
      
      
    }
  }
  
  
  #Time interval n.int.
  if (n.int>1){
    #Assume same censor rate as average over previous time intervals.
    n.censor[n.int]<- min(round(sum(n.censor[1:(n.int-1)])*(t.S[upper[n.int]]-
                                                              t.S[lower[n.int]])/(t.S[upper[(n.int-1)]]-t.S[lower[1]])), n.risk[n.int])
  }
  if (n.int==1){n.censor[n.int]<-0}
  if (n.censor[n.int] <= 0){
    cen[lower[n.int]:(upper[n.int]-1)]<-0
    n.censor[n.int]<-0
  }
  if (n.censor[n.int]>0){
    cen.t<-rep(0,n.censor[n.int])
    for (j in 1:n.censor[n.int]){
      cen.t[j]<- t.S[lower[n.int]] +
        j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
    }
    cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],
                                             plot=F)$counts
  }
  
  
  #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
  n.hat[lower[n.int]]<-n.risk[n.int]
  last<-last.i[n.int]
  for (k in lower[n.int]:upper[n.int]){
    if(KM.hat[last] !=0){
      d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))} else {d[k]<-0}
    KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
    n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
    #No. at risk cannot be negative
    if (n.hat[k+1] < 0) {
      n.hat[k+1]<-0
      cen[k]<-n.hat[k] - d[k]
    }
    if (d[k] != 0) last<-k
  }
  
  
  
  
  
  #If total no. of events reported, adjust no. censored so that total no. of events agrees.
  if (tot.events != "NA"){
    if (n.int>1){
      sumdL<-sum(d[1:upper[(n.int-1)]])
      #If total no. events already too big, then set events and censoring = 0 on all further time intervals
      if (sumdL >= tot.events){
        d[lower[n.int]:upper[n.int]]<- rep(0,(upper[n.int]-lower[n.int]+1))
        cen[lower[n.int]:(upper[n.int]-1)]<- rep(0,(upper[n.int]-lower[n.int]))
        n.hat[(lower[n.int]+1):(upper[n.int]+1)]<- rep(n.risk[n.int],(upper[n.int]+1-lower[n.int]))
      }
    }
    #Otherwise adjust no. censored to give correct total no. events
    if ((sumdL < tot.events)|| (n.int==1)){
      sumd<-sum(d[1:upper[n.int]])
      while ((sumd > tot.events)||((sumd< tot.events)&&(n.censor[n.int]>0))){
        n.censor[n.int]<- n.censor[n.int] + (sumd - tot.events)
        if (n.censor[n.int]<=0){
          cen[lower[n.int]:(upper[n.int]-1)]<-0
          n.censor[n.int]<-0
        }
        if (n.censor[n.int]>0){
          cen.t<-rep(0,n.censor[n.int])
          for (j in 1:n.censor[n.int]){
            cen.t[j]<- t.S[lower[n.int]] +
              j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
          }
          
          cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],
                                                   plot=F)$counts
        }
        n.hat[lower[n.int]]<-n.risk[n.int]
        last<-last.i[n.int]
        for (k in lower[n.int]:upper[n.int]){
          d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
          KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          if (k != upper[n.int]){
            n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
            #No. at risk cannot be negative
            if (n.hat[k+1] < 0) {
              n.hat[k+1]<-0
              cen[k]<-n.hat[k] - d[k]
            }
          }
          if (d[k] != 0) last<-k
        }
        sumd<- sum(d[1:upper[n.int]])
      }
    }
  }
  
  
  ### Now form IPD ###
  #Initialise vectors
  d[d<0] <- 0       ## ADD 
  t.IPD<-rep(t.S[n.t],n.risk[1])
  event.IPD<-rep(0,n.risk[1])
  #Write event time and event indicator (=1) for each event, as separate row in t.IPD and event.IPD
  k=1
  for (j in 1:n.t){
    if(d[j]!=0){
      t.IPD[k:(k+d[j]-1)]<- rep(t.S[j],d[j])
      event.IPD[k:(k+d[j]-1)]<- rep(1,d[j])
      k<-k+d[j]
    }
  }
  #Write censor time and event indicator (=0) for each censor, as separate row in t.IPD and event.IPD
  for (j in 1:(n.t-1)){
    if(cen[j]!=0){
      t.IPD[k:(k+cen[j]-1)]<- rep(((t.S[j]+t.S[j+1])/2),cen[j])
      event.IPD[k:(k+cen[j]-1)]<- rep(0,cen[j])
      k<-k+cen[j]
    }
  }
  
  #Output IPD for more than 2 groups
  IPD<-rbind(IPD, matrix(c(t.IPD,event.IPD,rep(gr.number,length(t.IPD))),ncol=3,byrow=F))
  
  
  # add column names
  colnames(IPD) <- c("Surv.Time","Censor","Group") 
  # Find Kaplan-Meier estimates
  IPD <- as.data.frame(IPD)
  
  IPD$Surv.Time <- as.numeric(as.character(IPD$Surv.Time))
  IPD$Censor <- as.numeric(as.character(IPD$Censor))
  
  return(IPD)
}



