#' @title Extract Data from KM Curve When Number at Risk Table Is Available
#' 
#' @description This function derives or extracts IPD when only number at risk table information are available
#' @param t.risk.T vector of time points from the number at risk table
#' @param n.risk.T vector of patients' number at risk at time points from the table
#' @param group.name name of treatment group
#' @param data data frame with 2 columns; 1st column has digitised time information and 2nd column has digitised survival rate information
#' @export
#' @return data frame with time to event, censoring, and group name information
#' @references 
#' \insertRef{guyot2012enhanced}{extractKM}
#' @examples
#' df <- ex.event[,1:2]
#' df <- df[!is.na(df$x1) & !is.na(df$y1),]  # only valid data
#' res1 <- km.table(t.risk.T = c(0,25,50,75,100,125),
#'                  n.risk.T = c(22,18,10,7,5,2),
#'                  group.name="group1",
#'                  data=df);res1
#' 
#' df <- ex.event[,3:4]
#' df <- df[!is.na(df$x2) & !is.na(df$y2),]
#' res2 <- km.table(t.risk.T = c(0,25,50,75),
#'                  n.risk.T = c(5,1,1,0),
#'                  group.name="group2",
#'                  data=df);res2
#' 
#' library(survival)
#' result <- rbind(res1,res2); 
#' km <- survfit(Surv(Surv.Time, Censor) ~ Group, data=result);km
#' plot(km)


km.table <- function(t.risk.T, n.risk.T, group.name="group", data){
  
  df <- extractKM::data.org(data)
  
  # To count the points between time intervals -------------------------------
  tm.point <- t.risk.T[-1]
  time <- df[,1]
  tm.point <- tm.point[tm.point < time[length(time)-1]]    # Select time points that has information

    num.point <- length(tm.point)                            # Number of time points
  
  point.min <- rep(NA, num.point)
  point.max <- rep(NA, num.point)
  for (i in 1:num.point){
    point.min[i] <- min(which(time>tm.point[i]))
    point.max[i] <- max(which(time<tm.point[i]))
  }
  
  lower.T <- c(1,point.min)               # The lower bound of clicks
  upper.T <- c(point.max,length(time))    # The upper bound of clicks
  
  t.risk.T <- t.risk.T[1:(num.point+1)]  # Adjust number at risk table according to information
  n.risk.T <- n.risk.T[1:(num.point+1)]
  
  n.t.T	  <- as.numeric(length(time))       # number of clicks for each group
  n.int 	<- length(n.risk.T)			          # number of intervals
  
  res  <- extractKM::km.cal.tab(t.points=df[,1], s.points=df[,2], 
                                t.risk.T=t.risk.T, 
                                n.risk.T=n.risk.T,
                                lower.T=lower.T, upper.T=upper.T,
                                gr.number=group.name)
  return(res)
}





#' @title Extract Data from KM Curve When Both Number at Risk Table and Number of Events Are Available
#' 
#' @description This function derives or extracts IPD when both number at risk table and total number of events information are available
#' @param t.risk.T vector of time points from the number at risk table
#' @param n.risk.T vector of patients' number at risk at time points from the table
#' @param t.event numeric value; total number of events
#' @param group.name name of treatment group
#' @param data data frame with 2 columns; 1st column has digitised time information and 2nd column has digitised survival rate information
#' @export
#' @return data frame with time to event, censoring, and group name information
#' @references 
#' \insertRef{guyot2012enhanced}{extractKM}
#' @examples
#' df <- ex.tab.event[,1:2]
#' df <- df[!is.na(df$x1) & !is.na(df$y1),]  # only valid data
#' res1 <- km.tab.event(t.risk.T = c(0,24,48,72,96,120), 
#'                      n.risk.T = c(23,11,8,4,2,0), 
#'                      t.event=5, 
#'                      group.name="group1",
#'                      data=df);res1
#' 
#' df <- ex.tab.event[,3:4]
#' df <- df[!is.na(df$x2) & !is.na(df$y2),]
#' res2 <- km.tab.event(t.risk.T = c(0,24,48,72,96,120),
#'                      n.risk.T = c(12,3,3,2,1,1),
#'                      t.event=8,
#'                      group.name="group2",
#'                      data=df);res2
#' 
#' library(survival)
#' result <- rbind(res1,res2); 
#' km <- survfit(Surv(Surv.Time, Censor) ~ Group, data=result);km
#' plot(km)


km.tab.event <- function(t.risk.T, n.risk.T, t.event="NA", group.name="group", data){
  
  df <- extractKM::data.org(data)
  
  # To count the points between time intervals -------------------------------
  tm.point <- t.risk.T[-1]
  time <- df[,1]
  tm.point <- tm.point[tm.point < time[length(time)-1]]    # Select time points that has information

    num.point <- length(tm.point)                            # Number of time points
  
  point.min <- rep(NA, num.point)
  point.max <- rep(NA, num.point)
  for (i in 1:num.point){
    point.min[i] <- min(which(time>tm.point[i]))
    point.max[i] <- max(which(time<tm.point[i]))
  }
  
  lower.T <- c(1,point.min)               # The lower bound of clicks
  upper.T <- c(point.max,length(time))    # The upper bound of clicks
  
  t.risk.T <- t.risk.T[1:(num.point+1)]  # Adjust number at risk table according to information
  n.risk.T <- n.risk.T[1:(num.point+1)]
  
  n.t.T	  <- as.numeric(length(time))       # number of clicks for each group
  n.int 	<- length(n.risk.T)			          # number of intervals
  
  res  <- extractKM::km.cal.tab(t.points=df[,1], s.points=df[,2], 
                                t.risk.T=t.risk.T, 
                                n.risk.T=n.risk.T,
                                t.event=t.event,
                                lower.T=lower.T, upper.T=upper.T,
                                gr.number=group.name)
  return(res)
}







#' @title Extract Data from KM Curve When Number of Events Information Is Available
#' 
#' @description This function derives or extracts IPD when only total number of events information are available
#' @param ssize vector of patients' number at risk at time points from the table
#' @param t.event numeric value; total number of events
#' @param group.name name of treatment group
#' @param data data frame with 2 columns; 1st column has digitised time information and 2nd column has digitised survival rate information
#' @export
#' @return data frame with time to event, censoring, and group name information 
#' @references 
#' \insertRef{guyot2012enhanced}{extractKM}
#' @examples
#' df <- ex.event[,1:2]
#' df <- df[!is.na(df$x1) & !is.na(df$y1),]  # only valid data
#' res1 <- km.event(ssize=16, t.event=13, group.name="group1", data=df); res1
#' 
#' df <- ex.event[,3:4]
#' df <- df[!is.na(df$x2) & !is.na(df$y2),]
#' res2 <- km.event(ssize=39, t.event=28, group.name="group2", data=df); res2
#' 
#' library(survival)
#' result <- rbind(res1,res2); 
#' km <- survfit(Surv(Surv.Time, Censor) ~ Group, data=result);km
#' plot(km)


km.event <- function(ssize, t.event="NA", group.name="group", data){
  
  df <- extractKM::data.org(data)

  t.S.T <- time <- df[,1]  # Argument
  S.T   <- df[,2]          # Argument
  
  t.risk.T <- 0        # inside the function              
  n.risk.T <- ssize    # Argument
  t.event  <- t.event  # Argument
  
  lower.T <- 1               # inside the function
  upper.T <- length(time)    # inside the function
  
  n.group <- 1	           # inside the function: number of groups 
  n.int 	<- 1			       # inside the function: number of intervals
  
  res  <- extractKM::km.cal.tab(t.points=df[,1], s.points=df[,2], 
                                t.risk.T=t.risk.T, 
                                n.risk.T=n.risk.T,
                                t.event=t.event,
                                lower.T=lower.T, upper.T=upper.T,
                                gr.number=group.name)
  return(res)
}





