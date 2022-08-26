rm(list=ls())

## load libraries
{
library(shiny)
library(shinythemes)
library(shinyjs)
library(shinyWidgets)
library(vroom)
library(ggplot2)
library(ggiraph)
library(stringr)
library(data.table)
library(fontawesome)
library(dplyr)
library(reactlog)
library(segmented)
library(grid)
library(gridExtra)
library(zoo)
library(foreach)
library(doParallel)
library(purrr)
library(pracma)
library(forecast)
library(RColorBrewer)
#
library(tictoc)
}
#

##RGC
setwd("C:/Users/Sarah/Documents/Research/Littorinid HR")
master<-fread("master (rgc2020).csv",data.table=T)

##declare global variables
input <- output <- rv <- list()

##step 1
{ #input
  input$file_input$datapath<-"C:\\Users\\Sarah\\Documents\\Research\\Littorinid HR\\Raw csv\\rgc2020\\20220604D.zip"
}
{ #backend
  setwd("~/")
  rv$dname<-strsplit(unzip(input$file_input$datapath, list=T)[1, "Name"], split="/")[[1]][1] #get the name of the zipped folder
  unzip(input$file_input$datapath) #extract files from the zipped folder
  rv$folds<-str_sort(list.dirs(paste0("./", rv$dname), full.names=F)[-1], numeric=T) #names of directories
  rv$nfold<-length(list.dirs(paste0("./", rv$dname))[-1]) #number of directories
  rv$data_str<-data.table(files=rv$folds,
                          pages=sapply(as.list(paste0("./", rv$dname, "/", rv$folds)),
                                       function(x) length(list.files(x, full.name=F)))) #create a data.table displaying the data structure
  for (tempn in list.files(paste0("./", rv$dname, "/", rv$folds[1]), full.name=T)){
    temp<-suppressWarnings(fread(tempn, skip=1L))
    ifelse(nrow(temp)==0,next,{rv$ch_names<-suppressWarnings(colnames(fread(tempn)))#get column names
                               rm(tempn)
                               break})
  }#get the first csv as template (loop until not empty)
}
{ #output
  print(list(
    paste("Zipped file name: ",rv$dname,sep=""),
    paste("Number of files: ",rv$nfold,sep=""),
    rv$data_str,
    paste("Total duration: ",round(temp[nrow(temp),1]*sum(rv$data_str$pages)/60,digits=0)," mins",sep=""),
    paste("Number of channels: ",length(rv$ch_names)-1,sep=""),
    paste("Names of channels: ",paste(rv$ch_names[-1],collapse=", "),sep="")
  )) #print data structure for data preview
}

##step 2
{ #input
  input$reduce_res<-0.01 #default 0.01, can be 0.1 and 1
  input$an_in<-1 #analysis interval (1 min/5 min)
  input$ch_selected<- rv$ch_names[-1] #allow users to choose which channel to be analyzed
  input$acf_thres<-0.5 #acf threshold for analysis (can remove?)
}
{ #backend
  rv$raw_res<-as.numeric(temp[2,1])
  closest_interval<-which.min(unlist(abs(temp[, 1] - input$reduce_res))) - 1L #interval that is the closest to the desired resolution
  rv$actual_new_res<-unlist(temp[closest_interval + 1L, 1]) #new resolution
  rv$rawmaster<-list() #initialize master data.table
  for(f in 1:rv$nfold){
    ##for each directory
    rv$rawmaster<-c(rv$rawmaster, lapply(as.list(1:rv$data_str$pages[f]),
                            function(p){
                              tryCatch({
                                fread(list.files(paste0("./", rv$dname, "/", rv$folds[f]), full.name=T)[p], skip=1)
                              },warning=function(w){
                                message(paste("Empty csv file was detected: ",rv$folds[f],sprintf("_%02.0f.csv",p),sep=""))
                                return(data.table(matrix(0,ncol=ncol(temp),nrow=nrow(temp))))
                              }
                              )
                            })) #read csvs and add to list
    print(sprintf("Reading data: %.0f%%", f/rv$nfold*100)) #update progress bar
  }
  
  print("Finalizing...")
  rv$rawmaster<-rbindlist(rv$rawmaster) #bind data.tables into one master data.table  #very slow***
  
  names(rv$rawmaster)<-rv$ch_names #assign column names
  rv$rawmaster<-rv$rawmaster[which(!is.na(rv$rawmaster[, 1])), ] #remove extra rows
  rv$rawmaster<-rv$rawmaster[, colnames(rv$rawmaster) := lapply(.SD, as.numeric)] #convert data into numeric class #slow*
  rv$rawmaster[, 1] <- seq(0, by=rv$raw_res, length.out=nrow(rv$rawmaster)) #re-assign time
  
  rv$mindex<-seq(1,nrow(rv$rawmaster),by=closest_interval)
  rv$master<-rv$rawmaster[rv$mindex,]
  unlink(paste0("./", rv$dname), recursive=T) #remove extracted files from the working directory
}


##step 3
{ #input
}
{ #global
  
  fitness<-function(pop, segment, thres){ #determine whether populations are viable in the GA
    return(apply(pop, 1, function(x){
      ifelse(sum(segment[x["s"]:x["e"]],na.rm=T)==0,1,{
        acfs<-acf(segment[x["s"]:x["e"]], lag.max=x["e"] - x["s"], plot=F, na.action=na.pass)$acf
        ifelse(max(acfs[-c(1:(which(sign(acfs)==-1)[1]))], na.rm=T) > thres,
             x["e"] - x["s"], 1)})
      })) #return window width if max acf >= threshold, else return 1 (minimal window width)
  } 
  
  mutate<-function(pop, start, end, step){ #mutate gene e/s for each population
    return(
      data.table(
        t(apply(pop, 1, function(x){
          direction<-runif(1) - 0.5 * (x["s"] == start) + 0.5 * (x["e"] == end) #determine which gene (which can still be mutated) to mutate
          if(direction > 0.5){ #if mutate gene s
            return(c(s=as.integer(round(runif(1, start, max(start, x["s"] - step)))), x["e"], x["p"], x["f"])) #mutate gene s
          } else { #if mutate gene e
            return(c(x["s"], e=as.integer(round(runif(1, min(end, x["e"] + step), end))), x["p"], x["f"])) #mutate gene e
          }
        }))
      )
    ) #return mutated e or s
  }
  
  acffunc<-function(final, segment, thres){  #to identify wavelength (in indices) of recurrent beats by autocorrelation
    return(apply(final,1,function(x){
      acfout<-acf(segment[x["s"]:x["e"]],lag.max=as.integer(x["f"]), plot=F, na.action=na.pass)
      acfdt<-data.table(lag=acfout[["lag"]],
                        ACF=acfout[["acf"]]) #pair lag with acf
      clim<-qnorm((1 + 0.95)/2)/sqrt(acfout$n.used) #upper limit of 95% confidence interval
      acfdt<-acfdt[-c(1:(which(sign(acfdt$ACF)==-1)[1])),]
      allpeaks<-arrange(data.frame(findpeaks(acfdt$ACF,minpeakheight=clim,minpeakdistance=5))[,1:2],X2)
      localmax<-allpeaks[1:which.max(allpeaks[,1]),]
      return(data.table(ACF=localmax[,1],
                         lag=acfdt$lag[localmax[,2]],
                         hr=60/rv$actual_new_res/acfdt$lag[localmax[,2]]))
      }))
    #return a list of local maxima
  }
  
  ga<-function(master, hrdf, initial_pop, pop_size, max_gen, mutation_step, thres, patience){
    #genetic algorithm for calculating the heart rate of each 
    segment<-master[hrdf["start_int"]:hrdf["end_int"]] #extract the segment of signal
    initial_pop$f<-fitness(initial_pop, segment, thres) #initialize population and calculate fitness
    pop<-initial_pop #initialize population
    hof<-pop[which(pop$f > 1), ] #initialize hall of fame as the fit-enough individuals from the initial population
    wait<-0 #patience counter
    for(g in 2:max_gen){
      ##for a number of generations
      last_pop<-pop #save population (to be compared with the subsequent generation)
      pop<-pop[sample.int(pop_size, replace=T, prob=pop$f), ] #select individuals for reproducing the next generation (chance proportional to fitness)
      pop<-mutate(pop, 1, length(segment), mutation_step) #mutate the new generation
      pop$f<-fitness(pop, segment, thres) #calculate the fitness of each individual
      hof<-rbind(hof, pop[which(pop$f > 1), ]) #store fit-enough individuals into hall of fame
      wait<-ifelse(any(sapply(as.list(unique(pop$p)), function(x){
        ##compare the new generation with the last generation (along each lineage)
        max(pop$f[which(pop$p == x)]) > max(last_pop$f[which(last_pop$p == x)])
      })), 0, wait + 1) #if there are increments in at least one lineage, reset patience counter, else patience counter + 1
      if(any(pop$f == (length(segment)- 1L)) | wait == patience) break #break loop if the patience counter has reached the tolerance level
      if(all(pop$f == 1)) pop<-last_pop #if all individuals are unfit, revert to the previous population
    }
    
    final<-data.table(s=integer(0), e=integer(0), p=integer(0), f=integer(0)) #initialize data.tables to store outputs
    while(nrow(hof) > 0){
      ##while there are still individuals in hall of fame (not stored/removed)
      final<-rbind(hof[which.max(hof$f), ], final) #store the longest individual in final
      hof<-hof[-which(hof$s >= final$s[1] & hof$e <= final$e[1]), ] #remove individuals that are completely overlapped by the stored individual
    }
    if(nrow(final)==0){
      final<-data.table(s=NA,e=NA, p=NA,f=NA)
      counts<-list(data.table(ACF=NA,lag=NA,hr=NA))
    }else{
      counts<-acffunc(final,segment,thres)
    }
    return(list(final=final,
                counts=counts))
  }
  
  weight<-function(dt){ #calculate hr counts by weighting
    if(nrow(na.omit(dt))>1){ #multiple windows
      seqtally<-data.table(table(unlist(apply(na.omit(dt),1,function(x) seq(x[["s"]],x[["e"]],1))))) #count the occurrence
      if (length(which(seqtally$N>1))==0){ #no overlapping
        whr<-sum(apply(dt,1,function(x) x[["hr"]]*x[["f"]]),na.rm=T)/sum(dt[["f"]],na.rm=T) #weighted hr
      }else{ #with overlapping
        diffdt<-rle(diff(seqtally$N))
        diffdt<-data.table(lengths=diffdt$lengths,
                           values=diffdt$values) #find consecutive segments
        segm<-data.table(s=numeric(0),
                         e=numeric(0),
                         f=numeric(0),
                         wins=numeric(0),
                         hr=numeric(0))
        sx<-min(dt[["s"]],na.rm=T)
        for (j in c(which(diffdt$values!=0),NA)){ #manually introduced NA for the last segment
          ex<-ifelse(is.na(j),
                     max(dt[["e"]],na.rm=T), #for the last segment
                     as.numeric(seqtally$V1[sum(diffdt$lengths[1:(j)])]))
          
          #to discriminate back-to-back / non-overlaping segments from continuous segments
          sw<-which(apply(dt,1,function(x) sx>=x[["s"]] & sx<=x[["e"]])) #windows encompassing the start
          ew<-which(apply(dt,1,function(x) ex>=x[["s"]] & ex<=x[["e"]])) #windows encompassing the end
          if (identical(sw,ew)){ #continuous segments
            segm<-rbind(segm,data.table(s=sx,
                                        e=ex,
                                        f=ex-sx,
                                        wins=list(sw),
                                        hr=mean(dt[["hr"]][sw],na.rm=T)))
            sx<-ex+1
          }else{ #back-to-back segments
            sx<-sort(dt[["s"]][which(dt[["s"]]<ex & dt[["s"]]>=sx)])
            ex<-c(sx[-1]-1,ex)
            wins<-lapply(sx,function (y) which(apply(dt,1,function(x) y>=x[["s"]] & y<=x[["e"]])))
            segm<-rbind(segm,data.table(s=sx,
                                        e=ex,
                                        f=ex-sx,
                                        wins=wins,
                                        hr=lapply(wins,function(x) mean(dt[["hr"]][x],na.rm=T))))
            sx<-max(ex)+1
          }
        }
        whr<-sum(apply(segm,1,function(x) x[["hr"]]*x[["f"]]),na.rm=T)/sum(segm[["f"]],na.rm=T) #weighted hr
      }
      return(data.table(wACF=mean(dt[["ACF"]],na.rm=T),
                        whr=whr))
    }else{ 
      if(nrow(na.omit(dt))==1){ #only one countable window
        return(data.table(wACF=na.omit(dt[["ACF"]]),
                          whr=na.omit(dt[["hr"]])))
      }else{ #no countable windows
        return(data.table(wACF=NA,
                          whr=NA))
      }
    }
  }
  
  prelim<-function(out){ #calculate HR based on ACF only (selecting the lag with max ACF)
    prelimdt<-cbind(ix=rep(1:length(out[["final"]]),lapply(out[["final"]],nrow)),
                    win=unlist(lapply(lapply(out[["final"]],nrow),function(x) seq(1,x,1))),
                    rbindlist(out[["final"]]))
    
    prelimdt<-cbind(prelimdt,
                    rbindlist(lapply(out[["counts"]], function(i) rbindlist(lapply(i, function(w) if(is.na(w[1,hr])){w[1,]}
                                                                                   else{w %>% top_n(1,ACF)})))))
    return(prelimdt)
  }
  
  screen<-function(out){ #screen through the counts with a tracking record (past 5 intervals) through time
    finaldt<-cbind(ix=rep(1:length(out[["final"]]),lapply(out[["final"]],nrow)),
                   win=unlist(lapply(lapply(out[["final"]],nrow),function(x) seq(1,x,1))),
                   rbindlist(out[["final"]]),
                   ACF=as.numeric(NA),lag=as.numeric(NA),hr=as.numeric(NA))

    for (i in 1:length(out[["counts"]])){
      prev<-finaldt[ix %in% (i-5):(i-1),.(ix,hr)]
      est<-ifelse(nrow(na.omit(prev))==0,
                  NA,
                  ifelse(length(unique(na.omit(prev)[["ix"]]))>2 & summary(lm(hr~ix,prev))$r.squared>=0.7,
                         predict(lm(hr~ix,prev),data.frame(ix=i)),
                         median(prev$hr,na.rm=T)))
      for (w in 1:length(out[["counts"]][[i]])){
        if(is.na(out[["counts"]][[i]][[w]][1,hr])){ #NA count
          finaldt[ix==i & win==w,names(finaldt)[3:9]:=as.list(rep(NA,7))]
        }else{
          ifelse(is.na(est),
                 finaldt[ix==i & win==w,c("ACF","lag","hr"):=as.list(out[["counts"]][[i]][[w]] %>%
                                                                       top_n(1,ACF))], #return counts with max acf if no reference from previous counts
                 ifelse(abs(as.numeric(out[["counts"]][[i]][[w]][,.(hr)] %>% filter(abs(hr-est)==min(abs(hr-est))))-est)/est>0.4,
                        finaldt[ix==i & win==w,names(finaldt)[3:9]:=as.list(rep(NA,7))],
                        finaldt[ix==i & win==w,c("ACF","lag","hr"):=as.list(out[["counts"]][[i]][[w]] %>%
                                                                              filter(abs(hr-est)==min(abs(hr-est))))])) #return counts closest to previous counts
          
        }
        
        
      }
    }
    return(finaldt)
  }

  findpeaks<-function (x, nups = 1, ndowns = nups, zero = "+", peakpat = NULL,
                       minpeakheight = -Inf, minpeakdistance = 1, threshold = 0, 
                       npeaks = 0, sortstr = FALSE) {
    stopifnot(is.vector(x, mode = "numeric") || length(is.na(x)) == 
                0)
    if (!zero %in% c("0", "+", "-")) 
      stop("Argument 'zero' can only be '0', '+', or '-'.")
    xc <- paste(as.character(sign(diff(x))), collapse = "")
    xc <- gsub("1", "+", gsub("-1", "-", 
                              xc))
    if (zero != "0") 
      xc <- gsub("0", zero, xc)
    if (is.null(peakpat)) {
      peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
    }
    rc <- gregexpr(peakpat, xc)[[1]]
    if (rc[1] < 0) 
      return(NULL)
    x1 <- rc
    x2 <- rc + attr(rc, "match.length")
    attributes(x1) <- NULL
    attributes(x2) <- NULL
    n <- length(x1)
    xv <- xp <- numeric(n)
    for (i in 1:n) {
      tryCatch({
        xp[i] <- which.max(x[x1[i]:x2[i]]) + x1[i] - 1
        xv[i] <- x[xp[i]]
      },
      error=function(e){})
    }
    inds <- which(xv >= minpeakheight & xv - pmax(x[x1], x[x2]) >= 
                    threshold)
    X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])
    if (minpeakdistance < 1) 
      warning("Handling 'minpeakdistance < 1' is logically not possible.")
    if (sortstr || minpeakdistance > 1) {
      sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
      X <- X[sl, , drop = FALSE]
    }
    if (length(X) == 0) 
      return(c())
    if (minpeakdistance > 1) {
      no_peaks <- nrow(X)
      badpeaks <- rep(FALSE, no_peaks)
      for (i in 1:no_peaks) {
        ipos <- X[i, 2]
        if (!badpeaks[i]) {
          dpos <- abs(ipos - X[, 2])
          badpeaks <- badpeaks | (dpos > 0 & dpos < minpeakdistance)
        }
      }
      X <- X[!badpeaks, , drop = FALSE]
    }
    if (npeaks > 0 && npeaks < nrow(X)) {
      X <- X[1:npeaks, , drop = FALSE]
    }
    return(X)
  }
  
  reacf<-function(dt,segment){  #re-run acf at finest resolution if resolution is coarser than 2% of the bpm count
    acfout<-acf(segment,lag.max=as.integer((dt$lag+5)*round(rv$actual_new_res/rv$raw_res)), plot=F, na.action=na.pass)
    acfdt<-data.table(lag=acfout[["lag"]],
                      ACF=acfout[["acf"]]) #pair lag with acf
    clim<-qnorm((1 + 0.95)/2)/sqrt(acfout$n.used) #upper limit of 95% confidence interval
    acfdt<-acfdt[-c(1:(which(sign(acfdt$ACF)==-1)[1])),]
    allpeaks<-arrange(data.frame(findpeaks(acfdt$ACF,minpeakheight=clim,minpeakdistance=5))[,1:2],X2)
    i<-which.min(abs(apply(allpeaks,1,function(x) 60/rv$raw_res/acfdt$lag[x[2]])-dt$hr))
    return(data.table(ACF=allpeaks[i,1],
                      lag=acfdt$lag[allpeaks[i,2]],
                      hr=60/rv$raw_res/acfdt$lag[allpeaks[i,2]]))
  }

  pp<-function(dt){ 
    palette<-brewer.pal(n=11,name="RdYlBu")
    names(palette)<-seq(0,10,1)/10
    theme_set(theme_linedraw())
    theme_update(text=element_text(size=12,colour="#3D405B"),
                 axis.title=element_text(size=12,colour="#3D405B"),
                 axis.text=element_text(size=12,colour="#3D405B"),
                 axis.ticks=element_line(color="#3D405B",size=0.5,lineend="square"),
                 panel.grid.minor=element_blank(),
                 panel.grid.major=element_line(color="#3D405B",size=0.5,linetype="dotted"),
                 panel.border=element_rect(fill=NA,colour="#3D405B",size=1),
                 #axis.line=element_line(color="#3D405B",size=0.5),
                 panel.background=element_blank(),
                 legend.background=element_blank(),
                 legend.key=element_blank(),
                 plot.background=element_blank())
    return(as.list(
      ggplot()+
        geom_point(data=dt,aes(x=ix,y=hr,fill=factor(floor(ACF*10)/10)),colour="black",shape=21,size=4,alpha=.8)+
        scale_fill_manual(name="ACF",values=palette,na.value=NA,na.translate=F,guide="legend")+
        scale_y_continuous(name="Heart rate (bpm)",
                           breaks=seq(0,ifelse(sum(!is.na(dt[,hr]))>0,ceiling(max(dt[,hr],na.rm=T)/50)*50,50),50),
                           limits=c(0,ifelse(sum(!is.na(dt[,hr]))>0,ceiling(max(dt[,hr],na.rm=T)/50)*50,50)),
                           expand=c(0,0))+
        scale_x_continuous(name="Time (min)",
                           breaks=seq(0,ceiling(max(dt[,ix],na.rm=T)/30)*30,30),
                           limits=c(0,ceiling(max(dt[,ix],na.rm=T)/30)*30),
                           expand=c(0,0))
    ))
  }
  
}
{ #backend
  rv$pop_size=10L #declare population size
  rv$max_gen=20L #declare maximum no. of generations
  rv$an_index_in<-which.min(abs(rv$master[1:(floor(60*input$an_in/rv$master[2,Time])+200),Time]-(60*input$an_in)))-1L
  rv$an_index_length<-floor(rv$master[nrow(rv$master),Time]/(60*input$an_in))
  rv$windex<-data.table(start_int=seq(1,by=rv$an_index_in+1L,length.out=rv$an_index_length),
                        end_int=seq(rv$an_index_in+1L,by=rv$an_index_in+1L,length.out=rv$an_index_length),
                        start_min=seq(0,by=input$an_in,length.out=rv$an_index_length),
                        end_min=seq(input$an_in,length.out=rv$an_index_length))
  rv$wbound <- as.integer(seq(1,rv$an_index_in,length.out=rv$pop_size+1L))
  rv$initial_pop<-data.table(
    s=rv$wbound[1:rv$pop_size],
    e=rv$wbound[2:(rv$pop_size + 1L)] - 1L,
    p=1:rv$pop_size
  ) #s(starting index) and e(ending index) are genes, p is the lineage
  rv$mutation_step<-as.integer(ceiling(rv$an_index_in*(1 - 1/rv$pop_size)/rv$max_gen)) #determine the minimal step for mutation
  rv$thres_lag<-as.numeric(which(sapply(seq(1,rv$an_index_in,1),function(x) ((60/rv$actual_new_res/x)-(60/rv$actual_new_res/(x+1)))/(60/rv$actual_new_res/x))<=0.02)[1]) #lag threshold below which the resolution exceeds 2% of the calculated bpm
  
  for(channel in input$ch_selected){
    ##calculate heart rate for each channel
    print(sprintf("Calculating heart rate: %s...", channel)) #indicate which channel is being analyzed
    progress_counter<-0 #reset progress counter
    ncore<-detectCores() - 1L
    cl<-makeCluster(ncore) #parallelize
    registerDoParallel(cl) #parallelize
    
    #//
    out<-transpose(foreach(ti=1:nrow(rv$windex), .packages=c("data.table","pracma","dplyr")) %dopar% {
      ga(master=rv$master[[channel]], hrdf=unlist(rv$windex[ti, ]),
         initial_pop=rv$initial_pop, pop_size=rv$pop_size, max_gen=rv$max_gen,
         mutation_step=rv$mutation_step, thres=input$acf_thres, patience=2L)
    }) #use genetic algorithm to exclude the noisy sections and calculate heart rate
    
    stopCluster(cl) #de-parallelize
    stopImplicitCluster() #de-parallelize
    
    prelimdt<-prelim(out) 
    presults<-cbind(data.table(ix=unique(prelimdt[["ix"]])),
                    rbindlist(lapply(unique(prelimdt[["ix"]]),function(x) weight(prelimdt[ix==x]))))
    
    finaldt<-cbind(screen(out),data.table(res=rv$actual_new_res))
    for(i in which(finaldt$lag<rv$thres_lag)) { #re-run acf at finest resolution if resolution is coarser than 2% of the bpm count
      finaldt[i,7:10]<-cbind(reacf(dt=finaldt[i,],
                                      segment=rv$rawmaster[[channel]][rv$mindex[rv$windex[finaldt[i,ix],start_int]+finaldt[i,s]-1]:rv$mindex[rv$windex[finaldt[i,ix],start_int]+finaldt[i,e]-1]]),
                             res=rv$raw_res)
    }
    results<-cbind(data.table(ix=unique(finaldt[["ix"]])),
                   rbindlist(lapply(unique(finaldt[["ix"]]),function(x) weight(finaldt[ix==x]))))
    
    ind<-as.integer(master[Date==as.integer(substr(rv$dname,1,8)) & 
                             Batch==substr(rv$dname,nchar(rv$dname),nchar(rv$dname)) &
                             #Batch=="A" &
                             Channel==substr(channel,nchar(channel),nchar(channel)),Ind])
    setwd("C:/Users/Sarah/Documents/Research/Littorinid HR/Algo counts/Results v2.4")
    
    write.table(prelimdt,paste(sprintf("pwindow_ind%03.0f",ind),".csv",sep=""),row.names=F,col.names=T,sep=",")
    write.table(presults,paste(sprintf("pcount_ind%03.0f",ind),".csv",sep=""),row.names=F,col.names=T,sep=",")
    write.table(finaldt,paste(sprintf("swindow_ind%03.0f",ind),".csv",sep=""),row.names=F,col.names=T,sep=",")
    write.table(results,paste(sprintf("scount_ind%03.0f",ind),".csv",sep=""),row.names=F,col.names=T,sep=",")
    
    if(sum(!is.na(prelimdt[,hr]))>0){
      rv$pp[[channel]]<-c(rv$pp[[channel]],
                          list(pwindow=pp(prelimdt),
                               ptime=pp(presults[,.(ix=ix,ACF=wACF,hr=whr)]),
                               swindow=pp(finaldt),
                               stime=pp(results[,.(ix=ix,ACF=wACF,hr=whr)])))
      png(file=sprintf("ind%03.0f.png",ind),width=4800,height=2400,units="px",res=300,bg="transparent")
      grid.arrange(grobs=rv$pp[[channel]],
                   nrow=2,ncol=2,
                   top=textGrob(c("all windows","per minute",sprintf("ind%03.0f",ind)),x=unit(c(0.25,0.75,0.99),"npc"),hjust=c(0.5,0.5,1),gp=gpar(fontsize=c(12,12,18))))
      dev.off()
    }else{
      next
    }
  }
}
