#------------------------------------------------------------------------------#
# Results_0_Functions.R                                                        #
# Define functions to analyse and present intervention results                 #
# Last updated 2025-07-25 by KCH                                               #
#------------------------------------------------------------------------------#

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

measplot <- function(meas,covlist,durlist,algnum,plotyrs){
  plotlist <- list()
  plotnum  <- 0
  for (cov in covlist){
    for (int in durlist){
      plotnum <- plotnum + 1
      fulldatameas <- loadRData(file=paste("./Results_",prev,"/res/fulldata",meas,"_",cov,"_",int,".Rdata",sep=""))
      meassum <- array(NA,dim=c(algnum,plotyrs,3))
      for (i in c(1:algnum)){ 
        for (j in c(1:plotyrs)){
          meassum[i,j,1:3] <- quantile(fulldatameas[,i,j,1],probs=c(0.5,0.025,0.975))
        }
      }
      
      noint               <- as.data.frame(meassum[1,,])
      colnames(noint)     <- c("median","lower","upper")
      noint$time          <- c(1:plotyrs)/12+2024
      
      sympxpert           <- as.data.frame(meassum[2,,])
      colnames(sympxpert) <- c("median","lower","upper")
      sympxpert$time      <- c(1:plotyrs)/12+2024
      
      #cxrxpert            <- as.data.frame(meassum[3,,])
      #colnames(cxrxpert)  <- c("median","lower","upper")
      #cxrxpert$time       <- c(1:plotyrs)/12+2024
      
      xpert               <- as.data.frame(meassum[4,,])
      colnames(xpert)     <- c("median","lower","upper")
      xpert$time          <- c(1:plotyrs)/12+2024
      
      cxr                 <- as.data.frame(meassum[5,,])
      colnames(cxr)       <- c("median","lower","upper")
      cxr$time            <- c(1:plotyrs)/12+2024
      
      if (meas=="inc")    ytitle = "sTB incidence per 100,000"
      if (meas=="mort")   ytitle = "TB mortality per 100,000"
      
      if (cov==min(covlist)){
        if (int==1){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous(paste(" \n",int," rounds",sep=""),expand=c(0,0)) +
            ggtitle(paste(cov*100,"% coverage",sep="")) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.title=element_text(size=18),axis.text.x=element_blank(),
                  plot.title=element_text(hjust=0.5,size=18),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
        if (int==2){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2026.0833,xmax=2027,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous(paste(" \n",int," rounds",sep=""),expand=c(0,0)) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.text.x=element_blank(),axis.title=element_text(size=18),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
        if (int==3){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2026.0833,xmax=2027,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2027.0833,xmax=2028,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous(paste(ytitle,"\n",int," rounds",sep=""),expand=c(0,0)) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.text.x=element_blank(),axis.title=element_text(size=18),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
        if (int==4){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) + #C7CCDB
            geom_rect(aes(xmin=2026.0833,xmax=2027,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2027.0833,xmax=2028,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2028.0833,xmax=2029,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous(paste(" \n",int," rounds",sep=""),expand=c(0,0)) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.text.x=element_blank(),axis.title=element_text(size=18),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
        if (int==5){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2026.0833,xmax=2027,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2027.0833,xmax=2028,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2028.0833,xmax=2029,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2029.0833,xmax=2030,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous(paste(" \n",int," rounds",sep=""),expand=c(0,0)) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.title=element_text(size=18),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
      }else{
        if (int==1){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous("",expand=c(0,0)) +
            ggtitle(paste(cov*100,"% coverage",sep="")) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.title.x=element_blank(),axis.text.x=element_blank(),axis.text.y = element_blank(),
                  plot.title=element_text(hjust=0.5,size=18),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
        if (int==2){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2026.0833,xmax=2027,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous("",expand=c(0,0)) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.title.x=element_blank(),axis.text.x=element_blank(),axis.text.y = element_blank(),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
        if (int==3){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2026.0833,xmax=2027,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2027.0833,xmax=2028,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous("",expand=c(0,0)) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.title.x=element_blank(),axis.text.x=element_blank(),axis.text.y = element_blank(),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
        if (int==4){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2026.0833,xmax=2027,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2027.0833,xmax=2028,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2028.0833,xmax=2029,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous("",expand=c(0,0)) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.title.x=element_blank(),axis.text.x=element_blank(),axis.text.y = element_blank(),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
        if (int==5){
          incplot <- ggplot(data=noint) +
            geom_rect(aes(xmin=2025.0833,xmax=2026,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2026.0833,xmax=2027,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2027.0833,xmax=2028,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2028.0833,xmax=2029,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_rect(aes(xmin=2029.0833,xmax=2030,ymin=0,ymax=max(upper)),fill='#cecece',alpha=.05) +
            geom_line(aes(x=time,y=median,color="1"),size=1) +
            geom_ribbon(aes(x=time,min=lower,max=upper),size=1,fill=colnoint,alpha=.15) +
            geom_line(data=xpert,aes(x=time,y=median,color="3"),size=1) +
            geom_ribbon(data=xpert,aes(x=time,min=lower,max=upper),size=1,fill=colxpert,alpha=.15) +
            #geom_line(data=cxrxpert,aes(x=time,y=median),size=1,color=colcxrxpert) +
            #geom_ribbon(data=cxrxpert,aes(x=time,min=lower,max=upper),size=1,fill=colcxrxpert,alpha=.15) +
            geom_line(data=sympxpert,aes(x=time,y=median,color="2"),size=1) +
            geom_ribbon(data=sympxpert,aes(x=time,min=lower,max=upper),size=1,fill=colsympxpert,alpha=.15) +
            geom_line(data=cxr,aes(x=time,y=median,color="4"),size=1) +
            geom_ribbon(data=cxr,aes(x=time,min=lower,max=upper),size=1,fill=colcxr,alpha=.15) +
            scale_x_continuous("Year",lim=c(2025,2035),expand=c(0,0),breaks=c(2025,2030,2035)) + 
            scale_y_continuous("",expand=c(0,0)) +
            theme_classic(base_size=18) +
            theme(plot.margin=margin(10,20,0,10)) +
            theme(axis.text.y = element_blank()) +
            theme(plot.margin=margin(10,20,0,10),
                  axis.title.x=element_blank(),axis.text.y = element_blank(),
                  legend.text=element_text(size=18),legend.title=element_text(size=18)) + 
            scale_color_manual(name = 'Algorithm', 
                               values =c("1"=colnoint,"2"=colsympxpert,"3"=colxpert,"4"=colcxr), 
                               labels = c('No intervention','Cough+Xpert','Xpert','CXR'))
        }
      }
      plotlist[[plotnum]] <- incplot
    }
  }
  plotlist
}

meascalc <- function(meas,covlist,durlist,algnum,plotyrs){
  measpct <- c()
  for (cov in covlist){
    for (int in durlist){
      meastemp <- array(0,dim=c(1000,9))
      fulldatameas <- loadRData(file=paste("./Results_",prev,"/res/fulldata",meas,"_",cov,"_",int,".Rdata",sep=""))
      meastemp[,c(1:5)] <- fulldatameas[,,144,]
      meastemp[,6] <- (meastemp[,1] - meastemp[,2])/meastemp[,1]
      #meastemp[,7] <- (meastemp[,1] - meastemp[,3])/meastemp[,1]
      meastemp[,8] <- (meastemp[,1] - meastemp[,4])/meastemp[,1]
      meastemp[,9] <- (meastemp[,1] - meastemp[,5])/meastemp[,1]
      sympxpert <- quantile(meastemp[,6],probs=c(0.5,0.025,0.975))
      #cxrxpert  <- quantile(meastemp[,7],probs=c(0.5,0.025,0.975))
      xpert     <- quantile(meastemp[,8],probs=c(0.5,0.025,0.975))
      cxr       <- quantile(meastemp[,9],probs=c(0.5,0.025,0.975))
      measpct <- rbind(measpct,c(cov,int,sympxpert,xpert,cxr))
    }
  }
  colnames(measpct) <- c("cov","int","sympxpert_med","sympxpert_low","sympxpert_high",
                         #"cxrxpert_med","cxrxpert_low","cxrxpert_high",
                         "xpert_med","xpert_low","xpert_high", "cxr_med","cxr_low","cxr_high")
  measpct
}
