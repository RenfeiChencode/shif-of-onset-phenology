library(ggplot2)
library(scales)
library(patchwork)
library(pracma)

set.seed(123) 
t1<- rnorm(50, mean = 100, sd = 10)
t2<- rnorm(50, mean = 150, sd = 10)
t0_bi=c(t1,t2)
t0_uni<- rnorm(100, mean = 150, sd = 10)
alpha=1;beta=1;tem=290;k=8.62*10^(-5);E=0.65;cycle=10
G=alpha*tem/(1+beta*tem)

#######################
### initial no grazing: bimodal distribution
df=data.frame(class=rep(1,length(t0_bi)),flower=t0_bi)
nograze_bi=ggplot(df,aes(x = t0_bi,fill = class)) + geom_density(alpha=0.5)+
  labs(title = "No grazed \n Bimodal distribution", x = "",y="Probability denstiy")  +
  annotate("text",x=0, y=0.019, label="A",angle = 0,size=15,family="serif")+
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.text.y=element_text(size=25,family="serif",angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
        axis.text.x=element_text(size=25,family="serif",angle=0,colour = "black",vjust = 0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "none", 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 23,color='black'),
        legend.title = element_text(size = 16,color=NA),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2.5,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))+
  scale_x_continuous(breaks=seq(0,200,by=50),limits=c(0,200))+
  scale_y_continuous(breaks=seq(0,0.02,by=0.02/4),limits=c(0,0.02))


### grazing at time tG; start from the initial stage with no grazing
t0=t0_bi;rh=30;tG=200 # grazing time
tf=matrix(nrow = cycle,ncol = 100);t_op=NA
for (j in 1:cycle){
for (i in 1:100){
  if(is.na(t0[i])){next #
  }else if(t0[i]<tG){
   RHO=-1*rh
   tf[j,i]=t0[i]-exp(-E/(k*tem))+RHO*G
 } else if (t0[i]==tG){RHO=0;tf[j,i]=t0[i]-exp(-E/(k*tem))+RHO*G
 }else{RHO=rh;tf[j,i]=t0[i]-exp(-E/(k*tem))+RHO*G}
  if (tf[j,i]<=20){tf[j,i]=NA}
}
t0=tf[j,];t_op=c(t_op,t0)

if (length(na.omit(t0))<=50) break
}

class=factor(rep(1:j,each=100))
df=data.frame(Generation=class,flower=t_op[-1])
df=na.omit(df)
result_sigbi=list();peakTbi <- data.frame()
for (i in 1:j){
  sub_df=subset(df,df$Generation==i)
  result_sigbi[[i]]=wilcox.test( sub_df$flower, t0_bi, exact = FALSE)
  print(median(sub_df$flower))
  pdf_obj=density(sub_df$flower)
  peak_indices=findpeaks(pdf_obj$y, npeaks=5,minpeakheight  = 0.0015, minpeakdistance = 1,sortstr=TRUE)
  peak_values <- pdf_obj$x[peak_indices[,2]]
  peak_densities <- peak_indices[,1]
  peak_df=data.frame(peak_values=peak_values,Generation=i)
  peakTbi=rbind(peakTbi,peak_df)
}
graze_bi=ggplot(df,aes(x = flower,fill = Generation)) + geom_density(alpha=0.5)+
  labs(title = "Grazed selection \n Initial bimodal distribution", x = "",y="") +
  geom_vline(xintercept = tG,size=2, col = "red", lty = "dashed")+
  annotate("text",x=0, y=0.059, label="B",angle = 0,size=15,family="serif")+
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.text.y=element_text(size=25,family="serif",angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
        axis.text.x=element_text(size=25,family="serif",angle=0,colour = "black",vjust = 0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.8,0.75), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 23,color='black'),
        legend.title = element_text(size = 20,color="black"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2.5,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))+
  scale_x_continuous(breaks=seq(0,200,by=50),limits=c(0,200))+
  scale_y_continuous(breaks=seq(0,0.06,by=0.02),limits=c(0,0.06))

#######################
### initial no grazing: unimodal distribution
df=data.frame(class=rep(1,length(t0_uni)),flower=t0_uni)
nograze_uni=ggplot(df,aes(x = t0_uni,fill = class)) + geom_density(alpha=0.5)+
  labs(title = "Unimodal distribution", x = "Onset of flowering",y="Probability denstiy")  +
  annotate("text",x=45, y=0.0448, label="C",angle = 0,size=15,family="serif")+
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.text.y=element_text(size=25,family="serif",angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
        axis.text.x=element_text(size=25,family="serif",angle=0,colour = "black",vjust = 0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "none", 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 23,color='black'),
        legend.title = element_text(size = 16,color=NA),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2.5,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))+
  scale_x_continuous(breaks=seq(45,245,by=50),limits=c(45,245))+
  scale_y_continuous(breaks=seq(0,0.045,by=0.045/3),limits=c(0,0.045))


### grazing at time tG; start from the initial stage with no grazing
t0=t0_uni;rh=10;tG=145;cycle=4 # grazing time
tf=matrix(nrow = cycle,ncol = 100);t_op=NA
for (j in 1:cycle){
  for (i in 1:100){
    if(is.na(t0[i])){next 
    }else if(t0[i]<tG){
      RHO=-1*rh
      tf[j,i]=t0[i]-exp(-E/(k*tem))+RHO*G
    } else if (t0[i]==tG){RHO=0;tf[j,i]=t0[i]-exp(-E/(k*tem))+RHO*G
    }else{RHO=rh;tf[j,i]=t0[i]-exp(-E/(k*tem))+RHO*G}
    if (tf[j,i]<=20){tf[j,i]=NA}
  }
  t0=tf[j,];t_op=c(t_op,t0)
  
  if (length(na.omit(t0))<=50) break
}

class=factor(rep(1:j,each=100))
df=data.frame(Generation=class,flower=t_op[-1])
df=na.omit(df)
result_siguni=list();peakTuni <- data.frame()
for (i in 1:j){
  sub_df=subset(df,df$Generation==i)
  result_siguni[[i]]=wilcox.test( sub_df$flower, t0_uni, exact = FALSE)
  print(median(sub_df$flower))
  pdf_obj=density(sub_df$flower)
  peak_indices=findpeaks(pdf_obj$y, npeaks=5,minpeakheight  = 0.0015, minpeakdistance = 1,sortstr=TRUE)
  peak_values <- pdf_obj$x[peak_indices[,2]]
  peak_densities <- peak_indices[,1]
  peak_df=data.frame(peak_values=peak_values,Generation=i)
  peakTuni=rbind(peakTuni,peak_df)
  }
graze_uni=ggplot(df,aes(x = flower,fill = Generation)) + geom_density(alpha=0.5)+
  labs(title = "Initial unimodal distribution", x = "Onset of flowering",y="") +
  geom_vline(xintercept = tG,size=2, col = "red", lty = "dashed")+
  annotate("text",x=45, y=0.039, label="D",angle = 0,size=15,family="serif")+
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.text.y=element_text(size=25,family="serif",angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
        axis.text.x=element_text(size=25,family="serif",angle=0,colour = "black",vjust = 0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.8,0.75), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 23,color='black'),
        legend.title = element_text(size = 20,color="black"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2.5,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))+
  scale_x_continuous(breaks=seq(45,245,by=50),limits=c(45,245))+
  scale_y_continuous(breaks=seq(0,0.04,by=0.01),limits=c(0,0.04))
OUTFIG=(nograze_bi|graze_bi)/(nograze_uni|graze_uni)
OUTFIG
#ggsave("simulation3.pdf",OUTFIG,width = 45, height = 40, units = "cm", dpi = 300) 
result_sigbi
result_siguni
peakTbi

peakTuni
