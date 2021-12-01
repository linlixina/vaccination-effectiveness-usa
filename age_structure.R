if(FALSE){
require('readxl')
y<-as.data.frame(read_excel('report_34_table_2.xlsx',col_name=F))
y[,2]<-gsub(","," ",y[,2])
y[,2]<-gsub("[(]"," ",y[,2])
y[,2]<-gsub("[)]"," ",y[,2])
y<-0.01*matrix(as.numeric(unlist(strsplit(y[,2],split=" "))),ncol=6,byrow=T)[,c(1,3,5)]

x<-read.csv('covid_hk.csv')
j<-which(x$Case=='L' & x$Date.of.onset!='Asymptomatic')  # L for local
i<-which(x[,6] %in% 'Deceased' & x$Case=='L')
# we consider all local cases
# we consider all deceased, most of imported are young
# presumbly deceased cases are local. Thus it is reasonable to use local cases as the baseline
h1<-hist(x$Age[j],breaks=c(seq(0,90,by=5),120),plot=F)
h2<-hist(x$Age[i],breaks=h1$breaks,plot=F)
age<-h2$breaks[-1]-5
cfr<-h2$counts/h1$counts
print(cfr)
cfr<-as.data.frame(cbind(cfr,y))
cfr$max<-apply(cfr[,1:2],1,max)
cfr$max[1]<-cfr$max[2]
}
x<-read.csv('INFLUD-11-10-2021.csv',sep=";")
x1<-read.csv('INFLUD21-11-10-2021.csv',sep=";")
#x<-rbind(x,x1)
x<-rbind(x,x1[,match(names(x),names(x1))])
print(range(as.Date(x[,3],format="%d/%m/%Y")))
x[,3]<-as.Date(x[,3],format="%d/%m/%Y")
x[which(x$TP_IDADE==1),'NU_IDADE_N']<-x[which(x$TP_IDADE==1),'NU_IDADE_N']/365
x[which(x$TP_IDADE==2),'NU_IDADE_N']<-x[which(x$TP_IDADE==2),'NU_IDADE_N']/12
x[which(x$NU_IDADE_N>100),'NU_IDADE_N']<-100
x[which(x$NU_IDADE_N<0),'NU_IDADE_N']<-0
top<-names(sort(table(x$ID_MUNICIP),dec=T))[1:18]
x1<-x[which(x$EVOLUCAO %in% 2:3),]
pdf('SARI_top_age_death6.pdf',width=12,height=10)
par(mfrow=c(3,4),mar=c(4,4,1,1),xaxs='i',yaxs='i',las=1)
result<-as.null()
p<-0
for(city in top[c(1:12)]){
p<-p+1
#b<-x[which(x$ID_MUNICIP==city&x$EVOLUCAO %in% 2:3),]
b<-x[which(x$ID_MUNICIP==city),]
d<-x1[which(x1$ID_MUNICIP==city),]
#d<-table(b[,1])
#plot(as.Date(names(d)),as.numeric(d),type='l',main=city)
h3<-hist(b$NU_IDADE_N[which(b[,3]<as.Date('2020-11-23')&b[,3]>as.Date('2020-3-13'))],breaks=c(seq(0,90,by=5),120),plot=F)
age<-h3$breaks[-1]-5
h4<-hist(b$NU_IDADE_N[which(b[,3]>as.Date('2020-11-23')&b[,3]<as.Date('2021-9-20'))],breaks=c(seq(0,90,by=5),120),plot=F)
h5<-hist(d$NU_IDADE_N[which(d[,3]<as.Date('2020-12-7')&d[,3]>as.Date('2020-3-27'))],breaks=c(seq(0,90,by=5),120),plot=F)
h6<-hist(d$NU_IDADE_N[which(d[,3]>as.Date('2020-12-7'))],breaks=c(seq(0,90,by=5),120),plot=F)
#cfr1<-h2$counts/h1$counts
cfr2<-h5$counts/h3$counts
cfr3<-h6$counts/h4$counts
#matplot(age,cbind(h3$density,h4$density),lty=1,type='l',main=city,xlab='age',ylab='density',lwd=2)
matplot(age,cbind(cfr2,cfr3),lty=1,type='l',main=city,xlab='age',ylab='SARI-fatality-ratio',lwd=2,ylim=c(0,1))
if(p==1)
legend(bty='n',"topleft",col=1:3,lty=1,lwd=2,legend=c(paste(c('first','second'),'wave')))
abline(v=seq(20,100,by=20),lty=2)
abline(h=seq(0.2,1,by=0.2),lty=2)
result<-rbind(result,summary(d$NU_IDADE_N))
}
#row.names(result)<-top
dev.off()
#write.csv(result,'result2.csv',quote=F)
