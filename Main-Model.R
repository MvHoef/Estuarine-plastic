setwd("C:/Users/...")
rm(list=ls())

###READ DATA

DAT0<-read.csv("NAME.csv", header=TRUE)
#Repeat the data X times for multiple years (standard 1)
DAT<-DAT0[rep(seq_len(nrow(DAT0)), 1), ]
DAT$X<-c(1:365) #Days numbered in data change for number of years repeated.

###MODEL PARAMETERS
X<-21

L1<--0 #Parameter for loss from moving plastic (M1) to trapped plastic (M2)
L2<--0 #Parameter for loss from moving plastic (M1) out of the system
L3<-0 #Parameter for loss from trapped plastic (M2) out of system
L4<-1 #Parameter correction for loss of M2 plastic as outflow !real parameter inside function!
L5<-0 #Parameter for loss of M2 plastic to micro plastic (M3)
L6<-NA #Reserve parameter

Tau3<-0 #Residence time for sediment (and M3)

###MODEL

#Function for incoming plastic residence of one day
FUNC<-function(d,i){
  ifelse(d>=DAT$X[i] & d<=(DAT$X[i]+DAT$Tau[i]),
         DAT$Mid[i]*exp((L1+L2)*(d-DAT$X[i])), 0)}
F_test<-function(d){FUNC(d,25)}
plot(F_test, from=0, to=365, n=100000, xlab="Day", ylab="Plastic (kg)")

#M1 the mass of the moving plastic
MCal1<-data.frame(matrix(NA, nrow = 365, ncol = 365))
for(i in 1:365){
  MCal1[,i]=FUNC(i)
}
M1<-data.frame(cbind(1:365, colSums(MCal1)))
colnames(M1)<-c("Day", "Mass")
###Plot of M1 in system
plot(M1$Day, M1$Mass, type="l", xlab="day", ylab="Plastic (kg)")
###

#Outflow of M1 plastic
Q1<-data.frame(matrix(NA, nrow = 365, ncol = 1))
Q1FUNC<-function(i){FUNC((DAT$X[i]+DAT$Tau[i]),i)}
for(i in 1:365){
  Q1[i,]<-Q1FUNC(i)
}
colnames(Q1)<-"Qout"
Q1$X<-DAT$X+DAT$Tau
Q1<-aggregate(Q1,by=list(Q1$X),FUN=sum)
Q1$X<-NULL
colnames(Q1)<-c("X","Qout")
###
plot(Q1$X, Q1$Qout, xlab="day", ylab="plastic output (kg/day)", ylim=c(0,25),xlim=c(0,365),type="l")
###

#M1 plastic getting trapped and becoming M2 plastic
#Function if no loss from M1 took place
FUNC0<-function(d,i){
  ifelse(d>=DAT$X[i] & d<=(DAT$X[i]+DAT$Tau[i]),
         DAT$Mid[i], 0)}
#Function that gives for each timestep how much has entered M2 from M1 for each day.
FUNC2<-function(d,i){
  ifelse(d<=(DAT$X[i]+DAT$Tau[i]),(FUNC0(d,i)-FUNC(d,i))*L1/(L1+L2)
         ,FUNC2((DAT$X[i]+DAT$Tau[i]),i))
}
###
F_test<-function(d){FUNC2(d,25)}
plot(F_test, from=0, to=365, n=10000, xlab="day", ylab="Plastic (kg)")
###
#Calculate masses of M2 for each day timepath (without loss) (may take some time depending on amount of years)
MCal2<-data.frame(matrix(NA, nrow = 365, ncol = 365))
for(d in 1:365){
  MCal2[d,]=FUNC2(d)
}
#Calculate mass of M2 for each day (without loss)
M2a<-data.frame(cbind(1:365, rowSums(MCal2)))
colnames(M2a)<-c("day", "mass")
#added mass per time step:
i=2
M2a$dif<-0
for(i in 2:365){
  M2a$dif[i]<-M2a$mass[i]-M2a$mass[i-1]
}
###
plot(M2a$day, M2a$mass , type="l", xlab="day", ylab="Plastic (kg)")
plot(M2a$day, M2a$dif , type="l", xlab="day", ylab="Plastic (kg)")
###

#M2 with loss
M2<-data.frame(DAT$X)
colnames(M2)<-"day"
M2$DOY<-DAT$DOY
#Calculate mass of M2 with loss:
M2$mass<-0
for(i in 2:365){
  M2$mass[i]<-M2$mass[i-1]*(1-L3-(max(DAT$Tau)*L4-DAT$Tau[i])*1/(max(DAT$Tau)*L4)-L5)+M2a$dif[i]
}
###
plot(M2$day, M2$mass , type="l", xlab="day", ylab="Plastic (kg)")
lines(M1$Day, M1$Mass)
###

#M3 (micro plastic)
M3<-data.frame(DAT$X)
colnames(M3)<-"day"
M3$DOY<-DAT$DOY

FUNC3<-function(d,i){
  ifelse(d>=DAT$X[i] & d<=(DAT$X[i]+Tau3),
         L5*M2$mass[i], 0)}
###
F_test<-function(d){FUNC3(d,25)}
plot(F_test, from=0, to=365, n=100000, xlab="day", ylab="Plastic (kg)")
###
MCal3<-data.frame(matrix(NA, nrow = 365, ncol = 365))
for(i in 1:365){
  MCal3[,i]=FUNC3(i)
}
M3a<-data.frame(cbind(1:365, colSums(MCal3)))
M3$mass<-M3a$X2

###
plot(M3$day, M3$mass , type="l", xlab="day", ylab="Plastic (kg)")
###