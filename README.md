# Semi-Objective-Weather-Classification-
This is a semi-objective weather classifiction R code
getwd()
library(dplyr)
library(ncdf4)
library(stringr)
library(progress)
setwd("C:\\Users\\Lenovo\\Desktop\\synopticclass")

siteall <- read.csv("location.csv")
lao <- siteall

# 进度条
pb$tick()
Sys.sleep(1/length(myfile))
nc_close(nc.file)
}

myarray <- myarraynew
#Zm <- matrix(Z, nrow = nl, ncol = lt,dimnames = c(,Z[[3]]))#将数组转换成矩阵，以便计算
#####从这里开始天气分型的计算
  myarray <- myarraynew
  lt <-  dim(myarray)[3] #研究时段天数
  lt1 <- lt-1
  nro <-dim(myarray)[1]#区域网格行数（经度）
  nco <- dim(myarray)[2]#区域网格列数（纬度）
  nl <- nro*nco#区域网格数
  nll <-nro+nco+1 
  
  nc <- myarray
  S <- 1.5*nro*nco#S Sr Sc可根据实际情况修改（选代表天时尽量用1.0 1.8 1.8），确定各分型的天数时放宽限值，尽可能分到各分型中
  Sr <- 2.8*nco
  Sc <- 2.8*nro
  
  Sm <- matrix(S,nrow = 1,ncol = 1)
  Srm <- matrix(Sr,nrow = nco,ncol = 1)#纬度
  Scm <- matrix(Sc,nrow = nro,ncol = 1)#经度
  Srb <- rbind(Sm,Srm,Scm)
  
  M <- apply(nc,3,mean)#每天各网格均值
  sd <- apply(nc, 3, sd)#每天各网格标准偏差
  MR <- rep(M,each = nl)#每天各网格均值重复，每天重复
  MRM <- array(MR,dim = c(nro,nco,lt))#与Z匹配的数组，以便计算
  sdr <- rep(sd,each = nl)#类似于均值
  sdrm <- array(sdr,dim = c(nro,nco,lt))
  
  Z <- (nc-MRM)/sdrm#slp数据归一化
  Zm <- matrix(Z, nrow = nl, ncol = lt)#将数组转换成矩阵，以便计算
  
  for(i in 1:lt) Scal <- function(i){#每天进行相同的计算，每天分别与剩下的天数进行计算
    Scal <- (Zm[,i]-Zm[,-i])^2
    return(Scal)
  }
  
  for (i in 1:lt) Scalm <- function(i){#矩阵再转换成数组
    Scalm <- array(Scal(i),dim = c(nro,nco,lt1) )
    return(Scalm)
  }
for (i in 1:lt) ScalmS <- function(i){#每天所有网格求和 ###st
  ScalmS <- apply(Scalm(i), 3, sum)
  return(ScalmS)
}


for (i in 1:lt) ScalmSr <- function(i){#每天列的网格求和（-1：除行以外）##sr
  ScalmSr <- apply(Scalm(i), -1, sum)
  return(ScalmSr)
}

for (i in 1:lt) ScalmSc <- function(i){#每天行的网格求和 ###sc
  ScalmSc <- apply(Scalm(i), -2, sum)
  return(ScalmSc)
}

for(i in 1:lt) ScalmSm <- function(i){
  ScalmSm <- matrix(ScalmS(i),nrow = 1,ncol = lt1)#重组矩阵，对比数值
  return(ScalmSm)
}

for(i in 1:lt) ScalmSrm <- function(i){
  ScalmSrm <- matrix(ScalmSr(i),nrow = nco,ncol = lt1)
  return(ScalmSrm)
}

for(i in 1:lt) ScalmScm <- function(i){
  ScalmScm <- matrix(ScalmSc(i),nrow = nro,ncol = lt1)
  return(ScalmScm)
}

for(i in 1:lt) ScalmSrb <- function(i){#合并每天所有网格、每天列网格、每天行网格的值
  ScalmSrb <- rbind(ScalmSm(i),ScalmSrm(i),ScalmScm(i))
  return(ScalmSrb)
}

for(i in 1:lt) ScalmScb <- function(i){#合并计算结果和限值，便于比较
  ScalmScb <- cbind(ScalmSrb(i),Srb)
  return(ScalmScb)
}

for(i in 1:lt) ScalmSi <- function(i){
  ScalmSi <- ifelse(ScalmScb(i)[,1:lt1] < ScalmScb(i)[,lt],1,2)
  return(ScalmSi)
}


for(i in 1:lt) ScalmSis <- function(i){
  ScalmSis <- apply(ScalmSi(i),2,sum)
  return(ScalmSis)
}


for(i in 1:lt) ScalmSiss <- function(i){
  ScalmSiss <- ifelse(ScalmSis(i)==nll,1,2)
  return(ScalmSiss)
}

for(i in 1:lt) ScalmSissc <- function(i){
  ScalmSissc <- length(which(ScalmSiss(i)==1)) 
  return(ScalmSissc)
}

###进度条
pb <- progress_bar$new(
  format = 'p1计算进度 [:bar] :percent 剩余时间: :eta',
  total = lt,clear = TRUE, width= 60)

rvv <- vector(mode="numeric",length = lt)
for (i in 1:lt) {
  rvv[i] <- ScalmSissc(i)
  # 进度条
  pb$tick()
  Sys.sleep(1/lt)
}

rvm <- max(rvv)
rvl <- which(rvv == rvm)#确定第一个代表天
rvr <- which(ScalmSiss(rvl[1])==1)
print(rvr)
nc1 <- nc[,,rvr]#第一个分型mslp值
num <- c(1:lt)
num1 <- num[rvr]

ncd1 <- nc1
ncdm1 <- matrix(ncd1,nrow = nl,ncol = rvm)
ncdmm1 <- matrix(rowMeans(ncdm1),nrow = nl,ncol = 1) ##行平均
lao1 <- cbind(lao,ncdmm1)
colnames(lao1) <- c("lon","lat","p1")
write.csv(lao1,file = "p1.csv",quote = TRUE, na = "NA",row.names = FALSE) #

rvr1date <- data.frame(rvr1date=dimnames(nc)[[3]][rvr])#第1个分型的代表天日期
write.csv(rvr1date,"p1_date.csv",row.names = F)

#第二个分型开始计算
Z2 <- (nc[,,-rvr]-MRM[,,-rvr])/sdrm[,,-rvr]
lt2 <- length(Z2)/nrow(Z2)/ncol(Z2)  #天数
ltl2 <- lt2-1

Zm2 <- matrix(Z2, nrow = nl, ncol = lt2)

for(i in 1:lt2) Scal2 <- function(i){
  Scal2 <- (Zm2[,i]-Zm2[,-i])^2
  return(Scal2)
}

for (i in 1:lt2) Scalm2 <- function(i){
  Scalm2 <- array(Scal2(i),dim = c(nro,nco,ltl2) )
  return(Scalm2)
}

for (i in 1:lt2) ScalmS2 <- function(i){
  ScalmS2 <- apply(Scalm2(i), 3, sum)
  return(ScalmS2)
}

for (i in 1:lt2) ScalmSr2 <- function(i){
  ScalmSr2 <- apply(Scalm2(i), -1, sum)#列求和
  return(ScalmSr2)
}

for (i in 1:lt2) ScalmSc2 <- function(i){
  ScalmSc2 <- apply(Scalm2(i), -2, sum)
  return(ScalmSc2)
}

for(i in 1:lt2) ScalmSm2 <- function(i){
  ScalmSm2 <- matrix(ScalmS2(i),nrow = 1,ncol = ltl2)
  return(ScalmSm2)
}

for(i in 1:lt2) ScalmSrm2 <- function(i){
  ScalmSrm2 <- matrix(ScalmSr2(i),nrow = nco,ncol = ltl2)
  return(ScalmSrm2)
}

for(i in 1:lt2) ScalmScm2 <- function(i){
  ScalmScm2 <- matrix(ScalmSc2(i),nrow = nro,ncol = ltl2)
  return(ScalmScm2)
}

for(i in 1:lt2) ScalmSrb2 <- function(i){
  ScalmSrb2 <- rbind(ScalmSm2(i),ScalmSrm2(i),ScalmScm2(i))
  return(ScalmSrb2)
}

for(i in 1:lt2) ScalmScb2 <- function(i){
  ScalmScb2 <- cbind(ScalmSrb2(i),Srb)
  return(ScalmScb2)
}

for(i in 1:lt2) ScalmSi2 <- function(i){
  ScalmSi2 <- ifelse(ScalmScb2(i)[,1:ltl2] < ScalmScb2(i)[,lt2],1,2)
  return(ScalmSi2)
}

for(i in 1:lt2) ScalmSis2 <- function(i){
  ScalmSis2 <- apply(ScalmSi2(i),2,sum)
  return(ScalmSis2)
}

for(i in 1:lt2) ScalmSiss2 <- function(i){
  ScalmSiss2 <- ifelse(ScalmSis2(i)==nll,1,2)
  return(ScalmSiss2)
}

for(i in 1:lt2) ScalmSissc2 <- function(i){
  ScalmSissc2 <- length(which(ScalmSiss2(i)==1)) 
  return(ScalmSissc2)
}

###进度条
pb <- progress_bar$new(
  format = 'p2计算进度 [:bar] :percent 剩余时间: :eta',
  total = lt2,clear = TRUE, width= 60)

rvv2 <- vector(mode="numeric",length = lt2)
for (i in 1:lt2) {
  rvv2[i] <- ScalmSissc2(i)
  # 进度条
  pb$tick()
  Sys.sleep(1/lt2)
}

rvm2 <- max(rvv2)
rvl2 <- which(rvv2 == rvm2)#确定第二个代表天
rvr2 <- which(ScalmSiss2(rvl2[1])==1)

num2 <- num[-rvr]
numm2 <- num2[rvr2]
np2 <- num2[rvl2] #第二个分型的代表天


nc2 <- nc[,,numm2]#第二个分型mslp值

ncd2 <- nc2
ncdm2 <- matrix(ncd2,nrow = nl,ncol = rvm2)
ncdmm2 <- matrix(rowMeans(ncdm2),nrow = nl,ncol = 1)
lao2 <- cbind(lao,ncdmm2)
colnames(lao2) <- c("lon","lat","p2")
#write.csv(lao2,file = "p2.csv",quote = TRUE, na = "NA",row.names = FALSE)

rvr2date <- data.frame(rvr2date=dimnames(nc)[[3]][rvr2])#第3个分型的代表天日期
#write.csv(rvr2date,"p2_date.csv",row.names = F)


#第三个分型开始计算
rvrb12 <- sort(c(rvr,numm2))

Z3 <- (nc[,,-rvrb12]-MRM[,,-rvrb12])/sdrm[,,-rvrb12]
lt3 <- length(Z3)/nrow(Z3)/ncol(Z3)  #天数
ltl3 <- lt3-1

Zm3 <- matrix(Z3, nrow = nl, ncol = lt3)

for(i in 1:lt3) Scal3 <- function(i){
  Scal3 <- (Zm3[,i]-Zm3[,-i])^2
  return(Scal3)
}

for (i in 1:lt3) Scalm3 <- function(i){
  Scalm3 <- array(Scal3(i),dim = c(nro,nco,ltl3) )
  return(Scalm3)
}

for (i in 1:lt3) ScalmS3 <- function(i){
  ScalmS3 <- apply(Scalm3(i), 3, sum)
  return(ScalmS3)
}

for (i in 1:lt3) ScalmSr3 <- function(i){
  ScalmSr3 <- apply(Scalm3(i), -1, sum)#列求和
  return(ScalmSr3)
}

for (i in 1:lt3) ScalmSc3 <- function(i){
  ScalmSc3 <- apply(Scalm3(i), -2, sum)
  return(ScalmSc3)
}

for(i in 1:lt3) ScalmSm3 <- function(i){
  ScalmSm3 <- matrix(ScalmS3(i),nrow = 1,ncol = ltl3)
  return(ScalmSm3)
}

for(i in 1:lt3) ScalmSrm3 <- function(i){
  ScalmSrm3 <- matrix(ScalmSr3(i),nrow = nco,ncol = ltl3)
  return(ScalmSrm3)
}

for(i in 1:lt3) ScalmScm3 <- function(i){
  ScalmScm3 <- matrix(ScalmSc3(i),nrow = nro,ncol = ltl3)
  return(ScalmScm3)
}

for(i in 1:lt3) ScalmSrb3 <- function(i){
  ScalmSrb3 <- rbind(ScalmSm3(i),ScalmSrm3(i),ScalmScm3(i))
  return(ScalmSrb3)
}

for(i in 1:lt3) ScalmScb3 <- function(i){
  ScalmScb3 <- cbind(ScalmSrb3(i),Srb)
  return(ScalmScb3)
}

for(i in 1:lt3) ScalmSi3 <- function(i){
  ScalmSi3 <- ifelse(ScalmScb3(i)[,1:ltl3] < ScalmScb3(i)[,lt3],1,2)
  return(ScalmSi3)
}

for(i in 1:lt3) ScalmSis3 <- function(i){
  ScalmSis3 <- apply(ScalmSi3(i),2,sum)
  return(ScalmSis3)
}

for(i in 1:lt3) ScalmSiss3 <- function(i){
  ScalmSiss3 <- ifelse(ScalmSis3(i)==nll,1,2)
  return(ScalmSiss3)
}

for(i in 1:lt3) ScalmSissc3 <- function(i){
  ScalmSissc3 <- length(which(ScalmSiss3(i)==1)) 
  return(ScalmSissc3)
}

pb <- progress_bar$new(
  format = 'p3计算进度 [:bar] :percent 剩余时间: :eta',
  total = lt3,clear = TRUE, width= 60)

rvv3 <- vector(mode="numeric",length = lt3)
for (i in 1:lt3) {
  rvv3[i] <- ScalmSissc3(i)
  # 进度条
  pb$tick()
  Sys.sleep(1/lt3)
}

rvm3 <- max(rvv3)

rvl3 <- which(rvv3 == rvm3)#确定第3个代表天

rvr3 <- which(ScalmSiss3(rvl3[1])==1)


num3 <- num[-rvrb12]
numm3 <- num3[rvr3]
np3 <- num3[rvl3] #第3个分型的代表天
nc3 <- nc[,,numm3]#第3个分型mslp值

ncd3 <- nc3
ncdm3 <- matrix(ncd3,nrow = nl,ncol = rvm3)
ncdmm3 <- matrix(rowMeans(ncdm3),nrow = nl,ncol = 1)
lao3 <- cbind(lao,ncdmm3)
colnames(lao3) <- c("lon","lat","p3")
#write.csv(lao3,file = "p3.csv",quote = TRUE, na = "NA",row.names = FALSE)
rvr3date <- data.frame(rvr3date=dimnames(nc)[[3]][rvr3])#第3个分型的代表天日期
#write.csv(rvr3date,"p3_date.csv",row.names = F)

#第四个分型开始计算
rvrb123 <- sort(c(rvr,numm2,numm3))

Z4 <- (nc[,,-rvrb123]-MRM[,,-rvrb123])/sdrm[,,-rvrb123]
lt4 <- length(Z4)/nrow(Z4)/ncol(Z4)  #天数
ltl4 <- lt4-1

Zm4 <- matrix(Z4, nrow = nl, ncol = lt4)

for(i in 1:lt4) Scal4 <- function(i){
  Scal4 <- (Zm4[,i]-Zm4[,-i])^2
  return(Scal4)
}

for (i in 1:lt4) Scalm4 <- function(i){
  Scalm4 <- array(Scal4(i),dim = c(nro,nco,ltl4) )
  return(Scalm4)
}

for (i in 1:lt4) ScalmS4 <- function(i){
  ScalmS4 <- apply(Scalm4(i), 3, sum)
  return(ScalmS4)
}

for (i in 1:lt4) ScalmSr4 <- function(i){
  ScalmSr4 <- apply(Scalm4(i), -1, sum)#列求和
  return(ScalmSr4)
}

for (i in 1:lt4) ScalmSc4 <- function(i){
  ScalmSc4 <- apply(Scalm4(i), -2, sum)
  return(ScalmSc4)
}

for(i in 1:lt4) ScalmSm4 <- function(i){
  ScalmSm4 <- matrix(ScalmS4(i),nrow = 1,ncol = ltl4)
  return(ScalmSm4)
}

for(i in 1:lt4) ScalmSrm4 <- function(i){
  ScalmSrm4 <- matrix(ScalmSr4(i),nrow = nco,ncol = ltl4)
  return(ScalmSrm4)
}

for(i in 1:lt4) ScalmScm4 <- function(i){
  ScalmScm4 <- matrix(ScalmSc4(i),nrow = nro,ncol = ltl4)
  return(ScalmScm4)
}

for(i in 1:lt4) ScalmSrb4 <- function(i){
  ScalmSrb4 <- rbind(ScalmSm4(i),ScalmSrm4(i),ScalmScm4(i))
  return(ScalmSrb4)
}

for(i in 1:lt4) ScalmScb4 <- function(i){
  ScalmScb4 <- cbind(ScalmSrb4(i),Srb)
  return(ScalmScb4)
}

for(i in 1:lt4) ScalmSi4 <- function(i){
  ScalmSi4 <- ifelse(ScalmScb4(i)[,1:ltl4] < ScalmScb4(i)[,lt4],1,2)
  return(ScalmSi4)
}

for(i in 1:lt4) ScalmSis4 <- function(i){
  ScalmSis4 <- apply(ScalmSi4(i),2,sum)
  return(ScalmSis4)
}

for(i in 1:lt4) ScalmSiss4 <- function(i){
  ScalmSiss4 <- ifelse(ScalmSis4(i)==nll,1,2)
  return(ScalmSiss4)
}

for(i in 1:lt4) ScalmSissc4 <- function(i){
  ScalmSissc4 <- length(which(ScalmSiss4(i)==1)) 
  return(ScalmSissc4)
}

pb <- progress_bar$new(
  format = 'p4计算进度 [:bar] :percent 剩余时间: :eta',
  total = lt4,clear = TRUE, width= 60)

rvv4 <- vector(mode="numeric",length = lt4)
for (i in 1:lt4) {
  rvv4[i] <- ScalmSissc4(i)
  # 进度条
  pb$tick()
  Sys.sleep(1/lt4)
}

rvm4 <- max(rvv4)

rvl4 <- which(rvv4 == rvm4)#确定第4个代表天
rvr4 <- which(ScalmSiss4(rvl4[1])==1)

num4 <- num[-rvrb123]
numm4 <- num4[rvr4]
np4 <- num4[rvl4] #第4个分型的代表天

nc4 <- nc[,,numm4]#第4个分型mslp值

ncd4 <- nc4
ncdm4 <- matrix(ncd4,nrow = nl,ncol = rvm4)
ncdmm4 <- matrix(rowMeans(ncdm4),nrow = nl,ncol = 1)
lao4 <- cbind(lao,ncdmm4)
colnames(lao4) <- c("lon","lat","p4")
write.csv(lao4,file = "p4.csv",quote = TRUE, na = "NA",row.names = FALSE)

rvr4date <- data.frame(rvr4date=dimnames(nc)[[3]][rvr4])#第4个分型的代表天日期
write.csv(rvr4date,"p4_date.csv",row.names = F)
