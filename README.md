# Semi-Objective-Weather-Classification
This is a semi-objective weather classifiction R code
getwd()
library(dplyr)
library(ncdf4)
library(stringr)
library(progress)
setwd("C:\\Users\\Lenovo\\Desktop\\synopticclass")

siteall <- read.csv("location.csv")
lao <- siteall

pb$tick()#进度条
Sys.sleep(1/length(myfile))
nc_close(nc.file)
}

myarray <- myarraynew
#Zm <- matrix(Z, nrow = nl, ncol = lt,dimnames = c(,Z[[3]]))#将数组转换成矩阵，以便计算
#####从这里开始天气分型的计算
 ##########xiaohua
  install.packages("prograss")
  library(progress)
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
write.csv(lao1,file = "G:\\zhangshu\\天气分型\\p1.csv",quote = TRUE, na = "NA",row.names = FALSE) #


rvr1date <- data.frame(rvr1date=dimnames(nc)[[3]][num1])#第1个分型的代表天日期

write.csv(rvr1date,"G:\\zhangshu\\天气分型\\p1_date.csv",row.names = F)
###代表天日期
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
write.csv(lao2,file = "G:\\zhangshu\\天气分型\\p2.csv",quote = TRUE, na = "NA",row.names = FALSE)

rvr2date <- data.frame(rvr2date=dimnames(nc)[[3]][numm2])#第2个分型的代表天日期
write.csv(rvr2date,"G:\\zhangshu\\天气分型\\p2_date.csv",row.names = F)


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
write.csv(lao3,file = "G:\\zhangshu\\天气分型\\p3.csv",quote = TRUE, na = "NA",row.names = FALSE)
rvr3date <- data.frame(rvr3date=dimnames(nc)[[3]][numm3])#第3个分型的代表天日期
write.csv(rvr3date,"G:\\zhangshu\\天气分型\\p3_date.csv",row.names = F)

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
write.csv(lao4,file = "G:\\zhangshu\\天气分型\\p4.csv",quote = TRUE, na = "NA",row.names = FALSE)

rvr4date <- data.frame(rvr4date=dimnames(nc)[[3]][numm4])#第4个分型的代表天日期
write.csv(rvr4date,"G:\\zhangshu\\天气分型\\p4_date.csv",row.names = F)



#第5个分型开始计算
rvrb5 <- sort(c(rvr,numm2,numm3,numm4))
rvrbl5 <- length(num[-rvrb5])
rvrbl15 <- rvrbl5-1

Z5 <- (nc[,,-rvrb5]-MRM[,,-rvrb5])/sdrm[,,-rvrb5]

Zm5 <- matrix(Z5, nrow = nl, ncol = rvrbl5)

for(i in 1:rvrbl5) Scal5 <- function(i){
  Scal5 <- (Zm5[,i]-Zm5[,-i])^2
  return(Scal5)
}

for (i in 1:rvrbl5) Scalm5 <- function(i){
  Scalm5 <- array(Scal5(i),dim = c(nro,nco,rvrbl15) )
  return(Scalm5)
}

for (i in 1:rvrbl5) ScalmS5 <- function(i){
  ScalmS5 <- apply(Scalm5(i), 3, sum)
  return(ScalmS5)
}

for (i in 1:rvrbl5) ScalmSr5 <- function(i){
  ScalmSr5 <- apply(Scalm5(i), -1, sum)#列求和
  return(ScalmSr5)
}

for (i in 1:rvrbl5) ScalmSc5 <- function(i){
  ScalmSc5 <- apply(Scalm5(i), -2, sum)
  return(ScalmSc5)
}

for(i in 1:rvrbl5) ScalmSm5 <- function(i){
  ScalmSm5 <- matrix(ScalmS5(i),nrow = 1,ncol = rvrbl15)
  return(ScalmSm5)
}

for(i in 1:rvrbl5) ScalmSrm5 <- function(i){
  ScalmSrm5 <- matrix(ScalmSr5(i),nrow = nco,ncol = rvrbl15)
  return(ScalmSrm5)
}

for(i in 1:rvrbl5) ScalmScm5 <- function(i){
  ScalmScm5 <- matrix(ScalmSc5(i),nrow = nro,ncol = rvrbl15)
  return(ScalmScm5)
}

for(i in 1:rvrbl5) ScalmSrb5 <- function(i){
  ScalmSrb5 <- rbind(ScalmSm5(i),ScalmSrm5(i),ScalmScm5(i))
  return(ScalmSrb5)
}

for(i in 1:rvrbl5) ScalmScb5 <- function(i){
  ScalmScb5 <- cbind(ScalmSrb5(i),Srb)
  return(ScalmScb5)
}

for(i in 1:rvrbl5) ScalmSi5 <- function(i){
  ScalmSi5 <- ifelse(ScalmScb5(i)[,1:rvrbl15] < ScalmScb5(i)[,rvrbl5],1,2)
  return(ScalmSi5)
}

for(i in 1:rvrbl5) ScalmSis5 <- function(i){
  ScalmSis5 <- apply(ScalmSi5(i),2,sum)
  return(ScalmSis5)
}

for(i in 1:rvrbl5) ScalmSiss5 <- function(i){
  ScalmSiss5 <- ifelse(ScalmSis5(i)==nll,1,2)
  return(ScalmSiss5)
}

for(i in 1:rvrbl5) ScalmSissc5 <- function(i){
  ScalmSissc5 <- length(which(ScalmSiss5(i)==1)) 
  return(ScalmSissc5)
}


pb <- progress_bar$new(
  format = 'p5计算进度 [:bar] :percent 剩余时间: :eta',
  total = rvrbl5,clear = TRUE, width= 60)

rvv5 <- vector(mode="numeric",length = rvrbl5)
for (i in 1:rvrbl5) {
  rvv5[i] <- ScalmSissc5(i)
  pb$tick()
  Sys.sleep(1/rvrbl5)
}

rvm5 <- max(rvv5)

rvl5 <- which(rvv5 == rvm5)#确定第5个代表天

rvr5 <- which(ScalmSiss5(rvl5[1])==1)

num5 <- num[-rvrb5]
numm5 <- num5[rvr5]
np5 <- num5[rvl5] #第5个分型的代表天

nc5 <- nc[,,numm5]#第5个分型mslp值

ncd5 <- nc5 #单位换算为hPa
ncdm5 <- matrix(ncd5,nrow = nl,ncol = rvm5)
ncdmm5 <- matrix(rowMeans(ncdm5),nrow = nl,ncol = 1)#均值
lao5 <- cbind(lao,ncdmm5)
colnames(lao5) <- c("lon","lat","p5")
write.csv(lao5,file = "G:\\zhangshu\\天气分型\\p5.csv",quote = TRUE, na = "NA",row.names = FALSE)

rvr5date <- data.frame(rvr5date=dimnames(nc)[[3]][numm5])#第5个分型的日期
write.csv(rvr5date,"G:\\zhangshu\\天气分型\\p5_date.csv",row.names = F)




#第6个分型开始计算
rvrb6 <- sort(c(rvr,numm2,numm3,numm4,numm5))
rvrbl6 <- length(num[-rvrb6])
rvrbl16 <- rvrbl6-1

Z6 <- (nc[,,-rvrb6]-MRM[,,-rvrb6])/sdrm[,,-rvrb6]

Zm6 <- matrix(Z6, nrow = nl, ncol = rvrbl6)

for(i in 1:rvrbl6) Scal6 <- function(i){
  Scal6 <- (Zm6[,i]-Zm6[,-i])^2
  return(Scal6)
}

for (i in 1:rvrbl6) Scalm6 <- function(i){
  Scalm6 <- array(Scal6(i),dim = c(nro,nco,rvrbl16) )
  return(Scalm6)
}

for (i in 1:rvrbl6) ScalmS6 <- function(i){
  ScalmS6 <- apply(Scalm6(i), 3, sum)
  return(ScalmS6)
}

for (i in 1:rvrbl6) ScalmSr6 <- function(i){
  ScalmSr6 <- apply(Scalm6(i), -1, sum)#列求和
  return(ScalmSr6)
}

for (i in 1:rvrbl6) ScalmSc6 <- function(i){
  ScalmSc6 <- apply(Scalm6(i), -2, sum)
  return(ScalmSc6)
}

for(i in 1:rvrbl6) ScalmSm6 <- function(i){
  ScalmSm6 <- matrix(ScalmS6(i),nrow = 1,ncol = rvrbl16)
  return(ScalmSm6)
}

for(i in 1:rvrbl6) ScalmSrm6 <- function(i){
  ScalmSrm6 <- matrix(ScalmSr6(i),nrow = nco,ncol = rvrbl16)
  return(ScalmSrm6)
}

for(i in 1:rvrbl6) ScalmScm6 <- function(i){
  ScalmScm6 <- matrix(ScalmSc6(i),nrow = nro,ncol = rvrbl16)
  return(ScalmScm6)
}

for(i in 1:rvrbl6) ScalmSrb6 <- function(i){
  ScalmSrb6 <- rbind(ScalmSm6(i),ScalmSrm6(i),ScalmScm6(i))
  return(ScalmSrb6)
}

for(i in 1:rvrbl6) ScalmScb6 <- function(i){
  ScalmScb6 <- cbind(ScalmSrb6(i),Srb)
  return(ScalmScb6)
}

for(i in 1:rvrbl6) ScalmSi6 <- function(i){
  ScalmSi6 <- ifelse(ScalmScb6(i)[,1:rvrbl16] < ScalmScb6(i)[,rvrbl6],1,2)
  return(ScalmSi6)
}

for(i in 1:rvrbl6) ScalmSis6 <- function(i){
  ScalmSis6 <- apply(ScalmSi6(i),2,sum)
  return(ScalmSis6)
}

for(i in 1:rvrbl6) ScalmSiss6 <- function(i){
  ScalmSiss6 <- ifelse(ScalmSis6(i)==nll,1,2)
  return(ScalmSiss6)
}

for(i in 1:rvrbl6) ScalmSissc6 <- function(i){
  ScalmSissc6 <- length(which(ScalmSiss6(i)==1)) 
  return(ScalmSissc6)
}


pb <- progress_bar$new(
  format = 'p6计算进度 [:bar] :percent 剩余时间: :eta',
  total = rvrbl6,clear = TRUE, width= 60)

rvv6 <- vector(mode="numeric",length = rvrbl6)
for (i in 1:rvrbl6) {
  rvv6[i] <- ScalmSissc6(i)
  pb$tick()
  Sys.sleep(1/rvrbl6)
}



rvm6 <- max(rvv6)

rvl6 <- which(rvv6 == rvm6)#确定第6个代表天

rvr6 <- which(ScalmSiss6(rvl6[1])==1)

num6 <- num[-rvrb6]
numm6 <- num6[rvr6]
np6 <- num6[rvl6] #第6个分型的代表天

nc6 <- nc[,,numm6]#第6个分型mslp值

ncd6 <- nc6#单位换算为hPa
ncdm6 <- matrix(ncd6,nrow = nl,ncol = rvm6)
ncdmm6 <- matrix(rowMeans(ncdm6),nrow = nl,ncol = 1)#均值
lao6 <- cbind(lao,ncdmm6)

colnames(lao6) <- c("lon","lat","p6")
write.csv(lao6,file = "G:\\zhangshu\\天气分型\\p6.csv",quote = TRUE, na = "NA",row.names = FALSE)

rvr6date <- data.frame(rvr6date=dimnames(nc)[[3]][numm6])#第6个分型的日期
write.csv(rvr6date,"G:\\zhangshu\\天气分型\\p6_date.csv",row.names = F)


#第7个分型开始计算
rvrb7 <- sort(c(rvr,numm2,numm3,numm4,numm5,numm6))
rvrbl7 <- length(num[-rvrb7])
rvrbl17 <- rvrbl7-1

Z7 <- (nc[,,-rvrb7]-MRM[,,-rvrb7])/sdrm[,,-rvrb7]

Zm7 <- matrix(Z7, nrow = nl, ncol = rvrbl7)

for(i in 1:rvrbl7) Scal7 <- function(i){
  Scal7 <- (Zm7[,i]-Zm7[,-i])^2
  return(Scal7)
}

for (i in 1:rvrbl7) Scalm7 <- function(i){
  Scalm7 <- array(Scal7(i),dim = c(nro,nco,rvrbl17) )
  return(Scalm7)
}

for (i in 1:rvrbl7) ScalmS7 <- function(i){
  ScalmS7 <- apply(Scalm7(i), 3, sum)
  return(ScalmS7)
}

for (i in 1:rvrbl7) ScalmSr7 <- function(i){
  ScalmSr7 <- apply(Scalm7(i), -1, sum)#列求和
  return(ScalmSr7)
}

for (i in 1:rvrbl7) ScalmSc7 <- function(i){
  ScalmSc7 <- apply(Scalm7(i), -2, sum)
  return(ScalmSc7)
}

for(i in 1:rvrbl7) ScalmSm7 <- function(i){
  ScalmSm7 <- matrix(ScalmS7(i),nrow = 1,ncol = rvrbl17)
  return(ScalmSm7)
}

for(i in 1:rvrbl7) ScalmSrm7 <- function(i){
  ScalmSrm7 <- matrix(ScalmSr7(i),nrow = nco,ncol = rvrbl17)
  return(ScalmSrm7)
}

for(i in 1:rvrbl7) ScalmScm7 <- function(i){
  ScalmScm7 <- matrix(ScalmSc7(i),nrow = nro,ncol = rvrbl17)
  return(ScalmScm7)
}

for(i in 1:rvrbl7) ScalmSrb7 <- function(i){
  ScalmSrb7 <- rbind(ScalmSm7(i),ScalmSrm7(i),ScalmScm7(i))
  return(ScalmSrb7)
}

for(i in 1:rvrbl7) ScalmScb7 <- function(i){
  ScalmScb7 <- cbind(ScalmSrb7(i),Srb)
  return(ScalmScb7)
}

for(i in 1:rvrbl7) ScalmSi7 <- function(i){
  ScalmSi7 <- ifelse(ScalmScb7(i)[,1:rvrbl17] < ScalmScb7(i)[,rvrbl7],1,2)
  return(ScalmSi7)
}

for(i in 1:rvrbl7) ScalmSis7 <- function(i){
  ScalmSis7 <- apply(ScalmSi7(i),2,sum)
  return(ScalmSis7)
}

for(i in 1:rvrbl7) ScalmSiss7 <- function(i){
  ScalmSiss7 <- ifelse(ScalmSis7(i)==nll,1,2)
  return(ScalmSiss7)
}

for(i in 1:rvrbl7) ScalmSissc7 <- function(i){
  ScalmSissc7 <- length(which(ScalmSiss7(i)==1)) 
  return(ScalmSissc7)
}

pb <- progress_bar$new(
  format = 'p7计算进度 [:bar] :percent 剩余时间: :eta',
  total = rvrbl7,clear = TRUE, width= 60)

rvv7 <- vector(mode="numeric",length = rvrbl7)
for (i in 1:rvrbl7) {
  rvv7[i] <- ScalmSissc7(i)
  pb$tick()
  Sys.sleep(1/rvrbl7)
}



rvm7 <- max(rvv7)

rvl7 <- which(rvv7 == rvm7)#确定第7个代表天

rvr7 <- which(ScalmSiss7(rvl7[1])==1)

num7 <- num[-rvrb7]
numm7 <- num7[rvr7]
np7 <- num7[rvl7] #第7个分型的代表天

nc7 <- nc[,,numm7]#第7个分型mslp值

ncd7 <- nc7 #单位换算为hPa
ncdm7 <- matrix(ncd7,nrow = nl,ncol = rvm7)
ncdmm7 <- matrix(rowMeans(ncdm7),nrow = nl,ncol = 1)#均值
lao7 <- cbind(lao,ncdmm7)
colnames(lao7) <- c("lon","lat","p7")
write.csv(lao7,file = "G:\\zhangshu\\天气分型\\p7.csv",quote = TRUE, na = "NA",row.names = FALSE)

rvr7date <- data.frame(rvr7date=dimnames(nc)[[3]][numm7])#第7个分型的代表天日期
write.csv(rvr7date,"G:\\zhangshu\\天气分型\\p7_date.csv",row.names = F)



#第8个分型开始计算
rvrb8 <- sort(c(rvr,numm2,numm3,numm4,numm5,numm6,numm7))
rvrbl8 <- length(num[-rvrb8])
rvrbl18 <- rvrbl8-1

Z8 <- (nc[,,-rvrb8]-MRM[,,-rvrb8])/sdrm[,,-rvrb8]

Zm8 <- matrix(Z8, nrow = nl, ncol = rvrbl8)

for(i in 1:rvrbl8) Scal8 <- function(i){
  Scal8 <- (Zm8[,i]-Zm8[,-i])^2
  return(Scal8)
}

for (i in 1:rvrbl8) Scalm8 <- function(i){
  Scalm8 <- array(Scal8(i),dim = c(nro,nco,rvrbl18) )
  return(Scalm8)
}

for (i in 1:rvrbl8) ScalmS8 <- function(i){
  ScalmS8 <- apply(Scalm8(i), 3, sum)
  return(ScalmS8)
}

for (i in 1:rvrbl8) ScalmSr8 <- function(i){
  ScalmSr8 <- apply(Scalm8(i), -1, sum)#列求和
  return(ScalmSr8)
}

for (i in 1:rvrbl8) ScalmSc8 <- function(i){
  ScalmSc8 <- apply(Scalm8(i), -2, sum)
  return(ScalmSc8)
}

for(i in 1:rvrbl8) ScalmSm8 <- function(i){
  ScalmSm8 <- matrix(ScalmS8(i),nrow = 1,ncol = rvrbl18)
  return(ScalmSm8)
}

for(i in 1:rvrbl8) ScalmSrm8 <- function(i){
  ScalmSrm8 <- matrix(ScalmSr8(i),nrow = nco,ncol = rvrbl18)
  return(ScalmSrm8)
}

for(i in 1:rvrbl8) ScalmScm8 <- function(i){
  ScalmScm8 <- matrix(ScalmSc8(i),nrow = nro,ncol = rvrbl18)
  return(ScalmScm8)
}

for(i in 1:rvrbl8) ScalmSrb8 <- function(i){
  ScalmSrb8 <- rbind(ScalmSm8(i),ScalmSrm8(i),ScalmScm8(i))
  return(ScalmSrb8)
}

for(i in 1:rvrbl8) ScalmScb8 <- function(i){
  ScalmScb8 <- cbind(ScalmSrb8(i),Srb)
  return(ScalmScb8)
}

for(i in 1:rvrbl8) ScalmSi8 <- function(i){
  ScalmSi8 <- ifelse(ScalmScb8(i)[,1:rvrbl18] < ScalmScb8(i)[,rvrbl8],1,2)
  return(ScalmSi8)
}

for(i in 1:rvrbl8) ScalmSis8 <- function(i){
  ScalmSis8 <- apply(ScalmSi8(i),2,sum)
  return(ScalmSis8)
}

for(i in 1:rvrbl8) ScalmSiss8 <- function(i){
  ScalmSiss8 <- ifelse(ScalmSis8(i)==nll,1,2)
  return(ScalmSiss8)
}

for(i in 1:rvrbl8) ScalmSissc8 <- function(i){
  ScalmSissc8 <- length(which(ScalmSiss8(i)==1)) 
  return(ScalmSissc8)
}


pb <- progress_bar$new(
  format = 'p8计算进度 [:bar] :percent 剩余时间: :eta',
  total = rvrbl8,clear = TRUE, width= 60)

rvv8 <- vector(mode="numeric",length = rvrbl8)
for (i in 1:rvrbl8) {
  rvv8[i] <- ScalmSissc8(i)
  pb$tick()
  Sys.sleep(1/rvrbl8)
}



rvm8 <- max(rvv8)

rvl8 <- which(rvv8 == rvm8)#确定第8个代表天

rvr8 <- which(ScalmSiss8(rvl8[1])==1)

num8 <- num[-rvrb8]
numm8 <- num8[rvr8]
np8 <- num8[rvl8] #第8个分型的代表天

nc8 <- nc[,,numm8]#第8个分型mslp值

ncd8 <- nc8/100 #单位换算为hPa
ncdm8 <- matrix(ncd8,nrow = nl,ncol = rvm8)
ncdmm8 <- matrix(rowMeans(ncdm8),nrow = nl,ncol = 1)#均值
lao8 <- cbind(lao,ncdmm8)
colnames(lao8) <- c("lon","lat","p8")
write.csv(lao8,file = "G:\\zhangshu\\天气分型\\p8.csv",quote = TRUE, na = "NA",row.names = FALSE)
rvr8date <- data.frame(rvr8date=dimnames(nc)[[3]][numm8])#第8个分型的代表天日期
write.csv(rvr8date,"G:\\zhangshu\\天气分型\\p8_date.csv",row.names = F)



#第9个分型开始计算
rvrb9 <- sort(c(rvr,numm2,numm3,numm4,numm5,numm6,numm7,numm8))
rvrbl9 <- length(num[-rvrb9])
rvrbl19 <- rvrbl9-1

Z9 <- (nc[,,-rvrb9]-MRM[,,-rvrb9])/sdrm[,,-rvrb9]

Zm9 <- matrix(Z9, nrow = nl, ncol = rvrbl9)

for(i in 1:rvrbl9) Scal9 <- function(i){
  Scal9 <- (Zm9[,i]-Zm9[,-i])^2
  return(Scal9)
}

for (i in 1:rvrbl9) Scalm9 <- function(i){
  Scalm9 <- array(Scal9(i),dim = c(nro,nco,rvrbl19))
  return(Scalm9)
}

for (i in 1:rvrbl9) ScalmS9 <- function(i){
  ScalmS9 <- apply(Scalm9(i), 3, sum)
  return(ScalmS9)
}

for (i in 1:rvrbl9) ScalmSr9 <- function(i){
  ScalmSr9 <- apply(Scalm9(i), -1, sum)#列求和
  return(ScalmSr9)
}

for (i in 1:rvrbl9) ScalmSc9 <- function(i){
  ScalmSc9 <- apply(Scalm9(i), -2, sum)
  return(ScalmSc9)
}

for(i in 1:rvrbl9) ScalmSm9 <- function(i){
  ScalmSm9 <- matrix(ScalmS9(i),nrow = 1,ncol = rvrbl19)
  return(ScalmSm9)
}

for(i in 1:rvrbl9) ScalmSrm9 <- function(i){
  ScalmSrm9 <- matrix(ScalmSr9(i),nrow = nco,ncol = rvrbl19)
  return(ScalmSrm9)
}

for(i in 1:rvrbl9) ScalmScm9 <- function(i){
  ScalmScm9 <- matrix(ScalmSc9(i),nrow = nro,ncol = rvrbl19)
  return(ScalmScm9)
}

for(i in 1:rvrbl9) ScalmSrb9 <- function(i){
  ScalmSrb9 <- rbind(ScalmSm9(i),ScalmSrm9(i),ScalmScm9(i))
  return(ScalmSrb9)
}

for(i in 1:rvrbl9) ScalmScb9 <- function(i){
  ScalmScb9 <- cbind(ScalmSrb9(i),Srb)
  return(ScalmScb9)
}

for(i in 1:rvrbl9) ScalmSi9 <- function(i){
  ScalmSi9 <- ifelse(ScalmScb9(i)[,1:rvrbl19] < ScalmScb9(i)[,rvrbl9],1,2)
  return(ScalmSi9)
}

for(i in 1:rvrbl9) ScalmSis9 <- function(i){
  ScalmSis9 <- apply(ScalmSi9(i),2,sum)
  return(ScalmSis9)
}

for(i in 1:rvrbl9) ScalmSiss9 <- function(i){
  ScalmSiss9 <- ifelse(ScalmSis9(i)==nll,1,2)
  return(ScalmSiss9)
}

for(i in 1:rvrbl9) ScalmSissc9 <- function(i){
  ScalmSissc9 <- length(which(ScalmSiss9(i)==1)) 
  return(ScalmSissc9)
}

pb <- progress_bar$new(
  format = 'p9计算进度 [:bar] :percent 剩余时间: :eta',
  total = rvrbl9,clear = TRUE, width= 60)

rvv9 <- vector(mode="numeric",length = rvrbl9)
for (i in 1:rvrbl9) {
  rvv9[i] <- ScalmSissc9(i)
  pb$tick()
  Sys.sleep(1/rvrbl9)
}

rvm9 <- max(rvv9)

rvl9 <- which(rvv9 == rvm9)#确定第9个代表天

rvr9 <- which(ScalmSiss9(rvl9[1])==1)

num9 <- num[-rvrb9]
numm9 <- num9[rvr9]
np9 <- num9[rvl9] #第9个分型的代表天

nc9 <- nc[,,numm9]#第9个分型mslp值

ncd9 <- nc9/100 #单位换算为hPa
ncdm9 <- matrix(ncd9,nrow = nl,ncol = rvm9)
ncdmm9 <- matrix(rowMeans(ncdm9),nrow = nl,ncol = 1)#均值
lao9 <- cbind(lao,ncdmm9)
colnames(lao9) <- c("lon","lat","p9")
write.csv(lao9,file = "G:\\zhangshu\\天气分型\\p9.csv",quote = TRUE, na = "NA",row.names = FALSE)

rvr9date <- data.frame(rvr9date=dimnames(nc)[[3]][numm9])#第9个分型的代表天日期
write.csv(rvr9date,"G:\\zhangshu\\天气分型\\p9_date.csv",row.names = F)

rvm+rvm2+rvm3+rvm4+rvm5+rvm6+rvm7+rvm8+rvm9
