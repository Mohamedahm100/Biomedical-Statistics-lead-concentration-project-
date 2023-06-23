# Biomedical-Statistics-lead-concentration-project-
############ 0. Data reading:
Data <-get(load("D:/Uni/Project H.RData") )

library(dunn.test)
library(car)
library(multcomp)
library(report)
library(skimr)
skim(Data)

################################################
############ 1.	Descriptive statistics
summary(Data)

#Id
mean(Data$Id)
median(Data$Id,na.rm=TRUE)
min(Data$Id,na.rm=TRUE)
max(Data$Id,na.rm=TRUE)
quantile(Data$Id,na.rm=TRUE,c(0.25,0.75))

#Area     
mean(Data$Area)
median(Data$Area,na.rm=TRUE)
min(Data$Area,na.rm=TRUE)
max(Data$Area,na.rm=TRUE)
quantile(Data$Area,na.rm=TRUE,c(0.25,0.75))

#Age
mean(Data$Age)
median(Data$Age,na.rm=TRUE)
min(Data$Age,na.rm=TRUE)
max(Data$Age,na.rm=TRUE)
quantile(Data$Age,na.rm=TRUE,c(0.25,0.75))

#Sex
mean(Data$Sex)
median(Data$Sex,na.rm=TRUE)
min(Data$Sex,na.rm=TRUE)
max(Data$Sex,na.rm=TRUE)
quantile(Data$Sex,na.rm=TRUE,c(0.25,0.75))

#Iqf
mean(Data$Iqf)
median(Data$Iqf,na.rm=TRUE)
min(Data$Iqf,na.rm=TRUE)
max(Data$Iqf,na.rm=TRUE)
quantile(Data$Iqf,na.rm=TRUE,c(0.25,0.75))

#Lead_type
mean(Data$Lead_type)
median(Data$Lead_type,na.rm=TRUE)
min(Data$Lead_type,na.rm=TRUE)
max(Data$Lead_type,na.rm=TRUE)
quantile(Data$Lead_type,na.rm=TRUE,c(0.25,0.75))

#Ld72
mean(Data$Ld72)
median(Data$Ld72,na.rm=TRUE)
min(Data$Ld72,na.rm=TRUE)
max(Data$Ld72,na.rm=TRUE)
quantile(Data$Ld72,na.rm=TRUE,c(0.25,0.75))

#Ld73
mean(Data$Ld73)
median(Data$Ld73,na.rm=TRUE)
min(Data$Ld73,na.rm=TRUE)
max(Data$Ld73,na.rm=TRUE)
quantile(Data$Ld73,na.rm=TRUE,c(0.25,0.75))

#Totyrs
mean(Data$Totyrs)
median(Data$Totyrs,na.rm=TRUE)
min(Data$Totyrs,na.rm=TRUE)
max(Data$Totyrs,na.rm=TRUE)
quantile(Data$Totyrs,na.rm=TRUE,c(0.25,0.75))

#MAXFWT
mean(Data$MAXFWT)
median(Data$MAXFWT,na.rm=TRUE)
min(Data$MAXFWT,na.rm=TRUE)
max(Data$MAXFWT,na.rm=TRUE)
quantile(Data$MAXFWT,na.rm=TRUE,c(0.25,0.75))

#Exposed
mean(Data$Exposed)
median(Data$Exposed,na.rm=TRUE)
min(Data$Exposed,na.rm=TRUE)
max(Data$Exposed,na.rm=TRUE)
quantile(Data$Exposed,na.rm=TRUE,c(0.25,0.75))

# caterogical variables 
Data$Sex<-factor(Data$Sex,labels=c('male','female')) #assume 1= male, 2= female
class(Data$Sex  )
Category<-table(factor(Data[,"Sex"]))


# default correlation = Pearson
cor(Data$MAXFWT,Data$Ld72, use="complete.obs",method="spearman")
cor(Data$MAXFWT,Data$Ld73, use="complete.obs",method="spearman")

################################################################
############ 2. Graphs:
library(ggplot2) 
barplot(Category,main='Distribution of gender',horiz=FALSE,names.arg=c("male", "female"), xlab="Sex",ylab="Frequency")
barplot(tapply(Data$MAXFWT ,list(sex=Data$Sex),mean,na.rm=T),
        xlab="Sex",ylab="Mean MAXfWT ")

hist(Data$MAXFWT,xlab="MAXFWT",main="Distribution of MAXFWT")
hist(Data$Age,xlab="Age",main="Distribution of Age")

plot(lead$MAXFWT,lead$Ld72, col = "blue",xlab="MAXFWT",ylab="Ld72", main="Scatterplot")
plot(MAXFWT~Ld72, data=lead)
plot(lead$MAXFWT[lead$Sex=="1"],lead$Ld72[lead$Sex=="1"],xlab="MAXFWT",ylab="Ld72",main="Scatterplot",col=3)
points(lead$MAXFWT[lead$Sex=="2"],lead$Ld72[lead$Sex=="2"],col=2)
abline(lm(lead$MAXFWT[lead$Sex=="1"]~lead$Ld72[lead$Sex=="1"]),col=3)
abline(lm(lead$MAXFWT[lead$Sex=="2"]~lead$Ld72[lead$Sex=="2"]),col=2)


boxplot(Data$Age~ Data$Ld72, main='Box plot sperated by LD72  group', xlab=' Ld72 ',ylab=' Age ' )
boxplot(Data$Age~ Data$Ld73, main='Box plot sperated by LD73  group', xlab=' Treatment ',ylab=' ADWG0021 ' )

################################################
############ 3. Outliers:
outlier <- boxplot(Data)$out
outlier

boxplot(Data$Ld72)
outlier <- boxplot(Data$Ld72)$out
outlier

boxplot(Data$MAXFWT)
outlier <- boxplot(Data$MAXFWT)$out
outlier


outlier <- boxplot(Data$Id)$out
outlier

outlier <- boxplot(Data$Iqf)$out
outlier

outlier <- boxplot(Data$Exposed)$out
outlier

outlier <- boxplot(Data$Lead_type)$out
outlier

outlier <- boxplot(Data$Totyrs)$out
outlier

outlier <- boxplot(Data$Sex)$out
outlier

outlier <- boxplot(Data$Age)$out
outlier

outlier <- boxplot(Data$Area)$out
outlier

############ 4.	Testing for normality/ homoscedasticity ############## 
#normality
#using histogram
hist(Data$MAXFWT,xlab="MAXFWT", main='MAXFWT histogram')
hist(Data$Sex,xlab="gender", main='gender histogram')
hist(Data$Age,xlab="Age", main='Age histogram')
hist(Data$Ld72,xlab="Ld72", main='Ld72 histogram')
hist(Data$Ld73,xlab="Ld73", main='Ld73 histogram')

#using QQblot
qqnorm(Data$MAXFWT, main='MAXFWT QQ')
qqnorm(Data$Sex, main='Sex QQ')
qqnorm(Data$Age, main='Age QQ')
qqnorm(Data$Ld72, main='Ld72 QQ')
qqnorm(Data$Ld73, main='Ld73 QQ')

#using shapirotest
shapiro.test(Data$MAXFWT)
shapiro.test(Data$Sex)
shapiro.test(Data$Age)
shapiro.test(Data$Ld72)
shapiro.test(Data$Ld73)

#homoscedasticity
bartlett.test(list(Data$MAXFWT,Data$Ld72))
bartlett.test(list(Data$Ld73,Data$Ld72))

m <- lm(Data$MAXFWT ~ Data$Ld72, data = Data)
m$residuals
m$fitted
plot(m$residual ~ m$fitted.values,xlab="fitted.values",ylab="residual")

m <- lm(Data$Ld73 ~ Data$Ld72, data = Data)
m$residuals
m$fitted
plot(m$residual ~ m$fitted.values,xlab="fitted.values",ylab="residual")



############ 5.	Statistical Inference ############## 

t.test(Data$MAXFWT, Data$Sex, conf.level = 0.9)
t.test(Data$MAXFWT , Data$Sex, conf.level = 0.95)
t.test(Data$MAXFWT , Data$Sex, conf.level = 0.99)


############ 6.	Hypothesis testing ############## 
#1 we divide sex to male - female

#6.1
t.test(MAXFWT~Sex, data = Data)
# no significance in MAXFWT between males and females

#6.2
#test normality with shapiro test, homoscedasticity with var or lavine test
shapiro.test(MAXFWT~Sex)
shapiro.test(Data$Sex)

var.test(MAXFWT~Sex, data = Data)

#6.3
lower = lead[lead$Ld72 <= 40,]
higher = lead[lead$Ld72 > 40,]
wilcox.test( higher$MAXFWT ,lower$MAXFWT, alternative = "less" )

#6.4
# test heteroscedasiticy  with var or lavine test
var.test(higher$MAXFWT ,lower$MAXFWT, data = Data)

#6.5
#using anova test 

anovatest2 = aov(MAXFWT ~ Sex + Lead_type, data = Data)
summary(anovatest2)
eval(anovatest2)

##post_hoc testing##
library(multcomp)

post_test <- glht(anovatest2,linfct = mcp(SES = "Tukey"))

TukeyHSD(anovatest2)######
glht(anovatest2, linfct=mcp(SES = "Dunnett"))

############## 7.Linear model ##############
#7.1
#plot before and after model
plot(MAXFWT ~ Sex, data = lead)
plot(MAXFWT ~ Lead_type, data = lead)

LinMod <-lm(MAXFWT ~ Sex + Lead_type, data = lead)
summary(LinMod)
plot(LinMod)
abline(LinMod)

#7.2
sexCof <- confint(LinMod, 'Sex', level=0.95)
leadCof <- confint(LinMod, 'Lead_type', level=0.95)


#7.3

LinMod2 <-lm(MAXFWT ~ Ld73 , data = lead)
summary(LinMod2)

df <- data.frame (Ld73  = c(100))
predict(LinMod2, df)

#Using shapirotest and homoscedasticity to asume normality
shapiro.test(Data$MAXFWT)
bartlett.test(list(Data$MAXFWT,Data$Ld72))
#using shapirotest and bartlett-test, pvalue of MAXFWT smaller than alpha(0.05)
#so we reject null hypothesis
#they are not normally distributed
#so shapirotest meet for the test, but bartlett don't meet for the test
t.test(Data$MAXFWT, Data$Ld72, var.equal= TRUE)

var.test()
#using t-test, pvalue of MAXFWT larger than alpha(0.05)
# so we don't have evidence to reject null hypothesis
#it's normally distributed
#so it meet for the test

###


#p-value = 0.001518 is < 0.05 so it's significant and doesn't reject the null hypothesis

boxplot(MAXFWT ~ Sex, data = lead)
anovatest = anova(lm(MAXFWT ~ Sex + Lead_type, data = lead))
anovatest2 = aov(MAXFWT ~ Sex + Lead_type, data = lead)
summary(anovatest2)

summary(anovatest)
