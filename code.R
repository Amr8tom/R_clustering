#20198060 amr alaa ali                        
#20198001 abdala ahmed fathy

library(Biobase)
library(corrplot)
library(devtools)
library(broom)
library(antiProfilesData)


pdata=pData(apColonData)
edata=exprs(apColonData)
fdata = fData(apColonData)
#------------
head(pdata)
head(edata)
head(fdata)

# 1

typeof(pdata)
typeof(edata)
typeof(fdata)

# type of each column

sapply(pdata, class)
sapply(edata, class)
sapply(fdata, class)

# column names and rows name
row.names(pdata)
colnames(pdata)

row.names(edata)
colnames(edata)

row.names(fdata)
colnames(fdata)

# summary of each column
summary(pdata)
summary(edata)
summary(fdata)

#frequency of categorical data and  NA values 
table(pdata$filename,exclude=NULL)
table(pdata$DB_ID,exclude=NULL)
table(pdata$ExperimentID,exclude=NULL)
table(pdata$Tissue,exclude=NULL)
table(pdata$SubType,exclude=NULL)
table(pdata$ClinicalGroup,exclude=NULL)
table(pdata$Status,exclude=NULL)

# correlation and covariance between the first 10 columns 
edata0 <- as.matrix(exprs(apColonData))
res <- cor(edata0[ , 1:10] , use = "pairwise" , method = "pearson")
corrplot(res , type = "upper", order = "hclust" , tl.col = "blue" , tl.srt = 50 )
col <- colorRampPalette(c("dark blue" , "light blue" , "red"))(20)
heatmap(x = res ,col = col , symm = TRUE)

# genes: GSM95478,GSM95473 show the plot with a line of their relation.
edata1 <- as.matrix((edata))
lm1 <- lm(edata1[ ,"GSM95478"]~edata1[ ,"GSM95473"])
tidy(lm1)
plot(edata[ ,"GSM95478"]~edata[,"GSM95473"],col = 9)
abline(lm1 , col = 6 ,lwd = 3)


#svd
edata_centered <- edata - rowMeans(edata)             
svd1 <- svd(edata_centered)
names(svd1)
plot(svd1$v[,1],col = 9)  

#pca
pc1 <- prcomp(edata)
names(pc1)
plot(pc1$rotation[,1])  


plot(pc1$rotation[,1],svd1$v[,1])         #without normalization

edata_cenered2 <- t(t(edata)-colMeans(edata))         #with normalization
svd2 <- svd(edata_cenered2)
plot(pc1$rotation[,1], svd2$v[,1] , col = 2)


# hypothesis that zodiac signs 

zodiac<-as.factor(c(rep("Aries",29),rep("Taurus",24),rep("Gemini",22),rep("Cancer",19),rep("Leo",21),rep("Virgo",18),
                    rep("Libra",19),rep("Scorpio",20),rep("Sagittarius",23),rep("Capricorn",18),
                    rep("Aquarius",20),rep("Pisces",23)))
p<-c(1,1,1,1,1,1,1,1,1,1,1,1)
p<-p/sum(p)
table(zodiac)
chisq.test(table(zodiac),p=p)
print("accept null hypothesis and reject h1")

#we get chisquare value is 5.09 & p-value is >0.05 significance level there we accept the null hypothesis and 
#conclude that zodiac signs are evenly distributed across visual artists.

dist1 = dist(t(edata[,1:10]))
hclust1 = hclust(dist1)
plot(hclust1,hang = -1)
dev.off()
kmeans1 = kmeans(edata,centers=3)
z<-table(kmeans1$cluster)
z
write.table(z, file = "Result.txt")
