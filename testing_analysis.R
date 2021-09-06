### COVID  testing - Sao Gabriel da Cachoeira

ct <- read.csv("covid_testing.csv", header = TRUE, sep= ";")
ct$TIME_ELAPSED <- as.numeric(ct$TIME_ELAPSED)
ct <- na.omit(ct)
ct$TYPE <- ifelse(ct$TYPE=="0_TR_antigen", "Ag",
                  ifelse(ct$TYPE=="1_PCR", "PCR",
                         ifelse(ct$TYPE=="2_TR_anticorpos", "Ab", NA)))

ct$RESULT <- ifelse(ct$RESULT=="0_neg", 0,
                    ifelse(ct$RESULT=="1_pos", 1, NA))
ct$RESULT <- as.numeric(ct$RESULT)
ct$TEST_DATE <- as.Date(ct$TEST_DATE)
ct$DATE_FIRST_SYMPTOM <- as.Date(ct$DATE_FIRST_SYMPTOM)

ct <- ct[ct$TIME_ELAPSED<=14,] ### subset to symptoms reported up to 14 days

### epi weeks
install.packages("MMWRweek")
library(MMWRweek)



ct$epi_week_test <- MMWRweek(ct$TEST_DATE)$MMWRweek
ct$epi_week_test <- as.integer(ct$epi_week_test)


list("Ag"=summary(ct$TIME_ELAPSED[ct$TYPE=="Ag"]),
                      "PCR"=summary(ct$TIME_ELAPSED[ct$TYPE=="PCR"]),
                       "Ab"=summary(ct$TIME_ELAPSED[ct$TYPE=="Ab"]))

### TESTS DONE AND POSITIVITY
ct.ag <- length(ct$TYPE[ct$TYPE=="Ag"])
ct.ag1 <- sum(ct$RESULT[ct$TYPE=="Ag"]==1)
ct.pcr <- length(ct$TYPE[ct$TYPE=="PCR"])
ct.pcr1 <- sum(ct$RESULT[ct$TYPE=="PCR"]==1)
ct.ab <- length(ct$TYPE[ct$TYPE=="Ab"])
ct.ab1 <- sum(ct$RESULT[ct$TYPE=="Ab"]==1)

#Summary - ALL            
data.frame("Tests"=c("Ag", "PCR", "Ab"), "N_done"=c(ct.ag, ct.pcr, ct.ab), "N_pos"=c(ct.ag1, ct.pcr1, ct.ab1), "Positive_per100"=round(c(ct.ag1/ct.ag, ct.pcr1/ct.pcr, ct.ab1/ct.ab)*100, 1))
 


# POSITIVITY TIME
list("Ag"=summary(ct$TIME_ELAPSED[ct$TYPE=="Ag" & ct$RESULT==1]),
                        "PCR"=summary(ct$TIME_ELAPSED[ct$TYPE=="PCR" & ct$RESULT==1]),
                        "Ab"=summary(ct$TIME_ELAPSED[ct$TYPE=="Ab" & ct$RESULT==1]))

# Time since symptoms for negative tests
list("Ag"=summary(ct$TIME_ELAPSED[ct$TYPE=="Ag" & ct$RESULT==0]),
                        "PCR"=summary(ct$TIME_ELAPSED[ct$TYPE=="PCR" & ct$RESULT==0]),
                        "Ab"=summary(ct$TIME_ELAPSED[ct$TYPE=="Ab" & ct$RESULT==0]))


# Number of tests done by time since first symptoms
time.cat3 <- ifelse(ct$TIME_ELAPSED<=4, "0-4d",
                   ifelse(ct$TIME_ELAPSED>=5 & ct$TIME_ELAPSED<=9, "5-9d",
                          ifelse(ct$TIME_ELAPSED>=10 & ct$TIME_ELAPSED<=14, "10-14d", NA)))
time.cat3 <- factor(time.cat3, levels = c("0-4d", "5-9d", "10-14d"))
table(time.cat3)

# Number of tests done by epidemiological week bracket

time.catw <- ifelse(ct$epi_week_test<=3, "W1-3",
                    ifelse(ct$epi_week_test>=4 & ct$epi_week_test<=6, "W4-6",
                           ifelse(ct$epi_week_test>=7 & ct$epi_week_test<=9, "W7-9", NA)))
time.catw <- factor(time.catw, levels = c("W1-3", "W4-6", "W7-9"))
table(time.catw)

# Number of AG tests done epi week x time since first symptoms
addmargins(table(time.catw[ct$TYPE=="Ag"], time.cat3[ct$TYPE=="Ag"]))

# AG Test positive rate epi week x time since first symptoms
round(table(time.catw[ct$TYPE=="Ag" & ct$RESULT==1], time.cat3[ct$TYPE=="Ag" & ct$RESULT==1])/table(time.catw[ct$TYPE=="Ag"], time.cat3[ct$TYPE=="Ag"])*100,1)

# Number of PCR tests done epi week x time since first symptoms
addmargins(table(time.catw[ct$TYPE=="PCR"], time.cat3[ct$TYPE=="PCR"]))
# PCR Test positive rate epi week x time since first symptoms
round(table(time.catw[ct$TYPE=="PCR" & ct$RESULT==1], time.cat3[ct$TYPE=="PCR" & ct$RESULT==1])/table(time.catw[ct$TYPE=="PCR"], time.cat3[ct$TYPE=="PCR"])*100,1)

# Number of AB tests done epi week x time since first symptoms
addmargins(table(time.catw[ct$TYPE=="Ab"], time.cat3[ct$TYPE=="Ab"]))
# AB Test positive rate epi week x time since first symptoms
round(table(time.catw[ct$TYPE=="Ab" & ct$RESULT==1], time.cat3[ct$TYPE=="Ab" & ct$RESULT==1])/table(time.catw[ct$TYPE=="Ab"], time.cat3[ct$TYPE=="Ab"])*100,1)

# Pooled AG/PCR positivity by symptom time bracket
ct.pos1 <- tapply(ct$RESULT[ct$TYPE=="Ag" | ct$TYPE=="PCR"], time.cat3[ct$TYPE=="Ag" | ct$TYPE=="PCR"], sum)
data.frame(round(ct.pos1/table(time.cat3[ct$TYPE=="Ag" | ct$TYPE=="PCR"])*100,1))


# AG positivity by symptom time bracket
ct.pos.week.ag <- table(ct$RESULT[ct$TYPE=="Ag"], time.cat3[ct$TYPE=="Ag"])["1",]
ct.tests.week.ag <- table(time.cat3[ct$TYPE=="Ag"])
ag.week <- data.frame(round(ct.pos.week.ag/ct.tests.week.ag*100,1))
ag.week

# PCR positivity by symptom time bracket
ct.pos.week.pcr <- table(ct$RESULT[ct$TYPE=="PCR"], time.cat3[ct$TYPE=="PCR"])["1",]
ct.tests.week.pcr <- table(time.cat3[ct$TYPE=="PCR"])
pcr.week <- data.frame(round(ct.pos.week.pcr/ct.tests.week.pcr*100,1))
pcr.week

# AB positivity by symptom time bracket
ct.pos.week.ab <- table(ct$RESULT[ct$TYPE=="Ab"], time.cat3[ct$TYPE=="Ab"])["1",]
ct.tests.week.ab <- table(time.cat3[ct$TYPE=="Ab"])
ab.week <- data.frame(round(ct.pos.week.ab/ct.tests.week.ab*100,1))
ab.week

library(dplyr)



### Time from first symptoms to testing by week

for (week in sort(unique(ct$epi_week_test))) {print(cbind(week, summary(ct$TIME_ELAPSED[ct$epi_week_test==week])))
}

#### Interestingly, it seems patients presented earlier and earlier as weeks went by - effect of programming???



# AG
ct.ag.pos <- tapply(ct$RESULT[ct$TYPE=="Ag"], time.cat3[ct$TYPE=="Ag"], sum)
ct.ag.tests <- table(time.cat3[ct$TYPE=="Ag"])
data.frame(round(ct.ag.pos/ct.ag.tests*100,1))


## confidence intervals by period - AG
ciag <- function(x){
prop.test(ct.ag.pos[x], ct.ag.tests[x])$conf.int
}
for (x in c(1,2,3)){print(ciag(x))}

# PCR
ct.pcr.pos <- tapply(ct$RESULT[ct$TYPE=="PCR"], time.cat3[ct$TYPE=="PCR"], sum)
ct.pcr.tests <- table(time.cat3[ct$TYPE=="PCR"])
data.frame(round(ct.pcr.pos/ct.pcr.tests*100, 1))


## confidence intervals by period - PCR
cipcr <- function(x){
prop.test(ct.pcr.pos[x], ct.pcr.tests[x])$conf.int
}
for (x in c(1,2,3)){print(cipcr(x))}

# AB
ct.ab.pos <- tapply(ct$RESULT[ct$TYPE=="Ab"], time.cat3[ct$TYPE=="Ab"], sum)
ct.ab.tests <- table(time.cat3[ct$TYPE=="Ab"])
data.frame(round(ct.ab.pos/ct.ab.tests*100,1)) ## POSITIVITY BY PERIOD

## confidence intervals by period - AB

ciab <- function(x){
  prop.test(ct.ab.pos[x], ct.ab.tests[x])$conf.int
}
for (x in c(1,2,3)){print(ciab(x))}


### MULTIPLE TESTS
mt <- read.csv("multiple.csv", header = TRUE, sep=";")
mt$TIME_ELAPSED <- as.integer(mt$TIME_ELAPSED)
mt <- mt[mt$TIME_ELAPSED<=14,]


mt$AG_RESULT <- ifelse(mt$AG_RESULT=="Negativo", 0,
                       ifelse(mt$AG_RESULT=="Positivo", 1, 999))


mt$AB_RESULT <- ifelse(mt$AB_RESULT=="Negativo", 0,
                       ifelse(mt$AB_RESULT=="Positivo", 1, 999))
mt$PCR_RESULT <- ifelse(mt$PCR_RESULT=="Negativo", 0,
                       ifelse(mt$PCR_RESULT=="Positivo",1, 999))



### AG
m.ag <- length(mt$AG_RESULT[mt$AG_RESULT!=999])
m.ag.1 <- length(mt$AG_RESULT[mt$AG_RESULT==1])


### PCR
m.pcr <- length(mt$PCR_RESULT[mt$PCR_RESULT!=999])
m.pcr.1 <- length(mt$PCR_RESULT[mt$PCR_RESULT==1])


### AB
m.ab <- length(mt$AB_RESULT[mt$AB_RESULT!=999])
m.ab.1 <- length(mt$AB_RESULT[mt$AB_RESULT==1])



## SUMMARY - MULTIPLE TESTS
data.frame("Test"=c("AG", "PCR", "AB"), "N_done"=c(m.ag, m.pcr, m.ab), "N_positive"=c(m.ag.1, m.pcr.1, m.ab.1), "Pos_rate"=round(c(m.ag.1/m.ag, m.pcr.1/m.pcr, m.ab.1/m.ab)*100,1))

### Number of patients who got all 3 tests
m3 <- length(mt$ID[mt$AG_RESULT!=999 & mt$AB_RESULT!=999 & mt$PCR_RESULT!=999])

m3p <- length(mt$ID[mt$AG_RESULT==1 & mt$AB_RESULT==1 & mt$PCR_RESULT==1])
m3p/m3

## AT LEAST ONE POSITIVE
m1p <- length(mt$ID[mt$AG_RESULT==1 | mt$AB_RESULT==1 | mt$PCR_RESULT==1])


## BOTH AG AND PCR DONE
m.ag.p <- length(mt$ID[mt$AG_RESULT!=999 & mt$PCR_RESULT!=999])

## BOTH AG AND AB DONE
m.ag.ab <- length(mt$ID[mt$AG_RESULT!=999 & mt$AB_RESULT!=999])

## BOTH PCR AND AB DONE
m.p.ab<- length(mt$ID[mt$PCR_RESULT!=999 & mt$AB_RESULT!=999])

## EITHER AG OR PCR POSITIVE
m.ag1.or.p1 <- length(mt$ID[mt$AG_RESULT==1 | mt$PCR_RESULT==1])
m.ag1.or.p1/m.ag.p ### positivity if either positive
## AG POSITIVE & PCR NEGATIVE
m.ag1.p0 <- length(mt$ID[mt$AG_RESULT==1 & mt$PCR_RESULT==0])
m.ag1.p0/m.ag.p
## AG NEGATIVE & PCR POSITIVE
m.ag0.p1 <- length(mt$ID[mt$AG_RESULT==0 & mt$PCR_RESULT==1])
m.ag0.p1/m.ag.p
## AG POSITIVE & PCR POSITIVE
m.ag1.p1 <- length(mt$ID[mt$AG_RESULT==1 & mt$PCR_RESULT==1])
m.ag1.p1/m.ag.p


