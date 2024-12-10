## low birth weight Pakistan
## data for figures code

## import data
lbw <- read.csv("LBW_clean.csv") # this file can only be accessed upon request and with
                                 # only with simultaneous permission from DHS     
head(lbw)

# cases overall by month
xtabs(lbw$LBW_count ~ lbw$Month)

# cases overall by month by province
table(lbw$province)
xtabs(lbw$LBW_count ~ lbw$province)

BAL.case.mo <- xtabs(lbw[lbw$province == "Baluchistan",]$LBW_count ~ lbw[lbw$province == "Baluchistan",]$Month) # Baluchistan
GB.case.mo <- xtabs(lbw[lbw$province == "GB",]$LBW_count ~ lbw[lbw$province == "GB",]$Month) # GB
KPK.case.mo <- xtabs(lbw[lbw$province == "KPK",]$LBW_count ~ lbw[lbw$province == "KPK",]$Month) # KPK
PUN.case.mo <- xtabs(lbw[lbw$province == "Punjab",]$LBW_count ~ lbw[lbw$province == "Punjab",]$Month) # Punjab
SIN.case.mo <- xtabs(lbw[lbw$province == "Sindh",]$LBW_count ~ lbw[lbw$province == "Sindh",]$Month) # Sindh

case.mo.dat <- data.frame(as.numeric(BAL.case.mo), as.numeric(GB.case.mo), as.numeric(KPK.case.mo), as.numeric(PUN.case.mo), as.numeric(SIN.case.mo))
colnames(case.mo.dat) <- c("BAL", "GB", "KPK", "PUN", "SIN")
case.mo.dat

case.mo.ord <- data.frame(as.numeric(GB.case.mo), as.numeric(BAL.case.mo), as.numeric(KPK.case.mo), as.numeric(SIN.case.mo), as.numeric(PUN.case.mo))
colnames(case.mo.ord) <- c("GB", "BAL", "KPK", "SIN", "PUN")
case.mo.ord

cas.mo.cum <- data.frame(case.mo.ord[,1], rowSums(case.mo.ord[,1:2]), rowSums(case.mo.ord[,1:3]), rowSums(case.mo.ord[,1:4]), rowSums(case.mo.ord[,1:5]))
colnames(cas.mo.cum) <- c("GB", "BAL", "KPK", "SIN", "PUN")
cas.mo.cum


# heat index
GB <- subset(lbw, province == "GB")
GB.heat <- xtabs(GB$heat.index ~ GB$Month) / as.numeric(table(GB$Month))

BAL <- subset(lbw, province == "Baluchistan")
BAL.heat <- xtabs(BAL$heat.index ~ BAL$Month) / as.numeric(table(BAL$Month))

KPK <- subset(lbw, province == "KPK")
KPK.heat <- xtabs(KPK$heat.index ~ KPK$Month) / as.numeric(table(KPK$Month))

SIN <- subset(lbw, province == "Sindh")
SIN.heat <- xtabs(SIN$heat.index ~ SIN$Month) / as.numeric(table(SIN$Month))

PUN <- subset(lbw, province == "Punjab")
PUN.heat <- xtabs(PUN$heat.index ~ PUN$Month) / as.numeric(table(PUN$Month))

heat.mo.ord <- data.frame(as.numeric(GB.heat), as.numeric(BAL.heat), as.numeric(KPK.heat), as.numeric(SIN.heat), as.numeric(PUN.heat))
colnames(heat.mo.ord) <- c("GB", "BAL", "KPK", "SIN", "PUN")
heat.mo.ord

str(lbw)

# precip
table(GB$Year)
GB.precip <- xtabs(GB$total_precipitation ~ GB$Month) / as.numeric(table(GB$Month))
BAL.precip <- xtabs(BAL$total_precipitation ~ BAL$Month) / as.numeric(table(BAL$Month))
KPK.precip <- xtabs(KPK$total_precipitation ~ KPK$Month) / as.numeric(table(KPK$Month))
SIN.precip <- xtabs(SIN$total_precipitation ~ SIN$Month) / as.numeric(table(SIN$Month))
PUN.precip <- xtabs(PUN$total_precipitation ~ PUN$Month) / as.numeric(table(PUN$Month))
tot.precip <- 1000*xtabs(lbw$total_precipitation ~ lbw$Month) / as.numeric(table(lbw$Month)) # OVERALL
tot.precip
sum(tot.precip)

precip.mo.ord <- 1000*data.frame(as.numeric(GB.precip), as.numeric(BAL.precip), as.numeric(KPK.precip), as.numeric(SIN.precip), as.numeric(PUN.precip))
colnames(precip.mo.ord) <- c("GB", "BAL", "KPK", "SIN", "PUN")
precip.mo.ord
colSums(precip.mo.ord)
min(apply(precip.mo.ord, 1, min))
max(apply(precip.mo.ord, 1, max))


# air pollution (PM2.5)
GB.PM25 <- xtabs(GB$PM25 ~ GB$Month) / as.numeric(table(GB$Month))
BAL.PM25 <- xtabs(BAL$PM25 ~ BAL$Month) / as.numeric(table(BAL$Month))
KPK.PM25 <- xtabs(KPK$PM25 ~ KPK$Month) / as.numeric(table(KPK$Month))
SIN.PM25 <- xtabs(SIN$PM25 ~ SIN$Month) / as.numeric(table(SIN$Month))
PUN.PM25 <- xtabs(PUN$PM25 ~ PUN$Month) / as.numeric(table(PUN$Month))

xtabs(lbw$PM25 ~ lbw$Month) / as.numeric(table(lbw$Month)) # OVERALL

PM25.mo.ord <- data.frame(as.numeric(GB.PM25), as.numeric(BAL.PM25), as.numeric(KPK.PM25), as.numeric(SIN.PM25), as.numeric(PUN.PM25))
colnames(PM25.mo.ord) <- c("GB", "BAL", "KPK", "SIN", "PUN")
PM25.mo.ord

