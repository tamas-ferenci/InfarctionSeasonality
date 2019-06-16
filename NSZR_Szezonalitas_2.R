library( data.table )
library( Hmisc )
library( lattice )
library( mgcv )

setwd( "~/Egyetem/Kutatas/GOKI/Infarktusregiszter/Szezonalitas/" )

RawDataAdat <- fread( "2014_2017 AMI események előzményi adatokkal-Panasz kezdet.csv", dec = ",",
                      na.strings = c( "Nem ismert", "Nincs kitöltve" ) )
RawData <- fread( "2014_2017 AMI Módosított Kórelőzményi adatokkal.csv", dec = ",",
                  na.strings = c( "Nem ismert", "Nincs kitöltve" ) )
RawData <- merge( RawData, RawDataAdat[ , c( 3:5, 10, 12, 15:17, 28 ) ], by = names( RawData )[ 1:3 ] )
RawData$`(0) Halál dátuma` <- as.Date( RawData$`(0) Halál dátuma`, format = "%Y.%m.%d %H:%M" )
RawData$`(4) Születési dátum` <- as.Date( RawData$`(4) Születési dátum`, format = "%Y-%m-%d %H:%M:%S" )
RawData$`Esemény kezdete` <- as.Date( RawData$`Esemény kezdete`,
                                      format = "%Y-%m-%d %H:%M:%S" )
RawData$`(10) Panaszok fellépésének időpontja` <- as.POSIXct( RawData$`(10) Panaszok fellépésének időpontja`,
                                                              format = "%Y-%m-%d %H:%M:%S" )

RawData$SEX <- ifelse( RawData$`(2) Neme`=="Férfi", 1, 2 )
RawData$Életkor <- as.numeric( difftime( RawData$`Esemény kezdete`, RawData$`(4) Születési dátum`, unit = "days" ) )/365.24
RawData$AGE <- cut( RawData$Életkor, breaks = c( seq( 0, 85, 5 ), Inf ), labels = seq( 0, 85, 5 ) )
RawData$AGE <- as.numeric( levels( RawData$AGE ) )[ RawData$AGE ]
RawData$DATE <- RawData$`Esemény kezdete`
RawData$TYPE <- RawData$`(92) Diagnózis`

for( i in names( RawData )[ c( 20, 6:12, 16  ) ] ) {
  RawData[[ i ]] <- relevel( as.factor( RawData[[ i ]] ), ref = "Nem" )
}

### Descriptive statistics

write.csv2( print( summary( as.formula( paste0(
  "~`", paste0( names( RawData )[ c( 18, 25, 20, 22, 6:13, 16, 17  ) ], collapse = "`+`" ), "`" ) ),
  data = RawData ), prmsd = TRUE, digits = 3, pctdig = 1 ), "Table1.csv" )

mean( is.na( RawData$`(10) Panaszok fellépésének időpontja` ) )
RawData$MOD <- as.numeric( format( RawData$`(10) Panaszok fellépésének időpontja`, "%H" ) )*60 +
  as.numeric( format( RawData$`(10) Panaszok fellépésének időpontja`, "%M" ) )

p <- densityplot( ~MOD/60, groups = TYPE, data = RawData, plot.points = FALSE, ref = TRUE, from = 0, to = 24,
                  ylab = "", xlab = "Hour of day", auto.key = list( columns = 2 ),
                  scales = list( x = list( at = seq( 0, 24, 4 ) ), y = list( at = NULL ) ) )
cairo_pdf( "Figure1.pdf" )
p
dev.off()
cairo_pdf( "Figure1BW.pdf" )
update( p, par.settings =  standard.theme( color = FALSE ) )
dev.off()

### Incidence model

TimeStratified <- setkey( RawData, SEX, AGE, DATE, TYPE )[ CJ( unique( SEX ), seq( 0, 85, 5 ),
                                                               seq( as.Date( "2014-01-01" ), as.Date( "2017-12-31" ),
                                                                    by = 1 ), unique( TYPE ) ), .N, by = .EACHI ]

PopPyramid <- data.table( KSHStatinfoScraper::GetPopulationPyramidKSH( Years = 2014:2018, AgeGroup = "FiveYear",
                                                                       GeographicArea = "Total" ) )
PopPyramid <- PopPyramid[ , 1:4 ]
PopPyramid$SEX <- ifelse( PopPyramid$SEX=="Male", 1, 2 )

PopPyramid <- PopPyramid[ , with( approx( as.Date( paste0( YEAR, "-01-01" ) ), POPULATION,
                                          seq( as.Date( "2004-01-01" ), as.Date( "2017-12-31" ), by = 1 ) ),
                                  list( DATE = x, POPULATION = y ) ), .( AGE, SEX ) ]
PopPyramid <- PopPyramid[ , .( AGE, SEX, DATE, POPULATION = POPULATION/yearDays( DATE ) ) ]

TimeStratified <- merge( TimeStratified, PopPyramid, by = c( "SEX", "AGE", "DATE" ) )

TimeStratified <- cbind( TimeStratified, TemporalIndicators::TemporalIndicators( TimeStratified$DATE,
                                                                                 Congresses = "HUNCardiology" ) )

TimeStratified$DOW <- relevel( as.factor( TimeStratified$DOW ), ref = "1" )
TimeStratified$NWD <- relevel( as.factor( TimeStratified$NWD ), ref = "None" )
TimeStratified$SWAP <- relevel( as.factor( TimeStratified$SWAP ), ref = "None" )
TimeStratified$HUNCardiology <- relevel( as.factor( TimeStratified$HUNCardiology ), ref = "None" )
TimeStratified$SDAY <- relevel( as.factor( TimeStratified$SDAY ), ref = "None" )
TimeStratified$SEX <- as.factor( TimeStratified$SEX )

saveRDS( TimeStratified, "TimeStratified_NSZR_Szezonalitas.dat" )
TimeStratified <- readRDS( "TimeStratified_NSZR_Szezonalitas.dat" )

write.csv2( TimeStratified[ , .( sum( N )/sum( POPULATION )*100000 ), .( TYPE, SEASON ) ][ order( TYPE ) ], "BySeason.csv",
            row.names = FALSE )

p <- xyplot( Inc*100000 ~ WEEK | TYPE, groups = YEAR, type = "l", ylab = "Incidence [/100 000/year]", xlab = "Week",
             data = TimeStratified[ order( YEAR, WEEK ), .( Inc = sum( N )/sum( POPULATION ) ) , by = .( YEAR, WEEK, TYPE ) ],
             auto.key = list( columns = 4, points = FALSE, lines = TRUE ) )
cairo_pdf( "FigureS1.pdf", width = 10 )
p
dev.off()
cairo_pdf( "FigureS1BW.pdf" )
update( p, par.settings =  standard.theme( color = FALSE ) )
dev.off()

fmSTEMI <- gam( N ~ DOW + NWD + SWAP + HUNCardiology + SDAY + s( DATEnum, k = 50 ) +
                  te( DOY, AGE, bs = c( "cc", "tp" ), k = c( 50, 15 ), by = SEX ) + SEX,
                offset = log( POPULATION ), data = TimeStratified[ TYPE=="STEMI" ], family = nb( link = log ),
                method = "ML", knots = list( DOY = c( -(1/365.25)/2, 1+(1/365.25)/2 ) ) )

fmNSTEMI <- gam( N ~ DOW + NWD + SWAP + HUNCardiology + SDAY + s( DATEnum, k = 50 ) +
                   te( DOY, AGE, bs = c( "cc", "tp" ), k = c( 50, 15 ), by = SEX ) + SEX,
                 offset = log( POPULATION ), data = TimeStratified[ TYPE=="NSTEMI" ], family = nb( link = log ),
                 method = "ML", knots = list( DOY = c( -(1/365.25)/2, 1+(1/365.25)/2 ) ) )

fmSTEMI <- readRDS( "fmSTEMI.dat" )
fmNSTEMI <- readRDS( "fmNSTEMI.dat" )

gam.check( fmSTEMI )
gam.check( fmNSTEMI )

source( "temp.R" )
temp <- rbind( data.frame( Type = "STEMI", termplot2( fmSTEMI, terms = "DOW", se = TRUE, plot = FALSE ) ),
               data.frame( Type = "NSTEMI", termplot2( fmNSTEMI, terms = "DOW", se = TRUE, plot = FALSE ) ) )
temp$lwr <- exp( temp$DOW.y-1.96*temp$DOW.se )
temp$upr <- exp( temp$DOW.y+1.96*temp$DOW.se )
temp$DOW.y <- exp( temp$DOW.y )
names( temp )[ 2:4 ] <- c( "x", "y", "se" )

predgrid <- expand.grid( DOY = seq( 0, 1, by = 0.01 ), NWD = "None", SWAP = "None", HUNCardiology = "None", SDAY = "None",
                         DATEnum = as.numeric( as.Date( "2017-01-01" ) ), DOW = "1", SEX = 1:2, AGE = seq( 50, 80, 10 ) )
predgrid2 <- rbind( cbind( predgrid, transform( do.call( cbind, predict( fmSTEMI, predgrid, type = "link", se.fit = TRUE ) ),
                                                pred = fmSTEMI$family$linkinv( fit ),
                                                cilwr = fmSTEMI$family$linkinv( fit - 1.96*se.fit ),
                                                ciupr = fmSTEMI$family$linkinv( fit + 1.96*se.fit ), DG = "STEMI" ) ),
                    cbind( predgrid, transform( do.call( cbind, predict( fmNSTEMI, predgrid, type = "link", se.fit = TRUE ) ),
                                                pred = fmNSTEMI$family$linkinv( fit ),
                                                cilwr = fmNSTEMI$family$linkinv( fit - 1.96*se.fit ),
                                                ciupr = fmNSTEMI$family$linkinv( fit + 1.96*se.fit ), DG =  "NSTEMI" ) ) )
predgrid2$SEX <- ifelse( predgrid2$SEX==1, "Male", "Female" )

p <- xYplot( Cbind( pred*1e5, cilwr*1e5, ciupr*1e5 ) ~ DOY*12 | SEX + DG, groups = AGE, data = predgrid2, type = "l",
             method = "filled bands", col.fill = scales::alpha( trellis.par.get()$superpose.line$col, 0.1 ),
             xlab = "Month of year", ylab = "Incidence [/100,000/year]", layout = c( 4, 1 ) )
cairo_pdf( "Figure3.pdf" )
p
dev.off()
cairo_pdf( "Figure3BW.pdf" )
update( p, par.settings =  standard.theme( color = FALSE ), col.fill = "lightgray", label.curves = FALSE,
        auto.key = list( columns = 4, points = FALSE, lines = TRUE ) )
dev.off()

predgrid2 <- data.table( rbind( cbind( DG = "STEMI", predgrid, pred = predict( fmSTEMI, predgrid, type = "response" )*1e5 ),
                                cbind( DG = "NSTEMI", predgrid, pred = predict( fmNSTEMI, predgrid, type = "response" )*1e5 ) ) )
predgrid2$SEX <- ifelse( predgrid2$SEX==1, "Male", "Female" )
p <- xyplot( relative ~ DOY*12 | SEX + DG, groups = AGE, data = predgrid2[ , .( relative = pred/pred[1], DOY ),
                                                                           .( DG, SEX, AGE ) ],
        type = "l", xlab = "Month of year", ylab = "Relative change",
        auto.key = list( columns = 4, points = FALSE, lines = TRUE ) )
cairo_pdf( "Figure4.pdf" )
p
dev.off()
cairo_pdf( "Figure4BW.pdf" )
update( p, par.settings =  standard.theme( color = FALSE ) )
dev.off()

write.csv2( dcast( AGE ~ DG + SEX, data = predgrid2[ , .( max( pred )/min( pred ) ), .( AGE, SEX, DG ) ], value.var = "V1" ),
            "TroughToPeak.csv" )

TabGAM <- function( fit, df = FALSE, chisq = FALSE, eps = 0.001 ) {
  temp <- data.frame( anova( fit )$pTerms.table )
  temp$var <- rownames( temp )
  temp2 <- anova( fit )$p.table
  temp2 <- temp2[ rownames( temp2 )!="(Intercept)", ]
  temp <- data.frame( temp2, temp[ rep( 1:nrow( temp ), temp[ , 1 ] ), ] )
  temp2 <- anova( fit )$s.table
  temp2 <- data.frame( NA, NA, NA, NA, df = temp2[ , 1 ], Chi.sq = temp2[ , 3 ], p.value = temp2[ , 4 ], var = rownames( temp2 ) )
  names( temp2 ) <- names( temp )
  temp <- rbind( temp, temp2 )
  temp$Pr...z.. <- format.pval( temp$Pr...z.., eps = eps )
  temp$p.value <- format.pval( temp$p.value, eps = eps )
  if( !df ) temp <- temp[ , names( temp ) != "df" ]
  if( !chisq ) temp <- temp[ , names( temp ) != "Chi.sq" ]
  temp
}

write.csv2( cbind( TabGAM( fmSTEMI ), TabGAM( fmNSTEMI ) ), "TableS1.csv" )

### Mortality model

RawData[ `(0) Beteg jelenleg él-e`=="Nem"&is.na(RawData$`(0) Halál dátuma`) ]

RawData$STATUS <- !is.na( RawData$`(0) Halál dátuma` )
RawData$TIME <- as.numeric( dplyr::if_else( RawData$STATUS, RawData$`(0) Halál dátuma`, as.Date( "2018-05-10" ) )-
                              RawData$`Esemény kezdete` )
RawData[ TIME<0 ]

apply( RawData, 2, function( x ) mean( is.na( x ) ) )

RawData <- cbind( RawData, TemporalIndicators::TemporalIndicators( RawData$DATE, Congresses = "HUNCardiology", NWD = "Global" ) )

RawData$DOW <- relevel( as.factor( RawData$DOW ), ref = "1" )
RawData$NWD <- relevel( as.factor( RawData$NWD ), ref = "None" )
RawData$SWAP <- relevel( as.factor( RawData$SWAP ), ref = "None" )
RawData$HUNCardiology <- relevel( as.factor( RawData$HUNCardiology ), ref = "None" )
RawData$SDAY <- relevel( as.factor( RawData$SDAY ), ref = "None" )
RawData$SEASON <- factor( RawData$SEASON, levels = 1:4, labels = c( "Winter", "Spring", "Summer", "Autumn" ) )

RawData$TYPE <- as.factor( RawData$TYPE )

RawData$AGE <- RawData$Életkor
RawData$PCI <- as.factor( RawData$`(0) Esemény szintű PCI történt-e` )
RawData$PrevAMI <- as.factor( RawData$`(28) Kórelőzményben myocardialis infarctus` )
RawData$PrevHF <- as.factor( RawData$`(29) Kórelőzményben szívelégtelenség` )
RawData$HT <- as.factor( RawData$`(30) Kórelőzményben vagy a kezelés során megállapított hypertonia` )
RawData$STROKE <- as.factor( RawData$`(31) Kórelőzményben stroke` )
RawData$DM <- as.factor( RawData$`(32) Kórelőzményben, vagy a kezelés során megállapított diabetes` )
RawData$PAD <- as.factor( RawData$`(33) Kórelőzményben perifériális érbetegség` )

p <- xYplot( Cbind( surv*100, lower*100, upper*100 ) ~ time, groups = rep( sapply( strsplit( names( strata ), "=" ), `[`, 2 ),
                                                                           strata ),
             data = survfit( Surv( TIME, STATUS ) ~ SEASON, data = RawData ), xlab = "Time [days]", ylab = "Survival [%]",
             method = "filled bands", col.fill = scales::alpha( trellis.par.get()$superpose.line$col, 0.1 ), type = "l",
             ylim = c( 65, 101 ), xlim = c( 0, 3*365 ) )
cairo_pdf( "FigureS2.pdf" )
p
dev.off()
cairo_pdf( "FigureS2BW.pdf" )
update( p, par.settings =  standard.theme( color = FALSE ), col.fill = "lightgray", label.curves = FALSE,
        auto.key = list( columns = 4, points = FALSE, lines = TRUE ) )
dev.off()

fit <- gam( TIME ~ s( DOY, k = 50, bs = "cc" ), data = RawData, 
            family = cox.ph(), weights = STATUS, knots = list( DOY = c( -(1/365.25)/2, 1+(1/365.25)/2 ) ) )
summary( fit )
plot( fit, trans = exp )
fit <- gam( TIME ~ s( DOY, k = 50, bs = "cc" ) + s( AGE, k = 50 ), data = RawData, 
            family = cox.ph(), weights = STATUS, knots = list( DOY = c( -(1/365.25)/2, 1+(1/365.25)/2 ) ) )
summary( fit )
plot( fit, trans = exp, select = 1, scale = 0 )

p <- xYplot( Cbind( surv*100, lower*100, upper*100 ) ~ time | factor(rep( rep( c( 0, 40, 50, 60, 70, 80 ), 4 ), strata ) ),
             groups = rep( rep( levels( as.factor( RawData$SEASON ) ), each = 6 ), strata ),
             xlab = "Time [days]", ylab = "Survival [%]",
             data = survfit( Surv( TIME, STATUS ) ~ SEASON + cut( AGE, c( -Inf, 40, 50, 60, 70, 80, Inf ) ), data = RawData ),
             as.table = TRUE, method = "filled bands", col.fill = "lightgray", type = "l", xlim = c( 0, 3*365 ),
             label.curves = FALSE, scales = list( relation = "free" ),
             auto.key = list( columns = 4, points = FALSE, lines = TRUE ) )
cairo_pdf( "FigureS3.pdf" )
p
dev.off()
cairo_pdf( "FigureS3BW.pdf" )
update( p, par.settings =  standard.theme( color = FALSE ), col.fill = "lightgray", label.curves = FALSE )
dev.off()

p <- densityplot( ~ AGE, groups = SEASON, data = RawData, plot.points = FALSE, auto.key = list( columns = 4 ),
                  ylab = "", scales = list( y = list( at = NULL ) ) )
cairo_pdf( "FigureS4.pdf" )
p
dev.off()
cairo_pdf( "FigureS4BW.pdf" )
update( p, par.settings =  standard.theme( color = FALSE ), col.fill = "lightgray", label.curves = FALSE )
dev.off()

fit <- gam( TIME ~ DOW + NWD + SWAP + HUNCardiology + SDAY + s( DATEnum, k = 50 ) +
              s( DOY, k = 50, bs = "cc" ) + s( AGE, k = 50, by = SEX ) + SEX +
              PCI + PrevAMI + PrevHF + HT + STROKE + DM, data = RawData, 
            family = cox.ph(), weights = STATUS, knots = list( DOY = c( -(1/365.25)/2, 1+(1/365.25)/2 ) ) )

gam.check( fit )

temp2 <- termplot2( fit, terms = "DOW", se = TRUE, plot = FALSE ) 
temp2 <- temp2$DOW
temp2$lwr <- exp( temp2$y-1.96*temp2$se )
temp2$upr <- exp( temp2$y+1.96*temp2$se )
temp2$y <- exp( temp2$y )

temp3 <- rbind( cbind( temp, Model = "Incidence" ), cbind( Type = "Mortality", temp2, Model = "Mortality" ) )

p <- xYplot( Cbind( y, lwr, upr ) ~ rep( 1:7, 3 ) | Model, groups = Type, data = temp3, type = "l", ylim = c( 0.5, 1.3 ),
        xlab = "", ylab = "Rate ratio",
        scales = list( x = list( at = 1:7, labels = c( "Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
                                                       "Saturday", "Sunday" ) ) ),
        panel = function( ... ) {
          panel.xYplot( ..., label.curves = if( panel.number()==1 ) TRUE else FALSE )
        } )
cairo_pdf( "Figure2.pdf", width = 10 )
p
dev.off()
cairo_pdf( "Figure2BW.pdf", width = 10 )
update( p, par.settings =  standard.theme( color = FALSE ) )
dev.off()

write.csv2( TabGAM( fit ), "TableS2.csv" )
