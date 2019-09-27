#' StalAge Age model 1.0
#'
#' StalAge is an age modelling software developed by Denis Scholz & Dirk Hoffmann
#' Contact: Denis Scholz, Institute for Geosciences, University of Mainz, Becherweg 21, 55128 Mainz, Germany
#' Email: scholzd@uni-mainz.de
#'
#' @param Daten Dating file; column names: age, error and depth
#' @param x_raw sample depths
#' @param Daten_orig
#' @return Returns the median age and 2.5 \% and 97.5 \% quantiles.
#'
#' @references
#' \insertRef{scholz2011stalage}{SISAL.AM}
#'
#' \insertRef{R}{bibtex}
age_model<-function(Daten, x_raw, Daten_orig) {

library(Hmisc)

#attach(Daten)
age <- Daten$age
error <- Daten$error
depth <- Daten$depth


iter<-300							#Anzahl der Iterationen

anz<-length(age)						#Bestimmung der Anzahl der Alter

l<-0								#counter

up<-0
low<-0

x_raw<-sort(x_raw)						#Sortiert Vektor mit Tiefen der Proxymessungen nach der Tiefe

if (x_raw[length(x_raw)]>depth[length(depth)]) up<-x_raw[length(x_raw)]+100 else up<-depth[length(depth)]+100	#Suche Maximum und Minimum des Tiefen- und Datierungsvektors
if (x_raw[1]<depth[1]) low<-x_raw[1]-100 else low<-depth[1]-100

x<-seq(from=low, to=up, by=2)					#Erzeuge Vektor mit geringerer Aufl?sung f?r Tiefen der Proxies

y<-array(NA, c(iter*(anz-2)*(anz-1)/2, length(x)))		#Entsprechende modellierte Alter

x_start<-0							#Punkt f?r ?berlapp am Anfang des Intervalls
x_end<-0							#Punkt f?r ?berlapp am Ende des Intervalls

dist1<-0							#Abstand zwischen Anfangspunkt des Intervalls und Punkt davor
dist2<-0							#Abstand zwischen Anfangspunkt des Intervalls und Punkt danach

inv<-0								#Variable zum Testen auf Inversionen

upper<-0							#Variable f?r lokales Maximum im Ergebnisvektor (Entfernen von Inversionen)
lower<-0							#Variable f?r lokales Minimum im Ergebnisvektor (Entfernen von Inversionen)

age_sim<-0							#Variable f?r simulierte Alter


results<-array(0, c((anz-2)*(anz-1)/2, 8))			#Definition des Ergebnis-Arrays (Array f?r "Screening der Daten")
results2<-array(NA, c(anz, (anz-2)*(anz-1)/2, iter))		#Definition des 2. Ergebnis-Arrays (Array f?r Ergebnisse der Simulationen)
results3<-array(0, c(anz, 3))					#Definition des 3. Ergebnis-Arrays (Array f?r Mittelwerte und Fehler)
results4<-array(0, c(length(x), 3))				#Definition des 4. Ergebnis-Arrays (Array f?r modellierte Mittelwerte und Fehler)
results5<-array(0, c(length(x)-89, 3))				#Definition des 5. Ergebnis-Arrays (verk?rztes Array f?r modellierte Mittelwerte und Fehler)
results6<-array(0, c(length(x)-89, 3))				#Definition des 6. Ergebnis-Arrays (verk?rztes Array f?r modellierte Mittelwerte und Fehler, Fehler werden als absolute Alter angegeben, nicht als +-)
results7<-array(0, c(length(x_raw), 3))				#Definition des 7. Ergebnis-Arrays (Array f?r Spline ?ber modellierte Mittelwerte und Fehler)
results8<-array(0, c(length(x_raw), 3))				#Definition des 8. Ergebnis-Arrays (Array f?r Spline ?ber modellierte Mittelwerte und Fehler, Inversionen korrigiert)


for (i in (anz:3)) {						#Iterationsschleife f?r Anzahl der ber?cksichtigten Punkte

for (j in 1:(anz-(i-1))) {					#Iterationsschleife f?r einzelne Fits

test<-0

fit<-lm(age[j:(j+(i-1))] ~ depth [j:(j+(i-1))], weights=1/error[j:(j+(i-1))])			#Fit ?ber Intervall

l<-l+1								#Counter eins hoch

age_fit<-0
age_fit[j:(j+(i-1))]<-fit$fitted.values

for (k in j:(j+(i-1))) {					#Test ob linearer Fit m?glich

if (age_fit[k]>age[k]+error[k] || age_fit[k]<age[k]-error[k]) {	#wenn Outlier, raus aus Schleife

test<-1
break

}

}

attr(fit$coefficients, "names")<-NULL				#Schreiben der Daten ins Ergebnis-Array

results[l, 1]<-l
results[l, 2]<-j
results[l, 3]<-j+(i-1)
results[l, 4]<-test
results[l, 5]<-depth[j]
results[l, 6]<-depth[j+(i-1)]
results[l, 7]<-fit$coefficients[1]
results[l, 8]<-fit$coefficients[2]


if (test == 0) {						#wenn kein Outlier, simuliere Alter und fitte diese

for (m in (1:iter)) {

for (k in j:(j+(i-1))) age_sim[k]<-rnorm(1, age[k], error[k]/2)	#Simulation der Alter

fit<-lm(age_sim[j:(j+(i-1))] ~ depth [j:(j+(i-1))], weights=1/error[j:(j+(i-1))])		#Fit ?ber simulierte Alter
attr(fit$coefficients, "names")<-NULL

if (fit$coefficients[2]>0) {					#nur positive Steigungen werden zugelassen

for (k in j:(j+(i-1))) results2[k, l, m]<-fit$coefficients[2]*depth[k]+fit$coefficients[1]	#Schreiben der Daten ins Feld


if (j==1) dist1<-depth[j]-x[1] else dist1<-depth[j]-depth[j-1] ### Anfangsalter (gemessen) muss anders behandelt werden
if (j+(i-1)==length(depth)) dist2<-x[length(x)]-depth[j+(i-1)] else dist2<-depth[j+(i-1)]-depth[j+(i-1)-1] ### letztes Alter auch - Endalter (gemessen)

x_start<-rnorm(1, depth[j]-dist1/2, dist1/6)
x_end<-rnorm(1, depth[j+(i-1)]+dist2/2, dist2/6)


if (j+(i-1)!=length(depth)) {

if (fit$coefficients[2]*x_end+fit$coefficients[1]>=min(age[(j+i):length(age)]+3*error[(j+i):length(age)]/2)) x_end<-(min(age[(j+i):length(age)]+3*error[(j+i):length(age)]/2)-fit$coefficients[1])/fit$coefficients[2]	#wenn oberes Ende des Fits > 3-sigma Obergrenze der n?chsten Punkte, schneide ab

if (x_end<x_start) next						#wenn Endpunkt durch obigen Test vor Startpunkt geschoben wird, ?berspringe den Fit

if (x_end>depth[j+i]) x_end<-depth[j+i]				#wenn Fit ?ber den n?chsten Punkt hinausgeht, schneide dort ab

if (x_end==depth[j+i] && fit$coefficients[2]*x_end+fit$coefficients[1]<=age[j+i]-3*error[j+i]/2) next	#wenn Fit n?chsten Punkt schneidet und < 3-sigma Untergrenze, ?berspringe den Fit

}

if (j!=1) {

if (fit$coefficients[2]*x_start+fit$coefficients[1]<=max(age[1:(j-1)]-3*error[1:(j-1)]/2)) x_start<-(max(age[1:(j-1)]-3*error[1:(j-1)]/2)-fit$coefficients[1])/fit$coefficients[2]	#wenn unteres Ende des Fits < 3-sigma Untergrenze der vorherigen Punkte, schneide ab

if (x_start>x_end) next						#wenn Startpunkt durch obigen Test vor Endpunkt geschoben wird, ?berspringe den Fit

if (x_start<depth[j-1]) x_start<-depth[j-1]			#wenn Fit ?ber den vorherigen Punkt hinausgeht, schneide dort ab

if (x_start==depth[j-1] && fit$coefficients[2]*x_start+fit$coefficients[1]>=age[j-1]+3*error[j-1]/2) next	#wenn Fit vorherigen Punkt schneidet und > 3-sigma Obergrenze, ?berspringe den Fit

}


#curve(fit$coefficients[2]*x+fit$coefficients[1], from=x_start, to=x_end, add=TRUE, col=i)	#Plotten des Fits

for (k in 1:length(x)) if (x[k]>=x_start && x[k]<=x_end) y[(l-1)*iter+m, k]<-fit$coefficients[2]*x[k]+fit$coefficients[1]

}

}

}

else next

}

}

#write.table(results, file = 'results.csv', col.names = F, row.names = F, sep = ',')

for (i in (1:anz)) {						#Berechnen der Mittelwerte und Konfidenzintervalle

results3[i,1]<-median(results2[i,,], na.rm=TRUE)
results3[i,2]<-quantile(results2[i,,], probs=0.975, na.rm=TRUE, names=FALSE)-results3[i,1]
results3[i,3]<-results3[i,1]-quantile(results2[i,,], probs=0.025, na.rm=TRUE, names=FALSE)

}

for (i in (1:length(x))) {

results4[i,1]<-median(y[,i], na.rm=TRUE)
results4[i,2]<-quantile(y[,i], probs=0.975, na.rm=TRUE, names=FALSE)-results4[i,1]
results4[i,3]<-results4[i,1]-quantile(y[,i], probs=0.025, na.rm=TRUE, names=FALSE)

}


x<-x[46:(length(x)-44)]


results5[, 1]<-results4[46:(length(results4[, 1])-44), 1]
results5[, 2]<-results4[46:(length(results4[, 2])-44), 2]
results5[, 3]<-results4[46:(length(results4[, 3])-44), 3]


results6[, 1]<-results5[, 1]
results6[, 2]<-results5[, 1]+results5[, 2]
results6[, 3]<-results5[, 1]-results5[, 3]

results7[,1]<-spline(x=x, y=results6[,1], xout=x_raw)$y					#Spline ?ber das Altersmodell
results7[,2]<-spline(x=x, y=results6[,2], xout=x_raw)$y					#Spline ?ber den oberen Fehler des Altersmodells
results7[,3]<-spline(x=x, y=results6[,3], xout=x_raw)$y					#Spline ?ber den unteren Fehler des Altersmodells


for (i in (1:(length(results7[,1])-1))) if (results7[i+1,1]<results7[i,1]) inv[i]<-i else inv[i]<-"NA"	#Testen auf Inversionen des Altersmodells und Schreiben in Vektor


upper[1]<-results7[1, 1]								#Start-Wert des Maximum-Vektors

lower[length(results7[,1])]<-results7[length(results7[,1]), 1]				#Start-Wert des Minimum-Vektors


for (i in 2:length(results7[,1])) {							#Durchgehen des Ergebnis-Vektors von unten

if (results7[i, 1]>upper[i-1]) upper[i]<-results7[i, 1] else upper[i]<-upper[i-1]	#Wenn nicht monoton, dann nutze lokales Maximum

}


for (i in (length(results7[,1])-1):1) {							#Durchgehen des Ergebnis-Vektors von oben

if (results7[i, 1]<lower[i+1]) lower[i]<-results7[i, 1] else lower[i]<-lower[i+1]	#Wenn nicht monoton, dann nutze lokales Minimum

}


results8[,1]<-(upper+lower)/2								#Bilde Mittelwert und schreibe in Results-Array


upper<-0
lower<-0

upper[1]<-results7[1, 2]								#Start-Wert des Maximum-Vektors

lower[length(results7[,2])]<-results7[length(results7[,2]), 2]				#Start-Wert des Minimum-Vektors


for (i in 2:length(results7[,2])) {							#Durchgehen des Ergebnis-Vektors von unten

if (results7[i, 2]>upper[i-1]) upper[i]<-results7[i, 2] else upper[i]<-upper[i-1]	#Wenn nicht monoton, dann nutze lokales Maximum

}


for (i in (length(results7[,2])-1):1) {							#Durchgehen des Ergebnis-Vektors von oben

if (results7[i, 2]<lower[i+1]) lower[i]<-results7[i, 2] else lower[i]<-lower[i+1]	#Wenn nicht monoton, dann nutze lokales Minimum

}


results8[,2]<-(upper+lower)/2								#Bilde Mittelwert und schreibe in Results-Array


upper<-0
lower<-0

upper[1]<-results7[1, 3]								#Start-Wert des Maximum-Vektors

lower[length(results7[,3])]<-results7[length(results7[,3]), 3]				#Start-Wert des Minimum-Vektors


for (i in 2:length(results7[,3])) {							#Durchgehen des Ergebnis-Vektors von unten

if (results7[i, 3]>upper[i-1]) upper[i]<-results7[i, 3] else upper[i]<-upper[i-1]	#Wenn nicht monoton, dann nutze lokales Maximum

}


for (i in (length(results7[,3])-1):1) {							#Durchgehen des Ergebnis-Vektors von oben

if (results7[i, 3]<lower[i+1]) lower[i]<-results7[i, 3] else lower[i]<-lower[i+1]		#Wenn nicht monoton, dann nutze lokales Minimum

}


results8[,3]<-(upper+lower)/2								#Bilde Mittelwert und schreibe in Results-Array


for (i in (1:(length(results8[,1])-1))) if (results8[i+1,1]<=results8[i,1]) inv[i]<-i else inv[i]<-"NA"	#Testen auf gleiche Werte Altersmodell und Schreiben in Vektor

i<-1
l<-0
start<-0
end<-0

m<-0
b<-0

while (i<length(inv)) {

if (inv[i]!="NA") {									#wenn gleicher Wert,

start<-i
l<-i+1

while (inv[l]!="NA") l<-l+1

end<-l+1										#suche Ende des Intervalls

i<-l-1


m<-(results8[end, 1]-results8[start, 1])/(x_raw[end]-x_raw[start])
b<-results8[end, 1]-m*x_raw[end]

results8[start:end, 1]<-m*x_raw[start:end]+b						#fitte linear ?ber problematisches Intervall

}

i<-i+1

}


for (i in (1:(length(results8[,1])-1))) if (results8[i+1,1]<=results8[i,1]) inv[i]<-i else inv[i]<-"NA"


final1<-results8[,1]
final2<-results8[,2]
final3<-results8[,3]


results<-data.frame(x=x_raw, y=final1, y_plus=final2, y_minus=final3)


#x11()
#errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
#title(main="Final age model with screened errors")
#lines(x_raw, final1, col="green")
#lines(x_raw, final2, col="red")
#lines(x_raw, final3, col="red")

#detach(Daten)

#attach(Daten_orig)
depth <- Daten_orig$depth
age <- Daten_orig$age
error <- Daten_orig$error


#x11()
#errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
#title(main="Final age model with original errors")
#lines(x_raw, final1, col="green")
#lines(x_raw, final2, col="red")
#lines(x_raw, final3, col="red")

#detach(Daten_orig)

#write.table(results2, file = 'results2.csv', col.names = F, row.names = F, sep = ',')
#write.table(y, file = 'y.csv', col.names = F, row.names = F, sep = ',')
#write.table(x, file = 'x.csv', col.names = F, row.names = F, sep = ',')

return(results)

}

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#delete<-function(panel) {

#choice<<-1

#panel

#}


#enlarge<-function(panel) {

#choice<<-2

#panel

#}


#close<-function(panel) {

#choice<<-3

#panel

#}


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Broadly scan dating file for reversals/outliers and increase error.
#'
#' @param depth Dating depths
#' @param age Dates
#' @param error Age uncertainties for each date
#' @return Modified dating file.
#' @references
#' \insertRef{scholz2011stalage}{SISAL.AM}
#'
#' \insertRef{R}{bibtex}
#'
#' @importFrom Rdpack reprompt

scan<-function(depth, age, error) {

library(Hmisc)
#library(rpanel)
library(plotrix)

#x11()
#errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
#title(main="Original age data")

hilf<-0

test<-array(NA, c(length(depth), length(depth)))	#Test-Array f?r Screening


age<-age[order(depth)]		#Sortiert die Vektoren aufsteigend nach der Tiefe
error<-error[order(depth)]
depth<-depth[order(depth)]


for (i in 1:length(depth)) {	#Schleife zum Testen auf Inversionen. Wenn Inversion, schreibe 1 in Matrix.

j<-i+1

while (j<=length(depth)) {

if (age[j]+error[j]<age[i]-error[i]) test[i, j]<-1 else test[i, j]<-0	#wenn 2-sigma Fehler nicht ?berlappen, schreibe 1 ins Feld

j<-j+1

}

}

#print(test)

max<-1		#Variable f?r das Maximum der Inversionen
max_pos<-0	#Variable f?r den Punkt des Maximums


repeat {

for (i in 1:length(depth)) {

if (sum(test[i,], na.rm=TRUE)>max) {	#Gehe Zeilen durch ... ## NA z?hlen nicht da na.rm = T

max<-sum(test[i,], na.rm=TRUE)
max_pos<-i

}

if (sum(test[,i], na.rm=TRUE)>max) {	#Gehe Spalten durch ...

max<-sum(test[,i], na.rm=TRUE)
max_pos<-i

}

}

if (max>1) {				#Wenn ein Punkt mit mehr als einem anderen nicht passt, Abfrage, was getan werden soll

#x11()
#errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
#title(main="Screening age data for major outliers")

#radius<-(depth[length(depth)]-depth[1])/(age[length(depth)]-age[1])*error[max_pos]
#draw.circle(x=depth[max_pos], y=age[max_pos], radius=radius, border="red")

#panel<-rp.control(size=c(500, 120), title=paste("Point", max_pos, "is an outlier! What do you want to do?"), max_pos=max_pos)

#rp.button(panel, action=delete, title="Delete!", quitbutton=TRUE, pos=c(0,0,500, 40))
#rp.button(panel, action=enlarge, title="Enlarge!", quitbutton=TRUE, pos=c(0,41,500, 40))
#rp.button(panel, action=close, title="Close", quitbutton=TRUE, pos=c(0,81,500, 40))

#rp.block(panel)

#if (choice==1) {			#Auswahl: L?sche Punkt!

#test[max_pos,]<-NA
#test[,max_pos]<-NA
#age[max_pos]<-NA

#errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
#title(main="Screening age data for major outliers")

#}

#if (choice==2) {			#Auswahl: Vergr??ere Fehler!

if (max(test[max_pos,], na.rm=TRUE)==1) { ## was soll es sonst sein????? es muss mindestens 2 einser haben, da max>1

	hilf<-age[max_pos]

	for (i in (max_pos+1):length(depth)) {

		if (is.na(test[max_pos, i])) next ### wann kann das auftreten ???? wir gehen hier reihen durch

		if (test[max_pos, i]==1 && age[i]<hilf) hilf<-age[i]	#; print(hilf)

		error[max_pos]<-age[max_pos]-hilf

	}

}

if (max(test[, max_pos], na.rm=TRUE)==1) {

	hilf<-age[max_pos]

	for (i in 1:(max_pos-1)) {

		if (is.na(test[i, max_pos])) next

		if (test[i, max_pos]==1 && age[i]>hilf) hilf<-age[i]	#; print(hilf)

		error[max_pos]<-hilf-age[max_pos]

	}

}

test[max_pos,]<-NA
test[,max_pos]<-NA

#errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
#title(main="Screening age data for major outliers")

#}

#if (choice==3) {			#Auswahl: Weder noch -> Hinweis, dass das nicht m?glich ist!

#rp.messagebox("This is not possible. Select another option, please!", title = "STOP!")

#}

}

if (max==1) break

max<-1

}


test[is.na(test)]<-0			#Ersetze NA's in Matrix durch Nullen.


#x11()
#errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")
#title(main="Age data screened for major outliers")


depth<-depth[!is.na(age)]		#L?sche rausgeworfene Punkte
error<-error[!is.na(age)]
age<-age[!is.na(age)]


Daten<-data.frame(depth=depth, age=age, error=error)	#Speichere Punkte in Dataframe

return(Daten)

}


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Determines slope of three adjacent dates to determine reverslas/outliers
#'
#' @inheritParams scan
#' @return Percentage of failed fits.
#'
#' @references
#' \insertRef{scholz2011stalage}{SISAL.AM}
#'
#' \insertRef{R}{bibtex}
slope<-function(depth, age, error) {

iter<-200
count<-0

age_sim<-0

for (i in 1:iter) {

for (k in 1:3) age_sim[k]<-rnorm(1, age[k], error[k]/2)	#Simulation der Alter

fit<-lm(age_sim[1:3] ~ depth [1:3])			#Fit ?ber simulierte Alter
attr(fit$coefficients, "names")<-NULL

if (fit$coefficients[2]>0) count<-count+1	#wenn Steigung positiv, z?hle eins hoch

}

#print(count/iter)

return(count/iter)

}


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Fine scan dating input for reversals/outliers and increase errors at tagged depths.
#'
#' @param dating_tb Table containing dates, errors and depths.
#' @return Modified dating file.
#'
#'  @references
#' \insertRef{scholz2011stalage}{SISAL.AM}
#'
#' \insertRef{R}{bibtex}
scan_fine<-function(Daten){

#attach(Daten)
age <- Daten$age
error <- Daten$error
depth <- Daten$depth


#x11()
#errbar(depth, age, age+error, age-error, xlab="Distance from top [mm]", ylab="Age [a]")			#Plotten der Daten
#title(main="Screening for minor outliers")

#print(error)

Anteil<-0							#Vektoren f?r das Testen der Fits
counter<-0
part<-0

part[1]<-1							#Definition des Beteiligungs-Vektors
part[length(depth)]<-1
part[2]<-2
part[length(depth)-1]<-2

for (i in 3:(length(depth)-2)) part[i]<-3

#print(part)


test<-0								#Pruefvariable f?r ?bergeordneten Test
pruef<-0							#Pruefvariable f?r andere Tests


while (test==0) {						#Schleife f?r Iteration bis Fehler alle passen

test<-1								#Setzen der allgemeinen Pr?fvariable auf 1

for (i in 1:length(depth)) Anteil[i]<-0				#Nullsetzen von Anteil
for (i in 1:length(depth)) counter[i]<-0			#Nullsetzen von counter


for (i in 1:(length(depth)-2)) {

#print(i)


fit<-lm(age[i:(i+2)] ~ depth [i:(i+2)], weights=1/error[i:(i+2)])	#Fit ?ber 3-Punkt-Intervall
attr(fit$coefficients, "names")<-NULL

#curve(fit$coefficients[2]*x+fit$coefficients[1], from=depth[i], to=depth[i+2], add=TRUE, col=i)	#Plotten des Fits

age_fit<-0							#Bestimmung des Alters des Fits an den jeweiligen Punkten
age_fit[i:(i+2)]<-fit$fitted.values


pruef<-0							#Null-Setzen der Pr?fvariable

for (k in i:(i+2)) {						#Test ob linearer Fit m?glich

if (age_fit[k]>age[k]+error[k] || age_fit[k]<age[k]-error[k]) pruef<-pruef+1	#wenn Punkt au?erhalb der Fehlergrenzen, z?hle Pruefvariable um 1 hoch

}

#print (pruef)

if (pruef>0) {

counter[i:(i+2)]<-counter[i:(i+2)]+1				#wenn der Fit nicht geht, setze counter um 1 hoch

} else {

if (fit$coefficients[2]<0) {					#wenn die Steigung des Fits negativ ist

if (slope(depth[i:(i+2)], age[i:(i+2)], error[i:(i+2)])<0.2) counter[i:(i+2)]<-counter[i:(i+2)]+1	#wenn weniger als 30% der Fits eine positive Steigung haben, setze counter um 1 hoch

}

}

}

#print(counter)

Anteil<-counter/part

#print(Anteil)

for (i in 1:length(depth)) {				#gehe Anteil nach 1ern durch

if (Anteil[i]==1) {					#wenn Anteil=1, dann vergr??ere Fehler um 10% und setze allgemeine Pr?fvariable auf 1

error[i]<-error[i]*1.1
test<-0

}

}

#errbar(depth, age, age+error, age-error, add=TRUE)		#Plotten der Daten

}

#x11()
#errbar(depth, age, age+error, age-error)			#Plotten der Daten
#title(main="Age data screened for minor outliers")

Daten<-data.frame(depth=depth, age=age, error=error)

#detach(Daten)

return(Daten)

}

