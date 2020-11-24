# Loading the R packages dtw and astrochron----
library(dtw)
library(astrochron)

# Set working directory ----
setwd("/Users/pdoc3/Dropbox/DAAD - DANA/DAAD log DATA/") 

# Reading the U1463 downhole logging record in the depth domain ----
U1463=read.delim("U1463B_logging_fullcounts.txt")
U1463[c(1:512),2]=(U1463[c(1:512),2]+0.77)*4
plot(U1463,type="l")

##################  
#    Finucane    #
################## -----
Finucane=read.delim("finucane_1_20200713.txt", header = F)
Finucane=Finucane[c(1:4800),]
plot(Finucane, type = "l")
## Find the best match with the recursion formula
alignment <- dtw(Finucane[,2],U1463[,2], step.pattern=symmetricP1, window.type = "sakoechiba", window.size = 2000, keep = TRUE)
## Display the warping curve, i.e. the alignment curve
plot(alignment,type="threeway")
agemodel=cbind(alignment$index1, alignment$index2)
agemodel=sortNave(agemodel)
agemodel[,1]=Finucane[,1]
agemodel[,2]=U1463[agemodel[,2],1]
Finucane_tuned=tune(Finucane, agemodel)

#Plotting Finucane ----
dev.off()
pdf(file = "Finucane.pdf", width = 8.27, height = 11, paper = "a4")
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), widths=c(1,3), heights=c(0.75,1.25,0.5))
par(mai=c(0.75,0.75,0.25,0.25))

plot(U1463, type="l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
lines(Finucane_tuned, col = "darkorange", lwd = 1.2)
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.75,0.25,0.25))
plot(Finucane[,2],Finucane[,1], col = "darkorange", type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(Finucane[,1])),floor(min(Finucane[,1]))), xlim = c(0,50), xaxs = "i", yaxs = "i")
axis(2,cex.axis=1.5, at = c(floor(min(Finucane[,1])),200,400,600, ceiling(max(Finucane[,1]))))
mtext("Finucane Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,50))
mtext("NGR (gAPI)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.25,0.25,0.25))
plot(agemodel[,2], agemodel[,1], type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(Finucane[,1])),floor(min(Finucane[,1]))), xlim = c(0,433), xaxs = "i", yaxs = "i", lwd = 2)
grid(nx = NULL, ny = NULL , col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
axis(2,cex.axis=1.5, at = c(floor(min(Finucane[,1])),200,400,600, ceiling(max(Finucane[,1]))))
#mtext("Finucane Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
#mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)
text(300,230,"Dynamical Time Warping", cex = 2, col = "darkorange")
text(300,260,"between U1463 and Finucane", cex = 2, col = "darkorange",)

par(mai=c(0.75,0.25,0,0.25))
plot.new()

plot(U1463, type = "l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

dev.off()

# Finucane output ----
Finucane_output=cbind(Finucane[,1],Finucane_tuned)
colnames(Finucane_output)<- c("depth Finucane", "depth U1463", "gamma")
write.csv(Finucane_output, "Finucane_output.csv", row.names = F)

labelsFinucane=approx(agemodel[,2],agemodel[,1],xout = c(100,160,270,380))$y
labelsFinucane=cbind(c(100,160,270,380),round(labelsFinucane))
write.csv(labelsFinucane, "labelsFinucane.csv", row.names = F)

##################  
#   Goodwyn6     #
##################
Goodwyn6=read.delim("goodwyn_6.txt")
Goodwyn6=Goodwyn6[c(1:4215),c(2,7)]
plot(Goodwyn6,type="l")
## Find the best match with the recursion formula
alignment <- dtw(Goodwyn6[,2],U1463[,2], step.pattern=symmetricP1, window.type = "sakoechiba", window.size = 1500, keep = TRUE)
## Display the warping curve, i.e. the alignment curve
plot(alignment,type="threeway")
agemodel=cbind(alignment$index1, alignment$index2)
agemodel=sortNave(agemodel)
agemodel[,1]=Goodwyn6[,1]
agemodel[,2]=U1463[agemodel[,2],1]
Goodwyn6_tuned=tune(Goodwyn6, agemodel)

#Plotting Goodwyn6 ----
dev.off()
pdf(file = "Goodwyn6.pdf", width = 8.27, height = 11, paper = "a4")
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), widths=c(1,3), heights=c(0.75,1.25,0.5))
par(mai=c(0.75,0.75,0.25,0.25))

plot(U1463, type="l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
lines(Goodwyn6_tuned, col = "darkorange", lwd = 1.2)
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.75,0.25,0.25))
plot(Goodwyn6[,2],Goodwyn6[,1], col = "darkorange", type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(Goodwyn6[,1])),floor(min(Goodwyn6[,1]))), xlim = c(0,50), xaxs = "i", yaxs = "i")
axis(2,cex.axis=1.5, at = c(floor(min(Goodwyn6[,1])),200,400,600, ceiling(max(Goodwyn6[,1]))))
mtext("Goodwyn6 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,50))
mtext("NGR (gAPI)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.25,0.25,0.25))
plot(agemodel[,2], agemodel[,1], type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(Goodwyn6[,1])),floor(min(Goodwyn6[,1]))), xlim = c(0,433), xaxs = "i", yaxs = "i", lwd = 2)
grid(nx = NULL, ny = NULL , col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
axis(2,cex.axis=1.5, at = c(floor(min(Goodwyn6[,1])),200,400,600, ceiling(max(Goodwyn6[,1]))))
#mtext("Goodwyn6 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
#mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)
text(300,230,"Dynamical Time Warping", cex = 2, col = "darkorange")
text(300,260,"between U1463 and Goodwyn6", cex = 2, col = "darkorange",)

par(mai=c(0.75,0.25,0,0.25))
plot.new()

plot(U1463, type = "l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

dev.off()

#Goodwyn6 output ----
Goodwyn6_output=cbind(Goodwyn6[,1],Goodwyn6_tuned)
colnames(Goodwyn6_output)<- c("depth Goodwyn6", "depth U1463", "gamma")
write.csv(Goodwyn6_output, "Goodwyn6_output.csv", row.names = F)

labelsGoodwyn6=approx(agemodel[,2],agemodel[,1],xout = c(100,160,270,380))$y
labelsGoodwyn6=cbind(c(100,160,270,380),round(labelsGoodwyn6))
write.csv(labelsGoodwyn6, "labelsGoodwyn6.csv", row.names = F)

##################  
#    Goodwyn2    #
##################
Goodwyn2=read.delim("goodwyn_2.txt")
Goodwyn2=Goodwyn2[c(1:4595),c(2,7)]
plot(Goodwyn2,type="l")
## Find the best match with the recursion formula
alignment <- dtw(Goodwyn2[,2],U1463[,2], step.pattern=symmetricP1, window.type = "sakoechiba", window.size = 1800, keep = TRUE)
## Display the warping curve, i.e. the alignment curve
plot(alignment,type="threeway")
agemodel=cbind(alignment$index1, alignment$index2)
agemodel=sortNave(agemodel)
agemodel[,1]=Goodwyn2[,1]
agemodel[,2]=U1463[agemodel[,2],1]
Goodwyn2_tuned=tune(Goodwyn2, agemodel)

#Plotting Goodwyn2 ----
dev.off()
pdf(file = "Goodwyn2.pdf", width = 8.27, height = 11, paper = "a4")
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), widths=c(1,3), heights=c(0.75,1.25,0.5))
par(mai=c(0.75,0.75,0.25,0.25))

plot(U1463, type="l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
lines(Goodwyn2_tuned, col = "darkorange", lwd = 1.2)
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.75,0.25,0.25))
plot(Goodwyn2[,2],Goodwyn2[,1], col = "darkorange", type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(Goodwyn2[,1])),floor(min(Goodwyn2[,1]))), xlim = c(0,50), xaxs = "i", yaxs = "i")
axis(2,cex.axis=1.5, at = c(floor(min(Goodwyn2[,1])),200,400,600,800, ceiling(max(Goodwyn2[,1]))))
mtext("Goodwyn2 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,50))
mtext("NGR (gAPI)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.25,0.25,0.25))
plot(agemodel[,2], agemodel[,1], type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(Goodwyn2[,1])),floor(min(Goodwyn2[,1]))), xlim = c(0,433), xaxs = "i", yaxs = "i", lwd = 2)
grid(nx = NULL, ny = NULL , col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
axis(2,cex.axis=1.5, at = c(floor(min(Goodwyn2[,1])),200,400,600,800, ceiling(max(Goodwyn2[,1]))))
#mtext("Goodwyn2 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
#mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)
text(300,230,"Dynamical Time Warping", cex = 2, col = "darkorange")
text(300,260,"between U1463 and Goodwyn2", cex = 2, col = "darkorange",)

par(mai=c(0.75,0.25,0,0.25))
plot.new()

plot(U1463, type = "l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

dev.off()

#Goodwyn2 output ----
Goodwyn2_output=cbind(Goodwyn2[,1],Goodwyn2_tuned)
colnames(Goodwyn2_output)<- c("depth Goodwyn2", "depth U1463", "gamma")
write.csv(Goodwyn2_output, "Goodwyn2_output.csv", row.names = F)

labelsGoodwyn2=approx(agemodel[,2],agemodel[,1],xout = c(100,160,270,380))$y
labelsGoodwyn2=cbind(c(100,160,270,380),round(labelsGoodwyn2))
write.csv(labelsGoodwyn2, "labelsGoodwyn2.csv", row.names = F)

##################  
#     Angel2     #
##################
Angel2=read.csv("/Users/pdoc3/Dropbox/DAAD - DANA/DAAD log DATA/angel2.csv")
Angel2=Angel2[c(1:3400),c(1,2)]
plot(Angel2,type="l")
## Find the best match with the recursion formula
alignment <- dtw(Angel2[,2],U1463[,2], step.pattern=symmetricP1, window.type = "sakoechiba", window.size = 1800, keep = TRUE)
## Display the warping curve, i.e. the alignment curve
plot(alignment,type="threeway")
agemodel=cbind(alignment$index1, alignment$index2)
agemodel=sortNave(agemodel)
agemodel[,1]=Angel2[,1]
agemodel[,2]=U1463[agemodel[,2],1]
Angel2_tuned=tune(Angel2, agemodel)

#Plotting Angel2 ----
dev.off()
pdf(file = "Angel2.pdf", width = 8.27, height = 11, paper = "a4")
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), widths=c(1,3), heights=c(0.75,1.25,0.5))
par(mai=c(0.75,0.75,0.25,0.25))

plot(U1463, type="l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
lines(Angel2_tuned, col = "darkorange", lwd = 1.2)
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.75,0.25,0.25))
plot(Angel2[,2],Angel2[,1], col = "darkorange", type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(Angel2[,1])),floor(min(Angel2[,1]))), xlim = c(0,50), xaxs = "i", yaxs = "i")
axis(2,cex.axis=1.5, at = c(floor(min(Angel2[,1])),200,400,600,800, ceiling(max(Angel2[,1]))))
mtext("Angel2 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,50))
mtext("NGR (gAPI)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.25,0.25,0.25))
plot(agemodel[,2], agemodel[,1], type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(Angel2[,1])),floor(min(Angel2[,1]))), xlim = c(0,433), xaxs = "i", yaxs = "i", lwd = 2)
grid(nx = NULL, ny = NULL , col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
axis(2,cex.axis=1.5, at = c(floor(min(Angel2[,1])),200,400,600,800, ceiling(max(Angel2[,1]))))
#mtext("Angel2 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
#mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)
text(300,230,"Dynamical Time Warping", cex = 2, col = "darkorange")
text(300,260,"between U1463 and Angel2", cex = 2, col = "darkorange",)

par(mai=c(0.75,0.25,0,0.25))
plot.new()

plot(U1463, type = "l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

dev.off()

#Angel2 output ----
Angel2_output=cbind(Angel2[,1],Angel2_tuned)
colnames(Angel2_output)<- c("depth Angel2", "depth U1463", "gamma")
write.csv(Angel2_output, "Angel2_output.csv", row.names = F)

labelsAngel2=approx(agemodel[,2],agemodel[,1],xout = c(100,160,270,380))$y
labelsAngel2=cbind(c(100,160,270,380),round(labelsAngel2))
write.csv(labelsAngel2, "labelsAngel2.csv", row.names = F)

##################  
#    U1464       #
##################
U1464=read.delim("U1464C.txt")
U1464[c(1:553),2]=(U1464[c(1:553),2]+1.2)*2
U1464=U1464[c(1:2000),]
# At U1464, we first do "rough" depth-matching between U1464 and U1463, based on a priori knowledge, to facilitate the dtw algorithm
depthmatch1=matrix(c(304.8,100,0,432.9684,280,0), nrow = 3, ncol = 2)
U1464_dm1=tune(U1464,depthmatch1)
U1464_dm1=linterp(U1464_dm1, dt = 0.1)
# Subsequently, we continue the "normal" dtw procedure with the depth-matched U1464 series
alignment <- dtw(U1464_dm1[,2],U1463[,2], step.pattern=symmetricP1, window.type = "sakoechiba", window.size = 1500, keep = TRUE)
plot(alignment,type="threeway")
agemodel=cbind(alignment$index1, alignment$index2)
agemodel=sortNave(agemodel)
agemodel[,1]=U1464_dm1[,1]
agemodel[,2]=U1463[agemodel[,2],1]
# We have to correct for the depthmatching to obtain the relationship between the original U1464 depths and the corresponding depths in U1463 
agemodel=tune(agemodel, depthmatch1[,c(2,1)])
U1464_tuned=tune(U1464, agemodel, extrapolate = T)

#Plotting U1464 ----
dev.off()
pdf(file = "U1464.pdf", width = 8.27, height = 11, paper = "a4")
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), widths=c(1,3), heights=c(0.75,1.25,0.5))
par(mai=c(0.75,0.75,0.25,0.25))

plot(U1463, type="l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
lines(U1464_tuned, col = "darkorange", lwd = 1.2)
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.75,0.25,0.25))
plot(U1464[,2],U1464[,1], col = "darkorange", type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(U1464[,1])),floor(min(U1464[,1]))), xlim = c(0,50), xaxs = "i", yaxs = "i")
axis(2,cex.axis=1.5, at = c(floor(min(U1464[,1])),100,200,300, ceiling(max(U1464[,1]))))
mtext("U1464 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,50))
mtext("NGR (gAPI)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.25,0.25,0.25))
plot(agemodel[,2], agemodel[,1], type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(U1464[,1])),floor(min(U1464[,1]))), xlim = c(0,433), xaxs = "i", yaxs = "i", lwd = 2)
grid(nx = NULL, ny = NULL , col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
axis(2,cex.axis=1.5, at = c(floor(min(U1464[,1])),100,200,300, ceiling(max(U1464[,1]))))
#mtext("U1464 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
#mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)
text(300,50,"Dynamical Time Warping", cex = 2, col = "darkorange")
text(300,75,"between U1463 and U1464", cex = 2, col = "darkorange",)

par(mai=c(0.75,0.25,0,0.25))
plot.new()

plot(U1463, type = "l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

dev.off()

#U1464 output ----
U1464_output=cbind(U1464[,1],U1464_tuned)
colnames(U1464_output)<- c("depth U1464", "depth U1463", "gamma")
write.csv(U1464_output, "U1464_output.csv", row.names = F)

labelsU1464=approx(agemodel[,2],agemodel[,1],xout = c(100,160,270,380))$y
labelsU1464=cbind(c(100,160,270,380),round(labelsU1464))
write.csv(labelsU1464, "labelsU1464.csv", row.names = F)

##################  
#   U1462       #
##################
U1462=read.delim("/Users/pdoc3/Dropbox/DAAD - DANA/DAAD log DATA/U1462C.txt")
plot(U1462, type = "l")
# At U1462, we first do "rough" depth-matching between U1462 and U1463, based on a priori knowledge, to facilitate the dtw algorithm
depthmatch2=matrix(c(857.25,703.81,539.59,270,0,432.9684,335.8,255.21,50,0),nrow = 5, ncol = 2)
U1462_dm2=tune(U1462, depthmatch2)
U1462_dm2=linterp(U1462_dm2, dt = 0.1)
# Subsequently, we continue the "normal" dtw procedure with the depth-matched U1464 series
alignment <- dtw(U1462_dm2[,2],U1463[,2], step.pattern=symmetricP05, window.type = "sakoechiba", window.size = 1500,  keep = TRUE)
plot(alignment,type="threeway")

agemodel=cbind(alignment$index1, alignment$index2)
agemodel=sortNave(agemodel)
agemodel[,1]=U1462_dm2[,1]
agemodel[,2]=U1463[agemodel[,2],1]

agemodel=tune(agemodel, depthmatch2[,c(2,1)])
U1462_tuned=tune(U1462, agemodel)

#Plotting U1462 ----
dev.off()
pdf(file = "U1462.pdf", width = 8.27, height = 11, paper = "a4")
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), widths=c(1,3), heights=c(0.75,1.25,0.5))
par(mai=c(0.75,0.75,0.25,0.25))

plot(U1463, type="l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
lines(U1462_tuned, col = "darkorange", lwd = 1.2)
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.75,0.25,0.25))
plot(U1462[,2],U1462[,1], col = "darkorange", type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(U1462[,1])),floor(min(U1462[,1]))), xlim = c(0,50), xaxs = "i", yaxs = "i")
axis(2,cex.axis=1.5, at = c(floor(min(U1462[,1])),200,400,600,800, ceiling(max(U1462[,1]))))
mtext("U1462 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,50))
mtext("NGR (gAPI)", side=1, line=3, cex=1.2)

par(mai=c(0.5,0.25,0.25,0.25))
plot(agemodel[,2], agemodel[,1], type = "l", yaxt="n", xaxt = "n", ylab = "" , xlab = "", ylim = c(ceiling(max(U1462[,1])),floor(min(U1462[,1]))), xlim = c(0,433), xaxs = "i", yaxs = "i", lwd = 2)
grid(nx = NULL, ny = NULL , col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
axis(2,cex.axis=1.5, at = c(floor(min(U1462[,1])),200,400,600,800, ceiling(max(U1462[,1]))))
#mtext("U1462 Wireline Depth (m)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
#mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)
text(300,50,"Dynamical Time Warping", cex = 2, col = "darkorange")
text(300,100,"between U1463 and U1462", cex = 2, col = "darkorange",)

par(mai=c(0.75,0.25,0,0.25))
plot.new()

plot(U1463, type = "l", lwd = 1.2, yaxt="n", xaxt = "n", ylab = "", xlab = "", xlim = c(0,433), ylim = c(0,50), xaxs = "i", yaxs = "i", col = "gray42")
axis(2,cex.axis=1.5,  at = c(0,50))
mtext("NGR (gAPI)", side=2, line=3, cex=1.2)
axis(1,cex.axis=1.5, at = c(0,100,200,300,400,433))
mtext("U1463 Depth below seafloor (m WMSF)", side=1, line=3, cex=1.2)

dev.off()

#U1462 output ----
U1462_output=cbind(U1462[c(1:5625),1],U1462_tuned)
colnames(U1462_output)<- c("depth U1462", "depth U1463", "gamma")
write.csv(U1462_output, "U1462_output.csv", row.names = F)

labelsU1462=approx(agemodel[,2],agemodel[,1],xout = c(100,160,270,380))$y
labelsU1462=cbind(c(100,160,270,380),round(labelsU1462))
write.csv(labelsU1462, "labelsU1462.csv", row.names = F)