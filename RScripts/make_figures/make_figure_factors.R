# ==============================================================================
#             PLOTS UNOBSERVED FACTORS
# ==============================================================================

#### Y_t -----------------------------------------------------------------------
## Parameters in order to diplay well the chart

color_supply <- "grey45"
color_demand <- "black"
line_width <- 2.8
dotted_line_width <- 5.0

path <- paste0(path_graph, "Y.t.pdf")
pdf(path, width = 8, height = 5,
    pointsize = 11)

par(mfrow=c(2,2))
par(plt=c(.1,.95,.15,.8))

for(i in 1:Model.final$m){
  indic_demand <- sign(Model.final$delta.c[i,1]) == sign(Model.final$delta.c[i,2])
  if(indic_demand){
    color <- color_demand
    line_type <- 1
    lwd <- line_width
  }else{
    color <- color_supply
    line_type <- 3
    lwd <- dotted_line_width
  }
  
  name_factor <- paste("Y[",i,"*','*t]",sep="")
  main.t <- expression(paste('Factor ',Y[1*','*t],sep=''))
  
  eval(parse(text = paste("main.t <- expression(paste('Factor ',Y[",i,"*','*t],sep=''))",sep="")))
  
  plot(as.Date(vec.dates$date),KF.result4$xi.tT[,i], type ="n",
       xlab="",ylab="",las=1,
       main=main.t)
  grid()
  make_recessions()
  lines(as.Date(vec.dates$date), KF.result4$xi.tT[,i],
        col=color, lwd=lwd, lty=line_type)
  
  if(i==1){
    legend("bottomright",
           c("Demand factor", "Supply factor"),bg="white",
           col=c(color_demand, color_supply), lty=c(1, 3),
           lwd=c(line_width, dotted_line_width),
           cex=0.9, horiz = TRUE)
  }
}

dev.off()



#### z_t -----------------------------------------------------------------------

path <- paste0(path_graph,"z.t.pdf")
pdf(path, width = 8, height = 5,
    pointsize = 11)

par(mfrow=c(2,2))
par(plt=c(.1,.95,.15,.8))

main.t <- expression(paste(z[p*','*t]^s," and ",z[n*','*t]^s,sep=''))

z_ps <- KF.result4$xi.tT[,Model.final$n+1]
z_ns <- KF.result4$xi.tT[,Model.final$n+2]
ylim <- c(0,max(z_ps,z_ns))

plot(as.Date(vec.dates$date), z_ps, type ="n",
     xlab="",ylab="",las=1,
     ylim=ylim,
     main=main.t)
grid()
make_recessions()
lines(as.Date(vec.dates$date), z_ps,
      col=color_supply, lwd=line_width)
lines(as.Date(vec.dates$date), z_ns,
      col=color_supply, lwd=dotted_line_width, lty=3)

# legend("topright",
#        c(expression(paste(z[p*','*t]^s,sep='')),
#          expression(paste(z[n*','*t]^s,sep=''))),
#        bg="white",
#        lty=c(1,3),
#        lwd=2,
#        col=c(color_supply, color_supply),
#        cex=0.9)


main.t <- expression(paste(z[p*','*t]^d," and ",z[n*','*t]^d,sep=''))

z_pd <- KF.result4$xi.tT[,Model.final$n+3]
z_nd <- KF.result4$xi.tT[,Model.final$n+4]
ylim <- c(0,max(z_ps,z_ns))

plot(as.Date(vec.dates$date), z_pd, type ="n",
     xlab="",ylab="",las=1,
     ylim=ylim,
     main=main.t)
grid()
make_recessions()
lines(as.Date(vec.dates$date), z_pd,
      col=color_demand, lwd=line_width)
lines(as.Date(vec.dates$date), z_nd,
      col=color_demand, lwd=dotted_line_width, lty=3)

main.t <- expression(paste(z[v*','*t],sep=''))

z_v <- KF.result4$xi.tT[,Model.final$n+5]

plot(as.Date(vec.dates$date), z_v, type ="n",
     xlab="",ylab="",las=1,
     main=main.t)
grid()
make_recessions()
lines(as.Date(vec.dates$date), z_v,
      col="black", lwd=line_width)

plot.new()
legend("topleft",
       c(expression(paste(z[p*','*t]^s,", positive component of the supply factor ",Y[1*','*t],sep='')),
         expression(paste(z[n*','*t]^s,", negative component of the supply factor ",Y[1*','*t],sep='')),
         expression(paste(z[p*','*t]^d,", positive component of the demand factor ",Y[2*','*t],sep='')),
         expression(paste(z[n*','*t]^d,", negative component of the demand factor ",Y[2*','*t],sep='')),
         expression(paste(z[v*','*t],", volatility factor.",sep=''))),
       bg="white",
       lty=c(1,3,1,3,1),
       lwd=c(line_width, dotted_line_width, line_width,
             dotted_line_width, line_width),
       col=c(color_supply, color_supply,
             color_demand, color_demand,
             "black"),
       cex=0.9)



dev.off()


# 
# 
# if(i==1){
#   legend("bottomright",
#          c("Demand factor", "Supply factor"),bg="white",
#          fill=c("#ff8243", "#56B4E9"), cex=0.9, horiz = TRUE)
# }
# 
# dev.off()
# 
# 
# stop()
# 
# par(mar=c(2, 3, 1, 1))
# b <- Model.final$m
# m.1 <- matrix(c(0,2,
#                 1,2,
#                 1,4,
#                 3,4,
#                 3,5,
#                 0,5,
#                 6,6), nrow = 7, ncol=2, byrow = TRUE)
# #m.1 <- matrix(seq(1,ceiling(b/2)+1),nrow = ceiling((b/2+1)),ncol = 2,byrow = TRUE)
# #layout(mat = m.1, heights = c(rep(0.875/b,b),0.125))
# layout(mat = m.1)
# 
# color <- c("#56B4E9", "#ff8243", "#56B4E9", "#ff8243","#ff8243")
# 
# for(i in 1:Model.final$m){
#   
#   plot(as.Date(vec.dates$date),KF.result4$xi.tT[,i], type ="l", col=color[i], lwd=2)
#   make_recessions()
#   
#   if(i==1){
#     legend("bottomleft",legend=bquote(1^{st} ~ "factor"), cex=1, bty = "n")
#   }
#   
#   if(i==2){
#     legend("bottomleft",legend=bquote(2^{nd} ~ "factor"), cex=1, bty = "n")
#   }
#   
#   if(i==3){
#     legend("bottomleft",legend=bquote(3^{rd} ~ "factor"), cex=1, bty = "n")
#   }
#   
#   if(i>3){
#     legend("bottomleft",legend=bquote(.(i)^{th} ~ "factor"), cex=1, bty = "n")
#   }
#   
# }
# 
# par(mar=c(0.1, 4, 1, 1))
# plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
# legend("bottom",inset = 0, title="Variables:",
#        c("Demand factor", "Supply factor"),
#        fill=c("#ff8243", "#56B4E9"), cex=0.9, horiz = TRUE)
# dev.off()
# 
# #### z_t
# ## Parameters in order to diplay well the chart
# path <- paste0("graphs/US_2024/z.t.pdf")
# 
# pdf(path, width = 5, height = 4)
# par(mar=c(2, 3, 1, 1))
# b <- Model.final$m
# m.1 <- matrix(c(1,0,
#                 1,3,
#                 2,3,
#                 2,0,
#                 4,4), nrow = 5, ncol=2, byrow = TRUE)
# layout(mat = m.1, heights = c(0.22,0.22,0.22,0.22,0.12))
# 
# color <- c("#56B4E9", "#104E8B", "#ff8243", "#CD6600", "black")
# 
# for(i in c(1,3,5)){
#   
#   if(i<4){
#     plot(vec.dates$date, KF.result4$xi.tT[,Model.final$n+i], type ="l", col=color[i], lwd=2,
#          ylim=c(min(KF.result4$xi.tT[,Model.final$n+i],KF.result4$xi.tT[,Model.final$n+i+1]),
#                 max(KF.result4$xi.tT[,Model.final$n+i],KF.result4$xi.tT[,Model.final$n+i+1])))
#     lines(vec.dates$date, KF.result4$xi.tT[,Model.final$n+i+1], col=color[i+1], lwd=2)
#   } else{
#     plot(vec.dates,KF.result4$xi.tT[,Model.final$n+i], type ="l", col="black", lwd=2)
#   }
#   
#   if(i==1){
#     legend("topleft",legend=c(expression( z[p*",t"]^s),
#                               expression( z[n*",t"]^s)),
#            col = c("#56B4E9", "#104E8B"), lty=c(1,1), lwd=2,
#            cex=1, bty = "n", y.intersp=1.2)
#   }
#   
#   if(i==3){
#     legend("topleft",legend=c(expression( z[p*",t"]^d),
#                               expression( z[n*",t"]^d)),
#            col=c("#ff8243", "#CD6600"), lty=c(1,1), lwd=2,
#            cex=1, bty = "n")
#   }
#   
#   if(i>3){
#     legend("topleft",legend=c(expression( z[v*",t"])),
#            col=c("black"), lty=c(1), cex=1, bty = "n", lwd=2,)
#   }
#   
# }
# 
# par(mar=c(0.1, 4, 1, 1))
# plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
# legend("bottom",inset = 0, title="Variables:",
#        c("Demand factor", "Supply factor", "Volatility factor"),
#        fill=c("#ff8243", "#56B4E9", "black"), cex=0.9, horiz = TRUE)
# dev.off()
# 
# 
