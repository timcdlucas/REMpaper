
####################################
### Plot the regions            ####
####################################

library(RColorBrewer)

polys <- list(
              gas=list(c(2*pi, 2*pi, 2.3*pi, 2.3*pi),c(pi, 2.3*pi, 2.3*pi, pi)),
              p131=list(c(pi/2, pi/2, pi, pi),c(2*pi, 2.3*pi, 2.3*pi, 2*pi)),
              p141=list(c(0, 0, pi/2, pi/2),c(2*pi, 2.3*pi, 2.3*pi, 2*pi)),
              p221=list(c(pi,pi,2*pi,2*pi),c(2*pi, 2.3*pi, 2.3*pi, pi)),
              p222=list(c(3*pi/2, pi, 2*pi),c(pi, 2*pi, pi)),
              p223=list(c(pi, pi, 3*pi/2),c(pi, 2*pi, pi)),
              p231=list(c(pi/2, pi, pi),c(2*pi, 2*pi, pi)),
              p232=list(c(pi/2, pi/2, pi),c(3*pi/2, 2*pi, pi)),
              p233=list(c(pi/2, pi/2, pi),c(pi, 3*pi/2, pi)),
              p241=list(c(0, pi/2, pi/2),c(2*pi, 2*pi, 3*pi/2)),
              p242=list(c(0, pi/2, pi/2),c(2*pi, 3*pi/2, pi)),
              p243=list(c(0, 0, pi/2),c(pi, 2*pi, pi)),
              p311=list(c(2*pi, 2*pi, 2.3*pi, 2.3*pi),c(0,pi,pi,0)),
              p321=list(c(3*pi/2, 2*pi, 2*pi),c(pi, pi, 0)),
              p322=list(c(pi, 3*pi/2, 2*pi),c(pi, pi, 0)),
              p323=list(c(pi, 2*pi, pi),c(pi, 0, 0)),
              p331=list(c(pi/2, pi/2, pi),c(pi/2, pi, pi)),
              p332=list(c(pi/2, pi/2, pi),c(0, pi/2, pi)),
              p333=list(c(pi/2, pi, pi),c(0, pi, 0)),
              p341=list(c(pi/3, pi/2, pi/2),c(pi/3, pi/2, 0)),
              p342=list(c(pi/4,pi/2, pi/2, pi/3),c(pi/2, pi, pi/2, pi/3)),
              p343=list(c(0, pi/2, pi/4),c(pi, pi, pi/2)),   
              p344=list(c(0, pi/3, pi/2),c(0, pi/3, 0)),   
              p345=list(c(0, pi/4, pi/3),c(0, pi/2, pi/3)),
              p346=list(c(0, 0, pi/4),c(0, pi, pi/2))
)    

cols <- c( brewer.pal(9, 'Pastel1'), brewer.pal(10, 'Set3'), brewer.pal(9, 'Paired'))


pdf('~/Dropbox/PhD/Analysis/REM-chapter/imgs/regions.pdf')

plot(c(0,0), c(0,0) , type='l', ylim=c(0,2.3*pi), xlim=c(0,2.3*pi),
      xlab=expression(theta[s]), ylab=expression(theta[a]), xaxt='n', yaxt='n')

axis(1, c(0, pi, 2*pi), c(0,expression(pi), expression(paste('2', pi))) )
axis(2, c(0, pi, 2*pi), c(0,expression(pi), expression(paste('2', pi))) )


#cols <- sample(colors()[-c(1, 24, 151:250)], length(polys))



for(i in 1:length(polys)){
     polygon(polys[[i]][[1]], polys[[i]][[2]], col=cols[i])
}

for( i in 1: length(polys)){
	xV <- c(polys[[i]][[1]], polys[[i]][[1]][1])
	yV <- c(polys[[i]][[2]], polys[[i]][[2]][1])

	A <- 0.5*sum(sapply(c(1:(length(xV)-1)), function(j){ xV[j]*yV[j+1] - xV[j+1]*yV[j] }))

	Cx = sum(sapply(c(1:(length(xV)-1)), function(j){ (xV[j] + xV[j+1])*(xV[j]*yV[j+1] - xV[j+1]*yV[j]) }))/(6*A)
	Cy = sum(sapply(c(1:(length(xV)-1)), function(j){ (yV[j] + yV[j+1])*(xV[j]*yV[j+1] - xV[j+1]*yV[j]) }))/(6*A)

	text(Cx, Cy, names(polys)[i], cex=0.8)
}



#for(i in 1:length(polys)){
#    text(mean(polys[[i]][[1]]), mean(polys[[i]][[2]]), names(polys)[i], cex=0.8)
#}

dev.off()



# new regions plot showing similarity.

type = c(1, 2, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8)


#pdf('~/Dropbox/PhD/Analysis/REM-chapter/imgs/equalRegions.pdf')
png('~/Dropbox/PhD/Analysis/REM-chapter/imgs/equalRegions.png')

plot(c(0,0), c(0,0) , type='l', ylim=c(0,2.3*pi), xlim=c(0,2.3*pi),
      xlab=expression(theta[s]), ylab=expression(theta[a]), xaxt='n', yaxt='n')

axis(1, c(0, pi, 2*pi), c(0,expression(pi), expression(paste('2', pi))) )
axis(2, c(0, pi, 2*pi), c(0,expression(pi), expression(paste('2', pi))) )


#cols <- sample(colors()[-c(1, 24, 151:250)], length(polys))

cols <- c( brewer.pal(9, 'Pastel1'), brewer.pal(10, 'Set3'), brewer.pal(9, 'Paired'))

for(i in 1:length(polys)){
     polygon(polys[[i]][[1]], polys[[i]][[2]], col=cols[type[i]])
}

for( i in 1: length(polys)){
	xV <- c(polys[[i]][[1]], polys[[i]][[1]][1])
	yV <- c(polys[[i]][[2]], polys[[i]][[2]][1])

	A <- 0.5*sum(sapply(c(1:(length(xV)-1)), function(j){ xV[j]*yV[j+1] - xV[j+1]*yV[j] }))

	Cx = sum(sapply(c(1:(length(xV)-1)), function(j){ (xV[j] + xV[j+1])*(xV[j]*yV[j+1] - xV[j+1]*yV[j]) }))/(6*A)
	Cy = sum(sapply(c(1:(length(xV)-1)), function(j){ (yV[j] + yV[j+1])*(xV[j]*yV[j+1] - xV[j+1]*yV[j]) }))/(6*A)

	text(Cx, Cy, names(polys)[i], cex=0.8)
}

dev.off()



