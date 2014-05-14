# R convenience functions
# 
# Author: mg14
###############################################################################

#' Function to indicate significance
sig2star = function(s, breaks=c(-Inf, 0.001,0.01,0.05,0.1,1), labels=c("***","**","*",".","")) {
	r <- s
	r[] <- as.character(cut(s, breaks=breaks, labels)) 
	r[is.na(r)] <- ""
	r
}

#' Plot a donut
donut <- function(x, r0=0.5, r1=.75,...) {
	pie(x, border=NA, labels = paste(names(x), " (",round(100*x),"%)", sep=""), radius=r1, ...)
	polygon(cos(seq(0,2*pi,l=100))*r0, sin(seq(0,2*pi,l=100))*r0, col="white", border=NA)
}

require(RColorBrewer)
set1 <- brewer.pal(9,"Set1")

plotcvnet = function(cvnet,X, main="", simple.annot=NULL, col=set1, col0="grey", cex=sqrt(abs(colMeans(X, na.rm=TRUE))[n]+0.02)*5,xlim=c(0.5,ncol(X)),...){
	first=function(x) which(x)[1]
	lollyplot <- function(cvnet, beta, ..., X, main, cex) {
		r <- rank(cvnet$nzero[apply(beta!=0,1,first)], ties.method="first")
		R <- 1 - cvnet$cvm/cvnet$cvm[1]
		Rup <- 1 - cvnet$cvlo/cvnet$cvm[1]
		Rlo <- 1 - cvnet$cvup/cvnet$cvm[1]
		w <- apply(beta!=0,1,first)
		beta <- beta[,which.min(cvnet$cvm)]
		#n = names(beta)
		#ix <- which(c(diff(cvnet$nzero)>0, TRUE) & cvnet$nzero > 0)
		#xx <- c(0,cvnet$nzero[ix]) + 0.5
		
		ix <- na.omit(w)
		xx <- rank(	ix, ties.method="first") + 0.5
		o <- order(xx)
		ix <- ix[o]
		ix <- c(ix[1],ix)
		xx <- c(0.5,xx[o])
		#ix = c(ix[1],ix)
		plot(xx, R[ix], type="S", xaxt="n", xlab="", ylim=c(0,max(Rup)), ylab=expression(paste("Explained variance ",R^2)), main = main, xlim=xlim, col=set1[1],...)
		polygon(rep(c(xx,rev(xx)), each=2)[-2*length(ix) + c(0,-1)], rep(c(Rup[ix],rev(Rlo[ix])), each=2)[-c(1,4*length(ix))], border=NA, col = paste(set1[1],"22", sep=""))		
		lines(xx, cvnet$lambda[ix]/max(cvnet$lambda)*max(R), lty=1, type="S", col="black")
		u = par("usr")[3:4]
		#mtext(colnames(oncocyto), at=r, side=1, las=2, font=ifelse(grepl("[A-Z]", colnames(oncocyto)),3,1), cex=.7, line=0.5)
		n <- colnames(X)#names(ix)[-1]
		colLab <- ifelse(n %in% names(simple.annot), col[simple.annot[n]], "black")
		colLab[w > which.max(R) | is.na(w)] <- "grey"
		rotatedLabel(r, rep(0,length(r)), colnames(X), font=ifelse(grepl("^[A-Z]", n) &! grepl("PC", n),3,1), cex=0.66, col=colLab)
		scale <- diff(range(beta, na.rm=TRUE))*1.05
		if(scale ==0) scale=1
		shift <- mean(range(beta, na.rm=TRUE))#median(c, na.rm=TRUE)
		#lines(r,(beta - shift)/scale*diff(u)*.7 + mean(u), col="grey", type="h", lty=1)
		abline(h=(0 - shift)/scale*diff(u)*.7 + mean(u), col="grey", lty=2)
		points(r,(beta - shift)/scale*diff(u)*.7 + mean(u), col=colLab, cex=cex, pch=ifelse(beta==0,1,19))
		axis(4, labels=pretty(beta),at = (pretty(beta) - shift)/scale*diff(u)*.7 + mean(u), col="darkgrey", col.axis="darkgrey")
		mtext(side=4, line=3, "Coefficient", col="darkgrey", at=mean(u), las=0)
	}
	if(class(cvnet$glmnet.fit)[1] != "multnet"){
		#r <- rank(cvnet$nzero[apply(cvnet$glmnet.fit$beta!=0,1,first)], ties.method="first")
		#w = which.min(cvnet$cvm)
		#beta = cvnet$glmnet.fit$beta[,w]
		#c = c[c!=0]
		lollyplot(cvnet = cvnet, beta=cvnet$glmnet.fit$beta, ... = ..., X = X,   main=main, cex=cex)
	}else {j = 1; for(beta in cvnet$glmnet.fit$beta){
			#w = which.min(cvnet$cvm)
			#beta = beta[,w]
			#c = c[c!=0]
			#r <- rank(cvnet$nzero[apply(beta!=0,1,first)], ties.method="first")
			lollyplot(cvnet = cvnet, beta=beta, ... = ..., X = X,  main=paste(main, ": ", names(cvnet$glmnet.fit$beta)[j], sep=""), cex=cex)
			j = j+1
		}
	}
}

rotatedLabel <- function(x0, y0, labels, pos = 1, cex=1, srt=45, ...) {
	w <- strwidth(labels, units="user", cex=cex)
	h <- strheight(labels, units="user",cex=cex)
	u <- par('usr')
	p <- par('plt')
	f <- par("fin")
	xpd <- par("xpd")
	par(xpd=NA)
	text(x=x0 + ifelse(pos==1, -1,1) * w/2*cos(srt/360*2*base::pi), y = y0 + ifelse(pos==1, -1,1) * w/2 *sin(srt/360*2*base::pi) * (u[4]-u[3])/(u[2]-u[1]) / (p[4]-p[3]) * (p[2]-p[1])* f[1]/f[2] , labels, las=2, cex=cex, pos=pos, adj=1, srt=srt,...)
	par(xpd=xpd)
}


myscale <- function(data){
	apply(data,
			2,
			function(x){
				if(all(x %in% c(0,1,NA))) # Don't scale categorical variables, only re-center
					#x
					scale(x, scale=FALSE)
				else
					scale(x)
				#rank(x) / sum(!is.na(x)) - 0.5
			}
	)
}

violinJitter <- function(x, magnitude=1){
	d <- density(x)
	data.frame(x=x, y=runif(length(x),-magnitude/2, magnitude/2) * approxfun(d$x, d$y)(x))
}


stars <- function (x, full = TRUE, scale = TRUE, radius = TRUE, labels = dimnames(x)[[1L]], 
		locations = NULL, nrow = NULL, ncol = NULL, len = 1, key.loc = NULL, 
		key.labels = dimnames(x)[[2L]], key.xpd = TRUE, xlim = NULL, 
		ylim = NULL, flip.labels = NULL, draw.segments = FALSE, col.segments = 1L:n.seg, 
		col.stars = NA, col.lines = NA, axes = FALSE, frame.plot = axes, 
		main = NULL, sub = NULL, xlab = "", ylab = "", cex = 0.8, 
		lwd = 0.25, lty = par("lty"), xpd = FALSE, mar = pmin(par("mar"), 
				1.1 + c(2 * axes + (xlab != ""), 2 * axes + (ylab != 
									""), 1, 0)), add = FALSE, plot = TRUE, density=NULL, init.angle=0, ...) 
{
	if (is.data.frame(x)) 
		x <- data.matrix(x)
	else if (!is.matrix(x)) 
		stop("'x' must be a matrix or a data frame")
	if (!is.numeric(x)) 
		stop("data in 'x' must be numeric")
	n.loc <- nrow(x)
	n.seg <- ncol(x)
	if (is.null(locations)) {
		if (is.null(nrow)) 
			nrow <- ceiling(if (!is.numeric(ncol)) sqrt(n.loc) else n.loc/ncol)
		if (is.null(ncol)) 
			ncol <- ceiling(n.loc/nrow)
		if (nrow * ncol < n.loc) 
			stop("'nrow * ncol' is less than the number of observations")
		ff <- if (!is.null(labels)) 
					2.3
				else 2.1
		locations <- expand.grid(ff * 1L:ncol, ff * nrow:1)[1L:n.loc, 
		]
		if (!is.null(labels) && (missing(flip.labels) || !is.logical(flip.labels))) 
			flip.labels <- ncol * mean(nchar(labels, type = "c")) > 
					30
	}
	else {
		if (is.numeric(locations) && length(locations) == 2) {
			locations <- cbind(rep.int(locations[1L], n.loc), 
					rep.int(locations[2L], n.loc))
			if (!missing(labels) && n.loc > 1) 
				warning("labels do not make sense for a single location")
			else labels <- NULL
		}
		else {
			if (is.data.frame(locations)) 
				locations <- data.matrix(locations)
			if (!is.matrix(locations) || ncol(locations) != 2) 
				stop("'locations' must be a 2-column matrix.")
			if (n.loc != nrow(locations)) 
				stop("number of rows of 'locations' and 'x' must be equal.")
		}
		if (missing(flip.labels) || !is.logical(flip.labels)) 
			flip.labels <- FALSE
	}
	xloc <- locations[, 1]
	yloc <- locations[, 2]
	angles <- if (full) 
				seq.int(0, 2 * pi, length.out = n.seg + 1)[-(n.seg + 
									1)] + init.angle
			else if (draw.segments) 
				seq.int(0, pi, length.out = n.seg + 1)[-(n.seg + 1)] + init.angle
			else seq.int(0, pi, length.out = n.seg)
	if (length(angles) != n.seg) 
		stop("length of 'angles' must equal 'ncol(x)'")
	if (scale) {
		x <- apply(x, 2L, function(x) (x - min(x, na.rm = TRUE))/diff(range(x, 
									na.rm = TRUE)))
	}
	x[is.na(x)] <- 0
	mx <- max(x <- x * len)
	if (is.null(xlim)) 
		xlim <- range(xloc) + c(-mx, mx)
	if (is.null(ylim)) 
		ylim <- range(yloc) + c(-mx, mx)
	deg <- pi/180
	op <- par(mar = mar, xpd = xpd)
	on.exit(par(op))
	dev.hold()
	on.exit(dev.flush(), add = TRUE)
	if (plot && !add) 
		plot(0, type = "n", ..., xlim = xlim, ylim = ylim, main = main, 
				sub = sub, xlab = xlab, ylab = ylab, asp = 1, axes = axes)
	if (!plot) 
		return(locations)
	s.x <- xloc + x * rep.int(cos(angles), rep.int(n.loc, n.seg))
	s.y <- yloc + x * rep.int(sin(angles), rep.int(n.loc, n.seg))
	if (draw.segments) {
		aangl <- c(angles, if (full) 2 * pi + init.angle else pi+init.angle)
		for (i in 1L:n.loc) {
			px <- py <- numeric()
			for (j in 1L:n.seg) {
				k <- seq.int(from = aangl[j], to = aangl[j + 
										1], by = 1 * deg)
				px <- c(px, xloc[i], s.x[i, j], x[i, j] * cos(k) + 
								xloc[i], NA)
				py <- c(py, yloc[i], s.y[i, j], x[i, j] * sin(k) + 
								yloc[i], NA)
			}
			polygon(px, py, col = col.segments, lwd = lwd, lty = lty, border=col.lines, density=density)
		}
	}
	else {
		for (i in 1L:n.loc) {
			polygon(s.x[i, ], s.y[i, ], lwd = lwd, lty = lty, 
					col = col.stars[i])
			polygon(s.x[i, ], s.y[i, ], lwd = lwd, lty = lty, 
					border = col.lines[i], col = col.stars[i])
			if (radius) 
				segments(rep.int(xloc[i], n.seg), rep.int(yloc[i], 
								n.seg), s.x[i, ], s.y[i, ], lwd = lwd, lty = lty)
		}
	}
	if (!is.null(labels)) {
		y.off <- mx * (if (full) 
						1
					else 0.1)
		if (flip.labels) 
			y.off <- y.off + cex * par("cxy")[2L] * ((1L:n.loc)%%2 - 
						if (full) 
							0.4
						else 0)
		text(xloc, yloc - y.off, labels, cex = cex, adj = c(0.5, 
						1))
	}
	if (!is.null(key.loc)) {
		par(xpd = key.xpd)
		key.x <- len * cos(angles) + key.loc[1L]
		key.y <- len * sin(angles) + key.loc[2L]
		if (draw.segments) {
			px <- py <- numeric()
			for (j in 1L:n.seg) {
				k <- seq.int(from = aangl[j], to = aangl[j + 
										1], by = 1 * deg)
				px <- c(px, key.loc[1L], key.x[j], len * cos(k) + 
								key.loc[1L], NA)
				py <- c(py, key.loc[2L], key.y[j], len * sin(k) + 
								key.loc[2L], NA)
			}
			polygon(px, py, col = col.segments, lwd = lwd, lty = lty, border=col.lines, density=density)
		}
		else {
			polygon(key.x, key.y, lwd = lwd, lty = lty)
			if (radius) 
				segments(rep.int(key.loc[1L], n.seg), rep.int(key.loc[2L], 
								n.seg), key.x, key.y, lwd = lwd, lty = lty)
		}
		lab.angl <- angles + if (draw.segments) 
					(angles[2L] - angles[1L])/2
				else 0
		label.x <- 1.1 * len * cos(lab.angl) + key.loc[1L]
		label.y <- 1.1 * len * sin(lab.angl) + key.loc[2L]
		for (k in 1L:n.seg) {
			text.adj <- c(if (lab.angl[k] < 90 * deg || lab.angl[k] > 
									270 * deg) 0 else if (lab.angl[k] > 90 * deg && 
									lab.angl[k] < 270 * deg) 1 else 0.5, if (lab.angl[k] <= 
									90 * deg) (1 - lab.angl[k]/(90 * deg))/2 else if (lab.angl[k] <= 
									270 * deg) (lab.angl[k] - 90 * deg)/(180 * deg) else 1 - 
										(lab.angl[k] - 270 * deg)/(180 * deg))
			text(label.x[k], label.y[k], labels = key.labels[k], 
					cex = cex, adj = text.adj)
		}
	}
	if (frame.plot) 
		box(...)
	invisible(locations)
}


GetzPlot <- function(x, g, width=.8, ...){
	s <- split(x, g)
	plot(unlist(sapply(s, function(x) (rank(x)-1)/(length(x)-1) - .5)) * width + sort(as.numeric(g)), unlist(s), xaxt="n", ...)
	m <- sapply(s, mean)
	segments(1:nlevels(g)-width/2, m, 1:nlevels(g)+width/2, m)
	abline(v=1:nlevels(g), lty=3)
	axis(side=1, at=1:nlevels(g), labels = levels(g))
}

MoscowPlot <- function(x, sd, col=1:length(x), at = seq_along(x), width = 0.9, ...){
	y <- seq(0,max(x + 3*sd), l=100)
	plot(NA,NA, xlim=c(0,length(x)+1), ylim=range(y), xaxt="n", ...)
	for(i in seq_along(x)){
		xx <- width/2 * pnorm(y, x[i], sd[i], lower.tail = FALSE)
		polygon(i +c(-xx, rev(xx)), c(y, rev(y)) , col=col[i], border=NA)
	}
}

lasso <- function(..., theta=1,  scale=TRUE) {
	x <- cbind(...)
	nvar <- ncol(x)
	xname <- as.character(parse(text=substitute(cbind(...))))[-1]
	vars <- apply(x, 2, function(z) var(z[!is.na(z)]))
	class(x) <- 'coxph.penalty'
	
	if (missing(theta))
		stop("Missing theta")
	
	if (scale) 
		pfun <- function(coef,theta, ndead, scale) {
			list(penalty= sum(abs(coef) *scale)*theta,
					first  = theta*sign(coef)*scale,
					second = 0,
					flag=FALSE)
		}
	else
		pfun <- function(coef,theta, ndead, scale) {
			list(penalty= sum(abs(coef))*theta,
					first  = theta*sign(coef),
					second = 0,
					flag=FALSE)
		}
	
	
	temp <- list(pfun=pfun,
			diag=TRUE,
			cfun=function(parms, iter, history) {
				list(theta=parms$theta, done=TRUE) }, 
			cparm=list(theta= theta),
			pparm= vars,
			varname=paste('lasso(', xname, ')', sep=''))
	
	attributes(x) <- c(attributes(x), temp)
	x
}

colTrans <- function(col, f=2){
	hsv <- apply(col2rgb(col), 2, rgb2hsv)
	tmp <- car::logit(hsv[2,]*.9) - f
	hsv[2,] <- 1/(exp(-tmp) +1)
	tmp <- car::logit(hsv[3,]*.9) + f
	hsv[3,] <- 1/(exp(-tmp) +1)
	apply(hsv,2, function(x) hsv(x[1],x[2],x[3]))
}

jitterBox <- function(x,y, xlim = c(1,nlevels(x)) + c(-.5,.5),...){
	#boxplot(y ~ x, ..., notch=TRUE)
	plot(jitter(as.numeric(x)),y, xlim = xlim, ..., xaxt="n", pch=16)
	s <- split(y,x)
	m <- sapply(s, mean)
	n <- 1:length(m)
	l <- sapply(s, length)
	c <- qnorm(0.975,0,sapply(s, sd)/sqrt(l))
	rect(n-.2,m-c,n+.2,m+c,border=NA, col="#88888844")
	segments(n-.2,m,n+.2,m, lwd=2)
	axis(side = 1, at=n, labels=levels(x))
	if(length(unique(x))>1){
		a <- anova(lm(y~x))
		mtext(side=3, line=0, text= paste("P =", format(a$`Pr(>F)`[1],digits=2)))
	}
}


survdiff <- function (formula, data, subset, na.action, rho = 0, continuity = FALSE) 
{
	call <- match.call()
	m <- match.call(expand.dots = FALSE)
	m$rho <- NULL
	m$continuity <- NULL
	if (!inherits(formula, "formula")) 
		stop("The 'formula' argument is not a formula")
	Terms <- if (missing(data)) 
				terms(formula, "strata")
			else terms(formula, "strata", data = data)
	m$formula <- Terms
	m[[1]] <- as.name("model.frame")
	if (is.R()) 
		m <- eval(m, parent.frame())
	else m <- eval(m, sys.parent())
	y <- model.extract(m, "response")
	if (!inherits(y, "Surv")) 
		stop("Response must be a survival object")
	if (attr(y, "type") != "right") 
		stop("Right censored data only")
	ny <- ncol(y)
	n <- nrow(y)
	offset <- attr(Terms, "offset")
	if (!is.null(offset)) {
		offset <- as.numeric(m[[offset]])
		if (length(attr(Terms, "factors")) > 0) 
			stop("Cannot have both an offset and groups")
		if (any(offset < 0 | offset > 1)) 
			stop("The offset must be a survival probability")
		expected <- sum(-log(offset))
		observed <- sum(y[, ny])
		if (rho != 0) {
			num <- sum(1/rho - ((1/rho + y[, ny]) * offset^rho))
			var <- sum(1 - offset^(2 * rho))/(2 * rho)
		}
		else {
			var <- sum(-log(offset))
			num <- var - observed + if(continuity) 1/2 else 0
		}
		chi <- num * num/var
		rval <- list(n = n, obs = observed, exp = expected, var = var, 
				chisq = chi)
	}
	else {
		strats <- attr(Terms, "specials")$strata
		if (length(strats)) {
			temp <- untangle.specials(Terms, "strata", 1)
			dropx <- temp$terms
			if (length(temp$vars) == 1) 
				strata.keep <- m[[temp$vars]]
			else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
		}
		else strata.keep <- rep(1, nrow(m))
		if (length(strats)) 
			ll <- attr(Terms[-dropx], "term.labels")
		else ll <- attr(Terms, "term.labels")
		if (length(ll) == 0) 
			stop("No groups to test")
		else groups <- strata(m[ll])
		fit <- survdiff.fit(y, groups, strata.keep, rho)
		if (is.matrix(fit$observed)) {
			otmp <- apply(fit$observed, 1, sum)
			etmp <- apply(fit$expected, 1, sum)
		}
		else {
			otmp <- fit$observed
			etmp <- fit$expected
		}
		df <- (etmp > 0)
		if (sum(df) < 2) 
			chi <- 0
		else {
			temp2 <- ((otmp - etmp - if(continuity) 1/2 else 0)[df])[-1]
			vv <- (fit$var[df, df])[-1, -1, drop = FALSE]
			chi <- sum(solve(vv, temp2) * temp2)
		}
		rval <- list(n = table(groups), obs = fit$observed, exp = fit$expected, 
				var = fit$var, chisq = chi)
		if (length(strats)) 
			rval$strata <- table(strata.keep)
	}
	na.action <- attr(m, "na.action")
	if (length(na.action)) 
		rval$na.action <- na.action
	rval$call <- call
	if (is.R()) 
		class(rval) <- "survdiff"
	else oldClass(rval) <- "survdiff"
	rval
}