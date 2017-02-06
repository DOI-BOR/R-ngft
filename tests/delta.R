# test rgft with a delta function

# require(grDevices)

# make a delta function, centered in a window
#delta <- c(rep(0,13),1,rep(0,14))
#delta <- c(rep(0,15),1,rep(0,16))
#delta <- c(rep(0,23),1,rep(0,24))
delta <- c(rep(0,400),1,-1,rep(0,1600))

# choose time increment
dt <- 0.01
N <- length(delta)
df = 1 / (N * dt)

# transform the delta function
fst.delta <- ngft::fst(delta, dt)
if ( N < 50 )
  abs(fst.delta$image)

# plot the modulus. image() plots columns on horizontal axis,
# so also take transpose
img.plt <- t(abs(fst.delta$image))

xvals <- seq(from=1,to=fst.delta$wd,by=1)
yvals <- seq(from=0,to=fst.delta$ht,by=1)
nlev <- 2 * fst.delta$ht
breaks <- quantile(img.plt, probs=seq(0,1,1/nlev))
image(x=xvals, y=yvals, z=img.plt,
      col=rainbow(length(breaks)-1), breaks=breaks,
      xlab="time index", ylab="frequency index")

xvals <- dt * fst.delta$t.centers
yvals <- df * fst.delta$f.centers
#image2D(img.plt, x=xvals, y=yvals, xlab="time, sec.", ylab="frequency, Hz.")


