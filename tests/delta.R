# test rgft with a delta function

# require(grDevices)

# make a delta function, centered in a window
#delta <- c(rep(0,13),1,rep(0,14))
#delta <- c(rep(0,15),1,rep(0,16))
#delta <- c(rep(0,23),1,rep(0,24))
delta <- c(rep(0,101),1,rep(0,154))
delta <- c(rep(0,411),1,rep(0,612))
#delta <- c(rep(0,400),1,-1,rep(0,1600))
delta <- c(rep(0,108),rep(1,40),rep(0,108))

# choose time increment
dt <- 0.01
N <- length(delta)
df = 1 / (N * dt)

# transform the delta function
dec=1
fst.delta <- ngft::fst(delta, dt, eps=(N/dec-1), win.type='g', part.type="f", fw.width=dec, by.part=TRUE)
#fst.delta <- ngft::fst(delta, dt, eps=(N-1), win.type='g', part.type="d", fw.width=1, by.part=TRUE)
#fst.delta <- ngft::fst(delta, dt, eps=7, win.type='g', part.type="d", fw.width=1, by.part=TRUE)
#fst.delta <- ngft::fst(delta, dt, eps=(N-1), win.type='g', part.type="e", fw.width=1, by.part=TRUE)

# plot the modulus. image() plots columns on horizontal axis,
# so also take transpose
img.plt <- t(abs(fst.delta$image))
xvals <- seq(from=0,to=fst.delta$wd - 1,by=1)
yvals <- seq(from=0,to=fst.delta$ht - 1,by=1)
nlev <- fst.delta$wd # 0.25 * fst.delta$ht * fst.delta$wd
breaks <- quantile(img.plt, probs=seq(0,1,1/nlev))
image(x=xvals, y=yvals, z=img.plt,
      col=rainbow(length(breaks)-1), breaks=breaks,
      xlab="time index", ylab="frequency index")
# gray.colors, rainbow, heat.colors, topo.colors, terrain.colors

xvals <- dt * fst.delta$t.centers
yvals <- df * fst.delta$f.centers
img.plt <- t(abs(fst.delta$image))
plot3D::image2D(img.plt, x=xvals, y=yvals, xlab="time, sec.", ylab="frequency, Hz.")

