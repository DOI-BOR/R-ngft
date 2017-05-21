# test rgft with a delta function

# require(grDevices)

# make a delta function, centered in a window
#sig <- c(rep(0,13),1,rep(0,14))
#sig <- c(rep(0,15),1,rep(0,16))
#sig <- c(rep(0,23),1,rep(0,24))
sig <- c(rep(0,101),1,rep(0,154))
sig <- c(rep(0,411),1,rep(0,612))
#sig <- c(rep(0,400),1,-1,rep(0,1600))
sig <- c(rep(0,108),rep(1,40),rep(0,108))

# choose time increment
dt <- 0.01
N <- length(sig)
df = 1 / (N * dt)

# transform the delta function
dec=1
fst.fwd <- ngft::fst(sig, dt, eps=(N/dec-1), win.type='g', part.type="f", fw.width=dec, by.part=FALSE)
#fst.fwd <- ngft::fst(sig, dt, eps=(N-1), win.type='g', part.type="d", by.part=FALSE)
#fst.fwd <- ngft::fst(sig, dt, eps=0, win.type='g', part.type="d", by.part=FALSE)
#fst.fwd <- ngft::fst(sig, dt, eps=3, win.type='g', part.type="d", by.part=FALSE)
#fst.fwd <- ngft::fst(sig, dt, eps=(N-1), win.type='g', part.type="e", by.part=FALSE)
#fst.fwd <- ngft::fst(sig, dt, eps=(N/8-1), win.type='g', part.type="e", by.part=FALSE)
#fst.fwd <- ngft::fst(sig, dt, eps=(N-1), win.type='b', part.type="e", by.part=FALSE)

# plot the modulus. image() plots columns on horizontal axis,
# so also take transpose
img.plt <- t(abs(fst.fwd$image))
xvals <- seq(from=0,to=fst.fwd$wd - 1,by=1)
yvals <- seq(from=0,to=fst.fwd$ht - 1,by=1)
nlev <- fst.fwd$wd # 0.25 * fst.delta$ht * fst.delta$wd
breaks <- quantile(img.plt, probs=seq(0,1,1/nlev))
image(x=xvals, y=yvals, z=img.plt,
      col=rainbow(length(breaks)-1), breaks=breaks,
      xlab="time index", ylab="frequency index")
# gray.colors, rainbow, heat.colors, topo.colors, terrain.colors

xvals <- dt * fst.fwd$t.centers
yvals <- df * fst.fwd$f.centers
img.plt <- t(abs(fst.fwd$image))
plot3D::image2D(img.plt, x=xvals, y=yvals, xlab="time, sec.", ylab="frequency, Hz.")



### Test inverse

sig <- c(rep(0,7),1,rep(0,7))
sig <- c(rep(0,7),1,rep(0,8))
sig <- c(rep(0,108),rep(1,40),rep(0,108))
dt <- 0.01
N <- length(sig)

fst.fwd <- ngft::fst(sig, dt, eps=0, part.type="d", win.type="g", by.part=FALSE)
fst.inv <- ngft::ifst(fst.fwd$gft, N, dt, eps=0, part.type="d", win.type="g", inv.type="d")
fst.inv$TS

fst.fwd <- ngft::fst(sig, dt, eps=N-1, fw.width=1, part.type="f", win.type="g", by.part=FALSE)
fst.inv <- ngft::ifst(fst.fwd$gft, N, dt, eps=N-1, fw.width=1, part.type="f", win.type="g", inv.type="f")
fst.inv$TS
