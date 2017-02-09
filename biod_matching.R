k <- seq(0, 2, by=0.001)
m <- 2

# This computes the biodiversity ratio of m versus 1 fragments of the same shape
# if area of the m fragments is such that the total risk is the same in both cases,
# where the dilution effect is d
biod <- function(z, m, d=1) {
  # log(m ^ (z*(z*d - 1/2) / (z*d + 1/2)) * (1 - (1 - 1/m^z)^m))
  
  # question: Can we expand (1 - (1-1/m^z)^m) for large z and show some sort of decay?
  # e.g. for m=1 it's just 1
  #      for m=2 it's 1-(1 - 1/2^z)^2 = 1 - (1 - 2/(2^z) + 1/(2^z^2))
  #                                   = 2/2^z - 1/4^z -> 0 as z -> inf with rate similar to e^x ??
  #                                   = 1/2^z * (2 - 1/2^z)
  # then we also have                2^z / 2^(1/(z + 1/2))
  logb = (z*(z*d - 1/2) / (z*d + 1/2)) * log(m) + log(1 - ((m^z - 1)/m^(z*m))^m)
  exp(logb)
  #  m ^ (z*(z*d - 1/2) / (z*d + 1/2)) * (1 - (1 - 1/m^z)^m)
}

#
# 1. Plot of perimeter vs solidity for deBroglie with differing number of petals
# 2. Plot of dilution effect - split the colours up.
# 3. deBroglie cirles plots (i.e. what they look like)
# 4. Read + critique manuscript
# 5. Play with the maths???
# 
# This computes the risk ratio of m versus 1 fragments of the same shape
# if area of the m fragments is such that the total biodiversity is the same in both cases,
# where the dilution effect is d
risk <- function(z, m, d=1) {
  1 / (m ^ (z*d - 1/2) * (1 - (1 - 1/m^z)^m)^((z*d+1/2)/z))
}

plot(NULL, xlim=c(0,100), ylim=c(0.8,1.6), xlab="Power", ylab="Biodiversity")
for (i in 1:8) {
  plot(function(x) { biod(x, i)}, xlim=c(0, 100), col=i, add=TRUE, n=1000)
}

plot(NULL, xlim=c(0,3), ylim=c(0,4), xlab="Power", ylab="Risk")
for (i in 1:8) {
  plot(function(x) { risk(x, i, 0.7)}, xlim=c(0, 4), col=i, add=TRUE, n=1000)
}

# alternate computation for total biodiversity
area <- c(1,1,1)/3
power <- 0.5

total_area = sum(area)
prop_area = area / total_area
prop_biod = prop_area^power
prop_outside_biod = 1 - prop_biod
prop_outside_all  = prod(prop_outside_biod)
prop_union_biod   = 1 - prop_outside_all
total_biod = total_area^power
total_union = total_biod * prop_union_biod
total_union

# In the special case of power=0.5, the above is same as:
biod(0.5, 3)

# the biod function otherwise uses different areas in the m versus 1 case so that risk is balanced.

# now, for risk, we'll need to take area/k to the power z then multiply by perimeter and add?
3 * ((1/3)^power * sqrt(1/3))
(risk(0.5,3)/3)

1 - (1 - (1/m)^z)^m

find_power <- function(m, d) {
  lwr = risk(0.01,m,d)-1
  upr = risk(2,m,d)-1
  if (lwr*upr > 0) {
    return(NA);
  }
  uniroot(function(x, m, d) { risk(x,m,d)-1 }, interval=c(0,2), m=m, d=d)$root
}



##### Some later figure: Dilution effect stuff here

## Figure 6: Fragments vs power for different dilution effects
dilute_val <- c(0.9, 1, 1.2)
dilute_col <- c('purple', 'black', 'red')
dilute_lty <- c('dashed', 'solid', 'dashed')
dilute_lwd <- c(1, 2, 1)

balance <- expand.grid(fragments = 2:10, dilute = dilute_val)
balance$power <- apply(balance, 1, function(x) { find_power(x[1], x[2]) })

pdf("figures/power_vs_fragments_by_dilution.pdf", width=7, height = 7)
plot(NULL, xlim=c(2,10), ylim=c(0,3), xlab="Number of patches", ylab="power (z)")
for (i in 1:length(dilute_val))
  lines(power ~ fragments, data=subset(balance, dilute==dilute_val[i]), col=dilute_col[i], lty=dilute_lty[i], lwd=dilute_lwd[i])
legend('topright', bty='n', legend=dilute_val, lty = dilute_lty, col = dilute_col, lwd = dilute_lwd, title = 'd')
dev.off()



plot(NULL, xlim=c(0,4), ylim=c(0,3), xlab="Power", ylab="Risk")
for (i in 1:8) {
  plot(function(x) { risk(x, i)}, xlim=c(0, 4), col=i, add=TRUE, n=1000)
}
legend('topright', paste(1:8, "fragments"), fill=1:8)

#' Test out some de Broglie plots
#'
#' TODO:
#'   1. Change amplitude for a given k and then adjust such that area is preserved.
#'   2. Compute solidity (area to minimal circle area outside) and core (maximal circle area inside to full area?)
#'   3. Plot for various k and amplitudes, keeping area fixed (r will change)
#'   4. Compute perimeters for each one as well via numeric integration (already have this below!)
#'   
R = 1
a = 0.2
k = 10

theta = seq(0, 2*pi, length.out=1000)
r = R - a*sin(k*theta)

plot(r*cos(theta), r*sin(theta), type="l")

#' Montecarlo estimate of area.
#' Sample randomly R*sqrt(runif())

n = 10000
theta_hat = runif(n)*2*pi
r_hat = R*(1+a) * sqrt(runif(n))

sum(r_hat < R - R*a*sin(k*theta_hat))/n

area_debrolie = pi*R^2*(1 + a^2/2)
area_circle  = pi*R^2*(1+a)^2
area_debrolie/area_circle

# compute perimeter of radius R, amplitude a, petals k
arc_len <- function(R, a, k) {
  arclen_f <- function(theta, R, a, k) {
    r = R*(1 - a*sin(k*theta))
    dr = -R*a*k*cos(k*theta)
    sqrt(r^2 + dr^2)
  }
  alen <- numeric(length(R))
  if (length(a) == 1 && length(R) > 1)
    a = rep(a, length(R))
  for (i in seq_along(R))
    alen[i] = integrate(arclen_f, 0, 2*pi, subdivisions=1000,R=R[i], a=a[i], k=k)$value
  alen
}

# Amplitude should be between 0..1
# Area should be 1 (it's scale-independent, dummy)
debroglie <- function(petals, amplitude) {
  # 1. figure out the radius from the area, petals + amplitude
  #   area = pi*R^2*(1 + amplitude^2/2)
  R = 1/sqrt(1 + amplitude^2/2)
  # 2. figure out solidity + core ratio things
  solidity = (pi*1^2) / (pi*R^2*(1+amplitude)^2)
  core     = (pi*R^2*(1-amplitude)^2) / (pi*1^2)
  # 3. computer perimeter
  perimeter = arc_len(R, amplitude, petals)
  return(list(R=R, s=solidity, c=core, p=perimeter, k=petals))
}

a = seq(0,1, by=0.01)
out <- lapply(1:12, function(x) { as.data.frame(debroglie(x, a)) })
out <- do.call(rbind, out)
out$k <- factor(out$k)

library(lattice)
pdf("solidity_vs_core.pdf")
xyplot(s ~ c | k, data=out, type="l", xlab="core ratio", ylab="solidity ratio", main="Solidity ratio vs Core ratio")
xyplot(p ~ s | k, data=out, type="l", xlab="solidity ratio", ylab="perimeter", main="Perimeter vs Solidity ratio")
xyplot(p ~ c | k, data=out, type="l", xlab="core ratio", ylab="perimeter", main="Perimeter vs Core ratio")
dev.off()

out <- as.data.frame(debroglie(5, a))
tapply()
plot(s ~ c, data=out, col = k, type="l")

out <- as.data.frame(debroglie(10, a))

plot(p ~ s, data=out, type="l")

plot(r_hat*cos(theta_hat), r_hat*sin(theta_hat), pch='.', col=ifelse(r_hat < R* - R*a*sin(k*theta_hat), "red", "blue"))


plot(function(x) { ans = numeric(length(x)); for (i in seq_along(x)) ans[i] = arc_len(1, x[i], 5); ans }, xlim=c(0,1))
plot(function(x) { ans = numeric(length(x)); for (i in seq_along(x)) ans[i] = arc_len(1, x[i], 10); ans }, xlim=c(0,1), add=TRUE)

risk_petal <- function(z, m, a, k) {
  # start by balancing the biodiversity
  # Bm = (m*pi*rm^2*(1+a^2/2))^z * (1 - (1 - 1/m^z)^m)
  # B1 = (pi*r1^2*(1+a^2/2))^z
  # Bm = B1 -> m^z rm^(2z) (1 - (1 - 1/m^z)^m) = r1^(2z) ->
  r1 = 1
  rm = r1/(sqrt(m) * (1 - (1 - 1/m^z)^m)^(1/(2*z)))
  
  # risk...
  Rm = m * (pi*rm^2*(1+a^2/2))^z * arc_len(rm, a, k)
  R1 = (pi*r1^2*(1+a^2/2))^z * arc_len(r1, a, k)
  Rm / R1
}

solidity <- function(a) {
  (1 + a^2/2) / (1 + a)^2
}

pdf("z.pdf", width=12, height=12)
par(mfrow=c(3,3))
for (num_clusters in 1:9) {
  grid = expand.grid(amplitude=seq(0,1,by=0.01), complexity=0:10, z=NA)
  for (i in 1:nrow(grid)) {
    grid$z[i] = uniroot(function(x) { risk_petal(x, num_clusters, grid$amplitude[i], grid$complexity[i]) - 1 }, c(0.5, 10))$root
    if (i %% 100 == 0)
      cat("Up to", i, "of", nrow(grid), "\n")
  }
  image(seq(0,1,by=0.01), 0:10, matrix(grid$z, 101, 11), zlim=c(0.5, 10), xlab="Amplitude", ylab="Complexity (nodes)", main=paste(num_clusters, "clusters"))
}
dev.off()

pdf("risk.pdf", width=8, height=6)
plot(NULL, xlim=c(0,4), ylim=c(0,3), xlab="Power", ylab="Risk")
for (i in 1:8) {
  plot(function(x) { risk_petal(x, i, 0.8, 5)}, xlim=c(0, 4), col=i, add=TRUE, n=1000)
}
legend('topright', paste(1:8, "fragments"), fill=1:8)
dev.off()
