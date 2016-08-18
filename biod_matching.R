k <- seq(0, 2, by=0.001)
m <- 2

biod <- function(k, m) {
  m ^ ((k^2 - 1/2*k) / (k + 1/2)) * (1 - (1 - 1/m^k)^m)
}

risk <- function(z, m) {
  1 / (m ^ (z - 1/2) * (1 - (1 - 1/m^z)^m)^((z+1/2)/z))
}


plot(NULL, xlim=c(0,3), ylim=c(0.8,1.6), xlab="Power", ylab="Biodiversity")
for (i in 1:8) {
  plot(function(x) { biod(x, i)}, xlim=c(0, 4), col=i, add=TRUE, n=1000)
}

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

plot(r_hat*cos(theta_hat), r_hat*sin(theta_hat), pch='.', col=ifelse(r_hat < R* - R*a*sin(k*theta_hat), "red", "blue"))

arc_len <- function(R, a, k) {
  arclen_f <- function(theta, R, a, k) {
    r = R*(1 - a*sin(k*theta))
    dr = -R*a*k*cos(k*theta)
    sqrt(r^2 + dr^2)
  }
  alen <- numeric(length(R))
  for (i in seq_along(R))
    alen[i] = integrate(arclen_f, 0, 2*pi, subdivisions=1000,R=R[i], a=a, k=k)$value
  alen
}

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
