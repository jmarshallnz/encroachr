k <- seq(0, 2, by=0.001)
m <- 2

biod <- function(k, m) {
  m ^ ((k^2 - 1/2*k) / (k + 1/2)) * (1 - (1 - 1/m^k)^m)
}

risk <- function(k, m) {
  1 / (m ^ (k - 1/2) * (1 - (1 - 1/m^k)^m)^((k+1/2)/k))
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
r_hat = (R+a) * sqrt(runif(n))

sum(r_hat < R - a*sin(k*theta_hat))

area_debrolie = pi*R^2 + pi*a^2/2
area_circle  = pi*(R+a)^2
area_debrolie/area_circle

plot(r_hat*cos(theta_hat), r_hat*sin(theta_hat), pch='.', col=ifelse(r_hat < R - a*sin(k*theta_hat), "red", "blue"))

arc_len <- function(R, a, k) {
  arclen_f <- function(theta, R, a, k) {
    r = R - a*sin(k*theta)
    dr = -a*k*cos(k*theta)
    sqrt(r^2 + dr^2)
  }
  alen <- numeric(length(R))
  for (i in seq_along(R))
    alen[i] = integrate(arclen_f, 0, 2*pi, subdivisions=1000,R=R[i], a=a, k=k)$value
  alen
}

biod_petal <- function(z, m) {
  m ^ ((k^2 - 1/2*k) / (z + 1/2)) * (1 - (1 - 1/m^k)^m)
}

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


plot(NULL, xlim=c(0,4), ylim=c(0,3), xlab="Power", ylab="Risk")
for (i in 1:8) {
  plot(function(x) { risk_petal(x, i, 0.2, 3)}, xlim=c(0, 4), col=i, add=TRUE, n=1000)
}
legend('topright', paste(1:8, "fragments"), fill=1:8)

solidity <- function(a) {
  (1 + a^2/2) / (1 + a)^2
}
