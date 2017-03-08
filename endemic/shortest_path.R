# create a quartic graph and see how much space it takes up

# use igraph to generate the shortest path tree

# then accumulate it up over the tree and repeat for each node

# 1. Create matrix representing the pixel image
n_size <- 2^(1:4)
mem_use <- matrix(0, length(n_size), 1)
time_use <- matrix(0, length(n_size), 5)

for (i in seq_along(n_size)) {
  n <- n_size[i]

  time_use[i,1] <- system.time({m <- matrix(abs(rnorm(n*n)), n, n)})[3]
  
  # 2. Turn this into a graph. This is massively ineffeictict
  vnum <- function(x, y, n) {
    (x-1)*n + y
  }

  square_graph <- function(m) {
    
    edges <- matrix(0, (nrow(m)-1)*ncol(m) + (ncol(m)-1)*nrow(m), 3)
    enum <- 1
    for (x in 1:nrow(m)) {
      for (y in 1:ncol(m)) {
        if (x > 1) {
          # add link between x-1 and x
          edges[enum,] <- c(vnum(x, y, ncol(m)), vnum(x-1, y, ncol(m)), 1/(m[x,y]*m[x-1,y]))
          enum <- enum + 1
        }
        if (y > 1) {
          # add link between y-1 and y
          edges[enum,] <- c(vnum(x, y, ncol(m)), vnum(x, y-1, ncol(m)), 1/(m[x,y]*m[x,y-1]))
          enum <- enum + 1
        }
      }
    }
    make_graph(t(edges[,1:2]), dir=FALSE) %>% set_edge_attr('weight', value = edges[,3])
  }

  time_use[i,3] <- system.time({
  g <- square_graph(m)
  })[3]
  
  # now let's see how large it is
  mem_use[i,1] <- object.size(g)

  # compute the endemic risk to each vertex
  endemic_risk <- numeric(n*n)

  time_use[i,4] <- system.time({
  for (v in 1:(n*n)) {
    #cat('processing', v, 'of', n*n, '\n')

    # find shortest paths from v to all other vertices. Note this isn't symmetric
    # TODO: This would be scaled by the risk at v
    risk_contrib <- distances(g, v=v, algorithm='dijkstra')

    # add on to our endemic risk
    endemic_risk <- endemic_risk + risk_contrib[1,]
  }})[3]
  endemic_risk1 <- endemic_risk

  my_dist <- function(m, x, y) {
    # JM version of Dijkstra's algorithm for images. Returns distance
    # from (x,y) to all values in the image.
    
    # nodes are coded with (x,y) to make things easier to track
    # (may need to replace this later for efficiency purposes)

    # initialise our vertex distances to infinity, except for our source node
    dist <- matrix(Inf, nrow(m), ncol(m)); dist[x,y] <- 0;
#    prev <- array(NA, dim=c(nrow(m), ncol(m), 2))

    Q <- 1:(nrow(m)*ncol(m))

    # run over and iteratively improve
    while(length(Q) > 0) {
      # find node with smallest dist in the array inside Q.
      cat("length(Q) =", length(Q), "\n")
      row <- which.min(dist[Q]) # row in Q
      cat("smallest row is", row, "with dist", dist[Q[row]], "\n")
      v <- Q[row] # vertex on the grid
      # find the neighbours of v...
      nb <- function(v) {
        x <- (v - 1) %% nrow(m) + 1
        y <- (v - 1) %/% nrow(m) + 1
        nb <- NULL
        if (x > 1)
          nb <- c(nb, vnum(x-1, y, nrow(m)))
        if (y > 1)
          nb <- c(nb, vnum(x, y-1, nrow(m))) #v - nrow(m))
        if (x < nrow(m))
          nb <- c(nb, vnum(x+1, y, nrow(m))) #v + 1)
        if (y < ncol(m))
          nb <- c(nb, vnum(x, y+1, nrow(m))) #v + nrow(m))
        nb
      }
#      cat("this corresponds to x=", x, ", y=", y, "with dist", dist[x,y], "\n")
      md <- dist[v]
      Q <- Q[-row]; # removes the vertex from Q (super inefficient)

      # find and update the neighbours
      n <- nb(v)
      for (i in n) { # TODO: Only really need to check within Q
        d = dist[v] + 1/(m[v]*m[i])
        if (d < dist[i]) {
          dist[i] <- d
        }
      }
    }
    dist
  }
  endemic_risk <- numeric(n*n)
  time_use[i,5] <- system.time({
    d <- 20
    for (x in 1:n) {
      for (y in 1:n) {
        
        # find shortest paths from v to all other vertices. Note this isn't symmetric
        # TODO: This would be scaled by the risk at v
        risk_contrib <- my_dist(m, x, y)
        
        # add on to our endemic risk
        endemic_risk <- endemic_risk + as.numeric(t(risk_contrib))
      }
    }})[3]
#  print(all.equal(endemic_risk, endemic_risk1))
}

# repeat for each vertex