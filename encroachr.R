# load in some data

library(R.matlab)

# structure of files...

#<circ_num>_<circ_size>_<method>_<rep>.mat


# 
# for each rep (of which there are 50 or so) are the same for the method/circ_num/circ_size.
#
# method directly comparable.
#
# hopefully the methods are easily compared and we'll be able to find sign. diff.
# 

###_###_<method>_##_epi.mat

nums = seq(1,191,by=10)
sizes = seq(1,191,by=10)
reps = 1:50
methods = c("Box", "Linear", "Ribbon", "Vertical")

for (num in nums) {
  pdf(paste0("output_", num, ".pdf"), width=15, height=8)
  par(mfrow=c(4,5), mar=c(1,4,2,2))
  for (size in sizes) {
    for (method in methods) {
      cat("Processing location", num, ",", size, "for method", method, "\n")
      p = list()
      s = list()
      a = list()
      b = list()
      r = list()
      for (rep in reps) {
        # biodiversity data is in the _epi file
        filename = paste(num, size, method, rep, sep="_")

        # main data is here...
        d2 = readMat(file.path("data", paste0(filename, ".mat")))$simulation.data
        d2_epi = readMat(file.path("data", paste0(filename, "_epi.mat")))$biodiversity.data

        p[[rep]] = unlist(lapply(d2["edge.patches",,], sum))
        s[[rep]] = unlist(lapply(d2["shape",,], mean))
        a[[rep]] = unlist(lapply(d2["size.patches",,], sum))

        total_biodiversity <- function(area, power) {
          # basic idea is we're dividing a total area (with total biodiversity)
          # up into a bunch of areas. The biodiversity in each will be less than the total
          # but greater than it would be if it was proportional to area. Thus, the proportion
          # of biodiversity in each bit uses the power law. But there'll be overlap between
          # areas, so it's kinda like sampling without replacement from original popn into
          # sub popns where the probability of inclusion is the prop of biodiversity.
          # Total biodiversity is then the union of the biodiversities across the areas,
          # as you'll get species in common.
          total_area = sum(area)
          if (total_area == 0) {
            return(0)
          }
          prop_area = area / total_area
          prop_biod = prop_area^power
          prop_outside_biod = 1 - prop_biod
          prop_outside_all  = prod(prop_outside_biod)
          prop_union_biod   = 1 - prop_outside_all
          total_biod = total_area^power
          total_union = total_biod * prop_union_biod
          total_union
        }

        b[[rep]] = unlist(lapply(d2["size.patches",,], total_biodiversity, power=0.3))

        # not sure if this is right?
        r[[rep]] = unlist(d2_epi["total.risk",,])
      }

      # write a function!
      plot_thing <- function(p, name) {
        p_max = max(unlist(lapply(p, max)))
        main = paste(method, num, size, sep="_")
        plot(NULL, ylim=c(0,p_max), xlim = c(1,length(p[[1]])), xlab="", xaxt="n", ylab=name, main=main)
        lapply(p, lines)
      }

      plot_thing(p, "Total perimeter")
      plot_thing(s, "Average solidity")
      plot_thing(a, "Total area")
      plot_thing(b, "Total biodiversity")
      plot_thing(r, "Total risk")

      save_data <- function(x, file) {
        x = as.data.frame(x)
        names(x) = paste0("rep",1:ncol(x))
        x = cbind(time=1:nrow(x), x)
        write.csv(x, file, row.names=FALSE)
      }

      save_data(p, paste0(paste("out/perim", num, size, method, sep="_"), ".csv"))
      save_data(s, paste0(paste("out/shape", num, size, method, sep="_"), ".csv"))
      save_data(a, paste0(paste("out/area", num, size, method, sep="_"), ".csv"))
      save_data(b, paste0(paste("out/biodiversity", num, size, method, sep="_"), ".csv"))
      save_data(r, paste0(paste("out/risk", num, size, method, sep="_"), ".csv"))
    }
  }
  dev.off()
}

