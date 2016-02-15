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

pdf("output.pdf", width=12, height=8)
par(mfrow=c(3,4), mar=c(1,4,2,2))
for (num in nums) {
  for (size in sizes) {
    for (method in methods) {
      cat("Processing location", num, ",", size, "for method", method, "\n")
      p = list()
      s = list()
      a = list()
      for (rep in reps) {
        # biodiversity data is in the _epi file
        filename = paste(num, size, method, rep, sep="_")

        # main data is here...
        d2 = readMat(file.path("data", paste0(filename, ".mat")))$simulation.data

        p[[rep]] = unlist(lapply(d2["edge.patches",,], mean))
        s[[rep]] = unlist(lapply(d2["shape",,], mean))
        a[[rep]] = unlist(lapply(d2["size.patches",,], mean))
      }

      p_max = max(unlist(lapply(p, max)))
      main = paste(method, num, size, sep="_")
      plot(NULL, ylim=c(0,p_max), xlim = c(1,length(p[[1]])), xlab="", xaxt="n", ylab="Average perimeter", main=main)
      lapply(p, lines)

      s_max = max(unlist(lapply(s, max)))
      main = paste(method, num, size, sep="_")
      plot(NULL, ylim=c(0,s_max), xlim = c(1,length(s[[1]])), xlab="", xaxt="n", ylab="Average solidity", main=main)
      lapply(s, lines)

      a_max = max(unlist(lapply(a, max)))
      main = paste(method, num, size, sep="_")
      plot(NULL, ylim=c(0,a_max), xlim = c(1,length(a[[1]])), xlab="", xaxt="n", ylab="Average area", main=main)
      lapply(a, lines)

      p = as.data.frame(p); names(p) <- paste0("rep",1:ncol(p)); p = cbind(time=1:nrow(p), p)
      s = as.data.frame(s); names(s) <- paste0("rep",1:ncol(s)); s = cbind(time=1:nrow(s), s)
      a = as.data.frame(a); names(a) <- paste0("rep",1:ncol(a)); a = cbind(time=1:nrow(a), a)

      write.csv(p, paste0(paste("out/perim", num, size, method, sep="_"), ".csv"), row.names=FALSE)
      write.csv(s, paste0(paste("out/shape", num, size, method, sep="_"), ".csv"), row.names=FALSE)
      write.csv(a, paste0(paste("out/area", num, size, method, sep="_"), ".csv"), row.names=FALSE)
    }
  }
}
dev.off()

