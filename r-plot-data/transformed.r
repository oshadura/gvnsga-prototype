data=read.csv("/Users/oksanashadura/CERN_sources/gvnsga-prototype/r-plot-data/transformedulpca.csv", header=TRUE, sep=",")
library(scatterplot3d)
attach(data)
scatterplot3d(data$X1,data$X2,data$X3, color = "blue", pch=19, type = "h", main="Transformed data UPCA")
