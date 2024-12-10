## double checking some weird cells
cell.test <- "FOV4_CR6387105db-61105"
data.cell <- data[data$cell == cell.test,]
head(data.cell)

par(mfrow=c(1,2))
plot(data.cell$y, data.cell$z)
points(data.cell["980842",]$y, data.cell["980842",]$z, col='red', pch=16)
plot(data.cell$transformedY, data.cell$transformedZ)
points(data.cell["980842",]$transformedY, data.cell["980842",]$transformedZ, col='red', pch=16)
