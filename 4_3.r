x <- read.csv("sample-input.pdb", as.is = T,header = F, sep = "")
x <- x[which(x[,1] == "ATOM"),]
x$V2 <- as.numeric(x$V2)
x$V7 <- as.numeric(x$V7)
x$V8 <- as.numeric(x$V8)
x$V9 <- as.numeric(x$V9)
#Test:
CalcDihedral <- function(p1,p2,p3,p4){
  b1 <- as.vector(t(-1*(p2-p1)))
  b2 <- as.vector(t(p3-p2))
  b3 <- as.vector(t(p4-p3))
  #Make b2 normalized
  b2 <- b2/(sqrt(sum(b2^2)))
  # vector rejections
  # v = projection of b0 onto plane perpendicular to b1
  #   = b0 minus component that aligns with b1
  # w = projection of b2 onto plane perpendicular to b1
  #   = b2 minus component that aligns with b1
  v <- b1 - (sum(b1*b2)) * b2#n1
  w <- b3 - (sum(b3*b2)) * b2#n2
  
  # angle between v and w in a plane is the torsion angle
  # v and w may not be normalized but that's fine since tan is y/x
  x <- t(v) %*% w
  y <- t(crossProduct(b2,v)) %*% w #(t(b2) %*% v) %*% w
  res <- atan2(y,x)*(180/pi)
  return (res)
}
ids <- read.csv("sample-ids.txt",sep = "\t", header = F)
print(paste("Atom Names","Atom Ids","Angles",sep = "    "))
print(paste("A1","A2","A3","A4","A1","A2","A3","A4","Angle (°)",sep = "    "))
calc2 <- c()
pr2 <- c()
for(i in 1:nrow(ids)){
  calc <- CalcDihedral(x[x$V2 == ids[i,1],c(7:9)],x[x$V2 == ids[i,2],c(7:9)],
               x[x$V2 == ids[i,3],c(7:9)],x[x$V2 == ids[i,4],c(7:9)])
  calc2 <- c(calc2,calc)
  pr <- paste(x$V3[(ids[i,1])],x$V3[(ids[i,2])],x$V3[(ids[i,3])],x$V3[(ids[i,4])],
        ids[i,1],ids[i,2],ids[i,3],ids[i,4],calc)
  pr2 <- c(pr2, pr)
  print(pr)
}
write.table(paste("Atom Names","Atom Ids","Angles",sep = "    "),file = "Sample-Output2.txt",col.names = F,row.names = F, append = F)
write.table(paste("A1","A2","A3","A4","A1","A2","A3","A4","Angle (°)",sep = "    "),file = "Sample-Output2.txt",col.names = F,row.names = F,append = T)
write.table(pr2,file = "Sample-Output2.txt",col.names = F,row.names = F,append = T)
# ###Works
# plane1 <- matrix(unlist(c(p1,p2,p3)), byrow = T, ncol = 3)
# plane2 <- matrix(unlist(c(p2,p3,p4)), byrow = T, ncol = 3)
# AB_plane1 <- c(plane1[2,1]-plane1[1,1], plane1[2,2]-plane1[1,2],plane1[2,3]-plane1[1,3])
# AC_plane1 <-c(plane1[3,1]-plane1[1,1], plane1[3,2]-plane1[1,2],plane1[3,3]-plane1[1,3])
# 
# AB_plane2 <- c(plane2[2,1]-plane2[1,1], plane2[2,2]-plane2[1,2],plane2[2,3]-plane2[1,3])
# AC_plane2 <-c(plane2[3,1]-plane2[1,1], plane2[3,2]-plane2[1,2],plane2[3,3]-plane2[1,3])
# 
# crossp1 <- crossProduct(AB_plane1,AC_plane1)
# crossp2 <- crossProduct(AB_plane2,AC_plane2)
# 
# res <- (crossp1[1]*crossp2[1]+crossp1[2]*crossp2[2]+crossp1[3]*crossp2[3])/
#   (sqrt(crossp1[1]^2+crossp1[2]^2+crossp1[3]^2)*sqrt(crossp2[1]^2+crossp2[2]^2+crossp2[3]^2))
# acos(res)*(180/pi)

crossProduct <- function(ab,ac){
  abci = ab[2] * ac[3] - ac[2] * ab[3];
  abcj = ac[1] * ab[3] - ab[1] * ac[3];
  abck = ab[1] * ac[2] - ac[1] * ab[2];
  return (c(abci, abcj, abck))
}

###Works
