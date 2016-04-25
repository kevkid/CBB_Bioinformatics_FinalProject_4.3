library(Rpdb)
setwd("~/CBB_BioInformatics_4.3")
#some problem with the way it is supposed to normalize
CalcDihedral <- function(p1,p2,p3,p4){
  b1 <- as.vector(t(-1*(p2-p1)))
  b2 <- as.vector(t(p3-p2))
  b3 <- as.vector(t(p4-p3))
  #Make b1 normalized
  b1 <- b1/(sum(sqrt(b1^2)))
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
CalcDihedral(p1,p2,p3,p4)

x <- read.pdb("1g8p.pdb")
dihedral(x, 1,2,3,4)
x$atoms[1,9:11]
#Test:
p1 <- x$atoms[1,9:11]
p2 <- x$atoms[2,9:11]
p3 <- x$atoms[3,9:11]
p4 <- x$atoms[4,9:11]

b1 <- as.vector(t(-1*(p2-p1)))
b2 <- as.vector(t(p3-p2))
b3 <- as.vector(t(p4-p3))
#Make b1 normalized
b1 <- b1/(sum(sqrt(b1^2)))

# vector rejections
# v = projection of b0 onto plane perpendicular to b1
#   = b0 minus component that aligns with b1
# w = projection of b2 onto plane perpendicular to b1
#   = b2 minus component that aligns with b1
n1 <- b1 - ((b1%*%b2)) * b2#n1
n2 <- b3 - ((b3%*%b2)) * b2#n2

# angle between v and w in a plane is the torsion angle
# v and w may not be normalized but that's fine since tan is y/x
x <- sum(n1 * n2)
y <- (crossProduct(b2, n1)%*%n2) #(t(b2) %*% v) %*% w
res <- atan2(y,x)*(180/pi)
res

###Works
plane1 <- matrix(unlist(c(p1,p2,p3)), byrow = T, ncol = 3)
plane2 <- matrix(unlist(c(p2,p3,p4)), byrow = T, ncol = 3)
AB_plane1 <- c(plane1[2,1]-plane1[1,1], plane1[2,2]-plane1[1,2],plane1[2,3]-plane1[1,3])
AC_plane1 <-c(plane1[3,1]-plane1[1,1], plane1[3,2]-plane1[1,2],plane1[3,3]-plane1[1,3])

AB_plane2 <- c(plane2[2,1]-plane2[1,1], plane2[2,2]-plane2[1,2],plane2[2,3]-plane2[1,3])
AC_plane2 <-c(plane2[3,1]-plane2[1,1], plane2[3,2]-plane2[1,2],plane2[3,3]-plane2[1,3])

crossp1 <- crossProduct(AB_plane1,AC_plane1)
crossp2 <- crossProduct(AB_plane2,AC_plane2)

res <- (crossp1[1]*crossp2[1]+crossp1[2]*crossp2[2]+crossp1[3]*crossp2[3])/
  (sqrt(crossp1[1]^2+crossp1[2]^2+crossp1[3]^2)*sqrt(crossp2[1]^2+crossp2[2]^2+crossp2[3]^2))
acos(res)*(180/pi)

crossProduct <- function(ab,ac){
  abci = ab[2] * ac[3] - ac[2] * ab[3];
  abcj = ac[1] * ab[3] - ab[1] * ac[3];
  abck = ab[1] * ac[2] - ac[1] * ab[2];
  return (c(abci, abcj, abck))
}

###Works

std <- function(x){if(length(which(is.na(x)))==0) (x-mean(x))/sd(x) else
  
  (x-mean(x,na.rm=T))/sd(x,na.rm=T)
}