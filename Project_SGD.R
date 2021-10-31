install.packages("knitr")
install.packages("rmarkdown")
install.packages('markdown')
library(knitr)
# library(rmarkdown)
# library(markdown)
#-------------------------------------------#
#Part 1 : Simulation of various models

# theoretical model 
g<-function(n){
b<-rbinom(n,1,0.5)
X1<-rnorm(n,1.4,sqrt(0.6))
X2<-rnorm(n,4.8,sqrt(0.5))
Y<-b*X1+(1-b)*X2
return(Y)}

g1<-function(n){
  b<-rbinom(n,1,0.5)
  X1<-rnorm(n,-1,sqrt(0.5))
  X2<-rnorm(n,2,1)
  Y<-b*X1+(1-b)*X2
  return(Y)}


Y<-g(1000)
Z<-g1(1000)
plot(density(Y))
lines(density(Petal.Length),col="red")

dens<-function(x){
  0.5*dnorm(x,-1,sqrt(0.5))+0.5*dnorm(x,2,1)
}

x<-seq(-6,6,0.01)
plot(x,dens(x))
lines(density(Z))

#----------------------------------------------------#
#Part 1
#Objectif plot des 3 modèles
#Ver2.0
#PARTIE 1
#modèle 1
# g1=function(n){
#   u<-rnorm(n,2,1)
#   v<-rnorm(n,-1,sqrt(0.5))
#   b<-(runif(n,0,1)<1/2)
#   X=b*u+(1-b)*v
#   return(X)
# }

#---------------------------------------------------#

#Model 1
g1<-function(x){
  b<-rbinom(n,1,0.5)
  u<-rnorm(n,2,1)
  v<-rnorm(n,-1,sqrt(0.5))
  X=(b)*u+(1-b)*v
  return(X)
}
#Using Cross-Validation in denisty()'s arguments
n=1000
X<-g1(n)
hist(X,prob=T,breaks = 40)
lines(density(X,bw='ucv'),col="red")
plot(density(X),col="blue")

x=runif(n,-2,2)
X<-1/2*dnorm(x,0,0.5)+1/2*dnorm(0,-1,0.5)
plot(density(X))
X1<-g1(50000)
hist(X1,freq=FALSE)

# lines(density(X1,bw=0.01))
#print(hist(X1,freq=FALSE))


#Plot du Hist et density pour le model 1
par(mfrow=c(1,2))
n=1000
X<-g1(n)
hist(X,prob=T,breaks = 40,main="Histogramme_Model_1")
lines(density(X,kernel="gaussian",bw="nrd0"),col="red")
plot(density(X),col="blue",main = "Density_model_1")

#--------------------------------------------------#

#model 2

g2=function(n){
  u<-runif(n,0,1)
  v<-rnorm(n,-1,0.5)
  b<-(runif(n,0,1)<1/2)
  X=(b)*u+(1-b)*v
  return(X)
}
par(mfrow=c(1,2))
n=1000
X2<-g2(n)
hist(X2,breaks = 40,freq=FALSE,main = 'histogramme_model 2')
lines(density(X2,bw="nrd0"),col='red')
plot(density(X2),col="blue",main = "Density_model_2")

#model 3
g3=function(n){
  u<-rgamma(n,2,4)
  v<-rgamma(n,2,1)
  b<-(runif(n,0,1)<1/2)
  X=(b)*u+(1-b)*v
  return(X)
}
X3<-g3(50000)
hist(X3,freq=FALSE)

par(mfrow=c(1,2))
n=1000
X3<-g3(n)
hist(X3,breaks = 40,freq=FALSE,main = 'histogramme_model_3')
lines(density(X3,bw="nrd0"),col='red')#bw egale à h, on réduit
plot(density(X3),col="blue",main = "Density_model_3")


lines(density(X3,bw=0.1))

par(mfrow=c(1,1))

hist(X,prob=T,breaks=40)

#graphical approximation  of the bandwidth  h  
bws=c(0.3,0.39,0.5,0.1,0.17)

plot(density(X),main = "Approximation_h_optimal")
for (i in 1:length(bws)){
  lines(density(X,kernel="gaussian",bw=bws[i]),col=i)
}
legend(4,0.2,15,legend=c('bw'=0.3,'bw'=0.39,'bw'=0.5,'bw'=0.1,'bw'=0.17),
       col=c("black","red","green","blue","lightblue","magenta")
       ,lty=1:5)

# lines(density(X1,kernel="gaussian"))

for(i in 1:100){
  lines(density(X1,kernel="gaussian",n=i*30),col="red")
}
#-------------------------------------------------#
#Theoritical approximation of h using MISE (Optimization)
#Part2

n=1000

g<-function(x){
  b<-rbinom(n,1,0.5)
  u<-rnorm(n,2,1)
  v<-rnorm(n,-1,sqrt(0.5))
  X=(b)*u+(1-b)*v
  return(X)
}
X<-g(n) #Simulating X  using model 1
h=seq(0.001,2,0.01)# n=length(h)
a=min(X)#Value border inf integral
b=max(X)#Value border sup  integral


#Estimateur noyau gaussian test 
estimateur_noyau_gauss1=function(X,x,h){
  return(1/(length(X)*h)*sum(noyau_gauss(X-x)/h))
}
#Construction Estimateur_noyau_gaussien 

#Noyau gaussien K()
noyau_gauss=function(x){
  return(1/sqrt(2*pi)*exp(-x**2/2))
}

estimateur_noyau_gauss=function(X,x,h){
  #h vector, X vector, n scalar
  res=0
  n=length(X)
  # res_final=c()
  # for (i in 1:length(h))
  for(j in 1:n){
    res=res+noyau_gauss((x-X[j])/h)
  }
  return(1/(n*h)*res)}


#Plot de l'estimateur gaussian kernel 
n=1000
m<-seq(-6,6,0.01)
p<-g(n)
d=c()
for(i in 1:length(m)){
  d[i]=estimateur_noyau_gauss(p,m[i],0.3957)
}
plot(m,d,main = "Estimateur_noyau_gauss")
lines(density(X,kernel="gaussian",bw="nrd0"),col="red")
legend(x=1.5,y=0.28,legend = c("Estimateur noyau","Density(X)"),col= c('black',"red")
       ,lwd = 2)
densite = function(x){
  return(0.5*(dnorm(x,2,1))+0.5*(dnorm(x,-1,sqrt(0.5))))
}


#f_int<-function(x) (estimateur_noyau_gauss(x,n,h,X)-densite(x))**2

#h=seq(0.001,2,0.01)

#Fct qui permet d'avoir comme sortie un vecteur représetant les intégrales prises selon
#chaque valeur de h
fct_int=function(X,h){
  int_mc=c()
  a=min(X)
  b=max(X)
  for (i in (1:length(h))){
    f_int<-function(x) (estimateur_noyau_gauss(X,x,h=h[i])-densite(x))**2
    #Each integral in linked to a h[i] of h
    int_mc[i]=integrate(f_int,lower=a,upper=b,subdivisions =2000)$value
  }
  return(int_mc)}


#Calulating Matrix of Monte-carlo to approximate the mean
Mat_mc<-function(nb_iteration,n,h){
  M=matrix(nrow = nb_iteration,ncol=length(h))
  for (i in 1:nb_iteration){
    print(i)
    x=g(nb_iteration)
    M[i,]=fct_int(X,h)
  }
  MISE_final=c()
  for (i in 1:length(h)){
    MISE_final[i]=mean(M[,i])
  }
  return(MISE_final)
}   

j=Mat_mc(30,n,h)

plot(h,j,main = 'Trajectoire MISE')
m=which.min(j)
h_hat=h[m]
h_hat


#Graphical Comparison with the bandwidth obtained MISE and with cross-validation
#Plot using Cross-validation
plot(density(X,bw="ucv"),col='blue',main = 'Estimation par Cross validation')
lines(density(X,bw=h_hat),col='green')
lines(density(X))
legend(x="topright",legend = c("Estimation par VC","Estimation par h_hat",'Density (X)'),lty=c(1,1,1),col= c('blue',"green",1)
       ,lwd = 3)
#----------------------------------------------------------------#
#Partie 3
#Using Real data world to approximate the density (using kernel density)
data("iris")
attach(iris)
data1=Petal.Length
plot(density(data1,bw="nrd0"),col="red",main = "Density Sepal.Length")
#Comparison between density and approximation using CV
#Plot using Cross-validation
plot(density(data1,bw="ucv",kernel = 'gaussian'),col="blue",main = "Estimation de la density par validation croisee")
lines(density(data1),col="red")
legend(x="topright",legend = c("Estimation par VC",'Density (Sepal.Length)'),lty=c(1,1),col= c('blue',"red")
       ,lwd = 3)



#------------------------------------------------------------#
#Taille de l'échantillon
n=1000
X=g(n)
#Programme pour le cross-validation
cv_gaussien = function(h){
  new_data = sort(X)
  #membre de droite : 
  k = 0
  for (i in (1:n)){
    for (j in (1:n))
      if(i!=j)
        k=k+noyau_gauss((new_data[j]-new_data[i])/h)
  }
  k=k*(2/(n*h*(n-1)))
  
  #estimation
  integrale=0
  for(i in (2:n)){
    print(i)
    integrale = integrale + (new_data[i]-new_data[i-1])*(estimateur_noyau_gauss(X,new_data[i],h))^2
  }
  return(integrale-k)
}




#hmin
h_cv_gaussien = optimise(cv_gaussien,lower=0.001,upper=2,maximum=FALSE)$minimum

#Comparaison entre h_cv_gaussien(algo CV) et density() par Cross-validation
plot(density(X,bw=h_cv_gaussien),col="blue",main = 'Estimation h_cv_gaussien')
lines(density(X,bw="ucv"),col="green")
legend(x="topright",legend = c("h_cv_gaussien",'Density-CV'),lty=c(1,1),col= c('blue',"green")
       ,lwd = 3)



#Comparaison entre la fenetre h_hat obtenue par MISE 
#et h_cv_gaussien par Cross-validation 
bws=c(0.39,h_hat,h_cv_gaussien)
col=c('red','black','green')

plot(density(X),main = "Comparaison entre les estimateurs de densité")
for (i in 1:length(bws)){
  lines(density(X,kernel="gaussian",bw=bws[i]),col=col[i])
  legend(x='topright',legend=c("bw"='0.39','bw'="h_hat",'bw'="h_cv_gaussien"),
         col=c("red","black","green"),lty=c(1,1,1))
}