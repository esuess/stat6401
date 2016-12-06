
# Statistics 6401, Fall 2017

## Project

This R program transforms independent BVN to (X1,X2) correlated BVN (W1,W2)

## Simulation

Set up the parameters for the simulations.

The goal is to simulate the the $BVN(10, 25, 2, 3, -0.4)$


```R
n = 2000  # number of values simulated

muw1 = 10
muw2 = 25

sigmaw1 = 2
sigmaw2 = 3

sigmasqw1 = sigmaw1^2
sigmasqw2 = sigmaw2^2

rhow1w2 = -0.4
covw1w2 = rhow1w2*sigmaw1*sigmaw2

```

*****
## Figure 1

Plot the BVN


```R
mu1<-muw1 # setting the expected value of x1
mu2<-muw2 # setting the expected value of x2
s11<-sigmasqw1 # setting the variance of x1
s12<-covw1w2 # setting the covariance between x1 and x2
s22<-sigmasqw2 # setting the variance of x2
rho<-rhow1w2 # setting the correlation coefficient between x1 and x2
x1<-seq(mu1-10,mu1+10,length=41) # generating the vector series x1 
x2<-seq(mu2-10,mu2+10,length=41) 

f<-function(x1,x2){
	term1 <- 1/(2*pi*sqrt(s11*s22*(1-rho^2))) 
	term2 <- -1/(2*(1-rho^2))
	term3 <- (x1-mu1)^2/s11
	term4 <- (x2-mu2)^2/s22
	term5 <- -2*rho*((x1-mu1)*(x2-mu2))/(sqrt(s11)*sqrt(s22)) 
	term1*exp(term2*(term3+term4-term5)) 
} # setting up the function of the multivariate normal density

z<-outer(x1,x2,f) # calculating the density values 

persp(x1, x2, z, 
      main="Two dimensional Normal Distribution",
      sub=expression(italic(f)~(bold(x))==frac(1,2~pi~sqrt(sigma[11]~ 
                     sigma[22]~(1-rho^2)))~phantom(0)~exp~bgroup("{", 
	             list(-frac(1,2(1-rho^2)), 
	             bgroup("[", frac((x[1]~-~mu[1])^2, sigma[11])~-~2~rho~frac(x[1]~-~mu[1],
	             sqrt(sigma[11]))~ frac(x[2]~-~mu[2],sqrt(sigma[22]))~+~ 
	             frac((x[2]~-~mu[2])^2, sigma[22]),"]")),"}")),
      col="lightgreen", 
      theta=30, phi=20,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.75,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot

# adding a text line to the graph
mtext(expression(list(mu[1]==10,mu[2]==25,sigma[11]==2,sigma[22]==3,sigma[12]==-2.4,rho==-0.4)), side=3) 
```


![png](output_4_0.png)


*****
## Figure 2 

Simulate independent random uniforms


```R
u1 = runif(n)
u2 = runif(n)

par(mfrow=c(2,2))
hist(u1)
hist(u2)
plot(u1,u2)
```


![png](output_6_0.png)


Check the means, standard deviation, and the correlation between the uniforms


```R
cat('means')
umean = c(mean(u1),mean(u2)); umean
cat('standard deviations')
usd = c(sd(u1),sd(u2)); usd
cat('correlation')
ucor = cor(u1,u2); ucor
```

    means


<ol class=list-inline>
	<li>0.504959828843363</li>
	<li>0.507045170075144</li>
</ol>



    standard deviations


<ol class=list-inline>
	<li>0.289973714873899</li>
	<li>0.290818497431585</li>
</ol>



    correlation


0.0208548405232158


*****
## Figure 3

Simulate independent standard normals $BVN(0,0,1,1,0)$ - Box-Muller Method


```R
x1 = sqrt(-2*log(u1))*cos(2*pi*u2)
x2 = sqrt(-2*log(u1))*sin(2*pi*u2)

par(mfrow=c(2,2))
hist(x1)
hist(x2)
plot(x1,x2)
```


![png](output_10_0.png)


Check the means, standard deviation, and the correlation between the independent standard normals


```R
cat('means')
xmean = c(mean(x1),mean(x2));xmean
cat('standard deviationss')
xsd = c(sd(x1),sd(x2));xsd
cat('correlations')
xcor = cor(x1,x2);xcor
```

    means


<ol class=list-inline>
	<li>0.0182945023485414</li>
	<li>-0.0191937856762289</li>
</ol>



    standard deviationss


<ol class=list-inline>
	<li>1.00846822635981</li>
	<li>0.984440694041881</li>
</ol>



    correlations


0.00591407686078222


*****
## Figure 4

Transform to correlated normals $BVN(0,0,\sigma_{W_1}^2,\sigma_{W_2}^2,\rho)$


```R
c11 = sigmaw1
c21 = rhow1w2*sigmaw2
c22 = sigmaw2*sqrt(1-rhow1w2^2)

w1 = c11*x1
w2 = c21*x1 + c22*x2

par(mfrow=c(2,2))
hist(w1)
hist(w2)
plot(w1,w2)
```


![png](output_14_0.png)


Check the means, standard deviation, and the correlation between the correlations normals with mean zero


```R
cat('means')
wmean = c(mean(w1),mean(w2))
wmean
cat('standard deviations')
wsd = c(sd(w1),sd(w2))
wsd
cat('correlation')
wcor = cor(w1,w2)
wcor
```

    means


<ol class=list-inline>
	<li>0.0365890046970827</li>
	<li>-0.0747275882586964</li>
</ol>



    standard deviations


<ol class=list-inline>
	<li>2.01693645271963</li>
	<li>2.95843213405095</li>
</ol>



    correlation


-0.403644161780846


*****
## Figure 5

Transform to add means $BVN(\mu_{W_1},\mu_{W_2},\sigma_{W_1}^2,\sigma_{W_2}^2,\rho)$


```R
y1 = muw1 + w1
y2 = muw2 + w2

par(mfrow=c(2,2))
hist(y1)
hist(y2)
plot(y1,y2)
```


![png](output_18_0.png)


Check the means, standard deviation, and the correlation between the correlations normals


```R
cat('means')
ymean = c(mean(y1),mean(y2))
ymean
cat('standard deviation')
ysd = c(sd(y1),sd(y2))
ysd
cat('correlation')
ycor = cor(y1,y2);ycor
```

    means


<ol class=list-inline>
	<li>10.0365890046971</li>
	<li>24.9252724117413</li>
</ol>



    standard deviation


<ol class=list-inline>
	<li>2.01693645271963</li>
	<li>2.95843213405095</li>
</ol>



    correlation


-0.403644161780846


*****
## Figure 6

Now rotate by $\theta$ to make independent $BVN(\mu_{W_1},\mu_{W_2},\sigma_{W_1}^2,\sigma_{W_2}^2,0)$ again


```R
theta = 0.5*atan((2*rhow1w2*sigmaw1*sigmaw2)/(sigmasqw1-sigmasqw2))

cat('theta')
theta

v1 = y1*cos(theta) + y2*sin(theta)
v2 = -y1*sin(theta) + y2*cos(theta)
par(mfrow=c(2,2))
hist(v1)
hist(v2)
plot(v1,v2)
```

    theta


0.382496416355455



![png](output_22_2.png)


Check the means, standard deviation, and the correlation between the bivariate normals with zero correlation


```R
cat('means')
vmean = c(mean(v1),mean(v2))
vmean
cat('standard deviation')
vsd = c(sd(v1),sd(v2))
vsd
cat('correlation')
vcor = cor(v1,v2)
vcor
```

    means


<ol class=list-inline>
	<li>18.6143526075603</li>
	<li>19.3780339727778</li>
</ol>



    standard deviation


<ol class=list-inline>
	<li>1.74716865077951</li>
	<li>3.12534078971043</li>
</ol>



    correlation


-0.0211474306714942


*****
## Figure 7

Using the mvrnorm( ) function from the library MASS


```R
library(MASS)

cat('means')
mu = c(muw1,muw2);mu

cat('Variance-Covariance Matrix')
Sigma = matrix(c(sigmasqw1,covw1w2,covw1w2,sigmasqw2),2,2);Sigma

aa = mvrnorm(n, mu, Sigma)

par(mfrow=c(2,2))
hist(aa[,1])
hist(aa[,2])
plot(aa)
```

    means


<ol class=list-inline>
	<li>10</li>
	<li>25</li>
</ol>



    Variance-Covariance Matrix


<table>
<tbody>
	<tr><td> 4.0</td><td>-2.4</td></tr>
	<tr><td>-2.4</td><td> 9.0</td></tr>
</tbody>
</table>




![png](output_26_4.png)


Check the means, standard deviation, and the correlation between the correlations normals


```R
cat('means')
aamean = c(mean(aa[,1]),mean(aa[,2]))
aamean
cat('standard deviations')
aasd = c(sd(aa[,1]),sd(aa[,2]))
aasd
cat('correlation')
aacor = cor(aa[,1],aa[,2]);aacor
```

    means


<ol class=list-inline>
	<li>10.0814184259811</li>
	<li>24.9680273388844</li>
</ol>



    standard deviations


<ol class=list-inline>
	<li>1.99246311209493</li>
	<li>3.04521262059864</li>
</ol>



    correlation


-0.380091248739976


*****
## Scatterplots for large data

If you do not have the following R packages installed, you should uncomment the following two lines in the next cell.


```R
#install.packages("IDPmisc")
#install.packages("gplots")

library(IDPmisc)
library(gplots) 
```

## Figure 8

See Figure 2.

Plot of independent uniforms


```R
ipairs(matrix(c(u1,u2),n,2))
```


![png](output_32_0.png)


*****
## Figure 9

See Figure 3.

Plot of independent standard normals


```R
ipairs(matrix(c(x1,x2),n,2))
```


![png](output_34_0.png)


*****
## Figure 10

See Figure 4.

Plot of correlated normals with mean zero


```R
ipairs(matrix(c(w1,w2),n,2))
```


![png](output_36_0.png)


*****
## Figure 11

See Figure 5.

Plot of correlated normals


```R
ipairs(matrix(c(y1,y2),n,2))
```


![png](output_38_0.png)


*****
## Figure 12

See Figure 6.

Plot of uncorrelated normals


```R
ipairs(matrix(c(v1,v2),n,2))
```


![png](output_40_0.png)


*****
## Figure 13

See Figure 7.

Plot of correlated normals


```R
ipairs(aa)
```


![png](output_42_0.png)


*****
## Figure 14

The same plots but using the hist2d() function


```R
par(mfrow=c(2,2))
hist2d(u1,u2, nbins=50, col = c("white",heat.colors(16))) 
box()
hist2d(x1,x2, nbins=50, col = c("white",heat.colors(16)))
box()
hist2d(w1,w2, nbins=50, col = c("white",heat.colors(16)))
box()
hist2d(y1,y2, nbins=50, col = c("white",heat.colors(16)))
box()
```


    
    ----------------------------
    2-D Histogram Object
    ----------------------------
    
    Call: hist2d(x = u1, y = u2, nbins = 50, col = c("white", heat.colors(16)))
    
    Number of data points:  2000 
    Number of grid bins:  50 x 50 
    X range: ( 0.0001644774 , 0.9988243 )
    Y range: ( 0.0004471866 , 0.9996592 )




    
    ----------------------------
    2-D Histogram Object
    ----------------------------
    
    Call: hist2d(x = x1, y = x2, nbins = 50, col = c("white", heat.colors(16)))
    
    Number of data points:  2000 
    Number of grid bins:  50 x 50 
    X range: ( -3.556451 , 3.645411 )
    Y range: ( -3.566148 , 3.325337 )




    
    ----------------------------
    2-D Histogram Object
    ----------------------------
    
    Call: hist2d(x = w1, y = w2, nbins = 50, col = c("white", heat.colors(16)))
    
    Number of data points:  2000 
    Number of grid bins:  50 x 50 
    X range: ( -7.112901 , 7.290823 )
    Y range: ( -10.1398 , 11.58378 )




    
    ----------------------------
    2-D Histogram Object
    ----------------------------
    
    Call: hist2d(x = y1, y = y2, nbins = 50, col = c("white", heat.colors(16)))
    
    Number of data points:  2000 
    Number of grid bins:  50 x 50 
    X range: ( 2.887099 , 17.29082 )
    Y range: ( 14.8602 , 36.58378 )




![png](output_44_4.png)

