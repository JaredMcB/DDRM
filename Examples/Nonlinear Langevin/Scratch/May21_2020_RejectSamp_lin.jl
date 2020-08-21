"""
Rejection Sampling
"""

## Load libraries
using PyPlot, PlotUtils

u(x) = (x^2 - 1)^2 ## double-well
v(x) = 0.5*x^2 - 1 ## shift harmonic potential down so v(x) <= u(x) for all x

xgrid=-2.5:0.1:2.5;
let x = xgrid
    plot(x,u.(x))
    plot(x,v.(x))
end
grid()

p(x) = exp(-u(x))
q(x) = exp(-v(x))

let x=xgrid
    plot(x,p.(x))
    plot(x,q.(x))
end
grid()

x = randn(40000)

xy = [[x,rand()*q(x)] for x in x]

plot([xy[1] for xy in xy],[xy[2] for xy in xy],".",markersize=1)
let x=xgrid
    #plot(x,p.(x))
    plot(x,q.(x),"-",linewidth=3)
end
grid()

z = filter(xy->xy[2]<=p(xy[1]),xy)

plot([xy[1] for xy in xy],[xy[2] for xy in xy],".",markersize=1)
plot([xy[1] for xy in z],[xy[2] for xy in z],".",markersize=1)
let x=xgrid
plot(x,p.(x),"-",linewidth=3,label="p(x)")
plot(x,q.(x),"-",linewidth=3,label="q(x)")
end
grid()
legend()

α =length(z)/length(xy)
1/(sqrt(2*pi)*exp(1))/α

#plotpdf([z[1] for z in zs],".-";norm=:pdf,n=100)
let x=xgrid
    plot(x,0.511*p.(x),"-",linewidth=3,label="target pdf")
end
plotpdf([z[1] for z in z],".-";n=100,norm=:pdf,label="empirical pdf")
grid()
legend()
