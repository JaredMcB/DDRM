using Plots
using StatsBase
using Distributions
pyplot()

include("..\\DataGen.jl")


xs = 1:100
μs = log.(xs)
σs = rand(length(xs))

plot(xs, μs, grid=false, ribbon=σs, fillalpha=.5)

x = rand(400000)

# my_hist(x,50)
P = emp_cdf(x)

P = emp_pdf(x)

E = Exponential()

maximum(E)
minimum(E)

mean(E)

x = 0:.001:10
Y = map(z -> pdf(E,z), x)
plot(x,Y)
plot!(x,exp.(-x),line = (3,:dot))

x = [rand(E) for i = 1:10000]

emp_pdf(x)

E =

Y = map(z -> pdf(E,z),x)
plot(x,Y)
