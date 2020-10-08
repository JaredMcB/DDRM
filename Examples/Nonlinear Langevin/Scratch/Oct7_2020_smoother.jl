n = 4
p = 150
ty = "bin"
μ = _smoother(n,p; ty)
round(sum(μ);digits = 5) == 1.0 || println("bad smoother")
plot(-n*p:n*p,μ)

function smoother_plot(n,p,ty)
    μ = _smoother(n,p; ty)
    round(sum(μ);digits = 5) == 1.0 || println("bad smoother")
    plot(-n*p:n*p,μ)
end

foo(3,500,"ave")
