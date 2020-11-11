@assert VERSION >= v"0.6"

module plotutil

using PyPlot,FFTW,StatsBase

export   hist
function hist(data, nbins)
    r=fit(Histogram, data, nbins=nbins, closed=:left)
    return r.edges[1],r.weights
end

export   plotcdf
function plotcdf(xl, args...; flags...)
    let n=length(xl)
        plot(sort(xl),(0:(n-1))/n,args...; flags...)
    end
end

export   plotpdf
function plotpdf(xl, args...; n=50, norm=:count, plt=plot,flags...)
    xl,cl = hist(xl,n)
    Z=1
    if norm === :pdf
        Z=sum(cl)*(xl[2]-xl[1])
    elseif norm === :prob
        Z=sum(cl)
    elseif norm !== :count
        error("unknown normalization: $(norm)")
    end 
    plt(vcat(xl[1:1],0.5*(xl[1:(end-1)]+xl[2:end]),xl[end:end]),vcat([0],cl/Z,[0]),args...; flags...)
end

export   plotcirc
function plotcirc(x, args...; plt=plot, circ=2*pi, shift=false, flags...)
    n = length(x)
    if shift
        plt([i/n*circ-circ/2 for i in 0:(n-1)], fftshift(x), args...; flags...)
    else
        plt([i/n*circ for i in 0:(n-1)], x, args...; flags...)
    end
end

export   tplot
function tplot(x, args...; dt=1.0, t0=0, plt=plot, shift=false, flags...)
    n = length(x)
    if shift
        plt([t0+dt*i-(n-1)*dt/2 for i in 0:(n-1)],
            fftshift(x), args...;
            flags...)
    else
        plt([t0+dt*i for i in 0:(n-1)], x, args...; flags...)
    end
end

export   txtspy
function txtspy(A)
    m,n=size(A)
    print("  \\col\t")
    c=0
    for j=1:n
        if rem(j,10) in [0,5] || j==1
            print(j)
            c = floor(Int,log10(j))
        elseif c > 0
            c -= 1
        else
            print(" ")
        end
    end
    println()
    print("row\\   \t")
    for j=1:n
        if rem(j,10) in [0,5] || j==1
            print("|")
        else
            print(" ")
        end
    end
    println()
    for i=1:m
        print("$(i)\t")
        for j=1:n
            if A[i,j] == 0
                print(".")
            else
                print("*")
            end
        end
        println()
    end
end

################################
## Useful for saving figures in a notebook.
mutable struct PrintFigOpt
  doprint::Bool
  force::Bool
  prefix::String
end

pfopt = PrintFigOpt(false, false, "")

function printfig_enable()
    pfopt.doprint = true
end

function printfig_disable()
    pfopt.doprint = false
end

function printfig_setprefix(prefix)
    pfopt.prefix = prefix
end

function printfig(name, args...;
                  prefix  = pfopt.prefix,
                  doprint = pfopt.doprint,
                  ext     = ".pdf",
                  force   = pfopt.force,
                  flags...
                  )
    
    fullname = prefix*name*ext

    if isfile(fullname)
        if force
            println("# printfig(\"$fullname\"): removing existing file")
            rm(fullname)
        else
            println("# printfig(\"$fullname\"): file exists, skipping")
            return
        end
    end
    
    if doprint || pfopt.doprint
        println("# printfig(\"$fullname\"): printing")
        savefig(fullname, args...; flags...)
    else
        println("# printfig(\"$fullname\"): disabled")
    end
end

end#module
