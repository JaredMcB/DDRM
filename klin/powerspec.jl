# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Julia 1.5.0
#     language: julia
#     name: julia-1.5
# ---

# **Compare power spectra estimators**

using PyPlot

j=include("jaredm.jl")

W=include("KLPowerSpec.jl")

pu=include("plotutil.jl")

# **Example: white noise**

w=randn(Int(1e6));

pu.plotcirc(j.z_crossspect_dm(w,w);label="dm")
pu.plotcirc(j.z_crossspect_scalar_ASP(w,w);label="asp")
pu.plotcirc(W.powerspec(w;blks=1600);label="kl")
axis([0,2pi,0,1.3])
grid()
legend()

# **Overdamped Langevin with double well potential**
#
# $$\dot{x}_t = -x_t(x_t^2-1) + \sigma\dot{w}_t$$
#
# Here we solve the SDE with $\sigma=0.3$ and timestep $h=0.01$.

kr=include("kramers.jl")

sde=include("SDE.jl")

sigma=0.3
h=0.01
nskip=100
nsteps=Int(1e8+1e5)
ndrop=Int(1e5)

x=sde.onepath(kr.Drift(), nsteps, h; sigma=sigma)[ndrop:nskip:end];

## sanity check: make sure we really do have a particle jumping back and forth
plot(x[1:div(end,8000):end],".";markersize=1)

# *Test 0: cross spectra*

sig=x[2:end];
pred=x[1:end-1];
dm = j.z_crossspect_dm(sig,pred; L=5000, Nex=2^17);
asp = j.z_crossspect_scalar_ASP(sig,pred; n=2,p=5,nfft=2^17);
aspa = j.z_crossspect_scalar_ASP(sig,pred; nfft=2^17);
kl = W.powerspec(sig,pred; blks=8);

let plt=semilogx
    for fname in [:real,:imag]
        figure()
        let f=eval(fname)
            pu.plotcirc(f.(dm), "b-"; plt=plt, linewidth=1, label="dm ($(length(dm)))")
            pu.plotcirc(f.(asp), "r-"; plt=plt, linewidth=1, label="asp ($(length(asp)))")
            pu.plotcirc(f.(aspa), "r:"; plt=plt, linewidth=1, label="aspa ($(length(aspa)))")
            pu.plotcirc(f.(kl), "g-"; plt=plt, linewidth=1, label="kl ($(length(kl)))")
            grid()
            legend()
            title(fname)
        end
    end
end

# *Test 1: what happens to cross power spectra estimates if we skip more samples?*

x1=x[1:10:end];

sig1=x1[2:end];
pred1=x1[1:end-1];
dm1 = j.z_crossspect_dm(sig1,pred1; L=5000, Nex=2^17);
asp1 = j.z_crossspect_scalar_ASP(sig1,pred1; nfft=2^17);
kl1 = W.powerspec(sig1,pred1; blks=64);

let plt=semilogx
    for fname in [:real,:imag]
        figure()
        let f=eval(fname)
            pu.plotcirc(f.(dm), "b-"; plt=plt, linewidth=1, label="dm ($(length(dm)))")
            pu.plotcirc(f.(dm1), "b:"; plt=plt, linewidth=1, label="dm1 ($(length(dm1)))")
            pu.plotcirc(f.(asp), "r-"; plt=plt, linewidth=1, label="asp ($(length(asp)))")
            pu.plotcirc(f.(asp1), "r:"; plt=plt, linewidth=1, label="asp1 ($(length(asp1)))")
            pu.plotcirc(f.(kl), "g-"; plt=plt, linewidth=1, label="kl ($(length(kl)))")
            pu.plotcirc(f.(kl1), "g:"; plt=plt, linewidth=1, label="kl1 ($(length(kl1)))")
            grid()
            legend()
            title(fname)
        end
    end
end

# *Test 2: what happens if we use fewer points in frequency domain?*

dm2 = j.z_crossspect_dm(sig,pred; L=5000, Nex=2^16);
asp2 = j.z_crossspect_scalar_ASP(sig,pred; nfft=2^16);
kl2 = W.powerspec(sig,pred; blks=16);

let plt=semilogx
    for fname in [:real,:imag]
        figure()
        let f=eval(fname)
            pu.plotcirc(f.(dm), "b-"; plt=plt, linewidth=1, label="dm ($(length(dm)))")
            pu.plotcirc(f.(dm2), "b:"; plt=plt, linewidth=1, label="dm2 ($(length(dm2)))")
            pu.plotcirc(f.(asp), "r-"; plt=plt, linewidth=1, label="asp ($(length(asp)))")
            pu.plotcirc(f.(asp2), "r:"; plt=plt, linewidth=1, label="asp2 ($(length(asp2)))")
            pu.plotcirc(f.(kl), "g-"; plt=plt, linewidth=1, label="kl ($(length(kl)))")
            pu.plotcirc(f.(kl2), "g:"; plt=plt, linewidth=1, label="kl2 ($(length(kl2)))")
            grid()
            legend()
            title(fname)
        end
    end
end

# *Test 3: shorter time series, same h*

x3=x[1:div(end,2)];

sig3=x3[2:end];
pred3=x3[1:end-1];
dm3 = j.z_crossspect_dm(sig3,pred3; L=5000, Nex=2^17);
asp3 = j.z_crossspect_scalar_ASP(sig3,pred3; nfft=2^17);
kl3 = W.powerspec(sig3,pred3; blks=8);

let plt=semilogx
    for fname in [:real,:imag]
        figure()
        let f=eval(fname)
            pu.plotcirc(f.(dm), "b-"; plt=plt, linewidth=1, label="dm ($(length(dm)))")
            pu.plotcirc(f.(dm3), "b:"; plt=plt, linewidth=1, label="dm3 ($(length(dm3)))")
            pu.plotcirc(f.(asp), "r-"; plt=plt, linewidth=1, label="asp ($(length(asp)))")
            pu.plotcirc(f.(asp3), "r:"; plt=plt, linewidth=1, label="asp3 ($(length(asp3)))")
            pu.plotcirc(f.(kl), "g-"; plt=plt, linewidth=1, label="kl ($(length(kl)))")
            pu.plotcirc(f.(kl3), "g:"; plt=plt, linewidth=1, label="kl3 ($(length(kl3)))")
            grid()
            legend()
            title(fname)
        end
    end
end


