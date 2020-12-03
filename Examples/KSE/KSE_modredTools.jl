
module KSE_modredTools

using JLD
using DSP: conv # For conv function in Psi
using Dates
using Distributions

mr = include("../../Tools/Model_Reduction_Dev.jl")

function InvBurgRK4_1step(x;h,obs_gap,P,N)
 lx = length(x)
 function F(x)
     𝑥 = [conj(reverse(x, dims = 1));0; x]
     -im/2*(2π/P*(1:lx)/N) .* conv(𝑥,𝑥)[2*lx+2:3*lx+1]
 end

 Δt = h*obs_gap

 k1 = F(x)
 k2 = F(x .+ Δt*k1/2)
 k3 = F(x .+ Δt*k2/2)
 k4 = F(x .+ Δt*k3)
 A =  @. x + Δt/6*(k1 + 2k2 + 2k3 + k4)
end

function Inertialman_part(x)
  lx = length(x)
  𝑥(j) = ( j <= lx ? x[j] : im*sum(x[l]*x[j-l] for l = j-lx:lx) )

  L = complex(zeros(lx^2))
  for j = 1:lx
     for k = 1:lx
        L[ (j-1)*lx+k] = 𝑥(j+lx)*𝑥(j+lx-k)
     end
  end
  L
end

function Inertialman_part_short(x)
  lx = length(x)
  𝑥(j) = ( j <= lx ? x[j] : im*sum(x[l]*x[j-l] for l = j-lx:lx) )

  L = complex(zeros(binomial(lx+1,2)))
  i = 1
  for j = 1:lx
     for k = j:lx # k should normaly go from 1 to lx but i changed it to go from j to lx.
        L[i] = 𝑥(j+lx)*𝑥(j+lx-k)
        i += 1
     end
  end
  L
end

function PSI(x; short = true, h,obs_gap,P,N)
   short ? [x; InvBurgRK4_1step(x;h,obs_gap,P,N); Inertialman_part_short(x)] :
           [x; InvBurgRK4_1step(x;h,obs_gap,P,N); Inertialman_part(x)]
end



end #module
