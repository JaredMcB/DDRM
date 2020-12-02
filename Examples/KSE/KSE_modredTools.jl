
module KSE_modredTools

using JLD
using DSP: conv # For conv function in Psi
using Dates

mr = include("Model_Reduction_Dev.jl")

function InvBurgRK4_1step(x)
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

PSI(x; short = true) = short ? [x; InvBurgRK4_1step(x); Inertialman_part_short(x)] :
                [x; InvBurgRK4_1step(x); Inertialman_part(x)]


function modredrun(;
   X,
   h_wf,
   Psi,
   steps,
   discard,
   noise = false,
   Nosie_dist
   )

   d, nu, M_out = size(h_wf)

   pred = mr.get_pred(X[:,1:M_out],psi)


   ## Run reduced model with no noise
   steps = size(V_obs,2)
   V_rm = [V_obs[:,1:M_out] complex(zeros(d,steps-M_out))]
   nu = size(Psi(V_obs[:,1]),1)

   # load presamples
   PSI_past = complex(zeros(nu,steps))
   for i=1:M_out
       PSI_past[:,i] = Psi(V_obs[:,i])
   end

   # Move forward without original data
   print("Reduced Model Run Time: ")
   @time for i = M_out+1:steps
       V_rm[:,i] = sum(h_wf[:,:,k]*PSI_past[:,i-k] for k = 1:M_out)
       isnan(V_rm[1,i]) && break
       PSI_past[:,i] = Psi(V_rm[:,i])
   end

   # Process reduced model run
   steps_rm = size(V_rm,2)
   V_rm_end = conj(reverse(V_rm,dims = 1))
   VV_rm = [zeros(1,steps_rm); V_rm; zeros(N-2d-1,steps_rm);V_rm_end]
   UU_rm = real(ifft(VV_rm,1))
   tt_rm = tt[1:end]

   dat = Dict(
      "dat_UU_rm" => UU_rm,
      "dat_VV_rm" => VV_rm,
      "dat_tt_rm" => tt_rm)
   Data = merge(paramaters,dat)
   save(rmrun_file,Data)
   println("Reduced Model Run saved")
end

end #module
