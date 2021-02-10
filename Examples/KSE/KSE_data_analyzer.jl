### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 61204a48-6710-11eb-2fb8-099c1dbe495d
begin
	using PyPlot
	using JLD
	using Statistics: mean
	import AnalysisToolbox.jl: my_autocor
end

# ╔═╡ da51dedc-670f-11eb-3b0b-e9a19ce28437
md"""
# Analyze Data
"""

# ╔═╡ bb25ac5a-6712-11eb-276a-8dd819a3000c
vv = load("/media/extradrive1/jaredm/KSE_Data/KSE_sol_lin.jld","dat_vv")[:,1:50000]

# ╔═╡ 89c28b86-672d-11eb-351b-f9897c5ff59f
md"""
Here are the autocorrelations. These match pretty well, at least by eye.
"""

# ╔═╡ 3acc056c-6715-11eb-3d61-952f5674914e
begin
	close("all")
	for i = 1:5
		A = at.my_autocor(vv[i+1,:],0:500)
		subplot(5,1,i)
		plot(A)
		ylabel("k = $i")
	end
	gcf()
end

# ╔═╡ 6c3dbd94-672d-11eb-08c9-5984c9ebf4b8
md"""
Some trejectories, they seem to be in the right range. 
"""

# ╔═╡ 75242506-6717-11eb-3db7-bddd19087912
begin
	close("all")
	for i = 1:5
		subplot(5,1,i)
		plot(vv[i+1,1:900])
		ylabel("k = $i")
	end
	gcf()
end

# ╔═╡ 49e72792-6720-11eb-1cf5-4971756583a4
md"""
Here are the emperical marigan distributions.
"""

# ╔═╡ a7b9b4fa-6731-11eb-2aed-0972149d0628
md"""
# Mean Energy

Now, we look at the mean energy spectrum, observe that there is a discrepency here. A while back you told me that you switched scales at one point sort of mid paper, So this may just be the old scaleing. 
"""

# ╔═╡ 7d9542f2-6731-11eb-2465-01c9ddb0a714
nu = 1/sqrt(0.085)

# ╔═╡ ee188468-671c-11eb-05ea-9f2f29dd42f0
begin
	close("all")
	E = mean(abs2.(vv), dims = 2)[:]
	semilogy((1:18)/nu ,E[2:19])
	xlabel("K/ν")
	ylabel("⟨|vₖ|²⟩")
	title("Mean energy")
	gcf()
end

# ╔═╡ Cell order:
# ╟─da51dedc-670f-11eb-3b0b-e9a19ce28437
# ╠═61204a48-6710-11eb-2fb8-099c1dbe495d
# ╠═bb25ac5a-6712-11eb-276a-8dd819a3000c
# ╟─89c28b86-672d-11eb-351b-f9897c5ff59f
# ╟─3acc056c-6715-11eb-3d61-952f5674914e
# ╟─6c3dbd94-672d-11eb-08c9-5984c9ebf4b8
# ╟─75242506-6717-11eb-3db7-bddd19087912
# ╟─49e72792-6720-11eb-1cf5-4971756583a4
# ╟─a7b9b4fa-6731-11eb-2aed-0972149d0628
# ╠═7d9542f2-6731-11eb-2465-01c9ddb0a714
# ╟─ee188468-671c-11eb-05ea-9f2f29dd42f0
