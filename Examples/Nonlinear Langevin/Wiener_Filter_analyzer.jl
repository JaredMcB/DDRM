using StatsPlots
using DataFrames

function analyse_h_ens(h_wf_ens; plt = true)

    h_wf_mean = mean(h_wf_ens,dims = 4)[:,:,:,1]
    h_wf_var = var(h_wf_ens,dims = 4)[:,:,:,1]

    if plt
        T = h_wf_ens[1,1,:,:]'
        T[:,1] .-= 1

        S = h_wf_ens[1,2,:,:]'

        dfT, dfS = DataFrame(T), DataFrame(S)
        p1 = @df dfT boxplot(T,
            marker=(0.3,:orange,stroke(.5)),
            alpha=0.75,
            leg = :none,
            title = "First coefficients")
        p2 = @df dfS boxplot(S,
            marker=(0.3,:orange,stroke(.5)),
            alpha=0.75,
            leg = :none,
            title = "Second coefficients")

        P = plot(p1,p2, layout = (2,1))
    end

    plt == true ?  [h_wf_mean, h_wf_var, P] :  [h_wf_mean, h_wf_var]
end
