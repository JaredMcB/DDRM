function my_hist(series, bin_num)
    series = reshape(series,length(series))
    sort!(series)
    bins = (series[end] - series[1])/bin_num.*(0:bin_num) .+ series[1]
    cu_hist = zeros(bin_num+1)
    for i = 1:bin_num
        cu_hist[i+1] = sum(series .< bins[i+1])
    end
    hist = zeros(bin_num)
    bin = zeros(bin_num)
    for i=1:bin_num
        hist[i] = cu_hist[i+1] - cu_hist[i]
        bin[i] = (bins[i] + bins[i+1])/2
    end
    plot(bin,hist,t=[:bar])
end
