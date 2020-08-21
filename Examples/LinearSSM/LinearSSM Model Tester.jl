F = [0.1 10; 0 -.5] # Must be stable
G = [1 0; 0 1]
R = [1 0; 0 1]
Xo = [1; 1]
t_disc = 1000
gap = 1

t_start = 0
t_stop = 1e6
h = 1

X = modgen_LSSM(t_start,t_stop,h,
    F = F,
    G = G,
    R = R,
    Xo = Xo,
    t_disc = t_disc,
    gap = gap)

spect_ana(z) = ((I - z^(-1)*F)\R)/(I - z*F')
