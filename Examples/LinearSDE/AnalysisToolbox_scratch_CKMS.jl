# function spectfact_matrix_CKMS(P; N_ckms = 1500)

P = R_pred_smoothed
N_ckms = 200
d = size(P)[1];
m = size(P)[3] - 1

NN = reverse(P[:,:,2:end],dims = 3)
Re = Rr = p0 = P[:,:,1]

F = sparse([[zeros(d,d*(m-1)); I] zeros(d*m,d)])
h = sparse([zeros(d,d*(m-1)) I])

K = complex(zeros(d*m,d))
for i = 0 : m-1
    K[d*i + 1: d*(i+1),:] = NN[:,:,i+1]
end
L = K

# spectfactLog = zeros(4,N_ckms)

for i = 1:N_ckms
    global L, K, Re, Rr
    hL = h*L; FL = F*L

    K_new = K - FL/Rr*hL'
    L_new = FL - K/Re*hL
    Re_new = Re - hL/Rr*hL'
    Rr_new = Rr - hL'/Re*hL

    # spectfactLog[:,i] = [cond(Rr),
    #                      cond(Re),
    #                      norm(K - K_new),
    #                      norm(L - L_new)]

    K = K_new
    L = L_new
    Re = Re_new
    Rr = Rr_new
end

k = K/Re
re = Re

sqrt_re = sqrt(re)

l = complex(zeros(d,d,m+1))
l[:,:,1] = sqrt_re;
for i = m-1:-1:0
    l[:,:,m-i+1] = k[d*i + 1: d*(i+1),:]*sqrt_re
end

# save("Data\\CKMS_dat.jld",
#     "spectfactLog",
#     spectfactLog)

l
