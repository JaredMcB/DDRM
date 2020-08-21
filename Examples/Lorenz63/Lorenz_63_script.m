% Lorezn_63_script.m
obj = Lorenz_63;
obj = Load_model(obj);

obj = FE(obj);

obj = run(obj);
obj = get_obs(obj);
obj = get_PSI(obj);

PSI = @(x) [x, obj.F(x)];
[x,a,b] = gen_mod_reduc(obj.Z,PSI);

subplot(3,1,3);
hold on; plot(obj.t,x);
%ylim([-20,20])
