% Lorenz '63 Model
classdef Lorenz_mod
    
    properties
        % Model parameters
        sig = 10;
        bet = 8/3;
        rho = 28;
        
        % Run parameters
        steps = 10^4;
        start = 0;
        stop = 100;
        t = linspace(start,stop,steps);
        h = t(2) - t(1);
        x_int = [1;1;1];
        
        % Vector feild
        F;
        
        % Model
        mod_F;
        
    end
    
    methods
        
        function obj = gen_F(obj)
            obj.F = @(x) [sig*(x(2) - x(1));...
                          x(1)*(rho - x(3)) - x(2);...
                          x(1)*x(2) - bet*x(3)];
        end
        
        function obj = RK4(obj)
            obj.mod_F = @(x) mod_F_gen(x,F,h);
            
            function y = mod_F_gen(x,F,h)
                Y1 = x;
                Y2 = x + h/2*F(Y1);
                Y3 = x + h/2*F(Y2);
                Y4 = x + h*F(Y3);
                y = x + h/6*F(Y1) + h/3*F(Y2) + ...
                                h/3*F(Y3) + h/6*F(Y4);
            end
        end
    end
end

plot3(y(1,:),y(2,:),y(3,:))