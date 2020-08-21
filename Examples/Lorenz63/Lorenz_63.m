% Lorenz '63 Model
classdef Lorenz_63
    
    properties
        % Model parameters
        sig = 10;
        bet = 8/3;
        rho = 28;
        
        % Run parameters
        start = 0;
        stop = 10;
        steps = 10^3;
        t;
        h;
        x_int = [0;0;0];
        x;
        
        % Vector feild
        F;
        
        % Model
        mod_F;
        
        % Observations
        obs_ind = [1];
        dr;
        start_ob = 0;
        stop_ob = 100;
        steps_ob = 10^4;
        t_ob;
        Z;
        
        PSI;
        
    end
    
    methods
        
        function obj = Load_model(obj)
            obj.t = linspace(obj.start,obj.stop,obj.steps);
            obj.h = obj.t(2) - obj.t(1); 
            obj.F = @(x) [obj.sig*(x(2) - x(1));...
                          x(1)*(obj.rho - x(3)) - x(2);...
                          x(1)*x(2) - obj.bet*x(3)];
        end
        
        function obj = FE(obj)
            obj.mod_F = @(x) x + obj.h*obj.F(x);
            
        end

        function obj = RK4(obj)
            obj.mod_F = @(x) mod_F_gen(x,obj.F,obj.h);
            
            function y = mod_F_gen(x,F,h)
                Y1 = x;
                Y2 = x + h/2*F(Y1);
                Y3 = x + h/2*F(Y2);
                Y4 = x + h*F(Y3);
                y = x + h/6*F(Y1) + h/3*F(Y2) + ...
                                h/3*F(Y3) + h/6*F(Y4);
            end
        end
        
        function obj = run(obj)
            obj.x = zeros(3,obj.steps);
            obj.x(:,1) = obj.x_int;
            
            for i = 2 : obj.steps
                obj.x(:,i) = obj.mod_F(obj.x(:,i - 1))+...
                    mvnrnd([0,0,0],.001*eye(3))';
            end
            subplot(3,1,1);
            hold on; plot(obj.t,obj.x);
        end
        
        function obj = get_obs(obj)
            obj.Z = obj.x(1:3,:);
            subplot(3,1,2);
            hold on; plot(obj.t,obj.Z);
        end
        
        function obj = get_PSI(obj)
            obj.PSI = @(z) [obj.h*obj.sig*z(2) + (1-obj.h*obj.sig)*z(1),0;...
                       obj.h*obj.rho*z(1) + (1- obj.h)*z(2),sin(z(1))];
        end
        
    end
end
