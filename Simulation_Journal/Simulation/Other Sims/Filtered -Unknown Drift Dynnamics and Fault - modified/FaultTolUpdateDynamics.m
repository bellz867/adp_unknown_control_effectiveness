function [xNew,mu] = FaultTolUpdateDynamics(dt,time,x,u)

Wf = [-1,-0.5 ;1,0;0,-0.5]; %Drift dynamics actual weight
Wg = [0,1;0,2]; %control effectiveness actual weight

phif = [x(1);x(2);x(2)*(1-(cos(2*x(1))+2)^2)];
phig = [cos(2*x(1));1];



mu = [0.4]; %Actual fault


dx = Wf'*phif+Wg'*phig*(diag(mu))*u; %dynamics
% dx = Wf'*phif+Wg'*phig*(u-mu); %dynamics
xNew = x+dt*dx; %Update state

    
    
end


