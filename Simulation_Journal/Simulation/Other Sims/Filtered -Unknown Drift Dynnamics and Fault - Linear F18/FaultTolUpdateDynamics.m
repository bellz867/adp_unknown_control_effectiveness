function [xNew,mu] = FaultTolUpdateDynamics(dt,time,x,u)

Wf = [-1.175,0.9871 ;-8.458,-0.8776]'; %Drift dynamics actual weight
phif = [x(1);x(2)]; %drift dynamic basis
Wg = [-0.194,-0.03593;
       -19.29,-3.803]'; %control effectiveness actual weight
phig = eye(2); %control effectiveness basis

g = Wg'*phig;


mu = [0.7,0]; %actual fault

dx = Wf'*phif+g*(diag(mu))*u; %dynamics
% dx = Wf'*phif+Wg'*phig*(u-mu); %dynamics
xNew = x+dt*dx; %Update state

    
    
end


