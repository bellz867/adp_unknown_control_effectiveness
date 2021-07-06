function [xNew,mu] = FaultTolUpdateDynamics(dt,time,x,u)

Wf = [-1,-0.5 ;1,0;0,-0.5]; %Drift dynamics actual weight
Wg = [0,1;0,2]; %control effectiveness actual weight

phif = [x(1);x(2);x(2)*(1-(cos(2*x(1))+2)^2)];
psig = [cos(2*x(1));1];




% mu = [0.5+0.1*cos(3*time)*(1-tanh(x(1))*tanh(x(2)))*sin(x(1))]; %Actual fault
% mu = [0.5+0.2*(1-tanh(x(1))*tanh(x(2)))*sin(x(1))]; %Actual fault
mu = [0.5+0.2*(1-max(tanh(x(1)^2),tanh(x(2))^2)^2)*sin(x(1))]; %Actual fault


% wGmu = [0.6,0.3]';
% phiG = [1,tanh(x(1)^2+x(2)^2)]';

% if norm(x)>=1
%     mu = 0.7;
% elseif norm(x)<1 abs(x(1))>0.35 & abs(x(2))>0.35 
% 
% elseif abs(x(1))<=0.35 & abs(x(2))<=0.35


% if norm(x)^2>(1)^2
%     mu = 1;
% elseif norm(x)^2>(0.5)^2 && norm(x)^2<=(1)^2
%     mu = 1-0.25*exp(-0.2*norm(x)^2/2);
% else
%     mu = exp(-0.2*norm(x)^2/2)-0.5*exp(-0.2*norm(x)^2/2);
% end
%    
    
    
    
    
% mu = 0.6*(1-0.1*tanh(10*(min([0,max([x(1)^2-(0.25)^2,x(2)^2-(0.5)^2])]))));
    

% mu = wGmu'*phiG;

% if norm(x)>0.75
%     mu = 
    


dx = Wf'*phif+Wg'*psig*(diag(mu))*u; %dynamics
% dx = Wf'*phif+Wg'*phig*(u-mu); %dynamics
xNew = x+dt*dx; %Update state

    
    
end


