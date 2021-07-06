function [varargout] = FaultTolADPContinuous(time,x0,xHat0,ThetaHat0,...
    GammaTheta0,WcHat0,WaHat0,Gamma0,auxdata)

% Allocated Memory
x = zeros(length(x0),length(time));
xHat = zeros(length(xHat0),length(time));
ThetaHat = zeros(size(ThetaHat0,1),size(ThetaHat0,1),length(time));
GammaTheta = zeros(size(GammaTheta0,1),size(GammaTheta0,1),length(time));
muHat = zeros(auxdata.nodes,length(time));
WcHat = zeros(auxdata.nodes,length(time));
WaHat = zeros(auxdata.nodes,length(time));
Gamma = zeros(auxdata.nodes,auxdata.nodes,length(time));
uOUT = zeros(length(auxdata.Rvec),length(time));
WfHat = zeros(auxdata.pf,auxdata.dim,length(time));
WgHat = zeros(auxdata.pf+1:pf+1+pg,auxdata.dim,length(time));
WmuHat = zeros(auxdata.pf+2+auxdata.pg:auxdata.pf+2+auxdata.pg+...
    auxdata.pg*length(auxdata.Rvec),auxdata.dim,length(time));
PHI = zeros(auxdata.pf+auxdata.pg+auxdata.pg*length(auxdata.Rvec),length(time));


Gamma_MAT = zeros(length(WcHat0),length(WcHat0),length(time));


% Initialize first entry
x(:,1) = x0;
xHat(:,1) = xHat0;
ThetaHat(:,:,1) = ThetaHat0;
GammaTheta(:,:,1) = GammaTheta0; 
muHat(:,1) = muHat0;
WcHat(:,1) = WcHat0;
WaHat(:,1) = WaHat0;
Gamma(:,:,1) = Gamma0;
WfHat(:,:,1) = ThetaHat(1:auxdata.pf,1:auxdata.dim,1)
WgHat(:,:,1) = ThetaHat(auxdata.pf+1:auxdata.pf+1+auxdata.pg,1:auxdata.dim,1)
WmuHat(:,:,1) = ThetaHat(auxdata.pf+2+auxdata.pg:end,1:auxdata.dim,1);
PHI(:,1)

%Get initial muHat
[I,J] = (find(WgHat(:,:,1)~=0));
for k=1:auxdata.Rvec
    muHat(k,1) = min(max(1/WgHat(I(1),J(1),1)*WmuHat((k-1)*auxdata.pg+I(1),J(1),1),0),1);
end
    

%Get time steps
h = diff(time);

for i=1:length(time)-1
    %Get timestep
    dt = h(i);  
    

    %Get Basis 
    %Get basis functions for estimator
    [phiF,phiG] = FaultTolGetBasisDyn(x(:,i));    
    %Get StaF basis functions and their gradients
    [sig,sigPrime] = FaultTolGetBasisADP(x(:,i),auxdata);
    
    %Get Approximate Drift dynamics and G
    GHat = WgHat(:,:,i)*phiG*(eye(length(auxdata.Rvec))-diag(muHat(:,i)))
    FHat = WfHat(:,:,i)*phiF;
    
    %Get input input cost and musat(sigPrime(sgn()-Tanh())) term for WaHat
    [ u, Wu,GD ] = FaultTolGetInput(GHat,WaHat,sigPrime,auxdata);
    
    PHI(:,i) = [phiF;
           phiG*u;
           reshape(phiG*diag(u),auxdata.pg*length(auxdata.Rvec),1)];
    
    
    %Approximate Dynamics and Regressor
    FGHat = FHat+GHat*u;
    omega = sigPrime*FGHat;
    Inst_Cost = x(:,i)'*auxdata.Qx*x(:,i)+Wu;
    delta = omega'*Wc_hat(:,i)+Inst_Cost;
    rho = 1+auxdata.nu*(omega'*omega);
    
    
    
    
    %-------------------------------------------------------------------------%
    %------------------ Build Bellman Extrapolation Library ------------------%
    %-------------------------------------------------------------------------%



    nmx = (x(:,i)'*x(:,i)+auxdata.dbar)/(1+auxdata.scale2*x(:,i)'*x(:,i));
    % nmx=1;
    pts = unifrnd(-auxdata.scale*nmx,auxdata.scale*nmx,params.dim,auxdata.numpoints);
    % pts = auxdata.scale*nmx*randn(params.dim,auxdata.numpoints);
    samples = zeros(params.dim,auxdata.numpoints);
     
    for k = 1:auxdata.numpoints
        samples(:,k) = pts(:,k)+x(:,i);
        
        %Get extrapolated Basis 
        [phiF_k,phiG_k] = FaultTolGetBasisDyn(samples(:,k));    
        %Get StaF basis functions and their gradients
        [~,sigPrime_k] = FaultTolGetBasisADP(samples(:,k),auxdata);

        %Get Approximate Drift dynamics and G
        GHat_k = WgHat(:,:,i)*phiG_k*(eye(length(auxdata.Rvec))-diag(muHat(:,i)))
        FHat_k = WfHat(:,:,i)*phiF_k;

        %Get input input cost and musat(sigPrime(sgn()-Tanh())) term for WaHat
        [ u_k, Wu_k,GD_k ] = FaultTolGetInput(GHat_k,WaHat,sigPrime_k,auxdata);

         %Approximate Dynamics and Regressor
        omega_k = sigPrime*(FHat_k+GHat_k*u_k);
        delta_k = omega'*Wc_hat(:,i)+samples(:,k)'*auxdata.Qx*samples(:,k)+Wu_k;
        rho_k = 1+auxdata.nu*(omega_k'*omega_k);

        critic_sum = critic_sum + omega_k*delta_k/rho_k;
        actor_sum = actor_sum + sigPrime_k*GD_k*omega_k'/rho_k;
        Gamma_sum = Gamma_sum + (omega_i*omega_i')/(rho_i^2);
    end

   Gamma_MAT(:,:,i) = Gamma_sum; 

    
    xTilde = x(:,i)-xHat(:,i);
    
    
    
    
    
   
     
    
    
    
    
    
    
    %% Integrate Dynamics 
    [F,G] = FaultTolGetDynamics(x(:,i));
    FG = F+G(eye(length(auxdata.Rvec))-diag(mu))*u;
    x(:,i+1) = x(:,i)+dt*FG;
    
    
    uOUT(:,i) = u;
    
    %%
    % Get estimated Dynamics
    WfHat(:,:,i+1) = ThetaHat(1:auxdata.pf,1:auxdata.dim,i+1);
 nn   WgHat(:,:,i+1) = ThetaHat(auxdata.pf+1:pf+1+pg,1:auxdata.dim,i+1);
    WmuHat(:,:,i+1) = ThetaHat(pf+2+pg:end,1:auxdata.dim,i+1);
   
    %%
    %%Estimate muHat
    [I,J] = (find(WgHat(:,:,i+1)~=0));
    if isempty([I,J])==1
        muHat(:,i+1) = muHat(:,i);
    else
        for k = 1:auxdata.Rvec
            muHat(k+1,i+1) = min(max(1/WgHat(I(1),J(1),i+1)*WmuHat((k-1)*auxdata.pg+I(1),J(1),i+1),0),1);
        end
    end

 
end

xOUT = x;
xHatOUT = xHat;
ThetaHatOUT = reshape(ThetaHat);
GammaThetaOUT = rehsape();
muHatOUT = muHat;
WcHatOUT = WcHat;
WaHatOUT = WaHat;
GammaOUT = reshape(Gamma,auxdata.nodes.^2,length(time));
Gamma_MATOUT = reshape(Gamma_MAT,auxdata.nodes.^2,length(time));
%Calculate last
uOUT(:,length(time)) = getinput(

    
end
    
    
    