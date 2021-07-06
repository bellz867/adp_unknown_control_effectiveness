clear;
clc;
clear;
% cd('C:\Users\pdeptula\Desktop\Repo stuff\[ACC 2018] StaF ADP _ Unknown Drift and Control Effectiveness Fault\Simulation\Filtered -Unknown Drift Dynnamics and Fault')
%-------------------------------------------------------------------------%
%--------------------------Setup------------------------------------------%
%-------------------------------------------------------------------------%
plotSoln = 1;
savePlotsAndVars =0;
includeTitle = 0;
savePlotsInPdf = 0;
savePlotsInJpeg = 0;
savePlotsInPng = 0;
savePlotsInEps = 0;
includeComp = 0;
compareToComplementaryPaper = 0;

%Set Times
auxdata.finalTime = 15; %Simulation end time in seconds
auxdata.dTime = 1/1000; %Simulatons time step size
time = 0:auxdata.dTime:auxdata.finalTime; %simulation time
tLength = length(time);
avgtoc = 0;
numtocs = 0;
auxdata.Dt = 0.01; %Integration window

%-------------------------------------------------------------------------%
%------------------------Initial Conditions-------------------------------%
%-------------------------------------------------------------------------%
% Agent Initial Conditions
auxdata.x0 = [-1 1]';

n = length(auxdata.x0); %Number of State Dimensions
m = 1; %Number of Inputs
auxdata.NormDyn = norm(auxdata.x0); %Agent dynamics normalizations constant
auxdata.sysID = 1; %Sys identificiation (model free = 1) or model-based (0)
auxdata.Noisy = 0; %Add noise (=1) or No noise (=0);
auxdata.NoiseStd = 0.001; %Noise variance;

%-------------------------------------------------------------------------%
%--------------------------Agent ADP Parameters---------------------------%
%-------------------------------------------------------------------------%
%Number of basis functions
auxdata.nodes = n+1; %Number of basis functions
auxdata.scale = 1; %Normalization gain for ADP basis
auxdata.dbar = 0.00001; %Offset for ADP basis
auxdata.CenterOffsetsScale = 0.7; %Scaling for basis center offsets
auxdata.CenterOffsets =  SimplexVert(n); %Get StaF center offsets using vertices of n-simplex
auxdata.BasisType = 2; %Basis Type 1 for polynomial 2 for exponential
auxdata.ShrinkingCenterOffset = 1; %Shrinking basis offset or no
auxdata.numpoints = 1; %Number of BE extrapolation points
auxdata.MovePoints = 1; %move points via frequency (=1) or keep grid around current state (=0)
auxdata.BEscale = 0.5; %Scale for BE extraplation
auxdata.scale2 = 0.5; %Normalization gain for BE extrapolation
auxdata.dbar2 = 0.0001; %Offset for BE  extrapolation
auxdata.freq = unifrnd(500,1000,n,100); %frequency for BE extrapolation
auxdata.phase = unifrnd(0,2*pi,n,100); %phase for BE extrapolation

%ADP Gains
auxdata.kc1 = 0.001; %Critic gain 1
auxdata.kc2 = 0.25; %Critic gain 2
auxdata.ka1 = 0.5; %Actor gain 1
auxdata.ka2 = 0.02; %Actor gain 2
auxdata.nu1 = 0.05; %Normalization gain
auxdata.nu2 = auxdata.nu1; %Normalization gain fo BE extrapolation
auxdata.beta = 0.005; %Least Squares forgetting factor
auxdata.WaHatLower = -1000;
auxdata.WaHatUpper = 1000;
auxdata.Ka = 1*eye(auxdata.nodes);

if compareToComplementaryPaper == 0;
    auxdata.InputKp1 = 0.005*eye(2).*(auxdata.NormDyn); %Gain for input if adding Kx inside saturation
else
    auxdata.InputKp1 = 0*eye(2); %Gain for input if adding Kx inside saturation
    auxdata.KComp = [-0.5 -1];
end
auxdata.InputKp2 = 2;

% Weight Initial Conditions
auxdata.WcHat0 =  [1.557733614848178;1.344189320140740;1.586443459358885]'; %unifrnd(2,3,auxdata.nodes,1);%[0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4]'; ;[3.4, 3.2, 2.8]'; [3.2, 5.5, 3.1]';% [5.05, 5.9 3.48]';%
auxdata.WaHat0 = 0.7*auxdata.WcHat0;
auxdata.GammaADP0 = 10*eye(auxdata.nodes);

%Cost Parameters
auxdata.Qx = 5*eye(n);
auxdata.Rvec = 1*ones(1,m);
auxdata.R = diag(auxdata.Rvec);

%Actuator Saturation
auxdata.sat = 3;

%-------------------------------------------------------------------------%
%--------------------------Agent SYS ID Parameters------------------------%
%-------------------------------------------------------------------------%
auxdata.N = 25; %Stack Size
auxdata.pf = 3; %Length of drift dynamics weight;
auxdata.pGamma = n*auxdata.pf+m; %Dimension of CL Gain matrix
auxdata.L = 1; %Gain for Least Square Matrices P and Q
auxdata.kf = 1/10; % Gain for low pass filter
auxdata.kP1 = 0.0001; %Gain for Least Square Matrices P and Q
auxdata.kCL = 10; %CL gain
auxdata.betaCL = 10; %CL forgetting factor
auxdata.thetaLower = -10; %thetaHat lower bound
auxdata.thetaUpper = 10; %thetaHat upper bound
auxdata.CLNorm = 0;  %Update law normalization constant

auxdata.xf0 = auxdata.x0; %initial state estimate
auxdata.PHIf0 = zeros(n,auxdata.pGamma,1); %initial filtered basis
auxdata.Gf0 = zeros(n,1); %initial filtered G*u vector
auxdata.gf0 = zeros(n,m); %initial filtered control effectiveness matrix
auxdata.thetaHat0 = [-3.992912253745918;-4.852633581859568;-3.414365220772520;-1.607234481711097;3.016673527544185;0.901079981158318;3.453683303378552]'; %unifrnd(-5,5,auxdata.pGamma,1); %Initial sys id weight estimate 
auxdata.GammaCL0 = 50*eye(auxdata.pGamma); %Initial CL gain matrix
auxdata.P0 = zeros(auxdata.pGamma,auxdata.pGamma);
auxdata.Q0 = zeros(auxdata.pGamma,1);
stackIndex = 1; %initialize stack index


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Solve Using ADP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocated Memory
x = zeros(n,tLength); %agent state
PHI = zeros(n,auxdata.pGamma,tLength); %Regressor matrix
xf = zeros(n,tLength); %filtered state
Phif = zeros(n,auxdata.pGamma,tLength); %Filtered Basis 
Ff = zeros(n,tLength); % Filtered G*u
gf = zeros(n,m,tLength); %Filtered control effectiveness
mu = zeros(m,tLength);

% observer regressors
bufferSize = auxdata.Dt/auxdata.dTime;
scriptTimeBuffer = zeros(1,bufferSize);%time buffer for integrating
scriptYBuffer = zeros(n,auxdata.pGamma,bufferSize);%script y integration buffer
scriptFBuffer = zeros(n,1,bufferSize); %script G integration buffer
scriptYTYBuffer = zeros(auxdata.pGamma,auxdata.pGamma,bufferSize);%script y'y integration buffer
xBuffer = zeros(n,1,bufferSize);% delta x integration buffer
scriptYTYStack = zeros(auxdata.pGamma,auxdata.pGamma,auxdata.N);%stack of saved script y'y
scriptYTxfStack = zeros(auxdata.pGamma,1,auxdata.N);%stack of saved script y'[x(t)-x(t-T)-y*thetaHat]
minEigHistCL = zeros(1,tLength); %min eigenvalue for CL
minEigHistP = zeros(1,tLength); %min eigenvalue for P matrix

thetaHat = zeros(auxdata.pGamma,1,tLength); %estimate stack
GammaCL = zeros(auxdata.pGamma,auxdata.pGamma,tLength);%estimator update gain
P = zeros(auxdata.pGamma,auxdata.pGamma,tLength);
Q = zeros(auxdata.pGamma,1,tLength);
muHat = zeros(m,tLength);  %fault estimate
WfHat = zeros(n*auxdata.pf,tLength); %drift dynamics weight estimate


% ADP values
WcHat = zeros(auxdata.nodes,tLength); %critic weights
WaHat = zeros(auxdata.nodes,tLength); %actor weights
GammaADP = zeros(auxdata.nodes,auxdata.nodes,tLength); %Least squares adaptive gain   
ADPStackSum = zeros(auxdata.nodes,auxdata.nodes,tLength); %ADP CL stack
Cost = zeros(1,tLength); %total cost
CostBuffer = zeros(1,1,2);
integrationTimeBuffer = zeros(1,2);
delta = zeros(1,tLength); %Bellman Error
VHat = zeros(1,tLength); %Approximate Value function
minEigHistADP = zeros(1,tLength); %min eigenvalue for CL


% Initialize 
x(:,1) = auxdata.x0; %initialize state
xf(:,1) = auxdata.xf0; %initialize filtered state
PHIf(:,:,1) = auxdata.PHIf0; %intialize filtered basis;
Ff(:,1) = auxdata.Gf0; %intialize filtered G*u;
gf(:,:,1) = auxdata.gf0; %initialize filtered control effectiveness
thetaHat(:,:,1) = auxdata.thetaHat0; %initialize sys weight estimate
GammaCL(:,:,1) = auxdata.GammaCL0; %initialize CL learning gain matrix
P(:,:,1) = auxdata.P0; %initialize LS P matrix
Q(:,:,1) = auxdata.Q0; %initialze Q matrix
WcHat(:,1) = auxdata.WcHat0; %initialize critic weight
WaHat(:,1) = auxdata.WaHat0; %initialize actor weight
GammaADP(:,:,1) = auxdata.GammaADP0; %initialzie ADP learning gains
WfHat(:,1) = thetaHat(1:n*auxdata.pf,:,1); %initialize drift dynamics weight estimate
muHat(:,1) = thetaHat(n*auxdata.pf+1:auxdata.pGamma,:,1); %initialize control effective + fault weight estimate



%Get time steps
h = diff(time);
% h = waitbar(0,['time ',num2str(time(1)),', timeF ',num2str(finalTime),', minEig ',num2str(minEigHist(1)),', loopT ', num2str(0)]);

for t=1:length(time)-1
    tic
    dt = h(t); %Get timestep 
    if (auxdata.Noisy == 1)
        x(:,t) = x(:,t)+auxdata.NoiseStd*randn(n,1); %Add noise to state to see how system reacts
    end
    
        
    if (auxdata.sysID == 0)
        Wf = [-1,-0.5 ;1,0;0,-0.5];
        
        thetaHat(:,:,t) = [reshape(Wf',n*auxdata.pf,1);mu(:,t)];
        WfHat(:,t) = thetaHat(1:n*auxdata.pf,:,t); 
        muHat(:,t) = mu(:,t);
    end
    
   
    %Get Basis 
    [phif,g,fo] = FaultTolGetBasisDyn(x(:,t)); %Get basis functions for estimator   
    [sig,sigPrime] = FaultTolGetBasisADP(x(:,t),auxdata); %Get StaF basis functions and their gradients
    
    %Get Approximate Drift dynamics and G
    gHat = g*(diag(muHat(:,t)));

    [u1, uCost,Ga,Gc] = FaultTolGetInput(x(:,t),gHat,WaHat(:,t),sigPrime,auxdata); %Get input input cost and musat(sigPrime(sgn()-Tanh())) term for WaHat
    if compareToComplementaryPaper == 1;
        u0 = auxdata.KComp*x(:,t);
        u =u0+u1;
    else
        u0 = 0;
        u = u1;
    end
    PHI(:,:,t) = [kron(phif',eye(n)),g*diag(u)]; %Get Combined CL basis
    G = 0; %Get combined g*u
    

    %ADP
    [fgHat,omega,instCost,delta(:,t),rho] = FaultTolGetADPVar(x(:,t),fo,sigPrime,uCost,WcHat(:,t),thetaHat(:,:,t),PHI(:,:,t),auxdata,1);
    [criticSum,actorSum,ADPGammaSum,ADPStackSum(:,:,t)] = FaultTolGetExtrapBEVarStack(time(t),x(:,t),thetaHat(:,:,t),WcHat(:,t),WaHat(:,t),GammaADP(:,:,t),muHat(:,t),auxdata);
    [WcHat(:,t+1),WaHat(:,t+1),GammaADP(:,:,t+1)] = FaultTolUpdateADP(WcHat(:,t),WaHat(:,t),GammaADP(:,:,t),Ga,Gc,criticSum,actorSum,ADPGammaSum,omega,delta(:,t),rho,auxdata,dt);
    
    if (auxdata.sysID == 1)
        
        %ICL
        scriptTimeBuffer = [scriptTimeBuffer(2:end),time(t)]; %update integral buffer time]
        scriptYBuffer = cat(3,scriptYBuffer(:,:,2:end),PHI(:,:,t)); %update script Y Buffer
        scriptFBuffer = cat(3,scriptFBuffer(:,:,2:end),fo); %update script G Buffer
        xBuffer = cat(3,xBuffer(:,:,2:end),x(:,t)); %update delta x Buffer
        scriptYTYStackSum = zeros(size(scriptYBuffer,2),size(scriptYBuffer,2));
        scriptYTxStackSum = zeros(size(scriptYBuffer,2),size(x(:,t),2));
        %Update parameter estimates
        if time(t)>=auxdata.Dt
            deltax = xBuffer(:,:,end)-xBuffer(:,:,1);
            scriptY = FaultToIIntegrateBuffer(scriptTimeBuffer,scriptYBuffer);
            scriptG = FaultToIIntegrateBuffer(scriptTimeBuffer,scriptFBuffer);
            [scriptYTYStack,scriptYTxfStack,minEig,scriptYTYStackSum,scriptYTxStackSum,stackIndex] = FaultTolUpdateStacks(scriptYTYStack,scriptYTxfStack,scriptY,scriptG,deltax,thetaHat(:,:,t),stackIndex,auxdata); %Update History Stacks
            minEigHistCL(t) = minEig;      
        end

%
        %Estimator
        [xf(:,t+1),PHIf(:,:,t+1),Ff(:,t),gf(:,:,t+1),P(:,:,t+1),Q(:,:,t+1),thetaHat(:,:,t+1),GammaCL(:,:,t+1)] = FaultTolUpdateEstimator(x(:,t),xf(:,t),PHIf(:,:,t),Ff(:,t),gf(:,:,t),P(:,:,t),Q(:,:,t),thetaHat(:,:,t),GammaCL(:,:,t),PHI(:,:,t),G,g,scriptYTYStackSum,scriptYTxStackSum,auxdata,dt);
        minEigHistP(:,t+1) = min(eig(P(:,:,t+1)));
    else
        thetaHat(:,:,t+1) = thetaHat(:,:,t);
    end

    %
    %Get New Estimates for estimate matrix
    WfHat(:,t+1) = thetaHat(1:n*auxdata.pf,:,t+1); %initialize drift dynamics weight estimate
%     WmuHat(:,:,t+1) = thetaHat(auxdata.pf+1:auxdata.pf+n*m,:,t+1); %initialize control effective + fault weight estimate
    muHat(:,t+1)= thetaHat(n*auxdata.pf+1:auxdata.pGamma,:,t+1);


    
  
  
    
    
    % Integrate Dynamics 

    [x(:,t+1),mu(:,t+1)] = FaultTolUpdateDynamics(dt,time(t),x(:,t),u);
    Cost(t+1) = Cost(t)+dt*instCost;    
       
    p = min(isfinite(x(:,t)))&&min(isfinite(WcHat(:,t)))&&min(min(isfinite(ADPStackSum(:,:,t))));
    if p ~= 1
        numtocs = numtocs+1;
        avgtoc = avgtoc+toc;
        disp('A variable is not finite!!');
        break 
    end
    minEigHistADP(:,t) = min(eig(ADPStackSum(:,:,t)));
    CostBuffer = cat(3,CostBuffer(:,:,2:end),instCost);
    integrationTimeBuffer = [integrationTimeBuffer(2:end),time(t)];
    costInt = FaultToIIntegrateBuffer(integrationTimeBuffer,CostBuffer);
    Cost(:,t+1) = Cost(:,t)+costInt; %Total Cost
    

    avgtoc = avgtoc+toc;
    numtocs = numtocs+1;
%     waitbar(time(t)/time(end),h,['time ', num2str(time(t)), ', timeF ', num2str(finalTime), ', minEig ', num2str(minEigHist(t)), ', loopT ', num2str(avgtoc/numtocs)]);
 
end
% close(h)
avgtoc = avgtoc/numtocs


if plotSoln==1
    close all
    FaultTolPlot;
end




