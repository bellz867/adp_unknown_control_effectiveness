%% create directory to save plots and conditions
if (savePlotsAndVars==1)
    MainDir = pwd;
    cd('Plots and Conditions')
    
    foldername = inputdlg('Please enter output folder name:');
    dt = num2str(date);
    newSubFolder = sprintf('%s__%s', foldername{1}, dt);
    % Finally, create the folder if it doesn't exist already.
    if ~exist(newSubFolder, 'dir')
        mkdir(newSubFolder);
    else
        foldername2 = inputdlg('Please enter a different folder name:');
        newSubFolder2 = sprintf('%s__%s', foldername2{1}, dt);
        mkdir( newSubFolder2);
    end
    cd (MainDir)
    % Get current directory
    
    % Get desired directory
    FolderDest = sprintf('%s/Plots and Conditions/%s__%s', pwd, foldername{1}, dt);
end


%%

if p~=1
    tend = t-1;
else
    tend = tLength;
end

VHat = zeros(1,tend);
u0  = zeros(m,tend);  %complementary robust controller paper
u1 = zeros(m,tend);   %approximate saturated controller
u  = zeros(m,tend);   %total controller
gHat = zeros(n,m,tend);
gAct = zeros(n,m,tend);

S = zeros(1,tend);
for t=1:tend
    [phif,go] = FaultTolGetBasisDyn(x(:,t)); %Get basis functions for estimator   
    [sig,sigPrime] = FaultTolGetBasisADP(x(:,t)./auxdata.NormDyn,auxdata); %Get StaF basis functions and their gradients
    %Get Approximate Drift dynamics and G
    gHat(:,:,t) = go*(diag(muHat(:,t)));
     gAct(:,:,t) = go*(diag(mu(:,t)));
     y = x(:,t)./auxdata.NormDyn;
    S(:,t) = y'*auxdata.InputKp1*y/(1+auxdata.InputKp2*(y'*y))^2;
    VHat(t) = WcHat(:,t)'*sig+S(:,t);
    [u1(:,t), uCost,~,~,~] = FaultTolGetInput(x(:,t),gHat(:,:,t),WaHat(:,t),sigPrime,auxdata);
    if compareToComplementaryPaper == 1;
        u0(:,t) = auxdata.KComp*x(:,t);
        u(:,t) =u0(:,t)+u1(:,t);
    else
        u0(:,t) = 0;
        u(:,t) = u1(:,t);
    end
    
end

Time = time(1:tend);

fault = mu(:,end);

VStar = 1/2*x(1,:).^2+x(2,:).^2;
uStar = -fault*(cos(2*x(1,:))+2).*x(2,:);

TF = find(Time==15);


%%
f1 = figure(1);
box on
p = plot(Time(1:TF),x(:,1:TF),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);

if includeTitle==1
    title('States','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('x(t)','interpreter','latex');
l  = legend('$x_{1}(t)$','$x_{2}(t)$');
set(l,'interpreter','latex')

set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);



%%
f2 = figure(2);
box on
if includeComp==1
    p = plot(Time(1:TF),u(:,1:TF),Time(1:TF),uStar(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.5]);
    set(p(2),'Color',[0,0.5,0]);
    l  = legend('$\hat{u}(t)$: Constrained Approximate Input','$u^*(t)$: Unconstrained Optimal Input');
    set(l,'interpreter','latex','location','southeast');
    ylabel('$u(t)$','interpreter','latex');
else
    p = plot(Time(1:TF),u(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.5]);
    ylabel('$u(t)$','interpreter','latex');
end

if includeTitle==1
    title('Approximate Optimal Input','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');

ylim([min(min(u(:,1:tend)))-0.25,max(max(u(:,1:tend)))+0.25])



set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);


%%
f3 = figure(3);
box on
p = plot(Time(1:TF),WcHat(:,1:TF),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);
set(p(3),'Color',[0.5,0,0]);

if includeTitle==1
    title('Critic Weight Estimates','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{W}_{c}(t)$','interpreter','latex');
ylim([min(min(WcHat(:,1:tend)))-0.25,max(max(WcHat(:,1:tend)))+0.25])
l  = legend('$\hat{W}_{c1}(t)$','$\hat{W}_{c2}(t)$','$\hat{W}_{c3}(t)$');
set(l,'interpreter','latex');

set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);



%%
f4 = figure(4);
box on
p = plot(Time(1:TF),WaHat(:,1:TF),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);
set(p(3),'Color',[0.5,0,0]);

if includeTitle==1
    title('Actor Weight Estimates','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{W}_{a}(t)$','interpreter','latex');
ylim([min(min(WaHat(:,1:TF)))-0.25,max(max(WaHat(:,1:TF)))+0.25])
l  = legend('$\hat{W}_{a1}(t)$','$\hat{W}_{a2}(t)$','$\hat{W}_{a3}(t)$');
set(l,'interpreter','latex');

set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);


%%
f5 = figure(5);
box on
Wf =repmat(reshape([-1,-0.5 ;1,0;0,-0.5]',n*auxdata.pf,1),1,TF);
WfTilde = Wf(:,1:TF)-WfHat(:,1:TF);
p = plot(Time(1:TF),WfTilde,'LineWidth',1.25);

if includeTitle==1
    title('Drift Weight Estimate Error','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\tilde{W}_{1}(t)$','interpreter','latex');
ylim([min(min(WfTilde(:,1:TF)))-0.25,max(max(WfTilde(:,1:TF)))+0.25])
set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

%%
f6 = figure(6);
box on
p = plot(Time(1:TF),WfHat(:,1:TF),'LineWidth',1.25);

if includeTitle==1
    title('Drift Weight Estimates','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{W_{1}}(t)$','interpreter','latex');
ylim([min(min(WfHat(:,1:TF)))-0.25,max(max(WfHat(:,1:TF)))+0.25])
set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

%%
f7 = figure(7);
box on
muTilde = mu(:,1:TF)-muHat(:,1:TF);
p = plot(Time(1:TF),muTilde,'LineWidth',1.25);
set(p,'Color',[0,0,0.5]);

if includeTitle==1
    title('Control Fault Estimate Error','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\tilde{\mu}_{a}(t)$','interpreter','latex');
ylim([min(min(muTilde(:,1:TF)))-0.25,max(max(muTilde(:,1:TF)))+0.25])
set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

%%
f8 = figure(8);
box on
p = plot(Time(1:TF),muHat(:,1:TF),'LineWidth',1.25);
set(p,'Color',[0,0,0.5]);

if includeTitle==1
    title('Control Fault Estimate','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{\mu}_{a}(t)$','interpreter','latex');
ylim([min(min(muHat(:,1:TF)))-0.05,max(max(muHat(:,1:TF)))+0.05])
set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

%%
f9 = figure(9);
box on
p = plot(Time(1:TF),delta(:,1:TF),'LineWidth',1.25);
set(p,'Color',[0,0,0.5]);

if includeTitle==1
    title('Bellman Error','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\delta(t)$','interpreter','latex');

set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);


%%
f10 = figure(10);
box on
p = plot(Time(1:TF),Cost(:,1:TF),'LineWidth',1.25);
set(p,'Color',[0,0,0.5]);

if includeTitle==1
    title('Total Cost','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$J(x(t),u(t))$','interpreter','latex');

set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);


%%
f11 = figure(11);

box on
p = plot(Time(1:TF),minEigHistP(:,1:TF),Time(1:TF-1),minEigHistCL(:,1:TF-1),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);



if includeTitle==1
    title('Minimum Eigenvalues for System Identification','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\lambda_{min}$','interpreter','latex');
l = legend('$\lambda_{min} \left\{P(t)\right\}$','$\lambda_{min}\left\{\sum_{i}^{N}{\mathcal{Y}^T\mathcal{Y}}\right\}$');
set(l,'interpreter','latex');

set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6,5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

%%

%Save Plots
if (savePlotsAndVars==1)
    cd (FolderDest)
    if (savePlotsInPdf==1)
        print( f1, '-dpdf','FTStates');
        if includeComp==1
            print(f2, '-dpdf', 'FTInputComp');
        else
            print(f2, '-dpdf', 'FTInput');
        end
        print( f3, '-dpdf','FTCritEst');
        print( f4, '-dpdf','FTActEst');
        print( f5, '-dpdf','FTSysWErr');
        print( f6, '-dpdf','FTSysWEst');
        print( f7, '-dpdf','FaultErr');
        print( f8, '-dpdf','FaultEst');
        print( f9, '-dpdf','FTBE');
        print( f10, '-dpdf','FTTotCost');
        print( f11, '-dpdf','FTCLEigs');
    end
    if (savePlotsInJpeg==1)
        print( f1, '-djpeg', 'FTStates');
        if includeComp==1
            print( f2, '-djpeg', 'FTInputComp');
        else
            print( f2, '-djpeg', 'FTInput');
        end
        print( f3, '-djpeg' ,'FTCritEst')
        print( f4, '-djpeg' ,'FTActEst')
        print( f5, '-djpeg' ,'FTSysWeightErr')
        print( f6, '-djpeg' ,'FTSysWeightEst')
        print( f7, '-djpeg' ,'FaultErr')
        print( f8, '-djpeg' ,'FaultEst')
        print( f9, '-djpeg' ,'FTBE')
        print( f10, '-djpeg' ,'FTTotCost')
        print( f11, '-djpeg' ,'FTCLEigs')
    end
    if (savePlotsInPng==1)
        print( f1, '-dpng', 'FTStates');
        if includeComp==1
            print( f2, '-dpng', 'FTInputComp');
        else
            print( f2, '-dpng', 'FTInput');
        end
        print( f3, '-dpng' ,'FTCriticEst')
        print( f4, '-dpng' ,'FTActorEst')
        print( f5, '-dpng' ,'FTSysWeightErr')
        print( f6, '-dpng' ,'FTSysWeightEst')
        print( f7, '-dpng' ,'FaultErr')
        print( f8, '-dpng' ,'FaultEst')
        print( f9, '-dpng' ,'FTBE')
        print( f10, '-dpng' ,'FTTotCost')
        print( f11, '-dpng' ,'FTCLEigs')
    end
    if (savePlotsInEps==1)
        print( f1, '-deps', 'FTStates');
        print( f2, '-deps', 'FTInput');
        print( f3, '-deps' ,'FTCriticEst')
        print( f4, '-deps' ,'FTActorEst')
        print( f5, '-deps' ,'FTSysWeightErr')
        print( f6, '-deps' ,'FTSysWeightEst')
        print( f7, '-deps' ,'FTFaultEstErr')
        print( f8, '-deps' ,'FTFaultEst')
        print( f9, '-deps' ,'FTBE')
        print( f10, '-deps' ,'FTTotCost')
        print( f11, '-deps' ,'FTCLEigs')
    end
    cd(MainDir)
end 




%%
%Save Ouputs, Initial Conditions


if (savePlotsAndVars==1)
    cd(FolderDest)
    save('Conditions','auxdata')
    save('Outputs','time','tend','x','WcHat','WaHat','thetaHat','muHat','mu','WfHat','minEigHistP','minEigHistCL','u','u1','u0','delta','minEigHistADP','P','Q','GammaADP','GammaCL','Cost','TF')
    cd(MainDir)
end

