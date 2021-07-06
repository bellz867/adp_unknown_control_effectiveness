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
u = zeros(m,tend);
gHat = zeros(n,m,tend);
gAct = zeros(n,m,tend);
for t=1:tend
    [phif,go] = FaultTolGetBasisDyn(x(:,t)); %Get basis functions for estimator   
    [sig,sigPrime] = FaultTolGetBasisADP(x(:,t)./auxdata.NormDyn,auxdata); %Get StaF basis functions and their gradients
    %Get Approximate Drift dynamics and G
    gHat(:,:,t) = go*(diag(muHat(:,t)));
     gAct(:,:,t) = go*(diag(mu(:,t)));
    VHat(t) = WcHat(:,t)'*sig;
    [u(:,t), uCost,~ ] = FaultTolGetInput(x(:,t),gHat(:,:,t),WaHat(:,t),sigPrime,auxdata);
end

Time = time(1:tend);
   

%%
f1 = figure(1);
box on
p = plot(Time,x(:,1:tend),'LineWidth',1.25);
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
p = plot(Time,u(:,1:tend),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);

if includeTitle==1
    title('Approximate Optimal Input','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{u}(t)$','interpreter','latex');
ylim([min(min(u(:,1:tend)))-0.25,max(max(u(:,1:tend)))+0.25])
l  = legend('$\hat{u}(t)$');
set(l,'interpreter','latex');

set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);


%%
f3 = figure(3);
box on
p = plot(Time,WcHat(:,1:tend),'LineWidth',1.25);
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
p = plot(Time,WaHat(:,1:tend),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);
set(p(3),'Color',[0.5,0,0]);

if includeTitle==1
    title('Actor Weight Estimates','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{W}_{a}(t)$','interpreter','latex');
ylim([min(min(WaHat(:,1:tend)))-0.25,max(max(WaHat(:,1:tend)))+0.25])
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
Wf =repmat(reshape([-1,-0.5 ;1,0;0,-0.5]',n*auxdata.pf,1),1,tend);
WfTilde = Wf(:,1:tend)-WfHat(:,1:tend);
p = plot(Time,WfTilde,'LineWidth',1.25);

if includeTitle==1
    title('Drift Weight Estimate Error','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\tilde{W}_{1}(t)$','interpreter','latex');
ylim([min(min(WfTilde(:,1:tend)))-0.25,max(max(WfTilde(:,1:tend)))+0.25])
set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

%%
f6 = figure(6);
box on
p = plot(Time,WfHat(:,1:tend),'LineWidth',1.25);

if includeTitle==1
    title('Drift Weight Estimates','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{W_{1}}(t)$','interpreter','latex');
ylim([min(min(WfHat(:,1:tend)))-0.25,max(max(WfHat(:,1:tend)))+0.25])
set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

%%
f7 = figure(7);
box on
muTilde = mu(:,1:tend)-muHat(:,1:tend);
p = plot(Time,muTilde,'LineWidth',1.25);
set(p,'Color',[0,0,0.5]);

if includeTitle==1
    title('Control Fault Estimate Error','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\tilde{\mu}_{a}(t)$','interpreter','latex');
ylim([min(min(muTilde(:,1:tend)))-0.25,max(max(muTilde(:,1:tend)))+0.25])
set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

%%
f8 = figure(8);
box on
p = plot(Time,muHat(:,1:tend),'LineWidth',1.25);
set(p,'Color',[0,0,0.5]);

if includeTitle==1
    title('Control Fault Estimate','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{\mu}_{a}(t)$','interpreter','latex');
ylim([min(min(muHat(:,1:tend)))-0.05,max(max(muHat(:,1:tend)))+0.05])
set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

%%
f9 = figure(9);
box on
p = plot(Time,delta(:,1:tend),'LineWidth',1.25);
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
p = plot(Time,Cost(:,1:tend),'LineWidth',1.25);
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
p = plot(Time,minEigHistP(:,1:tend),Time(1:end-1),minEigHistCL(:,1:tend-1),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);



if includeTitle==1
    title('Minimum Eigenvalues for System Identification','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\lambda_{min}$','interpreter','latex');
l = legend('$\lambda_{min} \left\{P(t)\right\}$','$\lambda_{min}\left\{\sum_{i}^{M}{\mathcal{Y}^T\mathcal{Y}}\right\}$');
set(l,'interpreter','latex');

set(gca,'fontsize',12);
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
        print( f1, '-dpdf', 'FTStates', '-opengl');
        print( f2, '-dpdf', 'FTInput', '-opengl');
        print( f3, '-dpdf', 'FTlCriticWeightEstimates', '-opengl');
        print( f4, '-dpdf', 'FTActorWeightEstimates', '-opengl');
        print( f5, '-dpdf', 'FTSysWeightError', '-opengl');
        print( f6, '-dpdf', 'FTSysWeightEstimate', '-opengl');
        print( f7, '-dpdf', 'FTFaultEstimateError', '-opengl');
        print( f8, '-dpdf', 'FTFaultEstimate', '-opengl');
        print( f9, '-dpdf', 'FTBE', '-opengl');
        print( f10, '-dpdf', 'FTTotCost', '-opengl');
        print( f11, '-dpdf', 'FTSysIDCLMinEigs', '-opengl');
    end
    if (savePlotsInJpeg==1)
        print( f1, '-djpeg', 'FTStates');
        print( f2, '-djpeg', 'FTInput');
        print( f3, '-djpeg' ,'FTCriticWeightEstimates')
        print( f4, '-djpeg' ,'FTActorWeightEstimates')
        print( f5, '-djpeg' ,'FTSysWeightError')
        print( f6, '-djpeg' ,'FTSysWeightEstimate')
        print( f7, '-djpeg' ,'FTFaultEstimateError')
        print( f8, '-djpeg' ,'FTFaultEstimate')
        print( f9, '-djpeg' ,'FTBE')
        print( f10, '-djpeg' ,'FTTotCost')
        print( f11, '-djpeg' ,'FTSysIDCLMinEigs')
    end
    if (savePlotsInPng==1)
        print( f1, '-dpng', 'FTStates');
        print( f2, '-dpng', 'FTInput');
        print( f3, '-dpng' ,'FTCriticWeightEstimates')
        print( f4, '-dpng' ,'FTActorWeightEstimates')
        print( f5, '-dpng' ,'FTSysWeightError')
        print( f6, '-dpng' ,'FTSysWeightEstimate')
        print( f7, '-dpng' ,'FTFaultEstimateError')
        print( f8, '-dpng' ,'FTFaultEstimate')
        print( f9, '-dpng' ,'FTBE')
        print( f10, '-dpng' ,'FTTotCost')
        print( f11, '-dpng' ,'FTSysIDCLMinEigs')
    end
    if (savePlotsInEps==1)
        print( f1, '-deps', 'FTStates');
        print( f2, '-deps', 'FTInput');
        print( f3, '-deps' ,'FTCriticWeightEstimates')
        print( f4, '-deps' ,'FTActorWeightEstimates')
        print( f5, '-deps' ,'FTSysWeightError')
        print( f6, '-deps' ,'FTSysWeightEstimate')
        print( f7, '-deps' ,'FTFaultEstimateError')
        print( f8, '-deps' ,'FTFaultEstimate')
        print( f9, '-deps' ,'FTBE')
        print( f10, '-deps' ,'FTTotCost')
        print( f11, '-deps' ,'FTSysIDCLMinEigs')
    end
    cd(MainDir)
end 




%%
%Save Ouputs, Initial Conditions


if (savePlotsAndVars==1)
    cd(FolderDest)
    save('Conditions','auxdata')
    save('Outputs','time','tend','x','WcHat','WaHat','thetaHat','muHat','mu','WfHat','minEigHistP','minEigHistCL','u','delta','minEigHistADP','P','Q','GammaADP','GammaCL','Cost')
    cd(MainDir)
end

