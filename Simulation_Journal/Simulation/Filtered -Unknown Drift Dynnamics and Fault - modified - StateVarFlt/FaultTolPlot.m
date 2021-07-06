 savePlotsAndVars = 1;
 savePlotsInJpeg = 1;
 savePlotsInPdf = 0;
 savePlotsInPng = 0;
 includeTitle = 1;
 
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





VHat = zeros(1,tend);
u0  = zeros(m,tend);  %complementary robust controller paper
u1 = zeros(m,tend);   %approximate saturated controller
u  = zeros(m,tend);   %total controller
gHat = zeros(n,m,tend);
gActual = zeros(n,m,tend);
fActual = zeros(n,tend);
% fHat = zeros(n,tend); 
S = zeros(1,tend);
Wf = [-1,-0.5 ;1,0;0,-0.5];
auxInput=zeros(m,tend);
for t=1:tend
    phif = [x(1,t);x(2,t);x(2,t)*(1-(cos(2*x(1,t))+2)^2)];
    fActual(:,t) =  Wf'*phif;
    
    [phif,phig,go,fo] = FaultTolGetBasisDyn(x(:,t),auxdata); %Get basis functions for estimator 
%     fHat(:,t) = reshape(WfHat(:,t),auxdata.pf,n)'*phif;
    
    [sig,sigPrime] = FaultTolGetBasisADP(x(:,t)./auxdata.NormDyn,auxdata); %Get StaF basis functions and their gradients
    %Get Approximate Drift dynamics and G
    gHat(:,:,t) = go*(diag(WgHatTrans(:,:,t)*phig));
     gActual(:,:,t) = go*(diag(mu(:,t)));
     y = x(:,t)./auxdata.NormDyn;
    S(:,t) = y'*auxdata.InputKp1*y/(1+auxdata.InputKp2*(y'*y))^2;
    VHat(t) = WcHat(:,t)'*sig+S(:,t);
    [u1(:,t), uCost,~,~,auxInput(:,t)] = FaultTolGetInput(x(:,t),gHat(:,:,t),WaHat(:,t),sigPrime,auxdata);
    
    if compareToComplementaryPaper == 1;
        u0(:,t) = auxdata.KComp*x(:,t);
        u(:,t) =u0(:,t)+u1(:,t);
    else
        u0(:,t) = 0;
        u(:,t) = u1(:,t);
    end
    
end

%%

 
if ~exist('time')
    dtime = 0.001;
    time = 0:0.001:(length(x)-1)*dtime;
    if ~exist('tend')
        tend = 15000;
        tend = min(tend, length(x));
    end
end       
        
Time = time(1:tend);

fault = mu(:,end);

VStar = 1/2*x(1,:).^2+x(2,:).^2;
uStar = -fault*(cos(2*x(1,:))+2).*x(2,:);

if isempty(find(Time==15))~=1
    TF = min(find(Time==15),tend);
else
    TF = tend;
end


%%
f1 = figure(1);
box on

p = plot(Time(1:TF),x(:,1:TF),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);

if includeTitle==1
%     title('States','interpreter','latex');
    title('(a)','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('x(t)','interpreter','latex');
l  = legend('$x_{1}(t)$','$x_{2}(t)$');
set(l,'interpreter','latex')

set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on



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
    uSat = repmat([auxdata.sat;-auxdata.sat],1,length(u));
    hold on
    p = plot(Time(1:TF),uSat(:,1:TF),Time(1:TF),u(:,1:TF),'LineWidth',1.25);
    set(p(1:2),'Color',[0,0.5,0],'linestyle','--');
    set(p(end),'Color',[0,0,0.5],'linestyle','-');
    
    ylabel('$u(t)$','interpreter','latex');
    
end

if includeTitle==1
%     title('Approximate Optimal Input','interpreter','latex');
    title('(b)','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');

if compareToComplementaryPaper==1
    ylim([min(min([uSat(:,1:tend);u(:,1:tend)]))-0.25,max(max([uSat(:,1:tend);u(:,1:tend)]))+0.25])
else
    ylim([min(min(uSat(:,1:tend)))-0.25,max(max(uSat(:,1:tend)))+0.25])
end
l = legend([p(1),p(end)], '$\{ -\alpha, \alpha \}$', '$u(t)$');
set(l,'interpreter','latex');


set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on

%%
f3 = figure(3);
box on
p = plot(Time(1:TF),WcHat(:,1:TF),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);
set(p(3),'Color',[0.5,0,0]);

if includeTitle==1
%     title('Critic Weight Estimates','interpreter','latex');
    title('(e)','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{W}_{c}(t)$','interpreter','latex');
ylim([min(min(WcHat(:,1:tend)))-0.25,max(max(WcHat(:,1:tend)))+0.25])
l  = legend('$\hat{W}_{c1}(t)$','$\hat{W}_{c2}(t)$','$\hat{W}_{c3}(t)$');
set(l,'interpreter','latex');

set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on



%%
f4 = figure(4);
box on
p = plot(Time(1:TF),WaHat(:,1:TF),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);
set(p(3),'Color',[0.5,0,0]);

if includeTitle==1
%     title('Actor Weight Estimates','interpreter','latex');
    title('(f)','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{W}_{a}(t)$','interpreter','latex');
ylim([min(min(WaHat(:,1:TF)))-0.25,max(max(WaHat(:,1:TF)))+0.25])
l  = legend('$\hat{W}_{a1}(t)$','$\hat{W}_{a2}(t)$','$\hat{W}_{a3}(t)$');
set(l,'interpreter','latex');

set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on

%%

if auxdata.ExactBasis==1
    f5 = figure(5);
    box on
    WfTrans = repmat([-1,-0.5 ;1,0;0,-0.5]',1,1,TF);
    WfTilde = reshape(WfTrans - WfHatTrans(:,:,1:TF),n*auxdata.pf,TF);
    % WfTilde = Wf(:,1:TF)-reshape(WfHat(:,:,1:TF),n*auxdata.pf,TF);
    p = plot(Time(1:TF),WfTilde,'LineWidth',1.25);

    if includeTitle==1
        title('Drift Weight Estimate Error','interpreter','latex');
    end
    xlabel('time (sec)','interpreter','latex');
    ylabel('$\tilde{W}_{1}(t)$','interpreter','latex');
    ylim([min(min(min(WfTilde(:,1:TF))))-0.25,max(max(max(WfTilde(:,1:TF))))+0.25])
    set(gca,'fontsize',30);
    set(gca,'FontName','Times New Roman');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    grid on
end



%%
f6 = figure(6);
box on
p = plot(Time(1:TF),reshape(WfHatTrans(:,:,1:TF),n*auxdata.pf,TF),'LineWidth',1.25);

if includeTitle==1
    title('Drift Weight Estimates','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{W_{1}}(t)$','interpreter','latex');
ylim([min(min(min(WfHatTrans(:,:,1:TF))))-0.25,max(max(max(WfHatTrans(:,:,1:TF))))+0.25])
set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on
%%
f7 = figure(7);
box on
p = plot(Time(1:TF),reshape(WgHatTrans(:,:,1:TF),m*auxdata.pg,TF),'LineWidth',1.25);

if includeTitle==1
    title('Control Fault Weight Estimates','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{W_{2}}(t)$','interpreter','latex');
ylim([min(min(min(WgHatTrans(:,:,1:TF))))-0.25,max(max(max(WgHatTrans(:,:,1:TF))))+0.25])
set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on

%%
f8 = figure(8);
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
set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on
%%
f9 = figure(9);
box on
p = plot(Time(1:TF),muHat(:,1:TF),'LineWidth',1.25);
set(p,'Color',[0,0,0.5]);

if includeTitle==1
    title('Control Fault Estimate','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\hat{\mu}_{a}(t)$','interpreter','latex');
ylim([min(min(muHat(:,1:TF)))-0.05,max(max(muHat(:,1:TF)))+0.05])
set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on
%%
f10 = figure(10);
box on
p = plot(Time(1:TF),delta(:,1:TF),'LineWidth',1.25);
set(p,'Color',[0,0,0.5]);

if includeTitle==1
    title('Bellman Error','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\delta(t)$','interpreter','latex');

set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on

%%
f11 = figure(11);
box on
p = plot(Time(1:TF),Cost(:,1:TF),'LineWidth',1.25);
set(p,'Color',[0,0,0.5]);

if includeTitle==1
    title('Total Cost','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$J(x(t),u(t))$','interpreter','latex');

set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on

%%
f12 = figure(12);

box on
p = plot(Time(1:TF),minEigHistP(:,1:TF),Time(1:TF-1),minEigHistCL(:,1:TF-1),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);



if includeTitle==1
    title('System Identification Minimum Eigenvalues','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\lambda_{min}$','interpreter','latex');
l = legend('$\lambda_{min} \left\{P(t)\right\}$','$\lambda_{min}\left\{\sum_{i}^{N}{\mathcal{Y}^T\mathcal{Y}}\right\}$');
set(l,'interpreter','latex');

set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6,5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on


%%
f13 = figure(13);
box on

FTilde = fActual(:,1:TF)-fHat(:,1:TF);
p = plot(Time(1:TF),FTilde,'LineWidth',1.25);

if includeTitle==1
    title('Drift Dynamics Estimate Error','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$\tilde{f}(t)$','interpreter','latex');
ylim([min(min(FTilde(:,1:TF)))-0.25,max(max(FTilde(:,1:TF)))+0.25])
set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on
%%
f14 = figure(14);
box on

gTilde = reshape(gActual(:,:,1:TF)-gHat(:,:,1:TF),n*m,TF);
p = plot(Time(1:TF),gTilde,'LineWidth',1.25);

if includeTitle==1
    title('Control Effectiveness matrix $G(x(t))$ Error','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
ylabel('$g(x(t))\Lambda(\tilde{\mu}(t))$','interpreter','latex');
ylim([min(min(gTilde(:,1:TF)))-0.25,max(max(gTilde(:,1:TF)))+0.25])
set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on




%%
if size(mu,1)==1
    f15 = figure(15);
    box on
    p = plot(Time(2:TF),mu(:,2:TF),Time(1:TF),muHat(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0.5,0],'linestyle','--');
    set(p(2),'Color',[0,0,0.5]);
    

    if includeTitle==1
%         title('Control Fault Estimate','interpreter','latex');
        title('(d)','interpreter','latex');
    end
    xlabel('time (sec)','interpreter','latex');
    ylabel('$\hat{\mu}_{a}(t)$','interpreter','latex');
    l = legend([p(1),p(2)],'$\mu_{a}(x(t))$','$\hat{\mu}_{a}(x(t),\hat{W}_{2}(t))$');
    set(l,'interpreter','latex');
    ylim([min(min([muHat(:,2:TF);mu(:,2:TF)]))-0.05,max(max([muHat(:,1:TF);mu(:,1:TF)]))+0.05])
    set(gca,'fontsize',30);
    set(gca,'FontName','Times New Roman');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    grid on
end





%%

%%

if auxdata.ExactBasis==1
    f16 = figure(16);
    box on
    Wf = reshape(repmat([-1,-0.5 ;1,0;0,-0.5]',1,1,TF),n*auxdata.pf,TF);
    WfHat = reshape(WfHatTrans(:,:,1:TF),n*auxdata.pf,TF);
    
    
    % WfTilde = Wf(:,1:TF)-reshape(WfHat(:,:,1:TF),n*auxdata.pf,TF);
    p = plot(Time(1:TF),Wf,Time(1:TF),WfHat,'LineWidth',1.25);
    set(p(1:n*auxdata.pf),'color',[0,0.5,0],'linestyle','--');
    set(p(n*auxdata.pf+1:end),'color',[0,0,0.5]);
    
    if includeTitle==1
%         title('Drift Dynamics Weight Estimates','interpreter','latex');
        title('(c)','interpreter','latex');
    end
    xlabel('time (sec)','interpreter','latex');
    ylabel('$\hat{W}_{1}(t)$','interpreter','latex');
    l = legend([p(1),p(n*auxdata.pf+1)],'$W_{1}$','$\hat{W}_{1}(t)$');
    set(l,'interpreter','latex')
    ylim([min(min(min([Wf(:,1:TF),WfHat(:,1:TF)])))-0.25,max(max(max([Wf(:,1:TF),WfHat(:,1:TF)])))+0.25])
    set(gca,'fontsize',30);
    set(gca,'FontName','Times New Roman');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
    grid on

end

%%

if auxdata.ExactFaultBasis==1
    f17 = figure(17);
    box on
    WgHat = reshape(WgHatTrans(:,:,1:TF),m*auxdata.pg,TF);
    Wg = mu;
    p = plot(Time(2:TF),Wg(:,2:TF),Time(1:TF),WgHat(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0.5,0],'linestyle','--');
    set(p(2),'Color',[0,0,0.5]);
    

    if includeTitle==1
        title('Control Fault Weight Estimate','interpreter','latex');
    end
    xlabel('time (sec)','interpreter','latex');
    ylabel('$\hat{W}_{2}(t)$','interpreter','latex');
    l = legend([p(1),p(2)],'$W_{2}$','$\hat{W}_{2}(t)$');
    set(l,'interpreter','latex');
    ylim([min(min([muHat(:,2:TF);mu(:,2:TF)]))-0.05,max(max([muHat(:,1:TF);mu(:,1:TF)]))+0.05])
    set(gca,'fontsize',30);
    set(gca,'FontName','Times New Roman');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 5]);
    set(gcf, 'PaperPosition', [0 0 6 5]);
     grid on
end




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
        if auxdata.ExactBasis==1
            print( f5, '-dpdf','FTSysWErr');
        end
        print( f6, '-dpdf','FTSysWEst');
        print( f7, '-dpdf' ,'FTWeightEst')
        print( f8, '-dpdf','FaultErr');
        print( f9, '-dpdf','FaultEst');
        print( f10, '-dpdf','FTBE');
        print( f11, '-dpdf','FTTotCost');
        print( f12, '-dpdf','FTCLEigs');
        print( f13, '-dpdf','FTDrDynErr');
        print( f14, '-dpdf','FTCEffErr');
        print( f15, '-dpdf','FTFltandEst');
        print( f16, '-dpdf','FTWlW1Hat');
        if auxdata.ExactFaultBasis==1
            print( f17, '-dpdf','FTW2W2Hat');
        end
    end
    if (savePlotsInJpeg==1)
        print( f1, '-djpeg', 'FTStates');
        if includeComp==1
            print( f2, '-djpeg', 'FTInputComp');
        else
            print( f2, '-djpeg', 'FTInput');
        end
        print( f3, '-djpeg' ,'FTCritEst');
        print( f4, '-djpeg' ,'FTActEst');
        if auxdata.ExactBasis==1
            print( f5, '-djpeg' ,'FTSysWeightErr');
        end
        print( f6, '-djpeg' ,'FTSysWEst');
        print( f7, '-djpeg' ,'FTWeightEst');
        print( f8, '-djpeg' ,'FaultErr');
        print( f9, '-djpeg' ,'FaultEst');
        print( f10, '-djpeg' ,'FTBE');
        print( f11, '-djpeg' ,'FTTotCost');
        print( f12, '-djpeg' ,'FTCLEigs');
        print( f13, '-djpeg','FTDrDynErr');
        print( f14, '-djpeg','FTCEffErr');
        print( f15, '-djpeg','FTFltandEst');
        print( f16, '-djpeg','FTWlW1Hat');
        if auxdata.ExactFaultBasis==1
            print( f17, '-djpeg','FTW2W2Hat');
        end
    end
    if (savePlotsInPng==1)
        print( f1, '-dpng', 'FTStates');
        if includeComp==1
            print( f2, '-dpng', 'FTInputComp');
        else
            print( f2, '-dpng', 'FTInput');
        end
        print( f3, '-dpng' ,'FTCriticEst');
        print( f4, '-dpng' ,'FTActorEst');
        if auxdata.ExactBasis==1
            print( f5, '-dpng' ,'FTSysWeightErr');
        end
        print( f6, '-dpng' ,'FTSysWEst');
        print( f7, '-dpng' ,'FTWeightEst');
        print( f8, '-dpng' ,'FaultErr');
        print( f9, '-dpng' ,'FaultEst');
        print( f10, '-dpng' ,'FTBE');
        print( f11, '-dpng' ,'FTTotCost');
        print( f12, '-dpng' ,'FTCLEigs');
        print( f13, '-dpng','FTDrDynErr');
        print( f14, '-dpng','FTCEffErr');
        print( f15, '-dpng','FTFltandEst');
        print( f16, '-dpng','FTWlW1Hat');
        if auxdata.ExactFaultBasis==1
            print( f17, '-dpng','FTW2W2Hat');
        end
    end
    if (savePlotsInEps==1)
        print( f1, '-deps', 'FTStates');
        print( f2, '-deps', 'FTInput');
        print( f3, '-deps' ,'FTCriticEst');
        print( f4, '-deps' ,'FTActorEst');
        if auxdata.ExactBasis==1
            print( f5, '-deps' ,'FTSysWeightErr');
        end
        print( f6, '-deps' ,'FTSysWEst');
        print( f7, '-deps' ,'FTWeightEst');
        print( f8, '-deps' ,'FTFaultEstErr');
        print( f9, '-deps' ,'FTFaultEst');
        print( f10, '-deps' ,'FTBE');
        print( f11, '-deps' ,'FTTotCost');
        print( f12, '-deps' ,'FTCLEigs');
        print( f13, '-deps','FTDrDynErr');
        print( f14, '-deps','FTCEffErr');
        print( f15, '-deps','FTFltandEst');
        print( f16, '-deps','FTWlW1Hat');
        if auxdata.ExactFaultBasis==1
            print( f17, '-deps','FTW2W2Hat');
        end

    end
    cd(MainDir)
end 




%%
%Save Ouputs, Initial Conditions


if (savePlotsAndVars==1)
    cd(FolderDest)
    save('Conditions','auxdata')
    save('Outputs','time','tend','x','WcHat','WaHat','thetaHat','muHat','mu',...
        'WfHatTrans','WgHatTrans','minEigHistP','minEigHistCL','u','u1','u0','delta',...
        'minEigHistADP','P','Q','GammaADP','GammaCL','Cost','TF','fActual',...
        'fHat','gActual','gHat','mu','muHat','VHat', 'S','muTilde','gTilde',...
        'FTilde', 'm', 'n', 'TF', 'TFend','tend', 'uSat');
    cd(MainDir)
end

