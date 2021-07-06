clear all
close all


includeTitle = 1;
CompAll = 0;
CompNoS = 1;
%Get the current directory
MainDir = pwd;

if or(and(CompAll==1, CompNoS==1),or(and(CompAll==1, CompNoS==0),and(CompAll==0, CompNoS==0)))
    %Get input from complementary paper by Q. Fan, G.Yang
    cd('Plots and Conditions\StateVarying_SoftPlusGain\VarFT_CompPaper__27-Sep-2018')
    load Outputs.mat
    load Conditions.mat
    uComp = u;
    xComp = x;
    vHatComp = VHat;
    costComp = Cost;
    uSat = repmat([auxdata.sat;-auxdata.sat],1,length(u));
    cd(MainDir)
end

%Get input from my result
cd('Plots and Conditions\ConstFTComparisontoConf\ConstFTDesigned__27-Sep-2018')
load Outputs.mat
uDes = u;
xDes = x;
vHatDes = VHat;
costDes = Cost;
cd(MainDir)


if or(CompAll==1, CompNoS==1)
    %Get input from revious conference result
    cd('Plots and Conditions\ConstFTComparisontoConf\ConstFTPrev__27-Sep-2018')
    load Outputs.mat
    uNoS = u;
    xNoS = x;
    vHatNoS = VHat;
    costNoS = Cost;
    cd(MainDir)
end

if CompAll==1
    TF = min([length(uDes),length(uComp),length(uNoS)]);
elseif and(CompAll==0, CompNoS==1)
    TF = min([length(uDes),length(uNoS)]);
else 
    TF = min([length(uDes),length(uComp)]);
end

%%
f1 = figure(1);
box on
if CompAll==1
    p = plot(time(1:TF),uDes(:,1:TF),time(1:TF),uComp(:,1:TF),time(1:TF),uNoS(:,1:TF),time(1:TF),uSat(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.9]);
    set(p(2),'Color',[0.9,0.25,0],'linestyle','-.');
    set(p(3),'Color',[0,0.9,0],'linestyle','--');
    set(p(4:end),'Color',[0.25,0.25,0.25],'linestyle',':')
    hold on

    pMark = plot(time(1:750:TF),uDes(:,1:750:TF),time(1:1250:TF),uComp(:,1:1250:TF),time(1:1500:TF),uNoS(:,1:1500:TF),'Linestyle','none');
    set(pMark(1),'marker','s','markerfacecolor',get(p(1),'Color'),'color',get(p(1),'Color'),'markersize',10);
    set(pMark(2),'marker','v','markerfacecolor',get(p(2),'Color'),'color',get(p(2),'Color'),'markersize',10);
    set(pMark(3),'marker','o','markerfacecolor',get(p(3),'Color'),'color',get(p(3),'Color'),'markersize',10);

%     l  = legend([pMark(1),pMark(2),pMark(3),p(end)],'Developed controller',...
%         'Complementary controller',...
%         'Controller in preliminary result',...
%         'Controller saturations, $\{ -\alpha, \alpha \}$');
    
    MinY = min(min(min([uDes(:,1:TF);uComp(:,1:TF);uNoS(:,1:TF);uSat(:,1:TF)])));
    MaxY = max(max(max([uDes(:,1:TF);uComp(:,1:TF);uNoS(:,1:TF);uSat(:,1:TF)])));
    
elseif and(CompAll==0, CompNoS==1)
    p = plot(time(1:TF),uDes(:,1:TF),time(1:TF),uNoS(:,1:TF),time(1:TF),uSat(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.9]);
    set(p(2),'Color',[0,0.9,0],'linestyle','--');
    set(p(3:end),'Color',[0.25,0.25,0.25],'linestyle',':')
    hold on

    pMark = plot(time(1:750:TF),uDes(:,1:750:TF),time(1:1500:TF),uNoS(:,1:1500:TF),'Linestyle','none');
    set(pMark(1),'marker','s','markerfacecolor',get(p(1),'Color'),'color',get(p(1),'Color'),'markersize',10);
    set(pMark(2),'marker','o','markerfacecolor',get(p(2),'Color'),'color',get(p(2),'Color'),'markersize',10);

%     l  = legend([pMark(1),pMark(2),p(end)],'Developed controller',...
%         'Controller in preliminary result',...
%         'Controller saturations, $\{ -\alpha, \alpha \}$');
    
    MinY = min(min(min([uDes(:,1:TF);uNoS(:,1:TF);uSat(:,1:TF)])));
    MaxY = max(max(max([uDes(:,1:TF);uNoS(:,1:TF);uSat(:,1:TF)])));
    
else
    p = plot(time(1:TF),uDes(:,1:TF),time(1:TF),uComp(:,1:TF),time(1:TF),uSat(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.9]);
    set(p(2),'Color',[0.9,0.25,0],'linestyle','-.');
    set(p(3:end),'Color',[0.25,0.25,0.25],'linestyle',':')
    hold on

    pMark = plot(time(1:750:TF),uDes(:,1:750:TF),time(1:1250:TF),uComp(:,1:1250:TF),'Linestyle','none');
    set(pMark(1),'marker','s','markerfacecolor',get(p(1),'Color'),'color',get(p(1),'Color'),'markersize',10);
    set(pMark(2),'marker','v','markerfacecolor',get(p(2),'Color'),'color',get(p(2),'Color'),'markersize',10);

%     l  = legend([pMark(1),pMark(2),pMark(3),p(end)],'Developed controller',...
%         'Complementary controller',...
%         'Controller saturations, $\{ -\alpha, \alpha \}$');
    
    MinY = min(min(min([uDes(:,1:TF);uComp(:,1:TF);uSat(:,1:TF)])));
    MaxY = max(max(max([uDes(:,1:TF);uComp(:,1:TF);uSat(:,1:TF)])));
end

    
    
% set(l,'interpreter','latex','location','northeast');
ylabel('$u(t)$','interpreter','latex');

if includeTitle==1
%     title('Approximate Optimal Input','interpreter','latex');
    title('(a)','interpreter','latex');
%     title('(b)','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');

ylim([MinY-0.25,MaxY+0.25])



set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on
% set(l,'fontsize',15);

% print -dpdf ControllerComparison.pdf
print -djpeg Const_ControllerComparison.jpg


%%
f2 = figure(2);
box on
subplot(2,1,1)
if CompAll==1
    p1 = plot(time(1:TF),xDes(1,1:TF),time(1:TF),xComp(1,1:TF),time(1:TF),xNoS(1,1:TF),'LineWidth',1.25);
    set(p1(1),'Color',[0,0,0.9]);
    set(p1(2),'Color',[0.9,0.25,0],'linestyle','-.');
    set(p1(3),'Color',[0,0.9,0],'linestyle','--');
    hold on

    p1Mark = plot(time(1:750:TF),xDes(1,1:750:TF),time(1:1250:TF),xComp(1,1:1250:TF),time(1:1500:TF),xNoS(1,1:1500:TF),'Linestyle','none');
    set(p1Mark(1),'marker','s','markerfacecolor',get(p1(1),'Color'),'color',get(p1(1),'Color'),'markersize',10);
    set(p1Mark(2),'marker','v','markerfacecolor',get(p1(2),'Color'),'color',get(p1(2),'Color'),'markersize',10);
    set(p1Mark(3),'marker','o','markerfacecolor',get(p1(3),'Color'),'color',get(p1(3),'Color'),'markersize',10);
    
%     l  = legend([p1Mark(1),p1Mark(2),p1Mark(3)],'Developed controller', 'Complementary controller', 'Controller  in preliminary result');
    
    MinY = min(min(min([xDes(:,1:TF);xComp(:,1:TF);xNoS(:,1:TF)])));
    MaxY = max(max(max([xDes(:,1:TF);xComp(:,1:TF);xNoS(:,1:TF)])));
    
elseif and(CompAll==0, CompNoS==1)
    p1 = plot(time(1:TF),xDes(1,1:TF),time(1:TF),xNoS(1,1:TF),'LineWidth',1.25);
    set(p1(1),'Color',[0,0,0.9]);
    set(p1(2),'Color',[0,0.9,0],'linestyle','--');
    hold on

    p1Mark = plot(time(1:750:TF),xDes(1,1:750:TF),time(1:1500:TF),xNoS(1,1:1500:TF),'Linestyle','none');
    set(p1Mark(1),'marker','s','markerfacecolor',get(p1(1),'Color'),'color',get(p1(1),'Color'),'markersize',10);
    set(p1Mark(2),'marker','o','markerfacecolor',get(p1(2),'Color'),'color',get(p1(2),'Color'),'markersize',10);
    
%     l  = legend([p1Mark(1),p1Mark(2)],'Developed controller', 'Controller in preliminary result');
    
    MinY = min(min(min([xDes(1,1:TF);xNoS(1,1:TF)])));
    MaxY = max(max(max([xDes(1,1:TF);xNoS(1,1:TF)])));
    
else
    p1 = plot(time(1:TF),xDes(1,1:TF),time(1:TF),xComp(1,1:TF),'LineWidth',1.25);
    set(p1(1),'Color',[0,0,0.9]);
    set(p1(2),'Color',[0.9,0.25,0],'linestyle','-.');
    hold on

    p1Mark = plot(time(1:750:TF),xDes(1,1:750:TF),time(1:1250:TF),xComp(1,1:1250:TF),'Linestyle','none');
    set(p1Mark(1),'marker','s','markerfacecolor',get(p1(1),'Color'),'color',get(p1(1),'Color'),'markersize',10);
    set(p1Mark(2),'marker','v','markerfacecolor',get(p1(2),'Color'),'color',get(p1(2),'Color'),'markersize',10);
    
%     l  = legend([p1Mark(1),p1Mark(2),p1Mark(3)],'Developed controller', 'Complementary controller');
    
    MinY = min(min(min([xDes(:,1:TF);xComp(:,1:TF)])));
    MaxY = max(max(max([xDes(:,1:TF);xComp(:,1:TF)])));
end

if includeTitle==1
%         title('System State','interpreter','latex');
    H = title('(c)','interpreter','latex'); %Const Fault Plots
    set(H, 'Position',[7.5000 2.0 0],'fontsize',30);
%     H = title('(d)','interpreter','latex'); %Var Fault Plots
%     set(H, 'Position',[7.5000 1.2 0],'fontsize',30);
end


ylim([MinY-0.1,MaxY+0.1])
ylabel('$x_{1}(t)$','interpreter','latex');
POS1 = get(gca,'Position');
set(gca,'Position',[0.1700, 0.5575 0.7750, 0.3000]);
set(gca,'XTickLabel',[]);
set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
grid on

subplot(2,1,2)
if CompAll==1
    p2 = plot(time(1:TF),xDes(2,1:TF),time(1:TF),xComp(2,1:TF),time(1:TF),xNoS(2,1:TF),'LineWidth',1.25);
    set(p2(1),'Color',[0,0,0.9]);
    set(p2(2),'Color',[0.9,0.25,0],'linestyle','-.');
    set(p2(3),'Color',[0,0.9,0],'linestyle','--');
    hold on

    p2Mark = plot(time(1:750:TF),xDes(2,1:750:TF),time(1:1250:TF),xComp(2,1:1250:TF),time(1:1500:TF),xNoS(2,1:1500:TF),'Linestyle','none');
    set(p2Mark(1),'marker','s','markerfacecolor',get(p1(1),'Color'),'color',get(p1(1),'Color'),'markersize',10);
    set(p2Mark(2),'marker','v','markerfacecolor',get(p1(2),'Color'),'color',get(p1(2),'Color'),'markersize',10);
    set(p2Mark(3),'marker','o','markerfacecolor',get(p1(3),'Color'),'color',get(p1(3),'Color'),'markersize',10);
    
    MinY = min(min(min([xDes(:,1:TF);xComp(:,1:TF);xNoS(:,1:TF)])));
    MaxY = max(max(max([xDes(:,1:TF);xComp(:,1:TF);xNoS(:,1:TF)])));
    
    
elseif and(CompAll==0, CompNoS==1)
    p2 = plot(time(1:TF),xDes(2,1:TF),time(1:TF),xNoS(2,1:TF),'LineWidth',1.25);
    set(p2(1),'Color',[0,0,0.9]);
    set(p2(2),'Color',[0,0.9,0],'linestyle','--');
    hold on

    p2Mark = plot(time(1:750:TF),xDes(2,1:750:TF),time(1:1500:TF),xNoS(2,1:1500:TF),'Linestyle','none');
    set(p2Mark(1),'marker','s','markerfacecolor',get(p1(1),'Color'),'color',get(p1(1),'Color'),'markersize',10);
    set(p2Mark(2),'marker','o','markerfacecolor',get(p1(2),'Color'),'color',get(p1(2),'Color'),'markersize',10);

    MinY = min(min(min([xDes(2,1:TF);xNoS(2,1:TF)])));
    MaxY = max(max(max([xDes(2,1:TF);xNoS(2,1:TF)])));
    
else
    p2 = plot(time(1:TF),xDes(2,1:TF),time(1:TF),xComp(2,1:TF),'LineWidth',1.25);
    set(p2(1),'Color',[0,0,0.9]);
    set(p2(2),'Color',[0.9,0.25,0],'linestyle','-.');
    hold on

    p2Mark = plot(time(1:750:TF),xDes(2,1:750:TF),time(1:1250:TF),xComp(2,1:1250:TF),'Linestyle','none');
    set(p2Mark(1),'marker','s','markerfacecolor',get(p1(1),'Color'),'color',get(p1(1),'Color'),'markersize',10);
    set(p2Mark(2),'marker','v','markerfacecolor',get(p1(2),'Color'),'color',get(p1(2),'Color'),'markersize',10);
    
    MinY = min(min(min([xDes(:,1:TF);xComp(:,1:TF)])));
    MaxY = max(max(max([xDes(:,1:TF);xComp(:,1:TF)])));
    
end

ylim([MinY-0.1,MaxY+0.1])
ylabel('$x_{2}(t)$','interpreter','latex');
set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
grid on


POS2 = get(gca,'Position');
set(gca,'Position',[0.1700, 0.2400, 0.7750, 0.3000]);

% set(l,'interpreter','latex','Position',[0.1300, 0.8100, 0.7750, 0.1079], 'Units', 'normalized','fontsize',15);


xlabel('time (sec)','interpreter','latex');





set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

% print -dpdf StateComparison.pdf
print -djpeg Const_StateComparison.jpg

%%
f3 = figure(3);
box on
if CompAll==1
    p = plot(time(1:TF),vHatDes(:,1:TF),time(1:TF),vHatComp(:,1:TF),time(1:TF),vHatNoS(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.9]);
    set(p(2),'Color',[0.9,0.25,0],'linestyle','-.');
    set(p(3),'Color',[0,0.9,0],'linestyle','--');
    hold on

    pMark = plot(time(1:1000:TF),vHatDes(:,1:1000:TF),time(1:1250:TF),vHatComp(:,1:1250:TF),time(1:750:TF),vHatNoS(:,1:750:TF),'Linestyle','none');
    set(pMark(1),'marker','s','markerfacecolor',get(p(1),'Color'),'color',get(p(1),'Color'),'markersize',10);
    set(pMark(2),'marker','v','markerfacecolor',get(p(2),'Color'),'color',get(p(2),'Color'),'markersize',10);
    set(pMark(3),'marker','o','markerfacecolor',get(p(3),'Color'),'color',get(p(3),'Color'),'markersize',10);

    l  = legend([pMark(1),pMark(2),pMark(3)],'Developed controller',...
        'Complementary controller',...
        'Controller in preliminary result');
    
    MinY = min(min(min([vHatDes(:,1:TF);vHatComp(:,1:TF);vHatNoS(:,1:TF)])));
    MaxY = max(max(max([vHatDes(:,1:TF);vHatComp(:,1:TF);vHatNoS(:,1:TF)])));
    
elseif and(CompAll==0, CompNoS==1)
    p = plot(time(1:TF),vHatDes(:,1:TF),time(1:TF),vHatNoS(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.9]);
    set(p(2),'Color',[0,0.9,0],'linestyle','--');
    hold on

    pMark = plot(time(1:1000:TF),vHatDes(:,1:1000:TF),time(1:750:TF),vHatNoS(:,1:750:TF),'Linestyle','none');
    set(pMark(1),'marker','s','markerfacecolor',get(p(1),'Color'),'color',get(p(1),'Color'),'markersize',10);
    set(pMark(2),'marker','o','markerfacecolor',get(p(2),'Color'),'color',get(p(2),'Color'),'markersize',10);

    l  = legend([pMark(1),pMark(2)],'Developed controller',...
        'Controller  in preliminary result');
    
    MinY = min(min(min([vHatDes(:,1:TF);vHatNoS(:,1:TF)])));
    MaxY = max(max(max([vHatDes(:,1:TF);vHatNoS(:,1:TF)])));
    
else
    p = plot(time(1:TF),vHatDes(:,1:TF),time(1:TF),vHatComp(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.9]);
    set(p(2),'Color',[0.9,0.25,0],'linestyle','-.');
    hold on

    pMark = plot(time(1:1000:TF),vHatDes(:,1:1000:TF),time(1:1250:TF),vHatComp(:,1:1250:TF),'Linestyle','none');
    set(pMark(1),'marker','s','markerfacecolor',get(p(1),'Color'),'color',get(p(1),'Color'),'markersize',10);
    set(pMark(2),'marker','v','markerfacecolor',get(p(2),'Color'),'color',get(p(2),'Color'),'markersize',10);

    l  = legend([pMark(1),pMark(2),pMark(3)],'Developed controller',...
        'Complementary controller');
    
    MinY = min(min(min([vHatDes(:,1:TF);vHatComp(:,1:TF)])));
    MaxY = max(max(max([vHatDes(:,1:TF);vHatComp(:,1:TF)])));
end
set(l,'interpreter','latex','location','southeast');
ylabel('$\hat{V}(t)$','interpreter','latex');

if includeTitle==1
    title('System State','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');

ylim([MinY-0.25,MaxY+0.25])



set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on

% print -dpdf vHatComparison.pdf
% print -djpeg vHatComparison.jpg
%%
f4 = figure(4);
box on

if CompAll==1
    p = plot(time(1:TF),costDes(:,1:TF),time(1:TF),costComp(:,1:TF),time(1:TF),costNoS(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.9]);
    set(p(2),'Color',[0.9,0.25,0],'linestyle','-.');
    set(p(3),'Color',[0,0.9,0],'linestyle','--');
    hold on

    pMark = plot(time(1:1000:TF),costDes(:,1:1000:TF),time(1:1250:TF),costComp(:,1:1250:TF),time(1:750:TF),costNoS(:,1:750:TF),'Linestyle','none');
    set(pMark(1),'marker','s','markerfacecolor',get(p(1),'Color'),'color',get(p(1),'Color'),'markersize',10);
    set(pMark(2),'marker','v','markerfacecolor',get(p(2),'Color'),'color',get(p(2),'Color'),'markersize',10);
    set(pMark(3),'marker','o','markerfacecolor',get(p(3),'Color'),'color',get(p(3),'Color'),'markersize',10);

%     l  = legend([pMark(1),pMark(2),pMark(3)],'Developed controller',...
%         'Complementary controller',...
%         'Controller in preliminary result');
    
    MinY = min(min(min([costDes(:,1:TF);costComp(:,1:TF);costNoS(:,1:TF)])));
    MaxY = max(max(max([costDes(:,1:TF);costComp(:,1:TF);costNoS(:,1:TF)])));
    
elseif and(CompAll==0, CompNoS==1)
    p = plot(time(1:TF),costDes(:,1:TF),time(1:TF),costNoS(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.9]);
    set(p(2),'Color',[0,0.9,0],'linestyle','--');
    hold on

    pMark = plot(time(1:1000:TF),costDes(:,1:1000:TF),time(1:750:TF),costNoS(:,1:750:TF),'Linestyle','none');
    set(pMark(1),'marker','s','markerfacecolor',get(p(1),'Color'),'color',get(p(1),'Color'),'markersize',10);
    set(pMark(2),'marker','o','markerfacecolor',get(p(2),'Color'),'color',get(p(2),'Color'),'markersize',10);

%     l  = legend([pMark(1),pMark(2)],'Developed Controller','Controller in preliminary result');
    MinY = min(min(min([costDes(:,1:TF);costNoS(:,1:TF)])));
    MaxY = max(max(max([costDes(:,1:TF);costNoS(:,1:TF)])));
    
else
    p = plot(time(1:TF),costDes(:,1:TF),time(1:TF),costComp(:,1:TF),'LineWidth',1.25);
    set(p(1),'Color',[0,0,0.9]);
    set(p(2),'Color',[0.9,0.25,0],'linestyle','-.');
    hold on

    pMark = plot(time(1:1000:TF),costDes(:,1:1000:TF),time(1:1250:TF),costComp(:,1:1250:TF),'Linestyle','none');
    set(pMark(1),'marker','s','markerfacecolor',get(p(1),'Color'),'color',get(p(1),'Color'),'markersize',10);
    set(pMark(2),'marker','v','markerfacecolor',get(p(2),'Color'),'color',get(p(2),'Color'),'markersize',10);

%     l  = legend([pMark(1),pMark(2)],'Developed controller',...
%         'Complementary controller');
    
    MinY = min(min(min([costDes(:,1:TF);costComp(:,1:TF)])));
    MaxY = max(max(max([costDes(:,1:TF);costComp(:,1:TF)])));
end
% set(l,'interpreter','latex','location','southeast');
ylabel('Total Cost','interpreter','latex');

if includeTitle==1
  %     title('Total Cost','interpreter','latex');
    title('(e)','interpreter','latex'); %Const fault
%     title('(f)','interpreter','latex'); %Var Fault
end
xlabel('time (sec)','interpreter','latex');

ylim([MinY-0.25,MaxY+0.25])



set(gca,'fontsize',30);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);
grid on
% set(l,'fontsize',15);

% print -dpdf costComparison.pdf
print -djpeg Const_costComparison.jpg