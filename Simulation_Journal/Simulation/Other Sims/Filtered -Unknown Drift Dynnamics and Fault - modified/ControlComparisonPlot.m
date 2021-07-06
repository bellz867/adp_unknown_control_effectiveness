clear all
close all


includeTitle = 0;

%Get the current directory
MainDir = pwd;

%Get input from complementary paper by Q. Fan, G.Yang
cd('CompPpr__14-Sep-2017')
load Outputs.mat
uComp = u;
cd(MainDir)

%Get input from my result
cd('Pat__14-Sep-2017')
load Outputs.mat
uPat = u;
cd(MainDir)



f1 = figure(2);
box on
p = plot(time(1:TF),uPat(:,1:TF),time(1:TF),uComp(:,1:TF),'LineWidth',1.25);
set(p(1),'Color',[0,0,0.5]);
set(p(2),'Color',[0,0.5,0]);
l  = legend('Developed Controller','Controller from Q. Fan, G. Yang');
set(l,'interpreter','latex','location','southeast');
ylabel('$u(t)$','interpreter','latex');

if includeTitle==1
    title('Approximate Optimal Input','interpreter','latex');
end
xlabel('time (sec)','interpreter','latex');
MinY = min(min(min(uPat(:,1:tend))),min(min(uComp(:,1:tend))));
MaxY = max(max(max(uPat(:,1:tend))),max(max(uComp(:,1:tend))));
ylim([MinY-0.25,MaxY+0.25])



set(gca,'fontsize',13);
set(gca,'FontName','Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPosition', [0 0 6 5]);

print -dpdf ControllerComparison.pdf
print -djpeg ControllerComparison.jpg