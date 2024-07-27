clear all

%correlation time 10
load('..\r=0.1.mat')

col=1;
beta_{col}=beta;
sigb_{col}=sigb;
r(col)=r;

%correlation time 1
load('..\r=1.mat')

col=2;
beta_{col}=beta;
sigb_{col}=sigb;
r(col)=r;

%correlation time 0.1
load('..\r=10.mat')

col=3;
beta_{col}=beta;
sigb_{col}=sigb;
r(col)=r;

%only threshold noise
load('..\only threshold noise, r=0.1,eta=0.mat')

col=4;
beta_{col}=beta;
sigb_{col}=sigb;
r(col)=r;

%only partition noise
load('..\only partition noise, r=0,eta=0.mat')

col=5;
beta_{col}=beta;
sigb_{col}=sigb;
r(col)=r;

%only growth rate noise
load('..\only growth rate noise, r=0,eta=0.2.mat')

col=6;
beta_{col}=beta;
sigb_{col}=sigb;
sigepsilon_min(end)=100;
r(col)=r;


figure(1)
clf
subplot(1,2,1)
hold on
plot(beta_{5},sigb_{5},'.','color','r','markersize',10,'DisplayName','Partition noise')
plot(beta_{6},sigb_{6},'.','color','b','markersize',10,'DisplayName','Growth rate noise')
plot(beta_{4},sigb_{4},'.','color','c','markersize',10,'DisplayName','Division threshold noise, correlation time=10 hr.')
plot(betatot,sigepsilon_min,'k','linewidth',1.5,'DisplayName','Theoretical bound')
ylim([0.3 60])
pbaspect([1 1 1])
set(gca,'FontSize',35,'yscale','log')
xlabel('Homeostasis parameter, $\beta$','interpreter','latex')
ylabel('Rescaled size noise, $k \sigma^2_{b}/\bar{\alpha}\bar{b}^2$','interpreter','latex')
legend show

subplot(1,2,2)
hold on
plot(beta_{1},sigb_{1},'.','color','r','markersize',10,'DisplayName','Correlation time=10 hr.')
plot(beta_{2},sigb_{2},'.','color','b','markersize',10,'DisplayName','Correlation time=1 hr.')
plot(beta_{3},sigb_{3},'.','color','c','markersize',10,'DisplayName','Correlation time=0.1 hr.')
plot(betatot,sigepsilon_min,'k','linewidth',1.5,'DisplayName','Theoretical bound')
ylim([0.3 60])
xlim([0 1.1])
pbaspect([1 1 1])
set(gca,'FontSize',35,'yscale','log')
xlabel('Homeostasis parameter, $\beta$','interpreter','latex')
ylabel('Rescaled size noise, $k \sigma^2_{b}/\bar{\alpha}\bar{b}^2$','interpreter','latex')
title('All noise sources','interpreter','latex')
legend show