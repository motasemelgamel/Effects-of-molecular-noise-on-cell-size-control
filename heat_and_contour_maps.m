clear all

%Calculate Theoretical Noise
start=-1;lim=1;step=0.005;Z=(start:step:lim);
rho=ones(length(Z),1)*(10.^Z);
gamma=transpose(rho);P=1.75*60;ymax=60;

x0=gamma.*(2.^rho -1)./(rho.*(2.^(rho+1) -1)) + 1./(1+rho);
xstar=2*x0;
A=(1./(1+rho)).*(2 + 2.^(-rho) .*rho) + gamma.* 2.^(-rho) - x0.*rho.* 2.^(-rho);
B=(1./(1+rho)).*(2 - 2.^(-rho));
beta=B./A;
denom=1./(beta.*(2-beta));

sigeta=((1+gamma)/P).*(1./(A.^2)).*(xstar - x0.* 2.^(-2*rho));
sigepsilon=sigeta.*denom;

%Size noise as a function of beta
beta1=(0.0000001:0.0001:0.5);
beta2=(0.5:0.0001:1);betatot=[beta1 beta2];
sigepsilon_min=[((beta1-1).*(beta1+log(2)-2*beta1*log(2)) )./(P.*beta1.*(beta1-2))   (log(1./(1-beta2)).*(1+2.*beta2.*(beta2-2)) )./(log(4)*P.*beta2.*(beta2-2)) ];

%--------------------------------------------------------

figure(5)
contour(10.^Z,10.^Z,sigepsilon,'LevelList',0.0001:0.005:0.1,'linewidth',1.5)
set(gca,'YDir','normal','FontSize',55,'ColorScale','log','YScale','log','XScale','log')
ylim([min(10.^Z) max(10.^Z)])
xlim([min(10.^Z) max(10.^Z)])
xlabel('Degradation per growth, $\rho=\lambda/\alpha$','interpreter','latex')
ylabel('Production ratio, $\gamma=\nu/\mu\bar{b}$','interpreter','latex')
title('Size noise, $\sigma^2_b/\bar{b}^2$= homeostasis factor $\times$ timing noise','interpreter','latex')
pbaspect([1 1 1])
h=colorbar;
h.Ticks=[0.01 0.03 0.05 0.07];
% set(h,'ylim',[min(sigepsilon(:)) max(sigepsilon(:))])
%print('-painters','-dsvg','myVectorFile1')
%title(h,'$\frac{\sigma^2_b}{\bar{b}^2}$','interpreter','latex')

figure(4)
contour(10.^Z,10.^Z,denom,'LevelList',0.0001:0.1:6,'linewidth',1.5)
set(gca,'YDir','normal','FontSize',55,'ColorScale','log','YScale','log','XScale','log')
ylim([min(10.^Z) max(10.^Z)])
xlim([min(10.^Z) max(10.^Z)])
xlabel('Degradation per growth, $\rho=\lambda/\alpha$','interpreter','latex')
ylabel('Production ratio, $\gamma=\nu/\mu\bar{b}$','interpreter','latex')
title('Homeostasis factor, $\beta^{-1} (2-\beta)^{-1}$','interpreter','latex')
colormap winter
pbaspect([1 1 1])
h=colorbar;
% set(h,'ylim',[min(denom(:)) max(denom(:))])
%print('-painters','-dsvg','myVectorFile2')
%title(h,'$\frac{1}{\beta (2-\beta)}$','interpreter','latex')

figure(1)
contour(10.^Z,10.^Z,sigeta,'LevelList',0.0001:0.005:0.1,'linewidth',1.5)
set(gca,'YDir','normal','FontSize',55,'ColorScale','log','YScale','log','XScale','log')
ylim([min(10.^Z) max(10.^Z)])
xlim([min(10.^Z) max(10.^Z)])
xlabel('Degradation per growth, $\rho=\lambda/\alpha$','interpreter','latex')
ylabel('Production ratio, $\gamma=\nu/\mu\bar{b}$','interpreter','latex')
title('Timing noise, $\sigma^2_{\eta}$','interpreter','latex')
colormap copper
pbaspect([1 1 1])
h=colorbar;
h.Ticks=[0.01 0.03 0.05 0.07];
% set(h,'ylim',[min(sigeta(:)) max(sigeta(:))])
%print('-painters','-dsvg','myVectorFile3')
%title(h,'$\sigma^2_{\eta}$','interpreter','latex')

figure(3)
imagesc(10.^Z,10.^Z,beta)
set(gca,'YDir','normal','FontSize',55,'YScale','log','XScale','log')
ylim([min(10.^Z) max(10.^Z)])
xlim([min(10.^Z) max(10.^Z)])
xlabel('Degradation per growth, $\rho=\lambda/\alpha$','interpreter','latex')
ylabel('Production ratio, $\gamma=\nu/\mu\bar{b}$','interpreter','latex')
title('Homeostasis parameter, $\beta$','interpreter','latex')
pbaspect([1 1 1])
h=colorbar;
set(h,'ylim',[0 1])