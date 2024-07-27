clear all

%Size noise as a function of beta
beta1=(0.0000001:0.0001:0.5);
beta2=(0.5:0.0001:1);betatot=[beta1 beta2];

% %fitting alpha and an exponential function of beta
xstar=5000;
sigepsilon_min=[2*(beta1+log(2)*(1-2*beta1)).^2./(xstar.*beta1.*(2-beta1))   (4*beta2-2*beta2.^2 -1)./(xstar.*beta2.*(2-beta2)) ];

figure(2)
clf
hold on
plot(betatot,sigepsilon_min,'k','linewidth',1.5)
ylim([0 0.1])
xlim([0 1])
yticks([1e-4 1e-3 1e-2 1e-1 1])
pbaspect([1 1 1])
title('Fixed threshold, $x^{*}$','interpreter','latex')
set(gca,'FontSize',40,'yscale','log')
xlabel('Homeostasis parameter, $\beta$','interpreter','latex')
ylabel('Size noise, $\sigma^2_{b}/\bar{b}^2$','interpreter','latex')