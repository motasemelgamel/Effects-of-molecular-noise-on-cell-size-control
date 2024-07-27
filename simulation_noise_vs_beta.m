clear all

N_=1000;sigb=zeros(1,N_);
count=1;Z=logspace(-2,2,10*N_);
while count<=N_
    
idx1=randperm(length(Z),1);idx2=randperm(length(Z),1);idx3=randperm(length(Z),1);
    
rho=Z(idx1);gamma=Z(idx2);omega=Z(idx3);a_=omega/(1+gamma); N=1000;
beta(count)=(2^(rho+2) -4+2^(-rho))/(gamma*(rho+1) + 2^(rho+2) -2); %x0=x*/2
x_=ceil(((2*a_*gamma*(2^rho -1))/(rho*(2^(rho+1) -1))) + (2*a_/(1+rho))); %for x0=x*/2
v0=1;T=zeros(N,1);vb=zeros(N,1);vb(1)=v0;
x0=ceil(x_/2);

counter=1;
gen=1;t=0;t2=0;v=vb(1);x=x0;
time(counter)=t2;xT(counter)=x;
counter=counter+1;
a=[gamma*a_,rho*x,a_*v];
a0=sum(a);
while N>gen
    if x<x_
        r1=rand;r2=rand;
        mu2=1;
        tau=(1/a0) * log(1/r1);
        while ~((sum(a(1:mu2-1)) < r2*a0) && (r2*a0 <= sum(a(1:mu2))))
            mu2=mu2+1;
        end
        t=t+tau;
        t2=t2+tau;
        if mu2==1
            x=x+1;
        elseif mu2==2
            x=x-1;
        elseif mu2==3
            x=x+1;
        end
        time(counter)=t2;xT(counter)=x;
        v=v0*exp(t);
        counter=counter+1;
        a=[gamma*a_,rho*x,a_*v];
        a0=sum(a);
    else
        T(gen)=t;
        gen=gen+1;
        t=0;x0=ceil(x/2);x=x0;
        vb(gen)=v/2;v0=v/2;
        time(counter)=t2;xT(counter)=x;
        counter=counter+1;
        a=[gamma*a_,rho*x,a_*v];
        a0=sum(a);
    end
end
sigb(count)=omega*var(vb(5:end))/(mean(vb(5:end))^2);
count=count+1;
end

%Size noise as a function of beta
beta1=(0.0000001:0.0001:0.5);
beta2=(0.5:0.0001:1);betatot=[beta1 beta2];
sigepsilon_min=[((beta1-1).*(beta1+log(2)-2*beta1*log(2)) )./(beta1.*(beta1-2))   (log(1./(1-beta2)).*(1+2.*beta2.*(beta2-2)) )./(log(4)*beta2.*(beta2-2)) ];

figure(1)
hold on
plot(beta,sigb,'.','color',0.75*[1 1 1],'markersize',10,'DisplayName','Simulation')
plot(betatot,sigepsilon_min,'k','linewidth',1.5,'DisplayName','Theoretical bound')
ylim([0.3 60])
pbaspect([1 1 1])
set(gca,'FontSize',35,'yscale','log')
xlabel('Homeostasis parameter, $\beta$','interpreter','latex')
ylabel('Rescaled size noise, $k \sigma^2_{b}/\alpha \bar{b}^2$','interpreter','latex')
legend show
