clear all

N_=1000;sigb=zeros(1,N_);
count=1;Z=logspace(-2,2,10*N_);  P=0.5*60;
A=3.3980;B=1.6693; rng("shuffle")

while count<=N_    
idx1=randperm(length(Z),1);idx2=randperm(length(Z),1);idx3=randperm(length(Z),1);idx4=randperm(length(Z),1);idx5=randperm(length(Z),1);
    
gamma=Z(idx2);omega=Z(idx3);alpha=Z(idx4);a_mu=omega/(1+gamma);a_=a_mu/alpha;rho=Z(idx5);k=rho*alpha; 
N=1000;
x_=ceil(((2*a_*gamma*(2^rho -1))/(rho*(2^(rho+1) -1))) + (2*a_/(1+rho))); %for x0=x*/2
v0=1;T=zeros(N,1);alpha_=zeros(N,1);vb=zeros(N,1);vb(1)=v0;
x0=ceil(x_/2);

betath=(2^(rho+2) -4+2^(-rho))/(gamma*(rho+1) + 2^(rho+2) -2);
alpha=A*exp(-B*betath);

counter=1; eta=0;
gen=1;alpha_(gen)=alpha;t=0;t2=0;v=vb(1);x=x0;
%time(counter)=t2;xT(counter)=x;
y=x_; r=0;
%counter=counter+1;
a=[gamma*a_,rho*x,a_*v,r*x_,r*y];
a0=sum(a);
if x_>=10
while N>gen
    if x < y
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
        elseif mu2==4
            y=y+1;
        elseif mu2==5
            y=y-1;
        end
        %time(counter)=t2;
        %xT(counter)=x;
        v=v0*exp(t);
        %counter=counter+1;
        a=[gamma*a_,rho*x,a_*v,r*x_,r*y];
        a0=sum(a);
    else
        T(gen)=t;
        gen=gen+1;
        alpha_(gen)=alpha+normrnd(0,eta*alpha);
        a_=a_mu/(alpha_(gen)); rho=k/alpha_(gen);
        while alpha_(gen)<0
            alpha_(gen)=alpha+normrnd(0,eta*alpha);
            a_=a_mu/(alpha_(gen)); rho=k/alpha_(gen);
        end
        %time(counter)=t2;
        t=0;
        %x0=max(1,binornd(floor(x),0.5));
        x0=max(1,ceil(x/2));
        x=x0;
        vb(gen)=v/2;v0=v/2;
        %xT(counter)=x;
        %counter=counter+1;
        x_=ceil(((2*a_*gamma*(2^rho -1))/(rho*(2^(rho+1) -1))) + (2*a_/(1+rho)));
        a=[gamma*a_,rho*x,a_*v,r*x_,r*y];
        a0=sum(a);
    end
end

phi=T(100:end-1);epsilon=log(vb(100:end-1)./mean(vb(100:end-1)));delta=phi-mean(phi);
mdl=fitlm(epsilon,delta,'Intercept',false);
beta(count)=-mdl.Coefficients.Estimate(1);

betatheory(count)=(2^(rho+2) -4+2^(-rho))/(gamma*(rho+1) + 2^(rho+2) -2); %x0=x*/2
sigb(count)=omega*var(vb(100:end-1))/(mean(vb(100:end-1))^2)/P;
%sigb(count)=omega*A*exp(-B*beta(count)).*var(vb(100:end-1))/(mean(vb(100:end-1))^2)/alpha/ P;
alpha_2(count)=alpha;
count=count+1;
end
end

betatheory=betatheory(beta>=0);sigb=sigb(beta>=0);alpha_2=alpha_2(beta>=0);beta=beta(beta>=0);
%Size noise as a function of beta
beta1=(0.0000001:0.0001:0.5);
beta2=(0.5:0.0001:1);betatot=[beta1 beta2];
A=3.1432; B=1.6052;
alpha1=A*exp(-B*beta1);alpha2=A*exp(-B*beta2);
sigepsilon_min=[alpha1.*((beta1-1).*(beta1+log(2)-2*beta1*log(2)) )./(P.*beta1.*(beta1-2))   alpha2.*(log(1./(1-beta2)).*(1+2.*beta2.*(beta2-2)) )./(log(4)*P.*beta2.*(beta2-2)) ];
sig_old=0.7./(betatot.*(2-betatot));
%sigepsilon_min=[((beta1-1).*(beta1+log(2)-2*beta1*log(2)) )./(beta1.*(beta1-2))   (log(1./(1-beta2)).*(1+2.*beta2.*(beta2-2)) )./(log(4)*beta2.*(beta2-2)) ];
sigepsilon_min(end)=100;

figure(1)
clf
hold on
%plot(beta,sigb,'.','color',0.75*[1 1 1],'markersize',10,'DisplayName','Simulation')
plot(beta,sigb,'.','color',0.75*[1 1 1],'markersize',10,'DisplayName','Simulation')
plot(betatot,sigepsilon_min,'k','linewidth',1.5,'DisplayName','Theoretical bound')
ylim([0.3 60])
pbaspect([1 1 1])
set(gca,'FontSize',40,'yscale','log')
xlabel('Homeostasis parameter, $\beta$','interpreter','latex')
ylabel('Size noise, $\sigma^2_{b}/\bar{b}^2$','interpreter','latex')
% ylim([2e-1 60]./60)
%yticks([1e-2 1e-1 10])
pbaspect([1 1 1])
legend show

% subplot(1,2,2)
% hold on
% plot(betatheory,beta,'.','color','b','markersize',10)
% plot(betatheory,betatheory,'k','linewidth',1.5)
% pbaspect([1 1 1])
% set(gca,'FontSize',35)
% xlabel('Theoretical $\beta$','interpreter','latex')
% ylabel('Measured $\beta$ from simulation','interpreter','latex')