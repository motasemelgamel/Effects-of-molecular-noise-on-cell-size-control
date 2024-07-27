clear all

%Analysis of Taheri 2015 paper
[sebetaTaheri(1),betaTaheri(1),serateTaheri(1),rateTaheri(1)]=Taheri('..\Taheri 2015\glucose.txt');
[sebetaTaheri(2),betaTaheri(2),serateTaheri(2),rateTaheri(2)]=Taheri('..\Taheri 2015\glucose_6aa.txt');
[sebetaTaheri(3),betaTaheri(3),serateTaheri(3),rateTaheri(3)]=Taheri('..\Taheri 2015\glucose_12aa.txt');
[sebetaTaheri(4),betaTaheri(4),serateTaheri(4),rateTaheri(4)]=Taheri('..\Taheri 2015\glycerol.txt');
[sebetaTaheri(5),betaTaheri(5),serateTaheri(5),rateTaheri(5)]=Taheri('..\Taheri 2015\sorbitol.txt');
[sebetaTaheri(6),betaTaheri(6),serateTaheri(6),rateTaheri(6)]=Taheri('..\Taheri 2015\synthetic_rich.txt');
[sebetaTaheri(7),betaTaheri(7),serateTaheri(7),rateTaheri(7)]=Taheri('..\Taheri 2015\TSB.txt');
%plot(mean(rateTaheri),mean(betaTaheri),'sk','MarkerSize',10,'linewidth',1.5);


%Analysis of Wallden 2016
[sebetaWallden(1),betaWallden(1),serateWallden(1),rateWallden(1)]=Wallden('..\Wallden 2016\DnaQ_pooled_per_cell_cycle_fast_data.txt');
[sebetaWallden(2),betaWallden(2),serateWallden(2),rateWallden(2)]=Wallden('..\Wallden 2016\DnaQ_pooled_per_cell_cycle_intermediate_data.txt');
[sebetaWallden(3),betaWallden(3),serateWallden(3),rateWallden(3)]=Wallden('..\Wallden 2016\DnaQ_pooled_per_cell_cycle_slow_data.txt');
[sebetaWallden(4),betaWallden(4),serateWallden(4),rateWallden(4)]=Wallden('..\Wallden 2016\SeqA_per_cell_cycle_intermediate_data.txt');
[sebetaWallden(5),betaWallden(5),serateWallden(5),rateWallden(5)]=Wallden('..\Wallden 2016\SeqA_per_cell_cycle_slow_data.txt');
%plot(rateWallden,betaWallden,'Oc','MarkerSize',10,'linewidth',1.5);

%Analysis of Si 2019
[sebetaSi(1),betaSi(1),serateSi(1),rateSi(1)]=Si('..\Si2019.xlsx','MG1655 M9 acetate');
[sebetaSi(2),betaSi(2),serateSi(2),rateSi(2)]=Si('..\Si2019.xlsx','MG1655 MOPS glucose');
[sebetaSi(3),betaSi(3),serateSi(3),rateSi(3)]=Si('..\Si2019.xlsx','MG1655 MOPS glycerol 11aa');
[sebetaSi(4),betaSi(4),serateSi(4),rateSi(4)]=Si('..\Si2019.xlsx','NCM3722 MOPS arginine');
[sebetaSi(5),betaSi(5),serateSi(5),rateSi(5)]=Si('..\Si2019.xlsx','NCM3722 MOPS glucose 12aa');
[sebetaSi(6),betaSi(6),serateSi(6),rateSi(6)]=Si('..\Si2019.xlsx','NCM3722 MOPS glucose');
%plot(rateSi,betaSi,'Pk','MarkerSize',10,'linewidth',1.5);


%Analysis of Campos 2014
[sebetaCampos(1),betaCampos(1),serateCampos(1),rateCampos(1)]=Campos('..\Campos et al 2014_data.xlsx','E. coli LB');
[sebetaCampos(2),betaCampos(2),serateCampos(2),rateCampos(2)]=Campos('..\Campos et al 2014_data.xlsx','E. coli M9');
[sebetaCampos(3),betaCampos(3),serateCampos(3),rateCampos(3)]=Campos('..\Campos et al 2014_data.xlsx','E. coli M9 minC');
%plot(rateCampos,betaCampos,'Pb','MarkerSize',10,'linewidth',1.5);

%Analysis and plots of wang 2010 paper
[sebetaWang(1),betaWang(1),serateWang(1),rateWang(1)]=Wang('..\Wang 2010\E. coli MG1655 (CGSC 6300)\20090512\');
[sebetaWang(2),betaWang(2),serateWang(2),rateWang(2)]=Wang('..\Wang 2010\E. coli MG1655 (CGSC 6300)\20090129\');
[sebetaWang(3),betaWang(3),serateWang(3),rateWang(3)]=Wang('..\Wang 2010\E. coli MG1655 (CGSC 6300)\20090131\');
[sebetaWang(4),betaWang(4),serateWang(4),rateWang(4)]=Wang('..\Wang 2010\E. coli MG1655 (CGSC 6300)\20090210\');
[sebetaWang(5),betaWang(5),serateWang(5),rateWang(5)]=Wang('..\Wang 2010\E. coli MG1655 (CGSC 6300)\20090525\YFP0001\');
[sebetaWang(6),betaWang(6),serateWang(6),rateWang(6)]=Wang('..\Wang 2010\E. coli MG1655 (CGSC 6300)\20090525\YFP0002\');
[sebetaWang(7),betaWang(7),serateWang(7),rateWang(7)]=Wang('..\Wang 2010\E. coli MG1655 (CGSC 6300)\20090702\');
[sebetaWang(8),betaWang(8),serateWang(8),rateWang(8)]=Wang('..\Wang 2010\E. coli B_r SJ108 and SJ119\20090412 SJ108\');
[sebetaWang(9),betaWang(9),serateWang(9),rateWang(9)]=Wang('..\Wang 2010\E. coli B_r SJ108 and SJ119\20090529 SJ119\YFP0001\');
[sebetaWang(10),betaWang(10),serateWang(10),rateWang(10)]=Wang('..\Wang 2010\E. coli B_r SJ108 and SJ119\20090529 SJ119\YFP0002\');
[sebetaWang(11),betaWang(11),serateWang(11),rateWang(11)]=Wang('..\Wang 2010\E. coli MG1655 lexA3\20090922\');
[sebetaWang(12),betaWang(12),serateWang(12),rateWang(12)]=Wang('..\Wang 2010\E. coli MG1655 lexA3\20090923\');
[sebetaWang(13),betaWang(13),serateWang(13),rateWang(13)]=Wang('..\Wang 2010\E. coli MG1655 lexA3\20090930\');

%Vashistha 2021
[sebetaVashistha(1),betaVashistha(1),serateVashistha(1),rateVashistha(1)]=Vashistha('..\Vashistha 2021\ALL_LB37_Size.xls','Expt1_0529');
[sebetaVashistha(2),betaVashistha(2),serateVashistha(2),rateVashistha(2)]=Vashistha('..\Vashistha 2021\ALL_LB37_Size.xls','Expt2_050619');
[sebetaVashistha(3),betaVashistha(3),serateVashistha(3),rateVashistha(3)]=Vashistha('..\Vashistha 2021\ALL_LB32_Size.xls','Expt1_1012');
[sebetaVashistha(4),betaVashistha(4),serateVashistha(4),rateVashistha(4)]=Vashistha('..\Vashistha 2021\ALL_LB32_Size.xls','Expt1_1015');
[sebetaVashistha(5),betaVashistha(5),serateVashistha(5),rateVashistha(5)]=Vashistha('..\Vashistha 2021\ALL_LB32_Size.xls','Expt3_2807');
[sebetaVashistha(6),betaVashistha(6),serateVashistha(6),rateVashistha(6)]=Vashistha('..\Vashistha 2021\ALL_LB32_Size_nonsis.xls','Exp1-0728');
[sebetaVashistha(7),betaVashistha(7),serateVashistha(7),rateVashistha(7)]=Vashistha('..\Vashistha 2021\ALL_LB32_Size_nonsis.xls','1012');
[sebetaVashistha(8),betaVashistha(8),serateVashistha(8),rateVashistha(8)]=Vashistha('..\Vashistha 2021\ALL_LB32_Size_nonsis.xls','1015');
%plot(rateVashistha,betaVashistha,'mx','MarkerSize',10,'linewidth',1.5);

% %Fit exponential to data.
alpha_=[rateWang rateCampos rateSi rateWallden  rateTaheri rateVashistha];
beta_=[betaWang betaCampos betaSi betaWallden betaTaheri betaVashistha];

alpha=@(x,beta) x(1).*exp(x(2).*beta);
beta=@(x,beta) log(alpha(x,beta)./x(1))./x(2);
fun_=@(x) sum(sqrt(((beta_-beta(x,beta_))./mean(beta_)).^2  +  ((alpha_ - alpha(x,beta_))./mean(alpha_)).^2));
x=fminsearch(fun_,[1,1]);
betaTheory=(0.001:0.005:1);
rateTheory=x(1)*exp(x(2)*betaTheory);

figure(1)
clf
hold on
% %plot avergaed data points
plot(betaTaheri,rateTaheri,'s','color',[ 0.9100 0.4100 0.1700],'MarkerSize',25,'linewidth',1.5,'HandleVisibility','off')
errorbar(mean(betaTaheri),mean(rateTaheri),std(rateTaheri),std(rateTaheri),std(betaTaheri),std(betaTaheri),'s','color',[ 0.9100 0.4100 0.1700],'MarkerSize',25,'linewidth',1.5,'DisplayName','Taheri 2015')
plot(betaWallden,rateWallden,'O','color',[0.4 0 0],'MarkerSize',25,'linewidth',1.5,'HandleVisibility','off')
errorbar(mean(betaWallden),mean(rateWallden),std(rateWallden),std(rateWallden),std(betaWallden),std(betaWallden),'O','color',[0.4 0 0],'MarkerSize',25,'linewidth',1.5,'DisplayName','Wallden 2016')
plot(betaSi,rateSi,'vr','MarkerSize',25,'linewidth',1.5,'HandleVisibility','off')
errorbar(mean(betaSi),mean(rateSi),std(rateSi),std(rateSi),std(betaSi),std(betaSi),'vr','MarkerSize',25,'linewidth',1.5,'DisplayName','Si 2019')
plot(betaCampos,rateCampos,'^b','MarkerSize',25,'linewidth',1.5,'HandleVisibility','off')
errorbar(mean(betaCampos),mean(rateCampos),std(rateCampos),std(rateCampos),std(betaCampos),std(betaCampos),'^b','MarkerSize',25,'linewidth',1.5,'DisplayName','Campos 2014')
plot(betaWang,rateWang,'d','color',[0.4660 0.6740 0.1880],'MarkerSize',25,'linewidth',1.5,'HandleVisibility','off')
errorbar(mean(betaWang),mean(rateWang),std(rateWang),std(rateWang),std(betaWang),std(betaWang),'d','color',[0.4660 0.6740 0.1880],'MarkerSize',25,'linewidth',1.5,'DisplayName','Wang 2010')
plot(betaVashistha,rateVashistha,'mx','MarkerSize',25,'linewidth',1.5,'HandleVisibility','off')
errorbar(mean(betaVashistha),mean(rateVashistha),std(rateVashistha),std(rateVashistha),std(betaVashistha),std(betaVashistha),'mx','MarkerSize',25,'linewidth',1.5,'DisplayName','Vashistha 2021')
plot(betaTheory,rateTheory,'k','MarkerSize',10,'linewidth',1.5,'DisplayName',"Best fit "+round(x(1),1)+"exp("+round(x(2),1)+"\beta)")

pbaspect([1 1 1])
set(gca,'FontSize',40)
ylabel('Growth rate, $\alpha$ (hr$^{-1}$)','interpreter','latex')
xlabel('Homeostasis parameter, $\beta$','interpreter','latex')
legend show

function [sebeta,beta,serate,rate]=Wang(pth)
folders=string({dir(pth+"xy*").name});
b=[];phi=[];alpha=[];
for j =1:numel(folders)
    pth2=pth+folders(j)+"\";
    A=dir(pth2+"*cell0.dat");
    files = {A.name};
for i =1:numel(files)
    data{i}=importdata(strcat(pth2,files{i})).data;
    data_=data{i}(:,2:4);Time=data{i}(:,1);
    ind=find(data_(:,1)==1);
    if(length(ind)>0)
        lb=data_(ind,2);ld=data_(ind(2:end)-1,2);lb(end)=[];
        Tb=Time(ind);Td=Time(ind(2:end)-1);Tb(end)=[];Tn=abs(Td-Tb);
        check=(lb<ld & lb~=0 &ld~=0);lb=lb(check);ld=ld(check);Tn=Tn(check);
        b=[b; lb/mean(lb)];phi=[phi; log(ld./lb)];alpha=[alpha; log(ld./lb)./Tn];
    end
end
end
alpha=alpha*60;
delta=phi-mean(phi);epsilon=log(b);
mdl=fitlm(epsilon,delta,'Intercept',false);
beta=-mdl.Coefficients.Estimate(1);sebeta=mdl.Coefficients.SE(1);
varb=var(epsilon)/mean(alpha);
rate=mean(alpha);serate=std(alpha)/sqrt(length(alpha));sebeta=mdl.Coefficients.SE(1);
%varb=var(epsilon);
end

function [sebeta,beta,serate,rate]=Vashistha(file,sheet_name)
%give the function file name with path
data=xlsread(file,sheet_name);

%extract lineage data into arrays
k=1;
for i=2:4:size(data,2)
    T(:,k)=data(:,i-1);
    T(:,k+1)=data(:,i-1);
    x(:,k)=data(:,i);
    x(:,k+1)=data(:,i+1);
    k=k+2;
end

%initialize birth size xn
xn(1,:)=x(1,:);
Tb(1,:)=T(1,:);

%Finding xn (birth size) of each generation and growth rate
for i=1:size(x,2)
    k=2;
    for j=3:size(x,1) 
        %Finding xn and Tn
        if ((x(j-1,i)>x(j-2,i) && x(j,i)<x(j-1,i)) && x(j,i)~=0 && x(j-1,i)~=0 && (abs(round(x(j-1,i)-x(j,i)))~=0))
            xn(k,i)=x(j,i);
            Tb(k,i)=T(j,i);Td(k-1,i)=T(j-1,i);
            xd(k-1,i)=x(j-1,i);
            k=k+1; %Update indecies
        end
    end
end

%Calculating phi_n and storing separate lineage data in cell arrays
for i=1:size(xn,2)
    xngen{:,i}=xn(xn(:,i)~=0,i);
    xdgen{:,i}=xd(xd(:,i)~=0,i);
    Tbgen{:,i}=Tb(Tb(2:end,i)~=0,i);
    Tdgen{:,i}=Td(Td(:,i)~=0,i);
    Tgen{:,i}=abs(Tdgen{:,i}-Tbgen{:,i});
    phigen{:,i}=log(xdgen{:,i}./xngen{:,i}(1:end-1,:));
    alpha_{:,i}=phigen{:,i}./Tgen{:,i};
end

%N is Lineage Length
count=1;
N=25;b=[];delta=[];epsilon=[];alpha=[];
for j=1:size(xngen,2)
    if length(xngen{j}(1:end))>N
    epsilon=[epsilon ; log(xngen{j}(1:end-1)/mean(xngen{j}))];
    delta=[delta ; (phigen{j}-mean(phigen{j}))];
    alpha=[alpha ; alpha_{j}];
    b=[b ; xngen{j}(1:end)/mean(xngen{j})];
    end
end

mdl=fitlm(epsilon,delta,'Intercept',false);
beta=-mdl.Coefficients.Estimate(1);
varb=var(epsilon)/mean(alpha);
rate=mean(alpha);serate=std(alpha)/sqrt(length(alpha));sebeta=mdl.Coefficients.SE(1);
%varb=var(epsilon);
end

function [sebeta,beta,serate,rate]=Taheri(file)
data=importdata(file);
data_=data.data;

lb=data_(:,5);ld=data_(:,6);alpha=data_(:,3);
b=lb/mean(lb);w=11.3*0.0645;
epsilon=log(b);phi=log(ld./lb);delta=phi-mean(phi);
s=12*lb./(w.*(3*lb-w));s=s/mean(s);
mdl=fitlm(epsilon,delta,'Intercept',false);
beta=-mdl.Coefficients.Estimate(1);
alpha=alpha*60;
rate=mean(alpha);serate=std(alpha)/sqrt(length(alpha));sebeta=mdl.Coefficients.SE(1);
end

function [sebeta,beta,serate,rate]=Wallden(file)
data=importdata(file);w=11.3*0.0645;

vb=data(:,3);vd=data(:,4);alpha=data(:,1);
b=vb/mean(vb);phi=log(vd./vb);
epsilon=log(b);delta=phi-mean(phi);
s=12*vb./(w.*(3*vb-w));s=s/mean(s);
mdl=fitlm(epsilon,delta,'Intercept',false);
beta=-mdl.Coefficients.Estimate(1);
alpha=alpha*60;
rate=mean(alpha);serate=std(alpha)/sqrt(length(alpha));sebeta=mdl.Coefficients.SE(1);
end

function [sebeta,beta,serate,rate]=Si(file,sheet_name)
data=xlsread(file,sheet_name);
w=11.3*0.0645;

lb=data(1:end-1,9);ld=data(2:end,8);alpha=data(2:end,2);
b=lb/mean(lb);phi=log(ld./lb);
epsilon=log(b);delta=phi-mean(phi);
s=12*lb./(w.*(3*lb-w));s=s/mean(s);
mdl=fitlm(epsilon,delta,'Intercept',false);
beta=-mdl.Coefficients.Estimate(1);
rate=mean(alpha);serate=std(alpha)/sqrt(length(alpha));sebeta=mdl.Coefficients.SE(1);
end

function [sebeta,beta,serate,rate]=Campos(file,sheet_name)
data=xlsread(file,sheet_name);
w=11.3*0.0645;

lb=data(:,1);ld=data(:,2);alpha=data(:,4);
b=lb/mean(lb);phi=log(ld./lb);
epsilon=log(b);delta=phi-mean(phi);
s=12*lb./(w.*(3*lb-w));s=s/mean(s);
mdl=fitlm(epsilon,delta,'Intercept',false);
beta=-mdl.Coefficients.Estimate(1);
alpha=alpha*60;
rate=mean(alpha);serate=std(alpha)/sqrt(length(alpha));sebeta=mdl.Coefficients.SE(1);
end