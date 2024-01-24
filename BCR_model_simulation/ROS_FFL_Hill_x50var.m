%{ 
All content of this script refers to "Dynamic hydrogen peroxide levels reveal a 
rate-dependent sensitivity in B-cell lymphoma signaling" by Witmond et al.

The enclosed code can be used to simulate the computational model of the
simplified BCR network motif. For details, refer to the Materials and
Methods of the manuscript.
%}


%% Define model parameters and simulation parameters
species={'H_2O_2','pCD79a','pSyk','pPLCy2'};
newcolors=[rgb('BlueViolet');rgb('DodgerBlue');rgb('LightSeaGreen')];

p=0.00; 
rp=1.*[1 1 1]; rdp=10.*[1 1 1]; 
Km=ones(1,6);
Km([1 3 5])=1.*[2 0.5 0.125];
Km([2 4 6])=1.*[1 1 1];

kros=0.5.*[1 1 1];
H=1.*[1 1 1];
kpos=1;
x0=[0 0.01 0.01 0.01]; % initial protein concentrations
xT=[1 1 1]; % total protein concentrations

%  Simulate models to construct Dose Response curves
dose=logspace(-1,3,51);
response=zeros(length(dose),3);

tend1=100;
tsim1=linspace(0,tend1,1001);

for i=1:length(dose)
    tG0=0.01; Gmax1=dose(i);
    [~,y] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km,kros,kpos,H,rp,rdp,tG0,Gmax1,'linear',xT),tsim1,x0);
    response(i,:)=y(end,2:end);
end
response_sc=rescale(response,"InputMin",min(response),"InputMax",max(response));

figure
hold on
grid on
colororder(newcolors)
plot(dose,response,'linewidth',2)
legend(species(2:end),'location','northwest','box','off')
set(gca, 'XScale', 'log')
xlabel('H_2O_2 dose (mM)')
ylabel('Response')
xlim([dose(1) dose(end)])
ylim([0 1])
set(gca,'box','off')

% Fit DR curve to simulation results
DRcurve=@(p,ROS) p(1)+p(2)*(ROS.^p(4)./(p(3)^p(4)+ROS.^p(4)));
pDR=[0 1 3 3]; %b,Rmax,kROS,H
lb=[0 0 0 0]; ub=[1 xT(1) 100 10]; opts = optimset('Display','off');

pfit=zeros(3,length(pDR));
for i=1:3
    pfit(i,:) = lsqcurvefit(DRcurve,pDR,dose,response(:,i).',lb,ub,opts);
end


%% Hill version - dynamic input profiles
labels={'Linear 20 min','Linear 60 min','Quadr. 60 min'};

tend1=60;
tsim1=linspace(0,tend1,1001);

tG0=20; 
tG1=60; 
Gmax1=10;
[t0,y1] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km,kros,kpos,H,rp,rdp,tG0,Gmax1,'linear',xT),tsim1,x0);
[t1,y1] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km,kros,kpos,H,rp,rdp,tG1,Gmax1,'linear',xT),tsim1,x0);
[t2,y2] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km,kros,kpos,H,rp,rdp,tG1,Gmax1,'quadratic',xT),tsim1,x0);

figure
colororder(newcolors)
tcl=tiledlayout(1,2);
nexttile(tcl)
hold on
plot(t0,y1(:,1),'--k','linewidth',2)
plot(t1,y1(:,1),'-.k','linewidth',2);
p1=plot(t2,y2(:,1),'k','linewidth',2);
title('Input')
xlabel('Time (min.)')
ylabel('Concentration')
%legend(labels,'location','best','box','off','fontsize',12)

nexttile(tcl)
hold on
plot(t0,y1(:,2:end),'--','linewidth',2)
set(gca,'ColorOrderIndex',1)
p2=plot(t1,y1(:,2:end),'-.','linewidth',2);
set(gca,'ColorOrderIndex',1)
p2=plot(t2,y2(:,2:end),'linewidth',2);
title('Output')
xlabel('Time (min.)')

hL=legend([p1;p2],species,'box','off','fontsize',12);
hL.Layout.Tile='East';
tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';

figure
colororder(newcolors)
hold on
p3=plot(y1(:,1),y1(:,2:end),'k--','linewidth',2);
set(gca,'ColorOrderIndex',1)
plot(y1(:,1),y1(:,2:end),'--','linewidth',2);
p4=plot(y1(:,1),y1(:,2:end),'k-.','linewidth',2);
set(gca,'ColorOrderIndex',1)
plot(y1(:,1),y1(:,2:end),'-.','linewidth',2);
p5=plot(y2(:,1),y2(:,2:end),'k','linewidth',2);
set(gca,'ColorOrderIndex',1)
plot(y2(:,1),y2(:,2:end),'linewidth',2);
xlabel('H_2O_2 conc. (mM)','fontsize',12)
ylabel('Response','fontsize',12)
legend([p3(1),p4(1),p5(1)],labels,'box','off','location','best','fontsize',12)
xlim([0.1,Gmax1])
set(gca,'XScale','log')
title('Dose-response dynamic profiles')


%% Threshold variability
mu = @(m,v) log((m^2)./sqrt(v+m^2));
sigma = @(m,v) sqrt(log(v/(m^2)+1));

n=1000;

% kros thresholds
mx = kros; % mean activation threshold
vx = kros'.*[0.1 0.7 1.25].^2; % variance in activation threshold
kros_1=[lognrnd(mu(mx(1),vx(1)),sigma(mx(1),vx(1)),[n 1]), lognrnd(mu(mx(2),vx(1)),sigma(mx(2),vx(1)),[n 1]), lognrnd(mu(mx(3),vx(1)),sigma(mx(3),vx(1)),[n 1])];
kros_2=[lognrnd(mu(mx(1),vx(2)),sigma(mx(1),vx(2)),[n 1]), lognrnd(mu(mx(2),vx(2)),sigma(mx(2),vx(2)),[n 1]), lognrnd(mu(mx(3),vx(2)),sigma(mx(3),vx(2)),[n 1])];
kros_3=[lognrnd(mu(mx(1),vx(3)),sigma(mx(1),vx(3)),[n 1]), lognrnd(mu(mx(2),vx(3)),sigma(mx(2),vx(3)),[n 1]), lognrnd(mu(mx(3),vx(3)),sigma(mx(3),vx(3)),[n 1])];

% Km thresholds
mxK = Km; 
vxK = Km'.*[0.1 0.7 1.25].^2;
Km_1=[lognrnd(mu(mxK(1),vxK(1,1)),sigma(mxK(1),vxK(1,1)),[n 1]), lognrnd(mu(mxK(2),vxK(2,1)),sigma(mxK(2),vxK(2,1)),[n 1]), lognrnd(mu(mxK(3),vxK(3,1)),sigma(mxK(3),vxK(3,1)),[n 1]),...
      lognrnd(mu(mxK(4),vxK(4,1)),sigma(mxK(4),vxK(4,1)),[n 1]), lognrnd(mu(mxK(5),vxK(5,1)),sigma(mxK(5),vxK(5,1)),[n 1]), lognrnd(mu(mxK(6),vxK(6,1)),sigma(mxK(6),vxK(6,1)),[n 1])];
Km_2=[lognrnd(mu(mxK(1),vxK(1,2)),sigma(mxK(1),vxK(1,2)),[n 1]), lognrnd(mu(mxK(2),vxK(2,2)),sigma(mxK(2),vxK(2,2)),[n 1]), lognrnd(mu(mxK(3),vxK(3,2)),sigma(mxK(3),vxK(3,2)),[n 1]),...
      lognrnd(mu(mxK(4),vxK(4,2)),sigma(mxK(4),vxK(4,2)),[n 1]), lognrnd(mu(mxK(5),vxK(5,2)),sigma(mxK(5),vxK(5,2)),[n 1]), lognrnd(mu(mxK(6),vxK(6,2)),sigma(mxK(6),vxK(6,2)),[n 1])];
Km_3=[lognrnd(mu(mxK(1),vxK(1,3)),sigma(mxK(1),vxK(1,3)),[n 1]), lognrnd(mu(mxK(2),vxK(2,3)),sigma(mxK(2),vxK(2,3)),[n 1]), lognrnd(mu(mxK(3),vxK(3,3)),sigma(mxK(3),vxK(3,3)),[n 1]),...
      lognrnd(mu(mxK(4),vxK(4,3)),sigma(mxK(4),vxK(4,3)),[n 1]), lognrnd(mu(mxK(5),vxK(5,3)),sigma(mxK(5),vxK(5,3)),[n 1]), lognrnd(mu(mxK(6),vxK(6,3)),sigma(mxK(6),vxK(6,3)),[n 1])];


tend1=10;
tsim1=linspace(0,tend1,11);
tG0=0.01; Gmax1=5;

cells1=zeros(n,3);
cells2=zeros(n,3);
cells3=zeros(n,3);
for i=1:n
    [~,y1] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km_1(i,:),kros_1(i,:),kpos,H,rp,rdp,tG0,Gmax1,'linear',xT),tsim1,x0);
    [~,y2] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km_2(i,:),kros_2(i,:),kpos,H,rp,rdp,tG0,Gmax1,'linear',xT),tsim1,x0);
    [~,y3] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km_3(i,:),kros_3(i,:),kpos,H,rp,rdp,tG0,Gmax1,'linear',xT),tsim1,x0);
    cells1(i,:)=y1(end,2:end);
    cells2(i,:)=y2(end,2:end);
    cells3(i,:)=y3(end,2:end);
end

bins=linspace(0,max(xT),31);
figure
tl=tiledlayout(1,3);
for i=1:3
    nexttile
    hold on
    histogram(cells1(:,i),bins,'normalization','probability')
    histogram(cells2(:,i),bins,'normalization','probability')
    histogram(cells3(:,i),bins,'normalization','probability')
    title(species(i+1))
end
lgd=legend('Low','Medium','High','box','off','orientation','horizontal');
lgd.Layout.Tile = 'south';
lgd.Title.String = 'Level of threshold variability \sigma^2_{x_{50}}: ';
lgd.Title.FontSize = 10;
sgtitle(sprintf('H_2O_2:  %d mM',Gmax1))


%% DR curves with threshold variability
%  Simulate models to construct DR curves
dose2=logspace(-2,3,21);
response_var=zeros(length(dose2),n,3);

tend1=10;
tsim1=linspace(0,tend1,11);

for i=1:length(dose2)
    tG0=0.01; Gmax1=dose2(i);
    for j=1:n
        [~,y] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km_3(j,:),kros_3(j,:),kpos,H,rp,rdp,tG0,Gmax1,'linear',xT),tsim1,x0);
        response_var(i,j,:)=y(end,2:end);
    end
end

%% compute thresholds at 2.5%
% rescale data
response_var_sc=rescale(response_var,"InputMin",min(response_var),"InputMax",max(response_var));

figure
hold on
plot(dose2,mean(response_var(:,:,1),2))
plot(dose2,mean(response_var_sc(:,:,1),2))
set(gca,'xscale','log')


threshold=[0.3 0.3 0.3];%[0.5 0.5 0.5];
perc_on=zeros(length(dose2),3);
for i=1:3
    perc_on(:,i)=100.*sum(response_var(:,:,i)>threshold(i),2)/n;
end

figure
colororder(newcolors)
hold on
grid on
plot(dose2,perc_on,'linewidth',2)
set(gca,'xscale','log')
legend(species(2:end),'location','northwest','box','off')
xlabel('H_2O_2 dose (mM)')
ylabel('Percentage on')
xlim([dose2(1) dose2(end)])
ylim([0 100])
set(gca,'box','off')

% Fit DR curve to simulation results
pDR=[0 100 3 3]; %b,Rmax,kROS,H
lb=[0 0 0 0]; ub=[100 100 100 10];

perconfit=zeros(3,length(pDR));
for i=1:3
    perconfit(i,:) = lsqcurvefit(DRcurve,pDR,dose2,perc_on(:,i).',lb,ub,opts);
end


%% Dynamic input with threshold variability
labels={'Linear 20 min','Linear 60 min','Quadr. 60 min'};

% gradient parameters (max, duration)
Gmax1=10;
tG0=20; 
tG1=60;

tsim0=linspace(0,tG0,101);
tsim1=linspace(0,tG1,101);

[t0,y1] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km,kros,kpos,H,rp,rdp,tG0,Gmax1,'linear',xT),tsim0,x0);
[t1,y1] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km,kros,kpos,H,rp,rdp,tG1,Gmax1,'linear',xT),tsim1,x0);
[t2,y2] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km,kros,kpos,H,rp,rdp,tG1,Gmax1,'quadratic',xT),tsim1,x0);

% Compute at which time point H2O2 concentrations are identical
yq=linspace(0,Gmax1,31); % query H2O2 concentration points
t0q=interp1(y1(:,1),t0,yq);
t1q=interp1(y1(:,1),t1,yq);
t2q=interp1(y2(:,1),t2,yq);
t0q(end)=20;
t1q(end)=60;
t2q(end)=60;

LG20_var=zeros(length(t0q),n,3);
LG60_var=zeros(length(t1q),n,3);
QG60_var=zeros(length(t2q),n,3);
for i=1:n
    [t0var,y0var] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km_3(i,:),kros_3(i,:),kpos,H,rp,rdp,tG0,Gmax1,'linear',xT),t0q,x0);
    [t1var,y1var] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km_3(i,:),kros_3(i,:),kpos,H,rp,rdp,tG1,Gmax1,'linear',xT),t1q,x0);
    [t2var,y2var] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km_3(i,:),kros_3(i,:),kpos,H,rp,rdp,tG1,Gmax1,'quadratic',xT),t2q,x0);
    LG20_var(:,i,:)=y0var(:,2:end);
    LG60_var(:,i,:)=y1var(:,2:end);
    QG60_var(:,i,:)=y2var(:,2:end);
end

%% Plot the results
LG20_perc_on=zeros(length(t0q),3);
LG60_perc_on=zeros(length(t0q),3);
QG60_perc_on=zeros(length(t0q),3);
for i=1:3
    LG20_perc_on(:,i)=100.*sum(LG20_var(:,:,i)>threshold(i),2)/n;
    LG60_perc_on(:,i)=100.*sum(LG60_var(:,:,i)>threshold(i),2)/n;
    QG60_perc_on(:,i)=100.*sum(QG60_var(:,:,i)>threshold(i),2)/n;
end

figure
tiledlayout(1,3)
for i=1:3
    nexttile
    hold on
    colororder(newcolors)
    plot(y0var(:,1),LG20_perc_on(:,i),'linewidth',2)
    plot(y1var(:,1),LG60_perc_on(:,i),'linewidth',2)
    plot(y2var(:,1),QG60_perc_on(:,i),'linewidth',2)
    if i==1
        legend(labels,'box','off')
    end
    title(species(i+1))
    xlabel('H_2O_2 dose (mM)')
    ylabel('Percentage on')
    ylim([0 100])
    set(gca,'box','off')
end


%% Rate dependence
Gmax1=5;
Gmax2=10;
rates1=logspace(log(5),log(0.05),21);
rates2=logspace(log(5),log(0.05),21);
tGG1=Gmax1./rates1;
tGG2=Gmax2./rates2;

allsim1=zeros(length(tGG1),1001,4);
allsim2=zeros(length(tGG1),1001,4);
for i=1:length(tGG1)
    tend1=tGG1(i);
    tend2=tGG2(i);
    tsim1=linspace(0,tend1,1001);
    tsim2=linspace(0,tend2,1001);

    [t1,y1] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km,kros,kpos,H,rp,rdp,tGG1(i),Gmax1,'linear',xT),tsim1,x0);
    [t2,y2] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km,kros,kpos,H,rp,rdp,tGG2(i),Gmax2,'linear',xT),tsim2,x0);
    allsim1(i,:,:)=y1;
    allsim2(i,:,:)=y2;
end

allEC501=nan(length(tGG1),3);
for i=1:length(tGG1)
    for j=1:3
        ind=find(allsim1(i,:,j)>=0.5,1);
        if any(ind)
            allEC501(i,j)=allsim1(i,ind,1);
        end
    end
end

allEC502=nan(length(tGG1),3);
for i=1:length(tGG1)
    for j=1:3
        ind=find(allsim2(i,:,j)>=0.5,1);
        if any(ind)
            allEC502(i,j)=allsim2(i,ind,1);
        end
    end
end

figure
hold on
plot(rates1,allEC501(:,3),'-o','linewidth',2)
plot(rates2,allEC502(:,3),'--x','linewidth',2)
set(gca,'xscale','log')
legend(sprintf('max. %d mM gradients',Gmax1),sprintf('max. %d mM gradients',Gmax2))

%% Rate dependence with threshold variability
Gmax1=10;
rates1=[10 5 2 1 0.5 0.2 0.1 0.05 0.025];
tGG1=Gmax1./rates1;

allsim1=zeros(length(tGG1),n,1001,4);
for i=1:length(tGG1)
    for j=1:n
        tend1=tGG1(i);
        tsim1=linspace(0,tend1,1001);
    
        [t1,y1] = ode23(@(t,y) BCR_FFL_Hill(t,y,p,Km_3(j,:),kros_3(j,:),kpos,H,rp,rdp,tGG1(i),Gmax1,'linear',xT),tsim1,x0);
        allsim1(i,j,:,:)=y1;
    end
end

rate_percon=zeros(length(tGG1),1001,3);
for i=1:length(tGG1)
    for j=1:3
        rate_percon(i,:,j)=100.*sum(allsim1(i,:,:,j+1)>threshold(j),2)/n;
    end
end

rate_EC50=nan(length(tGG1),3);
for i=1:length(tGG1)
    for j=1:3
        ind=find(rate_percon(i,:,j)>=50,1);
        if any(ind)
            rate_EC50(i,j)=allsim1(i,1,ind,1);
        end
    end
end

figure
colororder(newcolors)
hold on
plot(rates1,rate_EC50(:,1),'-','linewidth',2)
plot(rates1,rate_EC50(:,2),'-','linewidth',2)
plot(rates1,rate_EC50(:,3),'-','linewidth',2)
set(gca,'xscale','log')
legend(species(2:end),'box', 'off','location','nw')
xlabel('Gradient rate (mM min^{-1})')
ylabel('EC_{50%} (mM)')
xlim([min(rates1),max(rates1)])

%% functions
function xdot=BCR_FFL_Hill(t,x,p,Km,kros,kpos,H,rp,rdp,tG,Gmax,gradienttype,xT)

xdot=zeros(4,1);

if strcmp(gradienttype,'linear')
    a=Gmax/tG;
    if t<tG
        xdot(1)=a-p*x(1);
    else
        xdot(1)=-p*x(1);
    end
elseif strcmp(gradienttype,'quadratic')
    a=Gmax/tG^3;
    if t<tG
        xdot(1)=3*a*t^2-p*x(1);
    else
        xdot(1)=-p*x(1);
    end
else
    error('Invalid input type')
end

% Total protein concentrations
tCD79a=xT(1);
tSYK=xT(2);
tPLCy2=xT(3);

xtot=[tCD79a,tSYK,tPLCy2];
uX=[(tCD79a-x(2)), (tSYK-x(3)), (tPLCy2-x(4))]./xtot; % unphosphorylated protein conc.
pX=[x(2), x(3), x(4)]./xtot; % phosphorylated protein conc.
if any(uX<0)||any(pX<0)
    if any(uX<0)
        j=find(uX<0);
        uX(j)=0;
        pX(j)=xT(j);
    else
        j=find(pX<0);
        pX(j)=0;
        uX(j)=xT(j);
    end
end
if any(uX>1)||any(pX>1)
    if any(uX>1)
        j=find(uX>1);
        uX(j)=1;
        pX(j)=0;
    else
        j=find(pX>1);
        pX(j)=1;
        uX(j)=0;
    end
end

ROSinh = (1./(1+(x(1)./kros).^H));
SYKact = uX(2)*pX(2)^H(2)/(kpos.^H(2)+pX(2)^H(2));

xdot(2)=rp(1) * (uX(1)^H(1)/(Km(1)^H(1)+uX(1)^H(1))) - rdp(1) * ROSinh(1) * pX(1)^H(1)/(Km(2)^H(1)+pX(1)^H(1));
xdot(3)=rp(2)*(pX(1) * (uX(2)^H(2)/(Km(3)^H(2)+uX(2)^H(2))) + SYKact) - rdp(2) * ROSinh(2) * pX(2)^H(2)/(Km(4)^H(2)+pX(2)^H(2));
xdot(4)=rp(3)*pX(2)*uX(3)^H(3)/(Km(5)^H(3)+uX(3)^H(3)) - rdp(3) * ROSinh(3) * pX(3)^H(3)/(Km(6)^H(3)+pX(3)^H(3));

end
