%% this is for part a)

clear all
% endtime_value_J=[15,300,400,500];
endtime_value_J=[15,250,800,2000];
% endtime_value_h=[15,250,800,2000];
endtime_value_h=[15,250,800,2000];
% species_change=[2,34,45,97];
species_change=[2,34,45,97];
% endtime_value=[200,300,500,600];
%endtime_value=[15,20,20,50];

ph = uipanel('Parent',figure,'BorderType','none');  
ph.Title = 'Population vs. time (days), number of species=[2,34,45,97], .7 sparse';
ph.TitlePosition = 'centertop'; 
ph.FontSize = 12;
ph.FontWeight = 'bold';
pJ = uipanel('Parent',figure,'BorderType','none');
pJ.Title = 'No ex.-threshold, Population vs. time (days), number of species=[2,34,45,97], .7 sparse';
pJ.TitlePosition = 'centertop'; 
pJ.FontSize = 12;
pJ.FontWeight = 'bold';

for k = 1:4
    h(k) = subplot(2,2,k,'Parent',ph);
end

for k = 1:4
    J(k) = subplot(2,2,k,'Parent',pJ);
end

for k=1:4
    run_diff_NoS_sbplt_num(species_change(k),k,h(k),J(k),endtime_value_h(k),endtime_value_J(k))
end

load('forformula.mat')
load('Display_legend_arrayk1.mat')
load('Display_legend_arrayk2.mat')
load('Display_legend_arrayk3.mat')
load('Display_legend_arrayk4.mat')
load('Display_legend_array1k1.mat')
load('Display_legend_array1k2.mat')
load('Display_legend_array1k3.mat')
load('Display_legend_array1k4.mat')

ratio_mat=string(zeros(4));
for i=1:4
    ratio_mat(1,i)=sprintf('%d/%d',surv_arrayk1(i),total_amount_arrayk1(i));
    ratio_mat(2,i)=sprintf('%d/%d',surv_arrayk2(i),total_amount_arrayk2(i));
    ratio_mat(3,i)=sprintf('%d/%d',surv_arrayk3(i),total_amount_arrayk3(i));
    ratio_mat(4,i)=sprintf('%d/%d',surv_arrayk4(i),total_amount_arrayk4(i));
end

%%% CREATING THE TABLE
table=array2table(ratio_mat,'VariableNames',{'Prey','Predators',...
'PreyLeaningHybrids','PredatorLeaningHybrids'},'RowNames',...
{'2 species','34 species','45 species','97 species'})

%% this is for the first figure in part b) (this could take some time to run)

% [200,1500,900,400]

clear all
% endtime_value_J=[15,300,400,500];
% endtime_value_h=[200,1500,900,400];
endtime_value_h=[1500,1500,1500,1500];
% endtime_value_h=[15,250,800,2000];
endtime_value_J=[50,50,50,50];
% endtime_value=[200,300,500,600];
%endtime_value=[15,20,20,50];
% number_of_species=13;
number_of_species=97;
% density_array=[0,.33,.66,.99];
density_array1=linspace(0,0.428571428571429,4);
density_array2=linspace(0.571428571428571,1,4);

ph = uipanel('Parent',figure,'BorderType','none');  
ph.Title = 'Population vs. time (days), number of species=13, sparsity=linspace(0,1,8)';
ph.TitlePosition = 'centertop'; 
ph.FontSize = 12;
ph.FontWeight = 'bold';
pJ = uipanel('Parent',figure,'BorderType','none');
pJ.Title = 'No ex.-threshold, Population vs. time (days), number of species=13, sparsity=linspace(0,1,8)';
pJ.TitlePosition = 'centertop'; 
pJ.FontSize = 12;
pJ.FontWeight = 'bold';

for k = 1:4
    h(k) = subplot(2,2,k,'Parent',ph);
end

for k = 1:4
    J(k) = subplot(2,2,k,'Parent',pJ);
end

for k=1:4
    run_diff_spars_sbplt_num(number_of_species,k,h(k),J(k),endtime_value_h(k),endtime_value_J(k),density_array1(k))
end

load('forformula.mat')
load('Display_legend_arrayk1.mat')
load('Display_legend_arrayk2.mat')
load('Display_legend_arrayk3.mat')
load('Display_legend_arrayk4.mat')
ratio_mat=string(zeros(4));
for i=1:4
    ratio_mat(1,i)=sprintf('%d/%d=%f',surv_arrayk1(i),total_amount_arrayk1(i),surv_arrayk1(i)/total_amount_arrayk1(i));
    ratio_mat(2,i)=sprintf('%d/%d=%f',surv_arrayk2(i),total_amount_arrayk2(i),surv_arrayk2(i)/total_amount_arrayk2(i));
    ratio_mat(3,i)=sprintf('%d/%d=%f',surv_arrayk3(i),total_amount_arrayk3(i),surv_arrayk3(i)/total_amount_arrayk3(i));
    ratio_mat(4,i)=sprintf('%d/%d=%f',surv_arrayk4(i),total_amount_arrayk4(i),surv_arrayk4(i)/total_amount_arrayk4(i));
end
save('first_plot_data','A1','A2','A3','A4','surv_arrayk1','surv_arrayk2','surv_arrayk3',...
    'surv_arrayk4','total_amount_arrayk1','total_amount_arrayk2','total_amount_arrayk3',...
    'total_amount_arrayk4','ratio_mat')

%% this is for the second figure in part b)
% []
clear all
% endtime_value_J=[15,300,400,500];
endtime_value_h=[1500,1500,1500,1500];
% endtime_value_h=[15,250,800,2000];
endtime_value_J=[50,50,50,50];
% endtime_value=[200,300,500,600];
%endtime_value=[15,20,20,50];
% number_of_species=13;
number_of_species=97;
% density_array=[0,.33,.66,.99];
density_array1=linspace(0,0.428571428571429,4);
density_array2=linspace(0.571428571428571,1,4);

ph = uipanel('Parent',figure,'BorderType','none');  
ph.Title = 'Population vs. time (days), number of species=13, sparsity=linspace(0,1,8)';
ph.TitlePosition = 'centertop'; 
ph.FontSize = 12;
ph.FontWeight = 'bold';
pJ = uipanel('Parent',figure,'BorderType','none');
pJ.Title = 'No ex.-threshold, Population vs. time (days), number of species=13, sparsity=linspace(0,1,8)';
pJ.TitlePosition = 'centertop'; 
pJ.FontSize = 12;
pJ.FontWeight = 'bold';

for k = 1:4
    h(k) = subplot(2,2,k,'Parent',ph);
end

for k = 1:4
    J(k) = subplot(2,2,k,'Parent',pJ);
end

for k=1:4
    run_diff_spars_sbplt_num(number_of_species,k,h(k),J(k),endtime_value_h(k),endtime_value_J(k),density_array2(k))
end

load('forformula.mat')
load('Display_legend_arrayk1.mat')
load('Display_legend_arrayk2.mat')
load('Display_legend_arrayk3.mat')
load('Display_legend_arrayk4.mat')
ratio_mat=string(zeros(4));
for i=1:4
    ratio_mat(1,i)=sprintf('%d/%d=%f',surv_arrayk1(i),total_amount_arrayk1(i),surv_arrayk1(i)/total_amount_arrayk1(i));
    ratio_mat(2,i)=sprintf('%d/%d=%f',surv_arrayk2(i),total_amount_arrayk2(i),surv_arrayk2(i)/total_amount_arrayk2(i));
    ratio_mat(3,i)=sprintf('%d/%d=%f',surv_arrayk3(i),total_amount_arrayk3(i),surv_arrayk3(i)/total_amount_arrayk3(i));
    ratio_mat(4,i)=sprintf('%d/%d=%f',surv_arrayk4(i),total_amount_arrayk4(i),surv_arrayk4(i)/total_amount_arrayk4(i));
end
save('second_plot_data','A1','A2','A3','A4','surv_arrayk1','surv_arrayk2','surv_arrayk3',...
    'surv_arrayk4','total_amount_arrayk1','total_amount_arrayk2','total_amount_arrayk3',...
    'total_amount_arrayk4','ratio_mat')

%% also for part b)

load('first_plot_data')
super_ratio_mat=ratio_mat;
load('second_plot_data')
super_ratio_mat=[super_ratio_mat;ratio_mat];

%%% CREATING THE TABLE
table=array2table(super_ratio_mat,'VariableNames',{'Prey','Predators',...
'PreyLeaningHybrids','PredatorLeaningHybrids'},'RowNames',...
{'sparcity=0','sparcity=0.14','sparcity=0.28','sparcity=0.42','sparcity=0.57',...
'sparcity=0.71','sparcity=0.85','sparcity=1'})


%% also for part b)

clear all
number_of_species=2;

% A=1/200.*A; % NORMALIZING
A=[0 1/200; 1/300 0];y=[300 500];r=[1,-1];

%%% EXTINCTION RATE AND POINT AND BOOKKEEPING
extinction_point=1;
k_extinction=1;
extinct=zeros(size(y));
save('forformula.mat','extinction_point','number_of_species','A','k_extinction','r','extinct');

%%% DEFINING THE FUNCTIONS AND TIMES
dpop=@(t,y) r'.*(ones(number_of_species,1)-flipud(A*y)).*y;fff=@ggg5;
dpop1=@(t,y) r'.*(ones(number_of_species,1)-(A*y)).*y;
starttime=0;endtime=20;

%%% RUNNING, MODEL WITHOUT EXTINCTION THRESHOLD
[t,p]=ode45(fff,[starttime endtime],y');

%%% PLOTTING
figure
plot(t,p,'LineWidth',2)
title('two species model')
xlabel('time elapsed in days')
ylabel('specie population')

A=[0 1/200; 1/300 -1/500];

%%% EXTINCTION RATE AND POINT AND BOOKKEEPING
extinction_point=1;
k_extinction=1;
extinct=zeros(size(y));
save('forformula.mat','extinction_point','number_of_species','A','k_extinction','r','extinct');

%%% DEFINING THE FUNCTIONS AND TIMES
dpop=@(t,y) r'.*(ones(number_of_species,1)-flipud(A*y)).*y;fff=@ggg5;
dpop1=@(t,y) r'.*(ones(number_of_species,1)-(A*y)).*y;
starttime=0;endtime=20;

%%% RUNNING, MODEL WITHOUT EXTINCTION THRESHOLD
[t,p]=ode45(fff,[starttime endtime],y');

%%% PLOTTING
figure
plot(t,p,'LineWidth',2)
title('two species model')
xlabel('time elapsed in days')
ylabel('specie population')


%% also for part b)

%%% CREATING THE PHASE SPACE

A=[0 1/200; 1/300 0];y=[300 500];r=[1,-1];
period = 6.5357;

pred_prey_ode=@(t,y) r'.*(ones(size(A,1),1)-(A*y)).*y;

[S,F] = meshgrid(0:50:800,0:50:500);
dS = zeros(size(S));
dF = zeros(size(F));
for i = 1:size(S, 1)
    for r1 = 1:size(S, 2)
        deriv = pred_prey_ode(0, [S(i, r1); F(i, r1)]);
        dS(i,r1) = deriv(1);
        dF(i,r1) = deriv(2);
    end
end

figure
s(1)=subplot(3,1,1);
quiver(S,F,dS,dF,2);
xlim([0 inf])
ylim([0 inf])
title(s(1),'unstable fixed point, no carrying-capacity')
hold on;
len = 100;
for n = linspace(100, 400, 5)
sol2 = ode45(pred_prey_ode,[0 3*period], [400 n]);
sol2y = deval(sol2, linspace(0, 3*period, len * 5));
plot(sol2y(1,:),sol2y(2,:), '-g', 'LineWidth', 2);
end

syms a b
f = [1.*(1-a*A(1,1)-b*A(1,2)).*a == 0, -1.*(1-a*A(2,1)-b*A(2,2)).*b == 0];
[solx,soly] = solve(f,[a b]);
steady_states = [solx,soly]
scatter(steady_states(:,1), steady_states(:,2), 'or');



A=[0 1/200; 1/300 -1/500];

pred_prey_ode=@(t,y) r'.*(ones(size(A,1),1)-(A*y)).*y;

[S,F] = meshgrid(0:50:900,0:50:800);
dS = zeros(size(S));
dF = zeros(size(F));
for i = 1:size(S, 1)
    for r1 = 1:size(S, 2)
        deriv = pred_prey_ode(0, [S(i, r1); F(i, r1)]);
        dS(i,r1) = deriv(1);
        dF(i,r1) = deriv(2);
    end
end

s(2)=subplot(3,1,2);
quiver(S,F,dS,dF,2);
ylim([0,inf])
xlim([0,inf])
title(s(2),'stable fixed point, 1/500 carrying-capacity')
hold on;
len = 100;
for n = linspace(100, 800, 5)
sol2 = ode45(pred_prey_ode,[0 3*period], [n 400]);
sol2y = deval(sol2, linspace(0, 3*period, len * 5));
plot(sol2y(1,:),sol2y(2,:), '-g', 'LineWidth', 2);
end

syms x y
f = [1.*(1-a*A(1,1)-b*A(1,2)).*a == 0, -1.*(1-a*A(2,1)-b*A(2,2)).*b == 0];
[solx,soly] = solve(f,[a b]);
steady_states = [solx,soly]
scatter(steady_states(:,1), steady_states(:,2), 'or');



A=[0 1/200; 1/300 -1/100];

pred_prey_ode=@(t,y) r'.*(ones(size(A,1),1)-(A*y)).*y;

[S,F] = meshgrid(0:50:1500,0:50:600);
dS = zeros(size(S));
dF = zeros(size(F));
for i = 1:size(S, 1)
    for r1 = 1:size(S, 2)
        deriv = pred_prey_ode(0, [S(i, r1); F(i, r1)]);
        dS(i,r1) = deriv(1);
        dF(i,r1) = deriv(2);
    end
end

s(3)=subplot(3,1,3);
quiver(S,F,dS,dF,2);
ylim([0,inf])
xlim([0,inf])
title(s(3),'stable fixed point, 1/100 carrying-capacity')
hold on;
len = 100;

for n = linspace(100, 1500, 5)
sol2 = ode45(pred_prey_ode,[0 3*period], [n 400]);
sol2y = deval(sol2, linspace(0, 3*period, len * 5));
plot(sol2y(1,:),sol2y(2,:), '-g', 'LineWidth', 2);
end

syms x y
f = [1.*(1-a*A(1,1)-b*A(1,2)).*a == 0, -1.*(1-a*A(2,1)-b*A(2,2)).*b == 0];
[solx,soly] = solve(f,[a b]);
steady_states = [solx,soly]
scatter(steady_states(:,1), steady_states(:,2), 'or');

%%

%% for part c)

%%% EXTEND ODE TO INVASION OR LOSS

clear all  
format long
rng(5)
density=0.428571428571429;
blah=" sparse";
number_of_species=97;

% MAKING THE MATRIX A
%%% RANDOM MATRIX BETWEEN D AND C
c=2;d=.1;A=(c-d).*sprand(number_of_species,number_of_species,density);

%%% CAN MAKE ALMOST HALF THE ENTRIES NEGATIVE
Like_A=zeros(size(A));
[i_array,j_array,s_array] = find(A);% nonzero indeces
for i=1:length(i_array)
    Like_A(i_array(i),j_array(i))=(s_array(i)-.5.*(c-d)+d);
end
A=Like_A;

%%% DEFINING THE DIAGONAL ENTRIES
A=full(A);
diagg=2.*rand(size(A,1),1);
for i=1:size(A,1)
    A(i,i)=diagg(i);
end

%%% CHOOSE WHICH ARE PREDATORS AND WHICH ARE PREY
%%% ALSO DEFINING THE Rs
types=rand(1,number_of_species);
r=zeros(1,number_of_species);% PREALLOCATION
for i=1:number_of_species
    if types(i)<=(1/4)
        types(i)=1;% prey
        r(i)=rand();
    elseif types(i)>(1/4) && types(i)<=(2/4)
        types(i)=2;% predator
        r(i)=-rand();
    elseif types(i)>(2/4) && types(i)<=(3/4)
        types(i)=3;% hybrid with positive r
        r(i)=rand();
    elseif types(i)>(3/4) && types(i)<=1
        types(i)=4;% hybrid with negative r
        r(i)=-rand();
    end
end

%%% NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if types(i)==1% prey
        for j=1:number_of_species% who is affecting
            if types(j)==1% both >=0
                if A(i,j)==0 % ignore
                    A(j,i)=0;
                elseif A(i,j)>0 % competition
                    A(j,i)=rand();
                elseif A(i,j)<0% mistake in A
                    if rand()<.5% ignore
                        A(i,j)=0;
                        A(j,i)=0;
                    else% competition
                        A(i,j)=-A(i,j);
                        A(j,i)=rand();
                    end
                    
                end
            elseif types(j)==2% both >=0
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(j,i)=rand();
                elseif A(i,j)<0% mistake in A
                    if rand()<.5% ignore
                        A(i,j)=0;
                        A(j,i)=0;
                    else% competition
                        A(i,j)=-A(i,j);
                        A(j,i)=rand();
                    end
                end
            elseif types(j)==3% positive r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                    end
                elseif A(i,j)<0% mistake in A
                    A(i,j)=-A(i,j);
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                    end
                end
            elseif types(j)==4% negative r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                    end
                elseif A(i,j)<0
                    A(i,j)=-A(i,j);
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                    end
                end
            end
        end
    elseif types(i)==2% predator
        for j=1:number_of_species% who is affecting
            if types(j)==1% both >=0
                if A(i,j)==0 % ignore
                    A(j,i)=0;
                elseif A(i,j)>0 % eating
                    A(j,i)=rand();
                elseif A(i,j)<0% mistake in A
                    if rand()<.5% ignore
                        A(i,j)=0;
                        A(j,i)=0;
                    else% eating
                        A(i,j)=-A(i,j);
                        A(j,i)=rand();
                    end
                end
            elseif types(j)==2% opposite signs
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% eating
                    A(j,i)=-rand();
                elseif A(i,j)<0% being eaten
                    A(j,i)=rand();
                end
            elseif types(j)==3% positive r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% eating
                    A(j,i)=rand();
                elseif A(i,j)<0% being eaten
                    A(j,i)=-rand();
                end
            elseif types(j)==4% negative r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% eating
                    A(j,i)=-rand();
                elseif A(i,j)<0% being eaten
                    A(j,i)=rand();
                end
            end
        end
    elseif types(i)==3% positive r hybrid
        for j=1:number_of_species% who is affecting
            if types(j)==1
                if A(i,j)==0 % ignore
                    A(j,i)=0;
                elseif A(i,j)>0 % competing
                    A(j,i)=rand();
                elseif A(i,j)<0% eating
                    A(j,i)=rand();
                end
            elseif types(j)==2
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(j,i)=rand();
                elseif A(i,j)<0% eating
                    A(j,i)=-rand();
                end
            elseif types(j)==3% positive r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                    end
                elseif A(i,j)<0% eating
                    A(i,j)=rand();
                end
            elseif types(j)==4% negative r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                    end
                elseif A(i,j)<0% eating
                    A(j,i)=-rand();
                end
            end
        end
    elseif types(i)==4% negative r hybrid
        for j=1:number_of_species% who is affecting
            if types(j)==1
                if A(i,j)==0 % ignore
                    A(j,i)=0;
                elseif A(i,j)>0 % competing
                    A(j,i)=rand();
                elseif A(i,j)<0% eating
                    A(j,i)=rand();
                end
            elseif types(j)==2
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(j,i)=-rand();
                elseif A(i,j)<0% eating
                    A(j,i)=rand();
                end
            elseif types(j)==3% positive r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(i,j)=rand();
                elseif A(i,j)<0% eating
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                    end
                end
            elseif types(j)==4% negative r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(j,i)=-rand();
                elseif A(i,j)<0% eating
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                    end
                end
            end
        end
    end
end

%%% INITIAL POPULATION
l=50.0;ll=30;y=(l-ll).*rand(1,number_of_species)+ll;

A=1/200.*A; % NORMALIZING

%%% EXTINCTION RATE AND POINT AND BOOKKEEPING
extinction_point=1;
k_extinction=.1;
extinct=zeros(size(y));
save('forformula.mat','extinction_point','number_of_species','A','k_extinction','r','extinct');

%%% DEFINING THE FUNCTIONS AND TIMES
fff=@ggg5;
starttime=0;endtime=150;
endtime4=550;
endtime3=350;
endtime2=250;
% endtime4=50;
% endtime3=50;
% endtime2=50;
l=50;

%%%%%% PART NUMBER 1

%%% RUNNING, MODEL WITH EXTINCTION THRESHOLD
[t,p]=ode45(fff,[starttime endtime],y');
k=1;

%%% MAKING THE LEGEND ENTRIES
format short
load('forformula.mat')
%extinct_values{k}=extinct;
Display_legend=string(zeros(1,number_of_species));
ext_array=zeros(1,4);% prey,pred,PreyLH,PredLH
ext_ind_array=cell(1,4);% ROW ARRAY, ALL THE EXTINCT INDICES SPECIES
surv_array=zeros(1,4);% prey,pred,PreyLH,PredLH
surv_ind_array=cell(1,4);% ROW ARRAY, ALL THE SURVIVING INDICES SPECIES
%%% RECORDING HOW MANY SURVIVED AND WENT EXTINCT
for i=1:number_of_species
    if types(i)==1
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,1)=ext_array(1,1)+1;
            ext_ind_array{1}=[ext_ind_array{1} i];
        else
            Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p(length(p),i) );
            surv_array(1,1)=surv_array(1,1)+1;
            surv_ind_array{1}=[surv_ind_array{1} i];
        end
    elseif types(i)==2
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,2)=ext_array(1,2)+1;
            ext_ind_array{2}=[ext_ind_array{2} i];
        else
            Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p(length(p),i) );
            surv_array(1,2)=surv_array(1,2)+1;
            surv_ind_array{2}=[surv_ind_array{2} i];
        end
    elseif types(i)==3
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,3)=ext_array(1,3)+1;
            ext_ind_array{3}=[ext_ind_array{3} i];
        else
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p(length(p),i) );
            surv_array(1,3)=surv_array(1,3)+1;
            surv_ind_array{3}=[surv_ind_array{3} i];
        end
    elseif types(i)==4
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,4)=ext_array(1,4)+1;
            ext_ind_array{4}=[ext_ind_array{4} i];
        else
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p(length(p),i) );
            surv_array(1,4)=surv_array(1,4)+1;
            surv_ind_array{4}=[surv_ind_array{4} i];
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    A1=A;
    ext_arrayk1=ext_array;
    ext_ind_arrayk1=ext_ind_array;
    surv_arrayk1=surv_array;
    surv_ind_arrayk1=surv_ind_array;
    total_amount_arrayk1=ext_arrayk1+surv_arrayk1;
    Display_legendk1=Display_legend;
    save('Display_legend_arrayk1.mat','Display_legendk1','total_amount_arrayk1',...
        'ext_arrayk1','surv_arrayk1','A1','ext_ind_arrayk1','surv_ind_arrayk1')
elseif k==2
    A2=A;
    ext_arrayk2=ext_array;
    ext_ind_arrayk2=ext_ind_array;
    surv_arrayk2=surv_array;
    surv_ind_arrayk2=surv_ind_array;
    total_amount_arrayk2=ext_arrayk2+surv_arrayk2;
    Display_legendk2=Display_legend;
    save('Display_legend_arrayk2.mat','Display_legendk2','total_amount_arrayk2',...
        'ext_arrayk2','surv_arrayk2','A2','ext_ind_arrayk2','surv_ind_arrayk2')
elseif k==3
    A3=A;
    ext_arrayk3=ext_array;
    surv_arrayk3=surv_array;
    ext_ind_arrayk3=ext_ind_array;
    surv_ind_arrayk3=surv_ind_array;
    total_amount_arrayk3=ext_arrayk3+surv_arrayk3;
    Display_legendk3=Display_legend;
    save('Display_legend_arrayk3.mat','Display_legendk3','total_amount_arrayk3',...
        'ext_arrayk3','surv_arrayk3','A3','ext_ind_arrayk3','surv_ind_arrayk3')
elseif k==4
    A4=A;
    ext_arrayk4=ext_array;
    surv_arrayk4=surv_array;
    ext_ind_arrayk4=ext_ind_array;
    surv_ind_arrayk4=surv_ind_array;
    total_amount_arrayk4=ext_arrayk4+surv_arrayk4;
    Display_legendk4=Display_legend;
    save('Display_legend_arrayk4.mat','Display_legendk4','total_amount_arrayk4',...
        'ext_arrayk4','surv_arrayk4','A4','ext_ind_arrayk4','surv_ind_arrayk4')
end



%%% SEE WHICH ARE EXTINCT
extinction_index_array=[];% has all the indexes that correspond to the extinct species
for i=1:97
    if p(size(p,1),i)<extinction_point
        extinction_index_array=[extinction_index_array i];
    end
end

%%% REPLACE THEM WITH INVADERS (WE'LL KEEP IT THE SAME SIZE)
types_new=rand(1,length(extinction_index_array));
r_new=zeros(1,length(extinction_index_array));% PREALLOCATION
for i=1:length(extinction_index_array)
    %DECIDE WHAT TYPES THE INVADERS ARE AND DEFINE THEIR Rs
    if types_new(i)<=(1/4)
        types_new(i)=1;% prey
        r_new(i)=rand();
    elseif types_new(i)>(1/4) && types_new(i)<=(2/4)
        types_new(i)=2;% predator
        r_new(i)=-rand();
    elseif types_new(i)>(2/4) && types_new(i)<=(3/4)
        types_new(i)=3;% hybrid with positive r
        r_new(i)=rand();
    elseif types_new(i)>(3/4) && types_new(i)<=1
        types_new(i)=4;% hybrid with negative r
        r_new(i)=-rand();
    end
    % REPLACE TYPES AND Rs
    types(extinction_index_array(i))=types_new(i);
    r(extinction_index_array(i))=r_new(i);
end

%%% REDO A MATRIX, NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if ismember(i,extinction_index_array)% BUT ONLY IF THE SPECIES IS EXTINCT
        if types(i)==1% prey
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competition
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end

                    end
                elseif types(j)==2% both >=0
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    elseif A(i,j)<0% mistake in A
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    elseif A(i,j)<0
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    end
                end
            end
        elseif types(i)==2% predator
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% eating
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==2% opposite signs
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=-rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                end
            end
        elseif types(i)==3% positive r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                        end
                    elseif A(i,j)<0% eating
                        A(i,j)=rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                end
            end
        elseif types(i)==4% negative r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(i,j)=rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                end
            end
        end
    end
end

%%% UPDATE THE STARTING POPULATION NUMBERS
y=p(size(p,1),:);
l=50.0;ll=30;
for i=1:length(extinction_index_array)
    y(extinction_index_array(i))=(l-ll).*rand()+ll;
end


%%%%%%%% PART NUMBER 2

%%% RUNNING, MODEL WITH INVASION/LOSS
[t1,p1]=ode45(fff,[endtime endtime2],y');
k=2;


%%% MAKING THE LEGEND ENTRIES
format short
load('forformula.mat')
%extinct_values{k}=extinct;
Display_legend=string(zeros(1,number_of_species));
ext_array=zeros(1,4);% prey,pred,PreyLH,PredLH
ext_ind_array=cell(1,4);% ROW ARRAY, ALL THE EXTINCT INDICES SPECIES
surv_array=zeros(1,4);% prey,pred,PreyLH,PredLH
surv_ind_array=cell(1,4);% ROW ARRAY, ALL THE SURVIVING INDICES SPECIES
%%% RECORDING HOW MANY SURVIVED AND WENT EXTINCT
for i=1:number_of_species
    if types(i)==1
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p1(length(p1),i) );
            ext_array(1,1)=ext_array(1,1)+1;
            ext_ind_array{1}=[ext_ind_array{1} i];
        else
            Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p1(length(p1),i) );
            surv_array(1,1)=surv_array(1,1)+1;
            surv_ind_array{1}=[surv_ind_array{1} i];
        end
    elseif types(i)==2
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p1(length(p1),i) );
            ext_array(1,2)=ext_array(1,2)+1;
            ext_ind_array{2}=[ext_ind_array{2} i];
        else
            Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p1(length(p1),i) );
            surv_array(1,2)=surv_array(1,2)+1;
            surv_ind_array{2}=[surv_ind_array{2} i];
        end
    elseif types(i)==3
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p1(length(p1),i) );
            ext_array(1,3)=ext_array(1,3)+1;
            ext_ind_array{3}=[ext_ind_array{3} i];
        else
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p1(length(p1),i) );
            surv_array(1,3)=surv_array(1,3)+1;
            surv_ind_array{3}=[surv_ind_array{3} i];
        end
    elseif types(i)==4
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p1(length(p1),i) );
            ext_array(1,4)=ext_array(1,4)+1;
            ext_ind_array{4}=[ext_ind_array{4} i];
        else
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p1(length(p1),i) );
            surv_array(1,4)=surv_array(1,4)+1;
            surv_ind_array{4}=[surv_ind_array{4} i];
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    A1=A;
    ext_arrayk1=ext_array;
    ext_ind_arrayk1=ext_ind_array;
    surv_arrayk1=surv_array;
    surv_ind_arrayk1=surv_ind_array;
    total_amount_arrayk1=ext_arrayk1+surv_arrayk1;
    Display_legendk1=Display_legend;
    save('Display_legend_arrayk1.mat','Display_legendk1','total_amount_arrayk1',...
        'ext_arrayk1','surv_arrayk1','A1','ext_ind_arrayk1','surv_ind_arrayk1')
elseif k==2
    A2=A;
    ext_arrayk2=ext_array;
    ext_ind_arrayk2=ext_ind_array;
    surv_arrayk2=surv_array;
    surv_ind_arrayk2=surv_ind_array;
    total_amount_arrayk2=ext_arrayk2+surv_arrayk2;
    Display_legendk2=Display_legend;
    save('Display_legend_arrayk2.mat','Display_legendk2','total_amount_arrayk2',...
        'ext_arrayk2','surv_arrayk2','A2','ext_ind_arrayk2','surv_ind_arrayk2')
elseif k==3
    A3=A;
    ext_arrayk3=ext_array;
    surv_arrayk3=surv_array;
    ext_ind_arrayk3=ext_ind_array;
    surv_ind_arrayk3=surv_ind_array;
    total_amount_arrayk3=ext_arrayk3+surv_arrayk3;
    Display_legendk3=Display_legend;
    save('Display_legend_arrayk3.mat','Display_legendk3','total_amount_arrayk3',...
        'ext_arrayk3','surv_arrayk3','A3','ext_ind_arrayk3','surv_ind_arrayk3')
elseif k==4
    A4=A;
    ext_arrayk4=ext_array;
    surv_arrayk4=surv_array;
    ext_ind_arrayk4=ext_ind_array;
    surv_ind_arrayk4=surv_ind_array;
    total_amount_arrayk4=ext_arrayk4+surv_arrayk4;
    Display_legendk4=Display_legend;
    save('Display_legend_arrayk4.mat','Display_legendk4','total_amount_arrayk4',...
        'ext_arrayk4','surv_arrayk4','A4','ext_ind_arrayk4','surv_ind_arrayk4')
end



%%% SEE WHICH ARE EXTINCT
extinction_index_array=[];% has all the indexes that correspond to the extinct species
for i=1:97
    if p1(size(p1,1),i)<extinction_point
        extinction_index_array=[extinction_index_array i];
    end
end

%%% REPLACE THEM WITH INVADERS (WE'LL KEEP IT THE SAME SIZE)
types_new=rand(1,length(extinction_index_array));
r_new=zeros(1,length(extinction_index_array));% PREALLOCATION
for i=1:length(extinction_index_array)
    %DECIDE WHAT TYPES THE INVADERS ARE AND DEFINE THEIR Rs
    if types_new(i)<=(1/4)
        types_new(i)=1;% prey
        r_new(i)=rand();
    elseif types_new(i)>(1/4) && types_new(i)<=(2/4)
        types_new(i)=2;% predator
        r_new(i)=-rand();
    elseif types_new(i)>(2/4) && types_new(i)<=(3/4)
        types_new(i)=3;% hybrid with positive r
        r_new(i)=rand();
    elseif types_new(i)>(3/4) && types_new(i)<=1
        types_new(i)=4;% hybrid with negative r
        r_new(i)=-rand();
    end
    % REPLACE TYPES AND Rs
    types(extinction_index_array(i))=types_new(i);
    r(extinction_index_array(i))=r_new(i);
end

%%% REDO A MATRIX, NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if ismember(i,extinction_index_array)% BUT ONLY IF THE SPECIES IS EXTINCT
        if types(i)==1% prey
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competition
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end

                    end
                elseif types(j)==2% both >=0
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    elseif A(i,j)<0% mistake in A
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    elseif A(i,j)<0
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    end
                end
            end
        elseif types(i)==2% predator
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% eating
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==2% opposite signs
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=-rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                end
            end
        elseif types(i)==3% positive r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                        end
                    elseif A(i,j)<0% eating
                        A(i,j)=rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                end
            end
        elseif types(i)==4% negative r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(i,j)=rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                end
            end
        end
    end
end

%%% UPDATE THE STARTING POPULATION NUMBERS
y=p1(size(p1,1),:);
l1=200.0;ll=30;
for i=1:length(extinction_index_array)
    y(extinction_index_array(i))=(l1-ll).*rand()+ll;
end

%%%%%%%%% PART NUMBER 3


%%% RUNNING, MODEL WITH INVASION/LOSS

[t2,p2]=ode45(fff,[250 350],y');
k=3;

%%% MAKING THE LEGEND ENTRIES
format short
load('forformula.mat')
%extinct_values{k}=extinct;
Display_legend=string(zeros(1,number_of_species));
ext_array=zeros(1,4);% prey,pred,PreyLH,PredLH
ext_ind_array=cell(1,4);% ROW ARRAY, ALL THE EXTINCT INDICES SPECIES
surv_array=zeros(1,4);% prey,pred,PreyLH,PredLH
surv_ind_array=cell(1,4);% ROW ARRAY, ALL THE SURVIVING INDICES SPECIES
%%% RECORDING HOW MANY SURVIVED AND WENT EXTINCT
for i=1:number_of_species
    if types(i)==1
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p2(length(p2),i) );
            ext_array(1,1)=ext_array(1,1)+1;
            ext_ind_array{1}=[ext_ind_array{1} i];
        else
            Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p2(length(p2),i) );
            surv_array(1,1)=surv_array(1,1)+1;
            surv_ind_array{1}=[surv_ind_array{1} i];
        end
    elseif types(i)==2
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p2(length(p2),i) );
            ext_array(1,2)=ext_array(1,2)+1;
            ext_ind_array{2}=[ext_ind_array{2} i];
        else
            Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p2(length(p2),i) );
            surv_array(1,2)=surv_array(1,2)+1;
            surv_ind_array{2}=[surv_ind_array{2} i];
        end
    elseif types(i)==3
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p2(length(p2),i) );
            ext_array(1,3)=ext_array(1,3)+1;
            ext_ind_array{3}=[ext_ind_array{3} i];
        else
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p2(length(p2),i) );
            surv_array(1,3)=surv_array(1,3)+1;
            surv_ind_array{3}=[surv_ind_array{3} i];
        end
    elseif types(i)==4
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p2(length(p2),i) );
            ext_array(1,4)=ext_array(1,4)+1;
            ext_ind_array{4}=[ext_ind_array{4} i];
        else
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p2(length(p2),i) );
            surv_array(1,4)=surv_array(1,4)+1;
            surv_ind_array{4}=[surv_ind_array{4} i];
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    A1=A;
    ext_arrayk1=ext_array;
    ext_ind_arrayk1=ext_ind_array;
    surv_arrayk1=surv_array;
    surv_ind_arrayk1=surv_ind_array;
    total_amount_arrayk1=ext_arrayk1+surv_arrayk1;
    Display_legendk1=Display_legend;
    save('Display_legend_arrayk1.mat','Display_legendk1','total_amount_arrayk1',...
        'ext_arrayk1','surv_arrayk1','A1','ext_ind_arrayk1','surv_ind_arrayk1')
elseif k==2
    A2=A;
    ext_arrayk2=ext_array;
    ext_ind_arrayk2=ext_ind_array;
    surv_arrayk2=surv_array;
    surv_ind_arrayk2=surv_ind_array;
    total_amount_arrayk2=ext_arrayk2+surv_arrayk2;
    Display_legendk2=Display_legend;
    save('Display_legend_arrayk2.mat','Display_legendk2','total_amount_arrayk2',...
        'ext_arrayk2','surv_arrayk2','A2','ext_ind_arrayk2','surv_ind_arrayk2')
elseif k==3
    A3=A;
    ext_arrayk3=ext_array;
    surv_arrayk3=surv_array;
    ext_ind_arrayk3=ext_ind_array;
    surv_ind_arrayk3=surv_ind_array;
    total_amount_arrayk3=ext_arrayk3+surv_arrayk3;
    Display_legendk3=Display_legend;
    save('Display_legend_arrayk3.mat','Display_legendk3','total_amount_arrayk3',...
        'ext_arrayk3','surv_arrayk3','A3','ext_ind_arrayk3','surv_ind_arrayk3')
elseif k==4
    A4=A;
    ext_arrayk4=ext_array;
    surv_arrayk4=surv_array;
    ext_ind_arrayk4=ext_ind_array;
    surv_ind_arrayk4=surv_ind_array;
    total_amount_arrayk4=ext_arrayk4+surv_arrayk4;
    Display_legendk4=Display_legend;
    save('Display_legend_arrayk4.mat','Display_legendk4','total_amount_arrayk4',...
        'ext_arrayk4','surv_arrayk4','A4','ext_ind_arrayk4','surv_ind_arrayk4')
end



%%% SEE WHICH ARE EXTINCT
extinction_index_array=[];% has all the indexes that correspond to the extinct species
for i=1:97
    if p2(size(p2,1),i)<extinction_point
        extinction_index_array=[extinction_index_array i];
    end
end

%%% REPLACE THEM WITH INVADERS (WE'LL KEEP IT THE SAME SIZE)
types_new=rand(1,length(extinction_index_array));
r_new=zeros(1,length(extinction_index_array));% PREALLOCATION
for i=1:length(extinction_index_array)
    %DECIDE WHAT TYPES THE INVADERS ARE AND DEFINE THEIR Rs
    if types_new(i)<=(1/4)
        types_new(i)=1;% prey
        r_new(i)=rand();
    elseif types_new(i)>(1/4) && types_new(i)<=(2/4)
        types_new(i)=2;% predator
        r_new(i)=-rand();
    elseif types_new(i)>(2/4) && types_new(i)<=(3/4)
        types_new(i)=3;% hybrid with positive r
        r_new(i)=rand();
    elseif types_new(i)>(3/4) && types_new(i)<=1
        types_new(i)=4;% hybrid with negative r
        r_new(i)=-rand();
    end
    % REPLACE TYPES AND Rs
    types(extinction_index_array(i))=types_new(i);
    r(extinction_index_array(i))=r_new(i);
end

%%% REDO A MATRIX, NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if ismember(i,extinction_index_array)% BUT ONLY IF THE SPECIES IS EXTINCT
        if types(i)==1% prey
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competition
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end

                    end
                elseif types(j)==2% both >=0
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    elseif A(i,j)<0% mistake in A
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    elseif A(i,j)<0
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    end
                end
            end
        elseif types(i)==2% predator
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% eating
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==2% opposite signs
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=-rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                end
            end
        elseif types(i)==3% positive r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                        end
                    elseif A(i,j)<0% eating
                        A(i,j)=rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                end
            end
        elseif types(i)==4% negative r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(i,j)=rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                end
            end
        end
    end
end

%%% UPDATE THE STARTING POPULATION NUMBERS
y=p2(size(p2,1),:);
l=300.0;ll=30;
for i=1:length(extinction_index_array)
    y(extinction_index_array(i))=(l-ll).*rand()+ll;
end


%%%%%%%%% PART NUMBER 4


%%% RUNNING, MODEL WITH INVASION/LOSS

[t3,p3]=ode45(fff,[350 550],y');
k=4;

%%% MAKING THE LEGEND ENTRIES
format short
load('forformula.mat')
%extinct_values{k}=extinct;
Display_legend=string(zeros(1,number_of_species));
ext_array=zeros(1,4);% prey,pred,PreyLH,PredLH
ext_ind_array=cell(1,4);% ROW ARRAY, ALL THE EXTINCT INDICES SPECIES
surv_array=zeros(1,4);% prey,pred,PreyLH,PredLH
surv_ind_array=cell(1,4);% ROW ARRAY, ALL THE SURVIVING INDICES SPECIES
%%% RECORDING HOW MANY SURVIVED AND WENT EXTINCT
for i=1:number_of_species
    if types(i)==1
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p3(length(p3),i) );
            ext_array(1,1)=ext_array(1,1)+1;
            ext_ind_array{1}=[ext_ind_array{1} i];
        else
            Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p3(length(p3),i) );
            surv_array(1,1)=surv_array(1,1)+1;
            surv_ind_array{1}=[surv_ind_array{1} i];
        end
    elseif types(i)==2
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p3(length(p3),i) );
            ext_array(1,2)=ext_array(1,2)+1;
            ext_ind_array{2}=[ext_ind_array{2} i];
        else
            Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p3(length(p3),i) );
            surv_array(1,2)=surv_array(1,2)+1;
            surv_ind_array{2}=[surv_ind_array{2} i];
        end
    elseif types(i)==3
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p3(length(p3),i) );
            ext_array(1,3)=ext_array(1,3)+1;
            ext_ind_array{3}=[ext_ind_array{3} i];
        else
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p3(length(p3),i) );
            surv_array(1,3)=surv_array(1,3)+1;
            surv_ind_array{3}=[surv_ind_array{3} i];
        end
    elseif types(i)==4
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p3(length(p3),i) );
            ext_array(1,4)=ext_array(1,4)+1;
            ext_ind_array{4}=[ext_ind_array{4} i];
        else
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p3(length(p3),i) );
            surv_array(1,4)=surv_array(1,4)+1;
            surv_ind_array{4}=[surv_ind_array{4} i];
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    A1=A;
    ext_arrayk1=ext_array;
    ext_ind_arrayk1=ext_ind_array;
    surv_arrayk1=surv_array;
    surv_ind_arrayk1=surv_ind_array;
    total_amount_arrayk1=ext_arrayk1+surv_arrayk1;
    Display_legendk1=Display_legend;
    save('Display_legend_arrayk1.mat','Display_legendk1','total_amount_arrayk1',...
        'ext_arrayk1','surv_arrayk1','A1','ext_ind_arrayk1','surv_ind_arrayk1')
elseif k==2
    A2=A;
    ext_arrayk2=ext_array;
    ext_ind_arrayk2=ext_ind_array;
    surv_arrayk2=surv_array;
    surv_ind_arrayk2=surv_ind_array;
    total_amount_arrayk2=ext_arrayk2+surv_arrayk2;
    Display_legendk2=Display_legend;
    save('Display_legend_arrayk2.mat','Display_legendk2','total_amount_arrayk2',...
        'ext_arrayk2','surv_arrayk2','A2','ext_ind_arrayk2','surv_ind_arrayk2')
elseif k==3
    A3=A;
    ext_arrayk3=ext_array;
    surv_arrayk3=surv_array;
    ext_ind_arrayk3=ext_ind_array;
    surv_ind_arrayk3=surv_ind_array;
    total_amount_arrayk3=ext_arrayk3+surv_arrayk3;
    Display_legendk3=Display_legend;
    save('Display_legend_arrayk3.mat','Display_legendk3','total_amount_arrayk3',...
        'ext_arrayk3','surv_arrayk3','A3','ext_ind_arrayk3','surv_ind_arrayk3')
elseif k==4
    A4=A;
    ext_arrayk4=ext_array;
    surv_arrayk4=surv_array;
    ext_ind_arrayk4=ext_ind_array;
    surv_ind_arrayk4=surv_ind_array;
    total_amount_arrayk4=ext_arrayk4+surv_arrayk4;
    Display_legendk4=Display_legend;
    save('Display_legend_arrayk4.mat','Display_legendk4','total_amount_arrayk4',...
        'ext_arrayk4','surv_arrayk4','A4','ext_ind_arrayk4','surv_ind_arrayk4')
end



%%% SEE WHICH ARE EXTINCT
extinction_index_array=[];% has all the indexes that correspond to the extinct species
for i=1:97
    if p3(size(p3,1),i)<extinction_point
        extinction_index_array=[extinction_index_array i];
    end
end

%%% REPLACE THEM WITH INVADERS (WE'LL KEEP IT THE SAME SIZE)
types_new=rand(1,length(extinction_index_array));
r_new=zeros(1,length(extinction_index_array));% PREALLOCATION
for i=1:length(extinction_index_array)
    %DECIDE WHAT TYPES THE INVADERS ARE AND DEFINE THEIR Rs
    if types_new(i)<=(1/4)
        types_new(i)=1;% prey
        r_new(i)=rand();
    elseif types_new(i)>(1/4) && types_new(i)<=(2/4)
        types_new(i)=2;% predator
        r_new(i)=-rand();
    elseif types_new(i)>(2/4) && types_new(i)<=(3/4)
        types_new(i)=3;% hybrid with positive r
        r_new(i)=rand();
    elseif types_new(i)>(3/4) && types_new(i)<=1
        types_new(i)=4;% hybrid with negative r
        r_new(i)=-rand();
    end
    % REPLACE TYPES AND Rs
    types(extinction_index_array(i))=types_new(i);
    r(extinction_index_array(i))=r_new(i);
end

%%% REDO A MATRIX, NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if ismember(i,extinction_index_array)% BUT ONLY IF THE SPECIES IS EXTINCT
        if types(i)==1% prey
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competition
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end

                    end
                elseif types(j)==2% both >=0
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    elseif A(i,j)<0% mistake in A
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    elseif A(i,j)<0
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    end
                end
            end
        elseif types(i)==2% predator
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% eating
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==2% opposite signs
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=-rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                end
            end
        elseif types(i)==3% positive r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                        end
                    elseif A(i,j)<0% eating
                        A(i,j)=rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                end
            end
        elseif types(i)==4% negative r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(i,j)=rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                end
            end
        end
    end
end

%%% UPDATE THE STARTING POPULATION NUMBERS
y=p3(size(p3,1),:);
l=300.0;ll=30;
for i=1:length(extinction_index_array)
    y(extinction_index_array(i))=(l-ll).*rand()+ll;
end



%%% PLOTTING
figure
plot(t,p,'LineWidth',2)
hold on
plot(t1,p1,'LineWidth',2)
hold on
plot(t2,p2,'LineWidth',2)
hold on
plot(t3,p3,'LineWidth',2)
hold off
title('Model with invasion and loss')
xlabel('time elapsed in days')
ylabel('specie population')

load('Display_legend_arrayk1.mat')
load('Display_legend_arrayk2.mat')
load('Display_legend_arrayk3.mat')
load('Display_legend_arrayk4.mat')

%%% CHECKING IF A SPECIES i SURVIVED THEM ALL
survived_all=cell(1,4);
for i=1:number_of_species
    if ismember(i,surv_ind_arrayk1{1}) && ismember(i,surv_ind_arrayk2{1}) &&...
            ismember(i,surv_ind_arrayk3{1}) && ismember(i,surv_ind_arrayk4{1})
        survived_all{1}=[survived_all{1} i];
    elseif ismember(i,surv_ind_arrayk1{2}) && ismember(i,surv_ind_arrayk2{2}) &&...
            ismember(i,surv_ind_arrayk3{2}) && ismember(i,surv_ind_arrayk4{2})
        survived_all{2}=[survived_all{2} i];
    elseif ismember(i,surv_ind_arrayk1{3}) && ismember(i,surv_ind_arrayk2{3}) &&...
            ismember(i,surv_ind_arrayk3{3}) && ismember(i,surv_ind_arrayk4{3})
        survived_all{3}=[survived_all{3} i];
    elseif ismember(i,surv_ind_arrayk1{4}) && ismember(i,surv_ind_arrayk2{4}) &&...
            ismember(i,surv_ind_arrayk3{4}) && ismember(i,surv_ind_arrayk4{4})
        survived_all{4}=[survived_all{4} i];
    end
end

%% for part c)

ratio_mat=string(zeros(4));
for i=1:4
    ratio_mat(1,i)=sprintf('%d/%d=%f',surv_arrayk1(i),total_amount_arrayk1(i),surv_arrayk1(i)/total_amount_arrayk1(i));
    ratio_mat(2,i)=sprintf('%d/%d=%f',surv_arrayk2(i),total_amount_arrayk2(i),surv_arrayk2(i)/total_amount_arrayk2(i));
    ratio_mat(3,i)=sprintf('%d/%d=%f',surv_arrayk3(i),total_amount_arrayk3(i),surv_arrayk3(i)/total_amount_arrayk3(i));
    ratio_mat(4,i)=sprintf('%d/%d=%f',surv_arrayk4(i),total_amount_arrayk4(i),surv_arrayk4(i)/total_amount_arrayk4(i));
end

%%% CREATING THE TABLE
survival_rates=array2table(ratio_mat,'VariableNames',{'Prey','Predators',...
'PreyLeaningHybrids','PredatorLeaningHybrids'},'RowNames',...
{'first simulation','second simulation','third simulation','fourth simulation'})


%% for part c), second figure

%%% EXTEND ODE TO INVASION OR LOSS

clear all  
format long
rng(5)
density=0.428571428571429;
blah=" sparse";
number_of_species=97;

% MAKING THE MATRIX A
%%% RANDOM MATRIX BETWEEN D AND C
c=2;d=.1;A=(c-d).*sprand(number_of_species,number_of_species,density);

%%% CAN MAKE ALMOST HALF THE ENTRIES NEGATIVE
Like_A=zeros(size(A));
[i_array,j_array,s_array] = find(A);% nonzero indeces
for i=1:length(i_array)
    Like_A(i_array(i),j_array(i))=(s_array(i)-.5.*(c-d)+d);
end
A=Like_A;

%%% DEFINING THE DIAGONAL ENTRIES
A=full(A);
diagg=2.*rand(size(A,1),1);
for i=1:size(A,1)
    A(i,i)=diagg(i);
end

%%% CHOOSE WHICH ARE PREDATORS AND WHICH ARE PREY
%%% ALSO DEFINING THE Rs
types=rand(1,number_of_species);
r=zeros(1,number_of_species);% PREALLOCATION
for i=1:number_of_species
    if types(i)<=(1/4)
        types(i)=1;% prey
        r(i)=rand();
    elseif types(i)>(1/4) && types(i)<=(2/4)
        types(i)=2;% predator
        r(i)=-rand();
    elseif types(i)>(2/4) && types(i)<=(3/4)
        types(i)=3;% hybrid with positive r
        r(i)=rand();
    elseif types(i)>(3/4) && types(i)<=1
        types(i)=4;% hybrid with negative r
        r(i)=-rand();
    end
end

%%% NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if types(i)==1% prey
        for j=1:number_of_species% who is affecting
            if types(j)==1% both >=0
                if A(i,j)==0 % ignore
                    A(j,i)=0;
                elseif A(i,j)>0 % competition
                    A(j,i)=rand();
                elseif A(i,j)<0% mistake in A
                    if rand()<.5% ignore
                        A(i,j)=0;
                        A(j,i)=0;
                    else% competition
                        A(i,j)=-A(i,j);
                        A(j,i)=rand();
                    end
                    
                end
            elseif types(j)==2% both >=0
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(j,i)=rand();
                elseif A(i,j)<0% mistake in A
                    if rand()<.5% ignore
                        A(i,j)=0;
                        A(j,i)=0;
                    else% competition
                        A(i,j)=-A(i,j);
                        A(j,i)=rand();
                    end
                end
            elseif types(j)==3% positive r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                    end
                elseif A(i,j)<0% mistake in A
                    A(i,j)=-A(i,j);
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                    end
                end
            elseif types(j)==4% negative r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                    end
                elseif A(i,j)<0
                    A(i,j)=-A(i,j);
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                    end
                end
            end
        end
    elseif types(i)==2% predator
        for j=1:number_of_species% who is affecting
            if types(j)==1% both >=0
                if A(i,j)==0 % ignore
                    A(j,i)=0;
                elseif A(i,j)>0 % eating
                    A(j,i)=rand();
                elseif A(i,j)<0% mistake in A
                    if rand()<.5% ignore
                        A(i,j)=0;
                        A(j,i)=0;
                    else% eating
                        A(i,j)=-A(i,j);
                        A(j,i)=rand();
                    end
                end
            elseif types(j)==2% opposite signs
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% eating
                    A(j,i)=-rand();
                elseif A(i,j)<0% being eaten
                    A(j,i)=rand();
                end
            elseif types(j)==3% positive r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% eating
                    A(j,i)=rand();
                elseif A(i,j)<0% being eaten
                    A(j,i)=-rand();
                end
            elseif types(j)==4% negative r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% eating
                    A(j,i)=-rand();
                elseif A(i,j)<0% being eaten
                    A(j,i)=rand();
                end
            end
        end
    elseif types(i)==3% positive r hybrid
        for j=1:number_of_species% who is affecting
            if types(j)==1
                if A(i,j)==0 % ignore
                    A(j,i)=0;
                elseif A(i,j)>0 % competing
                    A(j,i)=rand();
                elseif A(i,j)<0% eating
                    A(j,i)=rand();
                end
            elseif types(j)==2
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(j,i)=rand();
                elseif A(i,j)<0% eating
                    A(j,i)=-rand();
                end
            elseif types(j)==3% positive r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                    end
                elseif A(i,j)<0% eating
                    A(i,j)=rand();
                end
            elseif types(j)==4% negative r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                    end
                elseif A(i,j)<0% eating
                    A(j,i)=-rand();
                end
            end
        end
    elseif types(i)==4% negative r hybrid
        for j=1:number_of_species% who is affecting
            if types(j)==1
                if A(i,j)==0 % ignore
                    A(j,i)=0;
                elseif A(i,j)>0 % competing
                    A(j,i)=rand();
                elseif A(i,j)<0% eating
                    A(j,i)=rand();
                end
            elseif types(j)==2
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(j,i)=-rand();
                elseif A(i,j)<0% eating
                    A(j,i)=rand();
                end
            elseif types(j)==3% positive r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(i,j)=rand();
                elseif A(i,j)<0% eating
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                    end
                end
            elseif types(j)==4% negative r hybrid
                if A(i,j)==0% ignore
                    A(j,i)=0;
                elseif A(i,j)>0% being eaten
                    A(j,i)=-rand();
                elseif A(i,j)<0% eating
                    if A(j,i)==0% if it isn't already picked
                        A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                    end
                end
            end
        end
    end
end

%%% INITIAL POPULATION
l=50.0;ll=30;y=(l-ll).*rand(1,number_of_species)+ll;

A=1/200.*A; % NORMALIZING

%%% EXTINCTION RATE AND POINT AND BOOKKEEPING
extinction_point=1;
k_extinction=.1;
extinct=zeros(size(y));
save('forformula.mat','extinction_point','number_of_species','A','k_extinction','r','extinct');

%%% DEFINING THE FUNCTIONS AND TIMES
fff=@ggg5;
starttime=0;endtime=150;
endtime4=550;
endtime3=350;
endtime2=250;
% endtime4=50;
% endtime3=50;
% endtime2=50;
l=50;

%%%%%% PART NUMBER 1

%%% RUNNING, MODEL WITH EXTINCTION THRESHOLD
[t,p]=ode45(fff,[starttime endtime],y');
k=1;

%%% MAKING THE LEGEND ENTRIES
format short
load('forformula.mat')
%extinct_values{k}=extinct;
Display_legend=string(zeros(1,number_of_species));
ext_array=zeros(1,4);% prey,pred,PreyLH,PredLH
ext_ind_array=cell(1,4);% ROW ARRAY, ALL THE EXTINCT INDICES SPECIES
surv_array=zeros(1,4);% prey,pred,PreyLH,PredLH
surv_ind_array=cell(1,4);% ROW ARRAY, ALL THE SURVIVING INDICES SPECIES
%%% RECORDING HOW MANY SURVIVED AND WENT EXTINCT
for i=1:number_of_species
    if types(i)==1
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,1)=ext_array(1,1)+1;
            ext_ind_array{1}=[ext_ind_array{1} i];
        else
            Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p(length(p),i) );
            surv_array(1,1)=surv_array(1,1)+1;
            surv_ind_array{1}=[surv_ind_array{1} i];
        end
    elseif types(i)==2
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,2)=ext_array(1,2)+1;
            ext_ind_array{2}=[ext_ind_array{2} i];
        else
            Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p(length(p),i) );
            surv_array(1,2)=surv_array(1,2)+1;
            surv_ind_array{2}=[surv_ind_array{2} i];
        end
    elseif types(i)==3
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,3)=ext_array(1,3)+1;
            ext_ind_array{3}=[ext_ind_array{3} i];
        else
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p(length(p),i) );
            surv_array(1,3)=surv_array(1,3)+1;
            surv_ind_array{3}=[surv_ind_array{3} i];
        end
    elseif types(i)==4
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,4)=ext_array(1,4)+1;
            ext_ind_array{4}=[ext_ind_array{4} i];
        else
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p(length(p),i) );
            surv_array(1,4)=surv_array(1,4)+1;
            surv_ind_array{4}=[surv_ind_array{4} i];
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    A1=A;
    ext_arrayk1=ext_array;
    ext_ind_arrayk1=ext_ind_array;
    surv_arrayk1=surv_array;
    surv_ind_arrayk1=surv_ind_array;
    total_amount_arrayk1=ext_arrayk1+surv_arrayk1;
    Display_legendk1=Display_legend;
    save('Display_legend_arrayk1.mat','Display_legendk1','total_amount_arrayk1',...
        'ext_arrayk1','surv_arrayk1','A1','ext_ind_arrayk1','surv_ind_arrayk1')
elseif k==2
    A2=A;
    ext_arrayk2=ext_array;
    ext_ind_arrayk2=ext_ind_array;
    surv_arrayk2=surv_array;
    surv_ind_arrayk2=surv_ind_array;
    total_amount_arrayk2=ext_arrayk2+surv_arrayk2;
    Display_legendk2=Display_legend;
    save('Display_legend_arrayk2.mat','Display_legendk2','total_amount_arrayk2',...
        'ext_arrayk2','surv_arrayk2','A2','ext_ind_arrayk2','surv_ind_arrayk2')
elseif k==3
    A3=A;
    ext_arrayk3=ext_array;
    surv_arrayk3=surv_array;
    ext_ind_arrayk3=ext_ind_array;
    surv_ind_arrayk3=surv_ind_array;
    total_amount_arrayk3=ext_arrayk3+surv_arrayk3;
    Display_legendk3=Display_legend;
    save('Display_legend_arrayk3.mat','Display_legendk3','total_amount_arrayk3',...
        'ext_arrayk3','surv_arrayk3','A3','ext_ind_arrayk3','surv_ind_arrayk3')
elseif k==4
    A4=A;
    ext_arrayk4=ext_array;
    surv_arrayk4=surv_array;
    ext_ind_arrayk4=ext_ind_array;
    surv_ind_arrayk4=surv_ind_array;
    total_amount_arrayk4=ext_arrayk4+surv_arrayk4;
    Display_legendk4=Display_legend;
    save('Display_legend_arrayk4.mat','Display_legendk4','total_amount_arrayk4',...
        'ext_arrayk4','surv_arrayk4','A4','ext_ind_arrayk4','surv_ind_arrayk4')
end



%%% SEE WHICH ARE EXTINCT
extinction_index_array=[];% has all the indexes that correspond to the extinct species
for i=1:97
    if p(size(p,1),i)<extinction_point
        extinction_index_array=[extinction_index_array i];
    end
end

%%% REPLACE THEM WITH INVADERS (WE'LL KEEP IT THE SAME SIZE)
types_new=rand(1,length(extinction_index_array));
r_new=zeros(1,length(extinction_index_array));% PREALLOCATION
for i=1:length(extinction_index_array)
    %DECIDE WHAT TYPES THE INVADERS ARE AND DEFINE THEIR Rs
    if types_new(i)<=(1/4)
        types_new(i)=1;% prey
        r_new(i)=rand();
    elseif types_new(i)>(1/4) && types_new(i)<=(2/4)
        types_new(i)=2;% predator
        r_new(i)=-rand();
    elseif types_new(i)>(2/4) && types_new(i)<=(3/4)
        types_new(i)=3;% hybrid with positive r
        r_new(i)=rand();
    elseif types_new(i)>(3/4) && types_new(i)<=1
        types_new(i)=4;% hybrid with negative r
        r_new(i)=-rand();
    end
    % REPLACE TYPES AND Rs
    types(extinction_index_array(i))=types_new(i);
    r(extinction_index_array(i))=r_new(i);
end

%%% REDO A MATRIX, NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if ismember(i,extinction_index_array)% BUT ONLY IF THE SPECIES IS EXTINCT
        if types(i)==1% prey
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competition
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end

                    end
                elseif types(j)==2% both >=0
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    elseif A(i,j)<0% mistake in A
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    elseif A(i,j)<0
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    end
                end
            end
        elseif types(i)==2% predator
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% eating
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==2% opposite signs
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=-rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                end
            end
        elseif types(i)==3% positive r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                        end
                    elseif A(i,j)<0% eating
                        A(i,j)=rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                end
            end
        elseif types(i)==4% negative r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(i,j)=rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                end
            end
        end
    end
end

%%% UPDATE THE STARTING POPULATION NUMBERS
y=p(size(p,1),:);
l=50.0;ll=30;
for i=1:length(extinction_index_array)
    y(extinction_index_array(i))=(l-ll).*rand()+ll;
end


%%%%%%%% PART NUMBER 2

%%% RUNNING, MODEL WITH INVASION/LOSS
[t1,p1]=ode45(fff,[endtime endtime2],y');
k=2;


%%% MAKING THE LEGEND ENTRIES
format short
load('forformula.mat')
%extinct_values{k}=extinct;
Display_legend=string(zeros(1,number_of_species));
ext_array=zeros(1,4);% prey,pred,PreyLH,PredLH
ext_ind_array=cell(1,4);% ROW ARRAY, ALL THE EXTINCT INDICES SPECIES
surv_array=zeros(1,4);% prey,pred,PreyLH,PredLH
surv_ind_array=cell(1,4);% ROW ARRAY, ALL THE SURVIVING INDICES SPECIES
%%% RECORDING HOW MANY SURVIVED AND WENT EXTINCT
for i=1:number_of_species
    if types(i)==1
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p1(length(p1),i) );
            ext_array(1,1)=ext_array(1,1)+1;
            ext_ind_array{1}=[ext_ind_array{1} i];
        else
            Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p1(length(p1),i) );
            surv_array(1,1)=surv_array(1,1)+1;
            surv_ind_array{1}=[surv_ind_array{1} i];
        end
    elseif types(i)==2
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p1(length(p1),i) );
            ext_array(1,2)=ext_array(1,2)+1;
            ext_ind_array{2}=[ext_ind_array{2} i];
        else
            Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p1(length(p1),i) );
            surv_array(1,2)=surv_array(1,2)+1;
            surv_ind_array{2}=[surv_ind_array{2} i];
        end
    elseif types(i)==3
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p1(length(p1),i) );
            ext_array(1,3)=ext_array(1,3)+1;
            ext_ind_array{3}=[ext_ind_array{3} i];
        else
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p1(length(p1),i) );
            surv_array(1,3)=surv_array(1,3)+1;
            surv_ind_array{3}=[surv_ind_array{3} i];
        end
    elseif types(i)==4
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p1(length(p1),i) );
            ext_array(1,4)=ext_array(1,4)+1;
            ext_ind_array{4}=[ext_ind_array{4} i];
        else
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p1(length(p1),i) );
            surv_array(1,4)=surv_array(1,4)+1;
            surv_ind_array{4}=[surv_ind_array{4} i];
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    A1=A;
    ext_arrayk1=ext_array;
    ext_ind_arrayk1=ext_ind_array;
    surv_arrayk1=surv_array;
    surv_ind_arrayk1=surv_ind_array;
    total_amount_arrayk1=ext_arrayk1+surv_arrayk1;
    Display_legendk1=Display_legend;
    save('Display_legend_arrayk1.mat','Display_legendk1','total_amount_arrayk1',...
        'ext_arrayk1','surv_arrayk1','A1','ext_ind_arrayk1','surv_ind_arrayk1')
elseif k==2
    A2=A;
    ext_arrayk2=ext_array;
    ext_ind_arrayk2=ext_ind_array;
    surv_arrayk2=surv_array;
    surv_ind_arrayk2=surv_ind_array;
    total_amount_arrayk2=ext_arrayk2+surv_arrayk2;
    Display_legendk2=Display_legend;
    save('Display_legend_arrayk2.mat','Display_legendk2','total_amount_arrayk2',...
        'ext_arrayk2','surv_arrayk2','A2','ext_ind_arrayk2','surv_ind_arrayk2')
elseif k==3
    A3=A;
    ext_arrayk3=ext_array;
    surv_arrayk3=surv_array;
    ext_ind_arrayk3=ext_ind_array;
    surv_ind_arrayk3=surv_ind_array;
    total_amount_arrayk3=ext_arrayk3+surv_arrayk3;
    Display_legendk3=Display_legend;
    save('Display_legend_arrayk3.mat','Display_legendk3','total_amount_arrayk3',...
        'ext_arrayk3','surv_arrayk3','A3','ext_ind_arrayk3','surv_ind_arrayk3')
elseif k==4
    A4=A;
    ext_arrayk4=ext_array;
    surv_arrayk4=surv_array;
    ext_ind_arrayk4=ext_ind_array;
    surv_ind_arrayk4=surv_ind_array;
    total_amount_arrayk4=ext_arrayk4+surv_arrayk4;
    Display_legendk4=Display_legend;
    save('Display_legend_arrayk4.mat','Display_legendk4','total_amount_arrayk4',...
        'ext_arrayk4','surv_arrayk4','A4','ext_ind_arrayk4','surv_ind_arrayk4')
end



%%% SEE WHICH ARE EXTINCT
extinction_index_array=[];% has all the indexes that correspond to the extinct species
for i=1:97
    if p1(size(p1,1),i)<extinction_point
        extinction_index_array=[extinction_index_array i];
    end
end

%%% REPLACE THEM WITH INVADERS (WE'LL KEEP IT THE SAME SIZE)
types_new=rand(1,length(extinction_index_array));
r_new=zeros(1,length(extinction_index_array));% PREALLOCATION
for i=1:length(extinction_index_array)
    %DECIDE WHAT TYPES THE INVADERS ARE AND DEFINE THEIR Rs
    if types_new(i)<=(1/4)
        types_new(i)=1;% prey
        r_new(i)=rand();
    elseif types_new(i)>(1/4) && types_new(i)<=(2/4)
        types_new(i)=2;% predator
        r_new(i)=-rand();
    elseif types_new(i)>(2/4) && types_new(i)<=(3/4)
        types_new(i)=3;% hybrid with positive r
        r_new(i)=rand();
    elseif types_new(i)>(3/4) && types_new(i)<=1
        types_new(i)=4;% hybrid with negative r
        r_new(i)=-rand();
    end
    % REPLACE TYPES AND Rs
    types(extinction_index_array(i))=types_new(i);
    r(extinction_index_array(i))=r_new(i);
end

%%% REDO A MATRIX, NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if ismember(i,extinction_index_array)% BUT ONLY IF THE SPECIES IS EXTINCT
        if types(i)==1% prey
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competition
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end

                    end
                elseif types(j)==2% both >=0
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    elseif A(i,j)<0% mistake in A
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    elseif A(i,j)<0
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    end
                end
            end
        elseif types(i)==2% predator
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% eating
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==2% opposite signs
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=-rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                end
            end
        elseif types(i)==3% positive r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                        end
                    elseif A(i,j)<0% eating
                        A(i,j)=rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                end
            end
        elseif types(i)==4% negative r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(i,j)=rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                end
            end
        end
    end
end

%%% UPDATE THE STARTING POPULATION NUMBERS
y=p1(size(p1,1),:);
l1=50.0;ll=30;
for i=1:length(extinction_index_array)
    y(extinction_index_array(i))=(l1-ll).*rand()+ll;
end

%%%%%%%%% PART NUMBER 3


%%% RUNNING, MODEL WITH INVASION/LOSS

[t2,p2]=ode45(fff,[250 350],y');
k=3;

%%% MAKING THE LEGEND ENTRIES
format short
load('forformula.mat')
%extinct_values{k}=extinct;
Display_legend=string(zeros(1,number_of_species));
ext_array=zeros(1,4);% prey,pred,PreyLH,PredLH
ext_ind_array=cell(1,4);% ROW ARRAY, ALL THE EXTINCT INDICES SPECIES
surv_array=zeros(1,4);% prey,pred,PreyLH,PredLH
surv_ind_array=cell(1,4);% ROW ARRAY, ALL THE SURVIVING INDICES SPECIES
%%% RECORDING HOW MANY SURVIVED AND WENT EXTINCT
for i=1:number_of_species
    if types(i)==1
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p2(length(p2),i) );
            ext_array(1,1)=ext_array(1,1)+1;
            ext_ind_array{1}=[ext_ind_array{1} i];
        else
            Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p2(length(p2),i) );
            surv_array(1,1)=surv_array(1,1)+1;
            surv_ind_array{1}=[surv_ind_array{1} i];
        end
    elseif types(i)==2
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p2(length(p2),i) );
            ext_array(1,2)=ext_array(1,2)+1;
            ext_ind_array{2}=[ext_ind_array{2} i];
        else
            Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p2(length(p2),i) );
            surv_array(1,2)=surv_array(1,2)+1;
            surv_ind_array{2}=[surv_ind_array{2} i];
        end
    elseif types(i)==3
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p2(length(p2),i) );
            ext_array(1,3)=ext_array(1,3)+1;
            ext_ind_array{3}=[ext_ind_array{3} i];
        else
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p2(length(p2),i) );
            surv_array(1,3)=surv_array(1,3)+1;
            surv_ind_array{3}=[surv_ind_array{3} i];
        end
    elseif types(i)==4
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p2(length(p2),i) );
            ext_array(1,4)=ext_array(1,4)+1;
            ext_ind_array{4}=[ext_ind_array{4} i];
        else
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p2(length(p2),i) );
            surv_array(1,4)=surv_array(1,4)+1;
            surv_ind_array{4}=[surv_ind_array{4} i];
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    A1=A;
    ext_arrayk1=ext_array;
    ext_ind_arrayk1=ext_ind_array;
    surv_arrayk1=surv_array;
    surv_ind_arrayk1=surv_ind_array;
    total_amount_arrayk1=ext_arrayk1+surv_arrayk1;
    Display_legendk1=Display_legend;
    save('Display_legend_arrayk1.mat','Display_legendk1','total_amount_arrayk1',...
        'ext_arrayk1','surv_arrayk1','A1','ext_ind_arrayk1','surv_ind_arrayk1')
elseif k==2
    A2=A;
    ext_arrayk2=ext_array;
    ext_ind_arrayk2=ext_ind_array;
    surv_arrayk2=surv_array;
    surv_ind_arrayk2=surv_ind_array;
    total_amount_arrayk2=ext_arrayk2+surv_arrayk2;
    Display_legendk2=Display_legend;
    save('Display_legend_arrayk2.mat','Display_legendk2','total_amount_arrayk2',...
        'ext_arrayk2','surv_arrayk2','A2','ext_ind_arrayk2','surv_ind_arrayk2')
elseif k==3
    A3=A;
    ext_arrayk3=ext_array;
    surv_arrayk3=surv_array;
    ext_ind_arrayk3=ext_ind_array;
    surv_ind_arrayk3=surv_ind_array;
    total_amount_arrayk3=ext_arrayk3+surv_arrayk3;
    Display_legendk3=Display_legend;
    save('Display_legend_arrayk3.mat','Display_legendk3','total_amount_arrayk3',...
        'ext_arrayk3','surv_arrayk3','A3','ext_ind_arrayk3','surv_ind_arrayk3')
elseif k==4
    A4=A;
    ext_arrayk4=ext_array;
    surv_arrayk4=surv_array;
    ext_ind_arrayk4=ext_ind_array;
    surv_ind_arrayk4=surv_ind_array;
    total_amount_arrayk4=ext_arrayk4+surv_arrayk4;
    Display_legendk4=Display_legend;
    save('Display_legend_arrayk4.mat','Display_legendk4','total_amount_arrayk4',...
        'ext_arrayk4','surv_arrayk4','A4','ext_ind_arrayk4','surv_ind_arrayk4')
end



%%% SEE WHICH ARE EXTINCT
extinction_index_array=[];% has all the indexes that correspond to the extinct species
for i=1:97
    if p2(size(p2,1),i)<extinction_point
        extinction_index_array=[extinction_index_array i];
    end
end

%%% REPLACE THEM WITH INVADERS (WE'LL KEEP IT THE SAME SIZE)
types_new=rand(1,length(extinction_index_array));
r_new=zeros(1,length(extinction_index_array));% PREALLOCATION
for i=1:length(extinction_index_array)
    %DECIDE WHAT TYPES THE INVADERS ARE AND DEFINE THEIR Rs
    if types_new(i)<=(1/4)
        types_new(i)=1;% prey
        r_new(i)=rand();
    elseif types_new(i)>(1/4) && types_new(i)<=(2/4)
        types_new(i)=2;% predator
        r_new(i)=-rand();
    elseif types_new(i)>(2/4) && types_new(i)<=(3/4)
        types_new(i)=3;% hybrid with positive r
        r_new(i)=rand();
    elseif types_new(i)>(3/4) && types_new(i)<=1
        types_new(i)=4;% hybrid with negative r
        r_new(i)=-rand();
    end
    % REPLACE TYPES AND Rs
    types(extinction_index_array(i))=types_new(i);
    r(extinction_index_array(i))=r_new(i);
end

%%% REDO A MATRIX, NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if ismember(i,extinction_index_array)% BUT ONLY IF THE SPECIES IS EXTINCT
        if types(i)==1% prey
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competition
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end

                    end
                elseif types(j)==2% both >=0
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    elseif A(i,j)<0% mistake in A
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    elseif A(i,j)<0
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    end
                end
            end
        elseif types(i)==2% predator
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% eating
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==2% opposite signs
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=-rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                end
            end
        elseif types(i)==3% positive r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                        end
                    elseif A(i,j)<0% eating
                        A(i,j)=rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                end
            end
        elseif types(i)==4% negative r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(i,j)=rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                end
            end
        end
    end
end

%%% UPDATE THE STARTING POPULATION NUMBERS
y=p2(size(p2,1),:);
l=50.0;ll=30;
for i=1:length(extinction_index_array)
    y(extinction_index_array(i))=(l-ll).*rand()+ll;
end


%%%%%%%%% PART NUMBER 4


%%% RUNNING, MODEL WITH INVASION/LOSS

[t3,p3]=ode45(fff,[350 550],y');
k=4;

%%% MAKING THE LEGEND ENTRIES
format short
load('forformula.mat')
%extinct_values{k}=extinct;
Display_legend=string(zeros(1,number_of_species));
ext_array=zeros(1,4);% prey,pred,PreyLH,PredLH
ext_ind_array=cell(1,4);% ROW ARRAY, ALL THE EXTINCT INDICES SPECIES
surv_array=zeros(1,4);% prey,pred,PreyLH,PredLH
surv_ind_array=cell(1,4);% ROW ARRAY, ALL THE SURVIVING INDICES SPECIES
%%% RECORDING HOW MANY SURVIVED AND WENT EXTINCT
for i=1:number_of_species
    if types(i)==1
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p3(length(p3),i) );
            ext_array(1,1)=ext_array(1,1)+1;
            ext_ind_array{1}=[ext_ind_array{1} i];
        else
            Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p3(length(p3),i) );
            surv_array(1,1)=surv_array(1,1)+1;
            surv_ind_array{1}=[surv_ind_array{1} i];
        end
    elseif types(i)==2
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p3(length(p3),i) );
            ext_array(1,2)=ext_array(1,2)+1;
            ext_ind_array{2}=[ext_ind_array{2} i];
        else
            Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p3(length(p3),i) );
            surv_array(1,2)=surv_array(1,2)+1;
            surv_ind_array{2}=[surv_ind_array{2} i];
        end
    elseif types(i)==3
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p3(length(p3),i) );
            ext_array(1,3)=ext_array(1,3)+1;
            ext_ind_array{3}=[ext_ind_array{3} i];
        else
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p3(length(p3),i) );
            surv_array(1,3)=surv_array(1,3)+1;
            surv_ind_array{3}=[surv_ind_array{3} i];
        end
    elseif types(i)==4
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p3(length(p3),i) );
            ext_array(1,4)=ext_array(1,4)+1;
            ext_ind_array{4}=[ext_ind_array{4} i];
        else
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p3(length(p3),i) );
            surv_array(1,4)=surv_array(1,4)+1;
            surv_ind_array{4}=[surv_ind_array{4} i];
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    A1=A;
    ext_arrayk1=ext_array;
    ext_ind_arrayk1=ext_ind_array;
    surv_arrayk1=surv_array;
    surv_ind_arrayk1=surv_ind_array;
    total_amount_arrayk1=ext_arrayk1+surv_arrayk1;
    Display_legendk1=Display_legend;
    save('Display_legend_arrayk1.mat','Display_legendk1','total_amount_arrayk1',...
        'ext_arrayk1','surv_arrayk1','A1','ext_ind_arrayk1','surv_ind_arrayk1')
elseif k==2
    A2=A;
    ext_arrayk2=ext_array;
    ext_ind_arrayk2=ext_ind_array;
    surv_arrayk2=surv_array;
    surv_ind_arrayk2=surv_ind_array;
    total_amount_arrayk2=ext_arrayk2+surv_arrayk2;
    Display_legendk2=Display_legend;
    save('Display_legend_arrayk2.mat','Display_legendk2','total_amount_arrayk2',...
        'ext_arrayk2','surv_arrayk2','A2','ext_ind_arrayk2','surv_ind_arrayk2')
elseif k==3
    A3=A;
    ext_arrayk3=ext_array;
    surv_arrayk3=surv_array;
    ext_ind_arrayk3=ext_ind_array;
    surv_ind_arrayk3=surv_ind_array;
    total_amount_arrayk3=ext_arrayk3+surv_arrayk3;
    Display_legendk3=Display_legend;
    save('Display_legend_arrayk3.mat','Display_legendk3','total_amount_arrayk3',...
        'ext_arrayk3','surv_arrayk3','A3','ext_ind_arrayk3','surv_ind_arrayk3')
elseif k==4
    A4=A;
    ext_arrayk4=ext_array;
    surv_arrayk4=surv_array;
    ext_ind_arrayk4=ext_ind_array;
    surv_ind_arrayk4=surv_ind_array;
    total_amount_arrayk4=ext_arrayk4+surv_arrayk4;
    Display_legendk4=Display_legend;
    save('Display_legend_arrayk4.mat','Display_legendk4','total_amount_arrayk4',...
        'ext_arrayk4','surv_arrayk4','A4','ext_ind_arrayk4','surv_ind_arrayk4')
end



%%% SEE WHICH ARE EXTINCT
extinction_index_array=[];% has all the indexes that correspond to the extinct species
for i=1:97
    if p3(size(p3,1),i)<extinction_point
        extinction_index_array=[extinction_index_array i];
    end
end

%%% REPLACE THEM WITH INVADERS (WE'LL KEEP IT THE SAME SIZE)
types_new=rand(1,length(extinction_index_array));
r_new=zeros(1,length(extinction_index_array));% PREALLOCATION
for i=1:length(extinction_index_array)
    %DECIDE WHAT TYPES THE INVADERS ARE AND DEFINE THEIR Rs
    if types_new(i)<=(1/4)
        types_new(i)=1;% prey
        r_new(i)=rand();
    elseif types_new(i)>(1/4) && types_new(i)<=(2/4)
        types_new(i)=2;% predator
        r_new(i)=-rand();
    elseif types_new(i)>(2/4) && types_new(i)<=(3/4)
        types_new(i)=3;% hybrid with positive r
        r_new(i)=rand();
    elseif types_new(i)>(3/4) && types_new(i)<=1
        types_new(i)=4;% hybrid with negative r
        r_new(i)=-rand();
    end
    % REPLACE TYPES AND Rs
    types(extinction_index_array(i))=types_new(i);
    r(extinction_index_array(i))=r_new(i);
end

%%% REDO A MATRIX, NO FEEDBACK LOOPS, INPUT LOGIC
for i=1:number_of_species% species affected
    if ismember(i,extinction_index_array)% BUT ONLY IF THE SPECIES IS EXTINCT
        if types(i)==1% prey
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competition
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end

                    end
                elseif types(j)==2% both >=0
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% competition
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    elseif A(i,j)<0% mistake in A
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means eaten
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    elseif A(i,j)<0
                        A(i,j)=-A(i,j);
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% <0 means competition, >0 means eaten
                        end
                    end
                end
            end
        elseif types(i)==2% predator
            for j=1:number_of_species% who is affecting
                if types(j)==1% both >=0
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% mistake in A
                        if rand()<.5% ignore
                            A(i,j)=0;
                            A(j,i)=0;
                        else% eating
                            A(i,j)=-A(i,j);
                            A(j,i)=rand();
                        end
                    end
                elseif types(j)==2% opposite signs
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=-rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% eating
                        A(j,i)=-rand();
                    elseif A(i,j)<0% being eaten
                        A(j,i)=rand();
                    end
                end
            end
        elseif types(i)==3% positive r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means competition, <0 means being eaten
                        end
                    elseif A(i,j)<0% eating
                        A(i,j)=rand();
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    elseif A(i,j)<0% eating
                        A(j,i)=-rand();
                    end
                end
            end
        elseif types(i)==4% negative r hybrid
            for j=1:number_of_species% who is affecting
                if types(j)==1
                    if A(i,j)==0 % ignore
                        A(j,i)=0;
                    elseif A(i,j)>0 % competing
                        A(j,i)=rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==2
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        A(j,i)=rand();
                    end
                elseif types(j)==3% positive r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(i,j)=rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                elseif types(j)==4% negative r hybrid
                    if A(i,j)==0% ignore
                        A(j,i)=0;
                    elseif A(i,j)>0% being eaten
                        A(j,i)=-rand();
                    elseif A(i,j)<0% eating
                        if A(j,i)==0% if it isn't already picked
                            A(j,i)=rand()-.5.*rand();% >0 means being eaten, <0 means competition
                        end
                    end
                end
            end
        end
    end
end

%%% UPDATE THE STARTING POPULATION NUMBERS
y=p3(size(p3,1),:);
l=50.0;ll=30;
for i=1:length(extinction_index_array)
    y(extinction_index_array(i))=(l-ll).*rand()+ll;
end



%%% PLOTTING
figure
plot(t,p,'LineWidth',2)
hold on
plot(t1,p1,'LineWidth',2)
hold on
plot(t2,p2,'LineWidth',2)
hold on
plot(t3,p3,'LineWidth',2)
hold off
title('Model with invasion and loss')
xlabel('time elapsed in days')
ylabel('specie population')

load('Display_legend_arrayk1.mat')
load('Display_legend_arrayk2.mat')
load('Display_legend_arrayk3.mat')
load('Display_legend_arrayk4.mat')

%%% CHECKING IF A SPECIES i SURVIVED THEM ALL
survived_all=cell(1,4);
for i=1:number_of_species
    if ismember(i,surv_ind_arrayk1{1}) && ismember(i,surv_ind_arrayk2{1}) &&...
            ismember(i,surv_ind_arrayk3{1}) && ismember(i,surv_ind_arrayk4{1})
        survived_all{1}=[survived_all{1} i];
    elseif ismember(i,surv_ind_arrayk1{2}) && ismember(i,surv_ind_arrayk2{2}) &&...
            ismember(i,surv_ind_arrayk3{2}) && ismember(i,surv_ind_arrayk4{2})
        survived_all{2}=[survived_all{2} i];
    elseif ismember(i,surv_ind_arrayk1{3}) && ismember(i,surv_ind_arrayk2{3}) &&...
            ismember(i,surv_ind_arrayk3{3}) && ismember(i,surv_ind_arrayk4{3})
        survived_all{3}=[survived_all{3} i];
    elseif ismember(i,surv_ind_arrayk1{4}) && ismember(i,surv_ind_arrayk2{4}) &&...
            ismember(i,surv_ind_arrayk3{4}) && ismember(i,surv_ind_arrayk4{4})
        survived_all{4}=[survived_all{4} i];
    end
end

%% for part c)

ratio_mat=string(zeros(4));
for i=1:4
    ratio_mat(1,i)=sprintf('%d/%d=%f',surv_arrayk1(i),total_amount_arrayk1(i),surv_arrayk1(i)/total_amount_arrayk1(i));
    ratio_mat(2,i)=sprintf('%d/%d=%f',surv_arrayk2(i),total_amount_arrayk2(i),surv_arrayk2(i)/total_amount_arrayk2(i));
    ratio_mat(3,i)=sprintf('%d/%d=%f',surv_arrayk3(i),total_amount_arrayk3(i),surv_arrayk3(i)/total_amount_arrayk3(i));
    ratio_mat(4,i)=sprintf('%d/%d=%f',surv_arrayk4(i),total_amount_arrayk4(i),surv_arrayk4(i)/total_amount_arrayk4(i));
end

%%% CREATING THE TABLE
survival_rates=array2table(ratio_mat,'VariableNames',{'Prey','Predators',...
'PreyLeaningHybrids','PredatorLeaningHybrids'},'RowNames',...
{'first simulation','second simulation','third simulation','fourth simulation'})


