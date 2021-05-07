function run_diff_NoS_sbplt_num(number_of_species,k,h,J,endtime_value_h,endtime_value_J)
% k is the plot index
format long
rng(5)
density=.3;
blah=" sparse";

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
%b=2;a=1;diagg = (b-a).*rand(size(A,1),1) + a;% positive, diagonal entries between a and b
diagg=2.*rand(size(A,1),1);
%diagg=zeros(size(A,1),1);
for i=1:size(A,1)
    A(i,i)=diagg(i);
end

%%% CHOOSE WHICH ARE PREDATORS AND WHICH ARE PREY
%%% ALSO DEFINING THE Rs
types=rand(1,number_of_species);
r=zeros(1,number_of_species);% PREALLOCATION
ttt=1;
for i=1:number_of_species
    if types(i)<=(1/4)
        types(i)=1;% prey
%         A(i,i)=ttt*rand();
        r(i)=rand();
    elseif types(i)>(1/4) && types(i)<=(2/4)
        types(i)=2;% predator
%         A(i,i)=ttt*rand();
        r(i)=-rand();
    elseif types(i)>(2/4) && types(i)<=(3/4)
        types(i)=3;% hybrid with positive r
%         A(i,i)=ttt*rand();
        r(i)=rand();
    elseif types(i)>(3/4) && types(i)<=1
        types(i)=4;% hybrid with negative r
%         A(i,i)=ttt*rand();
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
k_extinction=.1;%CHANGED THIS
extinct=zeros(size(y));

save('forformula.mat','extinction_point'...
    ,'number_of_species','A','k_extinction','r','extinct');

%%% DEFINING THE FUNCTIONS AND TIMES
% dpop=@(t,y) r'.*(ones(number_of_species,1)-flipud(A*y)).*y;
fff=@ggg5;
dpop1=@(t,y) r'.*(ones(number_of_species,1)-(A*y)).*y;
starttime=0;

%%% RUNNING, MODEL WITH EXTINCTION THRESHOLD
[t,p]=ode45(fff,[starttime endtime_value_h],y');
% stable_values{k}=p(length(p));

%%% MAKING THE LEGEND ENTRIES
format short
load('forformula.mat')
%extinct_values{k}=extinct;
Display_legend=string(zeros(1,number_of_species));
ext_array=zeros(1,4);% prey,pred,PreyLH,PredLH
surv_array=zeros(1,4);% prey,pred,PreyLH,PredLH
for i=1:number_of_species
    if types(i)==1
        %Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p(length(p),i) );
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,1)=ext_array(1,1)+1;
        else
            Display_legend(i)=sprintf( "#%d. prey, stable value=%d",i,p(length(p),i) );
            surv_array(1,1)=surv_array(1,1)+1;
        end
    elseif types(i)==2
        %Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p(length(p),i) );
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,2)=ext_array(1,2)+1;
        else
            Display_legend(i)=sprintf( "#%d. predator, stable value=%d",i,p(length(p),i) );
            surv_array(1,2)=surv_array(1,2)+1;
        end
    elseif types(i)==3
%         Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p(length(p),i) );
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,3)=ext_array(1,3)+1;
        else
            Display_legend(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p(length(p),i) );
            surv_array(1,3)=surv_array(1,3)+1;
        end
    elseif types(i)==4
%         Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p(length(p),i) );
        if extinct(i)==1
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p(length(p),i) );
            ext_array(1,4)=ext_array(1,4)+1;
        else
            Display_legend(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p(length(p),i) );
            surv_array(1,4)=surv_array(1,4)+1;
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    ext_arrayk1=ext_array;
    surv_arrayk1=surv_array;
    total_amount_arrayk1=ext_arrayk1+surv_arrayk1;
    Display_legendk1=Display_legend;
    save('Display_legend_arrayk1.mat','Display_legendk1','total_amount_arrayk1',...
        'ext_arrayk1','surv_arrayk1')
elseif k==2
    ext_arrayk2=ext_array;
    surv_arrayk2=surv_array;
    total_amount_arrayk2=ext_arrayk2+surv_arrayk2;
    Display_legendk2=Display_legend;
    save('Display_legend_arrayk2.mat','Display_legendk2','total_amount_arrayk2',...
        'ext_arrayk2','surv_arrayk2')
elseif k==3
    ext_arrayk3=ext_array;
    surv_arrayk3=surv_array;
    total_amount_arrayk3=ext_arrayk3+surv_arrayk3;
    Display_legendk3=Display_legend;
    save('Display_legend_arrayk3.mat','Display_legendk3','total_amount_arrayk3',...
        'ext_arrayk3','surv_arrayk3')
elseif k==4
    ext_arrayk4=ext_array;
    surv_arrayk4=surv_array;
    total_amount_arrayk4=ext_arrayk4+surv_arrayk4;
    Display_legendk4=Display_legend;
    save('Display_legend_arrayk4.mat','Display_legendk4','total_amount_arrayk4',...
        'ext_arrayk4','surv_arrayk4')
end
% save('forformula.mat','extinction_point'...
%     ,'number_of_species','A','k_extinction','r','extinct')

%%% PLOTTING

plot(h,t,p,'LineWidth',2);h.YLimMode='manual';
h.YLim=[0,inf];
% ax = gca;
% ax.YLim = [0 inf];
%xlabel('time elapsed in days')
%ylabel('specie population')
%legend(Display_legend,'Location','northeast','Orientation','vertical')
%lgd = legend('boxoff');
%lgd.FontWeight = 'bold';

%%% RUNNING, MODEL WITHOUT EXTINCTION THRESHOLD
% if k==4
%     endtime=200;
% end
[t1,p1]=ode45(dpop1,[starttime endtime_value_J],y');
%stable_values1{k}=p1(length(p1));

%%% MAKING THE LEGEND ENTRIES AND SAVING WHOS EXTINCT DATA
format short
load('forformula.mat')
%extinct_values1{k}=extinct;
Display_legend1=string(zeros(1,number_of_species));
for i=1:number_of_species
    if types(i)==1
        if extinct(i)==1
            Display_legend1(i)=sprintf( "#%d. prey(extinct), stable value=%d",i,p1(length(p1),i) );
        else
            Display_legend1(i)=sprintf( "#%d. prey, stable value=%d",i,p1(length(p1),i) );
        end
    elseif types(i)==2
        if extinct(i)==1
            Display_legend1(i)=sprintf( "#%d. predator(extinct), stable value=%d",i,p1(length(p1),i) );
        else
            Display_legend1(i)=sprintf( "#%d. predator, stable value=%d",i,p1(length(p1),i) );
        end
    elseif types(i)==3
        if extinct(i)==1
            Display_legend1(i)=sprintf( "#%d. prey-leaning hybrid(extinct), stable value=%d",i,p1(length(p1),i) );
        else
            Display_legend1(i)=sprintf( "#%d. prey-leaning hybrid, stable value=%d",i,p1(length(p1),i) );
        end
    elseif types(i)==4
        if extinct(i)==1
            Display_legend1(i)=sprintf( "#%d. predator-leaning hybrid(extinct), stable value=%d",i,p1(length(p1),i) );
        else
            Display_legend1(i)=sprintf( "#%d. predator-leaning hybrid, stable value=%d",i,p1(length(p1),i) );
        end
    end
end
format long
%%% SAVING HOW MANY SURVIVED AND HOW MANY WENT EXTINCT
if k==1
    Display_legend1k1=Display_legend1;
    save('Display_legend_array1k1.mat','Display_legend1k1')
elseif k==2
    Display_legend1k2=Display_legend1;
    save('Display_legend_array1k2.mat','Display_legend1k2')
elseif k==3
    Display_legend1k3=Display_legend1;
    save('Display_legend_array1k3.mat','Display_legend1k3')
elseif k==4
    Display_legend1k4=Display_legend1;
    save('Display_legend_array1k4.mat','Display_legend1k4')
end
% save('forformula.mat','extinction_point'...
%     ,'number_of_species','A','k_extinction','r','extinct')

%%% PLOTTING

plot(J,t1,p1,'LineWidth',2);
J.YLimMode='manual';
J.YLim=[0,inf];
% ax = gca;
% ax.YLim = [0 inf];
% xlabel('time elapsed in days')
% ylabel('specie population')
% legend(Display_legend1,'Location','northeast','Orientation','vertical')
% lgd = legend('boxoff');
% lgd.FontWeight = 'bold';

% lowest_pops=min(p)
% lowest_pops1=min(p1)
% extinction_point


end
