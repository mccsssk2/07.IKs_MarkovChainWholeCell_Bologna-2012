% time is in ms, voltage is in mV, current is whatever Gks x mV is.
function [cost] = hiksiv(p)

global filecounter;
global CaseCase;
global CostCost;
% global baseCost;
global printOutput;
global useSingle;
global fixedGKs;

global vcc;
global defGKs;

% multiplierIndex = [1 3 6 7 9 10 12 14];
% if sum(p(multiplierIndex)<0) || p(4)>0 || p(16)<1 || p(16)>200 
%     if sum(p(multiplierIndex)<0)
%         [ans b]=min(p(multiplierIndex));
%         p(multiplierIndex(b)) = 1e-10;
%     end
%     if p(4)>0
%         p(4)=-1e-4;
%     end
%    if p(16)<1 || p(16)>200 
%        p(16)=defGKs;
%    end
%     
%     cost = 1e15;
% else
if CaseCase ~= 3
    minGKs = 50; %nS
else
    minGKs = 40;
end

nonNegativeIndex = [1 2 3 5 6 7 9 10 12 13 ];%14];
nonPositiveIndex = [4 8 11 15];
exponents = [2 4 6 8 11 13 15];
%constrains
if sum(p(nonNegativeIndex)<0) || sum(p(nonPositiveIndex)>0) || p(16)< minGKs || p(16)>200 || ...
        sum(abs(p([nonNegativeIndex nonPositiveIndex]))<1e-5) ||sum(abs(p(exponents))<1e-3) || ...
        sum(abs(p(exponents))>10) || p(14)<1e-10 %nS
        cost = 1e15;
else
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fixedGKs == 1
p(16)=defGKs;
end

% get rid of the intial transient. IMPORTANT - the fminsearch result
% depends on the initial conditions of the system variables. 
x0=[1 zeros(1,16)]; % Initial conditions. for initial conditions.
yp = x0;
% this is the solution for initial conditions of y at -80 mV. It MAKES a
% big difference.
[y1 y2]=ode23s(@hksrates,[0 1e7],yp,[],vcc(1),p);

yp = [];
yp=y2(length(y1),:); % take the final values of c's and put them in here.

% least squared cost. Put in weights by looking into Stefano's code.

cost = zeros(1,3);
cost(1)=costIV(p,yp);

if useSingle == 1
cost(2) = cost5p5(p,yp);

if CaseCase ~= 3
    cost(3) = 0.3*cost12(p,yp);
end

if CostCost==1
cost(2) = 0.1*cost(2);
end
end

if printOutput == 1
    if useSingle == 1
    disp('IV traces cost =')
    disp(' ')
    disp(cost(1))
    disp('5.5s trace cost =')
    disp(' ')
    disp(cost(2))
    disp('12s trace cost =')
    disp(' ')
    disp(cost(3))
    end
    disp('total cost =')
    disp(' ')
    disp(sum(cost))
end

% if filecounter == 1
%     baseCost = cost;
%     if baseCost <1e-15
%         error('Base cost is too small')
%     end
% end


% final value of cost
cost = sum(cost);
% cost=cost/baseCost;


if(cost~=cost)
    error('The cost function is not defined any more.'); 
end;

if(abs(cost)<1e-6) 
    error('My limits are reached although I have low expectations from fminsearch as of now.'); 
end;

filecounter = filecounter + 1;
output = [p cost];

% for a restart from here...
save('restart.dat','p','-ASCII');
save('estimates.txt','output','-ASCII','-APPEND');
end

% time is in ms, voltage is in mV, current is whatever Gks x mV is.
function [cost] = costIV(p,yp)

global filecounter;
global CaseCase;
global printOutput;
global CostCost;
global tracesIV
global weightIV;
global weightDataIV;

global erev;
global vcc;

% the IV protocol.
ivv=[vcc(2):vcc(3):vcc(4)]';   % Potentials to depolarize to (mV)

%%%%%%%%%%%%%%%%%%%%%%% load in the data. %%%%%%%%%%%%%%%%%%%%%

activationStep = 21978;
deactivationStep=activationStep + 1;

% read in some experimental data
 testexptdata = tracesIV;
% expdata = testexptdata;
 expttime = testexptdata(:,1);
 expdata = testexptdata(:,2:16);
 expdata = fliplr(expdata);
  
 % plot(expttime,expdata/max(expdata)); pause(100);

 exptIKs = [];

 % time stamps in the experimental data need to be sorted as above.

 missingtime=[95 expttime(1)];
 activationexpttime = expttime([1:activationStep]);
 deactivationexpttime = expttime([deactivationStep:end]);
 exptIKstemp = expdata(:,(15-length(ivv)+1):15); % /max(expdata(:,15));
 exptIKs = exptIKstemp([1:end],:);

 
Iact=zeros(activationStep,length(ivv));
Itail=zeros(length(expttime) - activationStep,length(ivv));

for i=1:1:length(ivv)
% the following solution time matches the solution to the expt data.
%    [y1 y2]=ode23s(@hksrates,expttime,yp,[],ivv(i),p);

[y0 yi]=ode23s(@hksrates,missingtime,yp,[],ivv(i),p);   % instead of the capacitive current in the experiment
[y1 y2]=ode23s(@hksrates,activationexpttime,yi(length(y0),:),[],ivv(i),p);
[y00 yii]=ode23s(@hksrates,[activationexpttime(end) deactivationexpttime(1)],y2(length(y1),:),[],vcc(5),p); % instead of the capacitive current in the experiment
[y3 y4]=ode23s(@hksrates,deactivationexpttime,yii(length(y00),:),[],vcc(5),p);
%     GKs(i)=max(expdata(:,i))/((ivv(i)-erev)*sum(y2(end,16:17),2));
Iact(:,i)=sum(y2(:,16:17),2)*(ivv(i)-erev);
Itail(:,i)=sum(y4(:,16:17),2)*(vcc(5)-erev);


end % end of for loop.

% GKs=max(expdata(:,15))/((max(ivv)-erev)*sum(y2(end,16:17),2));
Iact=p(16)*Iact;
Itail=p(16)*Itail;

data(:,1)=[y1;y3];
data(:,2:length(ivv)+1)=[Iact;Itail]; %/max(max([Iact;Itail])); 

if(printOutput==1)
    h = figure;                                  % Plot
    f1 = plot(data(:,1),data(:,2:length(ivv)+1));
    hold on
    f2 = plot([activationexpttime;deactivationexpttime],exptIKs);

    if CaseCase == 0
        plotName = sprintf('Control I_{Ks} IV');
    end
    if CaseCase == 1
        plotName = sprintf('Acute DHA I_{Ks} IV');
    end
    if CaseCase == 2
        plotName = sprintf('Acute EPA I_{Ks} IV');
    end
    if CaseCase == 3
        plotName = sprintf('Chronic EPA I_{Ks} IV');
    end

    title(plotName,'FontSize',18);
    xlabel('Time (ms)','FontSize',16);
    ylabel('I_{Ks} ','FontSize',16);
    set(f1,'Color',[0.5 0 1],'LineWidth',2);
    set(f2,'Color',[1 0 1],'LineWidth',2);
    set(gca,'FontSize',12,'FontName','Roman');
    % set(gca,'XLim',[0 8000],'YLim',[-0.05 1.05]);
    set(gca,'XTick',[0:2000:8000],'FontName','Symbol','FontSize',16);
    % set(gca,'YTick',[0:0.5:1],'FontName','Symbol','FontSize',16);
    pause(0.1)

    filename = sprintf('png/IVfig%d.png',filecounter);
    saveas(h,filename,'png');
    pause
    % clf
end

cost=0;
if CostCost==0
    for i=1:1:length(ivv)
    cost = cost + ((data(:,i+1) - exptIKs(:,i)))'*((data(:,i+1) - exptIKs(:,i)));
    end
    cost=cost/length(ivv)/length(data);
elseif CostCost==1
    weightTrace=weightIV;
    weightData = weightDataIV;

    for i=1:1:length(ivv)
        simulation=[data(1:activationStep,i+1); data(activationStep+1:end,i+1)];
        target=[exptIKs(1:activationStep,i); exptIKs(activationStep+1:end,i)];    

        cost = cost + weightTrace(i)*(weightData.*(simulation - target))'*(weightData.*(simulation - target));
    end
    cost=cost*length(ivv)/(sum(weightDataIV)*sum(weightTrace(1:i)));
elseif CostCost==2
   DexptIKsdt=zeros(size(exptIKs));
   weightTrace=weightIV;
   for i=1:length(data)-1
       DexptIKsdt(i+1,:)=(data(i+1,2:length(ivv)+1)-data(i,2:length(ivv)+1))./(data(i+1,1)-data(i,1));
   end
   weightD=zeros(length(data),length(ivv));
   for i=1:1:length(ivv)
   weightD(:,i) = weightTrace(i)*(abs(DexptIKsdt(:,i))+2);
   simulation=[data(1:activationStep,i+1); data(activationStep+1:end,i+1)];
   target=[exptIKs(1:activationStep,i); exptIKs(activationStep+1:end,i)];

   cost = cost + ((simulation - target).*(weightD(:,i)))'*((simulation - target).*(weightD(:,i)));     

   end 
   cost=cost*length(ivv)/sum(sum(weightD));
   
   if printOutput ==1
        h = figure(2);
        f1 = plot(data(:,1),weightD(:,1:length(ivv)),'.');
        title('Derivative cost function','FontSize',18);
        xlabel('ms','FontSize',16);
        ylabel('Weight','FontSize',16);
        set(f1,'Color',[0.5 0 1],'MarkerSize',4);
        set(gca,'FontSize',12,'FontName','Roman');
        set(gca,'XLim',[0 6000],'YLim',[0 40]);
        set(gca,'XTick',[0:2000:8000],'FontName','Symbol','FontSize',16);
        % set(gca,'YTick',[0:0.5:1],'FontName','Symbol','FontSize',16);
        pause

        filename = sprintf('png/derivative%d.png',filecounter);
        saveas(h,filename,'png');
        % clf
   end

else
    error('Specify the cost function to use')
end


% time is in ms, voltage is in mV, current is whatever Gks x mV is.
function [cost] = cost5p5(p,yp)


global filecounter;
global CaseCase;
global printOutput;
global CostCost;
global trace5p5;

global erev;
global vcc;

% the IV protocol.
ivv=[vcc(2):vcc(3):vcc(4)]';   % Potentials to depolarize to (mV)
exptIKs = [];
%%%%%%%%%%%%%%%%%%%%%%% load in the data. %%%%%%%%%%%%%%%%%%%%%

activationStep = 21969;
deactivationStep=activationStep +1;

% read in some experimental data
 testexptdata = trace5p5;
% expdata = testexptdata;
 expttime = testexptdata(:,1);
 expdata = testexptdata(:,2);

missingtime = [123.5 expttime(1)];
activationexpttime = expttime([1:activationStep]);
deactivationexpttime = expttime([deactivationStep:end]);
exptIKstemp = expdata(:,1); % /max(expdata(:,1));
exptIKs = exptIKstemp([1:end],:);


% the following solution time matches the solution to the expt data.
%    [y1 y2]=ode23s(@hksrates,expttime,yp,[],ivv(1),p);
[y0 yi]=ode23s(@hksrates,missingtime,yp,[],vcc(4),p);
[y1 y2]=ode23s(@hksrates,activationexpttime,yi(length(y0),:),[],vcc(4),p);
[y00 yii]=ode23s(@hksrates,[activationexpttime(end) deactivationexpttime(1)],y2(length(y1),:),[],vcc(5),p); % instead of the capacitive current in the experiment
[y3 y4]=ode23s(@hksrates,deactivationexpttime,yii(length(y00),:),[],vcc(5),p);
%     GKs(i)=max(expdata(:,i))/((ivv(1)-erev)*sum(y2(end,16:17),2));
Iact=sum(y2(:,16:17),2)*(vcc(4)-erev);
Itail=sum(y4(:,16:17),2)*(vcc(5)-erev);

GKs=max(expdata)/((max(ivv)-erev)*sum(y2(end,16:17),2));
Iact=p(16)*Iact;
Itail=p(16)*Itail;

data(:,1)=[y1;y3];
data(:,2)=[Iact;Itail]; %/max(max([Iact;Itail])); 

if(printOutput==1)

    h = figure;%(1);                                  % Plot
    f1 = plot(data(:,1),data(:,2));
    hold on
    f2 = plot([activationexpttime;deactivationexpttime],exptIKs);

    if CaseCase == 0
        plotName = sprintf('Control 5.5s pulse');
    end
    if CaseCase == 1
        plotName = sprintf('Acute DHA 5.5s pulse');
    end
    if CaseCase == 2
        plotName = sprintf('Acute EPA 5.5s pulse');
    end
    if CaseCase == 3
        plotName = sprintf('Chronic EPA 5.5s pulse');
    end

    title(plotName,'FontSize',18);
    xlabel('Time (ms)','FontSize',16);
    ylabel('I_{Ks} (arbitrary units) ','FontSize',16);
    set(f1,'Color',[0.5 0 1],'LineWidth',2);
    set(f2,'Color',[1 0 1],'LineWidth',2);
    set(gca,'FontSize',12,'FontName','Roman');
    % set(gca,'XLim',[0 8000],'YLim',[-0.05 1.05]);
    set(gca,'XTick',[0:2000:8000],'FontName','Symbol','FontSize',16);
    % set(gca,'YTick',[0:0.5:1],'FontName','Symbol','FontSize',16);

    filename = sprintf('png/5p5fig_%d.png',filecounter);
    saveas(h,filename,'png');
    pause
    % clf
end

cost = 0;
if CostCost==0
    cost = (data(:,2) - exptIKs)'*(data(:,2) - exptIKs);
    cost = cost/length(data);
elseif CostCost==1
    activationWeight = 0.1;
    deactivationWeight = 10;
    
    
    simulation=[activationWeight*data(1:activationStep,2); deactivationWeight*data(activationStep+1:end,2)];
    target=[activationWeight*exptIKs(1:activationStep); deactivationWeight*exptIKs(activationStep+1:end)];    
    cost = (simulation - target)'*(simulation - target);
    
    cost=cost/sum([activationWeight*ones(1,length(exptIKs(1:activationStep))) deactivationWeight*ones(1,length(exptIKs(deactivationStep:end)))]);
elseif CostCost==2
   DexptIKsdt=zeros(size(exptIKs));
   for i=1:length(data)-1
       if i<activationStep
          wgh=0.1;
       else 
          wgh = 1;
       end
       DexptIKsdt(i+1,:)=wgh*(data(i+1,2)-data(i,2))./(data(i+1,1)-data(i,1));
   end
   weightD = abs(DexptIKsdt)+2;
   simulation=[.1*data(1:activationStep,2); data(activationStep+1:end,2)];
   target=[.1*exptIKs(1:activationStep); exptIKs(activationStep+1:end)];

   cost = ((simulation - target).*(weightD))'*((simulation - target).*(weightD));     
  
   cost=cost/sum(weightD);
else
    error('Specify the cost function to use')
end



% time is in ms, voltage is in mV, current is whatever Gks x mV is.
function [cost] = cost12(p,yp)


global filecounter;
global CaseCase;
global printOutput;
global CostCost;
global trace12s;

global erev;
global vcc;

% the IV protocol.
ivv=[vcc(2):vcc(3):vcc(4)]';   % Potentials to depolarize to (mV)
exptIKs = [];
%%%%%%%%%%%%%%%%%%%%%%% load in the data. %%%%%%%%%%%%%%%%%%%%%


activationStep = 49945;
deactivationStep=activationStep +1;

% read in some experimental data
 testexptdata = trace12s;
% expdata = testexptdata;
 expttime = testexptdata(:,1);
 expdata = testexptdata(:,2);

missingtime = [236.75 expttime(1)]; 
activationexpttime = expttime([1:activationStep]);
deactivationexpttime = expttime([deactivationStep:end]);
exptIKstemp = expdata(:,1); % /max(expdata(:,1));
exptIKs = exptIKstemp([1:end],:);


% the following solution time matches the solution to the expt data.
%    [y1 y2]=ode23s(@hksrates,expttime,yp,[],ivv(1),p);
[y0 yi]=ode23s(@hksrates,missingtime,yp,[],vcc(4),p);
[y1 y2]=ode23s(@hksrates,activationexpttime,yi(length(y0),:),[],vcc(4),p);
[y00 yii]=ode23s(@hksrates,[activationexpttime(end) deactivationexpttime(1)],y2(length(y1),:),[],vcc(5),p); % instead of the capacitive current in the experiment
[y3 y4]=ode23s(@hksrates,deactivationexpttime,yii(length(y00),:),[],vcc(5),p);
%     GKs(i)=max(expdata(:,i))/((ivv(1)-erev)*sum(y2(end,16:17),2));
Iact=sum(y2(:,16:17),2)*(vcc(4)-erev);
Itail=sum(y4(:,16:17),2)*(vcc(5)-erev);

% i5p5=21970;

% GKs=expdata(i5p5)/((max(ivv)-erev)*sum(y2(i5p5,16:17),2));
Iact=p(16)*Iact;
Itail=p(16)*Itail;

data(:,1)=[y1;y3];
data(:,2)=[Iact;Itail]; %/max(max([Iact;Itail])); 

if(printOutput==1)

    h = figure(1);                                  % Plot
    f1 = plot(data(:,1),data(:,2));
    hold on
    f2 = plot([activationexpttime;deactivationexpttime],exptIKs);

    if CaseCase == 0
        plotName = sprintf('Control 12s pulse');
    end
    if CaseCase == 1
        plotName = sprintf('Acute DHA 12s pulse');
    end
    if CaseCase == 2
        plotName = sprintf('Acute EPA 12s pulse');
    end



    title(plotName,'FontSize',18);
    xlabel('Time (ms)','FontSize',16);
    ylabel('I_{Ks} (arbitrary units) ','FontSize',16);
    set(f1,'Color',[0.5 0 1],'LineWidth',2);
    set(f2,'Color',[1 0 1],'LineWidth',2);
    set(gca,'FontSize',12,'FontName','Roman');
    % set(gca,'XLim',[0 8000],'YLim',[-0.05 1.05]);
    set(gca,'XTick',[0:5000:20000],'FontName','Symbol','FontSize',16);
    % set(gca,'YTick',[0:0.5:1],'FontName','Symbol','FontSize',16);

    filename = sprintf('png/12fig_%d.png',filecounter);
    saveas(h,filename,'png');
    pause
    % clf
end

cost=0;
if CostCost==0
    cost = (data(:,2) - exptIKs)'*(data(:,2) - exptIKs);
    cost = cost/length(data);
elseif CostCost==1
          
    activationWeight5p5 = 0.1;
    activationWeight12 = 1;
    deactivationWeight =1;
    activationStep5p5=21978;
    simulation=[activationWeight5p5*data(1:activationStep5p5,2); activationWeight12*data(activationStep5p5+1:activationStep,2); deactivationWeight*data(activationStep+1:end,2)];
    target=[activationWeight5p5*exptIKs(1:activationStep5p5); activationWeight12*exptIKs(activationStep5p5+1:activationStep); deactivationWeight*exptIKs(activationStep+1:end)];    
    cost = (simulation - target)'*(simulation - target);
    
    cost=cost/sum([activationWeight5p5*ones(1,length(exptIKs(1:activationStep5p5))) activationWeight12*ones(1,length(exptIKs(activationStep5p5+1:activationStep))) deactivationWeight*ones(1,length(exptIKs(deactivationStep:end)))]);
elseif CostCost==2
   DexptIKsdt=zeros(size(exptIKs));
   activationStep5p5=21978;
   for i=1:length(data)-1       
       if i<activationStep5p5
          wgh=0.1;
       else
          wgh = 10;
       end
       DexptIKsdt(i+1,:)=(data(i+1,2)-data(i,2))./(data(i+1,1)-data(i,1));
   end
   weightD = abs(DexptIKsdt)+10;
   simulation=[data(1:activationStep,2); data(activationStep+1:end,2)];
   target=[exptIKs(1:activationStep); exptIKs(activationStep+1:end)];

   cost =((simulation - target).*(weightD))'*((simulation - target).*(weightD));     
  
   cost=cost/sum(weightD);
else
    error('Specify the cost function to use')
end
