clear all;
close all;

load('AllData.db','-mat');

%%

cv7=(data.covid(8:end)-data.covid(1:end-7));
crs=[];
ps=[];
aux=[];
vars=fields(data);
vars(1:2)=[];
for j=1:numel(vars)
    aux(j,:)=data.(vars{j})(8:end)-data.(vars{j})(1:end-7);
    [cr,p]=corr(aux(j,:)',cv7(:),'rows','pairwise');
    
    lm=fitlm(aux(j,:)',cv7(:));
    for k=-20:0
        [crs(j,k+21),ps(j,k+21)]=corr(cv7(max(1,1-k):min(end,end-k))',aux(j,max(1,1+k):min(end,end+k))','rows','pairwise');       
    end
    
end
%%  Correlation with meteorological variables for all lags
figure;
subplot(211)
kbins=-20:0;
plot(kbins,crs','LineWidth',2);

xlabel('Days')
ylabel('Correlation')

subplot(212)
semilogy(kbins,(ps'),'LineWidth',2);
hold on
plot([-20 0],.01*[1 1],'k--')
set(gca,'YTick',10.^[-20:5:20],'XTick',-20:5:20);
ylim([10^-14, 1])
legend(vars,'Location','NW');
legend boxoff
xlabel('Days')
ylabel('p-value')
printps(gcf,'Figure1B');

%% Co-evolution between covid and RH

fn=845;
dis=9;
figure;
RH7=(data.RH(8:end)-data.RH(1:end-7));
RH7_dis=RH7(fn-dis:end-dis);
 
    
yyaxis left
px=plot(cv7(fn:end),'Linewidth',2);
ylim([-500 700])
ylabel(['Δ COVID [cases]'])
set(gca,'YTick',-2000:250:2000);
    

yyaxis right
py=plot(-RH7_dis,'Linewidth',2);
ylim([-70 90]);
set(gca,'YTick',-200:50:100);
    
xlabel('Time')
ylabel(['-Δ RH [%]'])
xlim([0 Inf])

    mths=datevec(data.time(fn+7:end));
[v,indx]=unique(mths(:,2));
set(gca,'XTick',indx,'XTickLabel',datestr(datetime(1,v,1),'mmm'))
legend('Δ COVID(t)','-Δ RH(t-9 days)','Location','nw')
legend boxoff
     

%% Correlation with RH across months
figure;

caux=cv7(1+dis:end);
Raux=RH7(1:end-dis);

t=data.time(8+dis:end);
taux=datevec(t);
Rbins=-70:70;
for mths=3:10
   subplot(2,4,mths-2);
   sel= taux(:,1)==2020 & taux(:,2)==mths;
   sel=sel & ~isnan(Raux') & ~isnan(caux');
   plot(Raux(sel),caux(sel),'.','MarkerSize',20);
   hold on;
   p1=polyfit(Raux(sel),caux(sel),1);
   plot(Rbins,polyval(p1,Rbins'),'r');
   
%
    [cor,p]=corr(Raux(sel)',caux(sel)','rows','pairwise');
    lm=fitlm(Raux(sel)',caux(sel)');
   title([datestr(datetime(1,mths,1),'mmmm') ': ' trunc(cor,2)]);% ' - p: ' num2str(p)]);
   axis([70*[-1 1]  550*[-1 1]])
  
   sl(mths-2)=p1(1);
   in(mths-2)=p1(2);
   cr(mths-2)=cor;
   pv(mths-2)=p;
   
   axis square
end

%% Summary of monthly results
figure
subplot(131)
plot(sl,'.-','MarkerSize',30);
set(gca,'XTick',1:8,'xlim',[0 9],'ylim',[-6 1],'XTickLabel',datestr(datetime(1,(1:8)+2,1),'mmm'));
xlabel('Month')
ylabel('slope [patients / %]');

subplot(132)
plot(cr,'r.-','MarkerSize',30);
set(gca,'XTick',1:8,'xlim',[0 9],'ylim',[-.8 .3],'XTickLabel',datestr(datetime(1,(1:8)+2,1),'mmm'));
xlabel('Month')
ylabel('Correlation');

subplot(133)
semilogy((pv),'.-','Color',[.3 .5 0],'MarkerSize',30);
hold on;
plot([0 9],(0.05)*[1 1],':k');
plot([0 9],.01*[1 1],'--k');
set(gca,'XTick',1:8,'xlim',[0 9],'YTick',10.^[-6:2:0],'XTickLabel',datestr(datetime(1,(1:8)+2,1),'m'));
xlabel('Month')
ylabel('Significance');
ylim(10.^[-6 0])



%%

figure
lm=20;
win=-lm*1.5:lm*.5;

pksh=fliplr([ 100  200 ]);    
legs=[];
cls={'r','c'};
aux7(2,:)=data.RH;
aux7(4,:)=[NaN(1,7) cv7];
numpks=[];
for kp=1%1:numel(pksh)
    [~,pks]=findpeaks(aux7(4,:),'MinPeakDistance',7,'MinPeakHeight',pksh(kp));

    for j=[2 ]%1:4
        subplot(1,2,j/2)
                hold on;
                
        aux=[];
        for k=1:numel(pks)
            if pks(k)+lm<=numel(cv7)
                aux(k,:)=aux7(j,pks(k)+win);
            end
        end
        
        mn=nanmean(aux);
        st=nanstd(aux);
        a=area(win,[(mn-st);2*st]','EdgeColor','none');
        a(1).FaceColor='none';%[1 1 1];
        a(2).FaceColor=[cls{kp}];
        %a(2).FaceAlpha=.5;
        plot(win,mn,'k')
        % draw peaks of 
        sel=~isnan(cv7);
        %sel=mths>=5;
        mn=nanmean(aux7(j,sel));
        st=nanstd(aux7(j,sel));

        numpks(kp)=numel(pks);     
        plot(win([1 end]),(mn+st)*[1 1],'--k');
        plot(win([1 end]),(mn-st)*[1 1],'--k')

         set(gca,'box','on')    
    %     if j==1
    %         legs(1)=a(2);
    %         legs(2)=p;
    %         legs(3)=p2;
    %     end

        xlim([-20 0]);
        %ylim([50 100]);
    end
end

%% Cross-validation

crs=[];
ps=[];
aux=[];
for j=1:numel(vars)-1 
    aux(j,:)=data.(vars{j})(8:end)-data.(vars{j})(1:end-7);

    [cr,p]=corr(aux(j,:)',cv7(:),'rows','pairwise');
    
    lm=fitlm(aux(j,:)',cv7(:));
    for k=-25:25
        [crs(j,k+26),ps(j,k+26)]=corr(cv7(max(1,1-k):min(end,end-k))',aux(j,max(1,1+k):min(end,end+k))','rows','pairwise');       

    end
    
end

sel=6:26;
%% CX This is X but with every predictor displaced its own optimal amount
%CX=aux(1:end,1:end-dis);
CX=[];
paux=fliplr(-log10(ps(:,sel)')); % to get optimal dis
for j=1:size(aux,1)
    [~ , disj]=max(paux(j,:)); % 1 is dis=0
    disj=disj-1;
    CX(:,j)=[NaN(1,disj) aux(j,1:end-disj)];
end
y=cv7';
% CX(3,:)=[NaN(1,5) aux(3,1:end-dis-5)];
% CX(4,:)=[NaN(1,6) aux(4,1:end-dis-6)];
%CX=CX';
%y=cv7(10:end)';
    sel=all(~isnan([CX y]),2);
    CX=CX(sel,:); y=y(sel);

    %lm=fitlm(CX,y)
%     for j=[1 3:5]
%         lm=fitlm(CX(:,[j 2]),y,'linear');
%     end
    %lm=fitlm(CX(:,[2:4]),y,'linear')
%%

figure
warning('off','stats:LinearModel:RankDefDesignMat')
dis=9;
    if 1 % wether 9 days is used or the optimal displacement for each predictor
        X=aux(:,1:end-dis)';
        y=cv7(10:end)';
        sel=all(~isnan([X y]),2);  
        X=X(sel,:); y=y(sel);    
         %X=X(100:150,:);y=y(100:150);
    else
        X=CX;
    end
    
    indx=randperm(size(X,1));
    X=X(indx,:); y=y(indx);
    %X(:,2)=[];
    
cv_cum=NaN(1,size(X,2));
se_cum=cv_cum;
pars={'linear','interactions','poly3'};
var_cum=cell(1,size(X,2));
for depth=0:size(X,2)
    disp(depth)
    C = nchoosek(1:size(X,2),depth);
    var_cum{depth+1}=C(1,:);
    [cv_cum(depth+1),se_cum(depth+1)]=cross_validate(X(:,C(1,:)),y,10,pars{3});
    for jc=2:size(C,1)        
        [cv,se,~,lms]=cross_validate(X(:,C(jc,:)),y,10,pars{2});
        if cv<cv_cum(depth+1)
            cv_cum(depth+1)=cv;
            se_cum(depth+1)=se;
            var_cum{depth+1}=C(jc,:);
        end        
    end

end
%cv_cum=cv_cum.^.5
%plot(cv_cum)
errorbar(0:depth,cv_cum,se_cum);
vars(var_cum{2})
[mn,ind]=min(cv_cum);
hold on;
plot([-1 depth+1],(mn+se_cum(ind))*[1 1],'r--');
ylim([1.6 2.9]*10000)
xlim([-.5 8.5])
xlabel('Number of variables in model')
ylabel('best sum of squared errors')

% Best model
best_nr_var=find(cv_cum<mn+se_cum(ind),1,'first');
disp({['Cross-validation best model: ' num2str(best_nr_var-1) ' variables' ],char(vars(var_cum{best_nr_var}))});
%% Lasso

[B,FitInfo] = lasso(X,y,'Standardize',1,'CV',10,'Alpha',1);
lassoPlot(B,FitInfo,'PlotType','CV');
legend('show') % Show legend
%B(:,min(FitInfo.Index1SE,size(B,2)-1))'
disp('Lasso');
array2table(B(:,min(FitInfo.Index1SE,size(B,2)-1))','VariableNames',vars(1:8))

%%
function [cv,se,crs,lm]=cross_validate(X,y,k_fold,varargin)
    % varargin: to pass fitlm
    if ~exist('k_fold','var')
        k_fold=10;
    end
    
    % allow for poly555,etc
    if numel(varargin)>0
        if strcmp(varargin{1}(1:4),'poly')
            if size(X,2)==0                
                varargin{1}='linear';
            else
                varargin{1}=['poly' arrayfun(@(x) varargin{1}(5),1:size(X,2))];
            end
        end
    end
                
        
       
        
    
    N=size(X,1);
    x_st=round(linspace(1,N+1,k_fold+1));
    
    cv=NaN(1,k_fold);
    crs=cv;
    lm={};
    for j=1:k_fold
        test=x_st(j):x_st(j+1)-1;
        train=[1:test(1)-1 test(end)+1:N];
        
        lm{j}=fitlm(X(train,:),y(train),varargin{:});
        aux=predict(lm{j},X(test,:));
        cv(j)=nanmean((aux-y(test)).^2); 
        crs(j)=corr(aux,y(test),'rows','pairwise');
    end
    se=(nanstd(cv))/sqrt(k_fold);
    cv=nanmean(cv);
    crs=nanmean(crs);
end