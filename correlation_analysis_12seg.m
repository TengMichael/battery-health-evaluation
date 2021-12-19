clear;clc;close all;
run Init_HNEI.m
% run Init_NCA.m
% run Init_NMC.m
% run Init_LFP.m
%% extract features
Cyc=Cell(Train_cell).Cyc;
[Sample_seg,Label_seg]=feature(Cyc,para);
index=find(Label_seg(1,:)~=0);
Label_seg=Label_seg(:,index);
Sample_seg=Sample_seg(:,:,index);
L=size(Sample_seg,4);
pink=[255 179 204];lightblue=[0 115 189];
for i=1:3
    delta=(pink(i)-lightblue(i))/L;
    colordata(:,i)=[pink(i):-delta:lightblue(i)]/255;
end
close all;clc;cnt=1;
for i=1:size(Sample_seg,3)
    if (Label_seg(1,i)~=0)
        deltaQ_mean(cnt,:)=Sample_seg(2,:,i);% mean
        deltaQ_std(cnt,:)=Sample_seg(3,:,i);% std
        soh(cnt,1)=Label_seg(1,i);
        cnt=cnt+1;
    end
end
corr_mean=corrcoef([soh,deltaQ_mean]);
corr_std=corrcoef([soh,deltaQ_std]);
%%
close all;
figure;set(gcf,'Position',[100,100,900,800], 'color','w');Lw=1;fs=14;
for i=1:size(Sample_seg,2)
    x=deltaQ_mean(:,i);y=soh;p=polyfit(x,y,1);f=polyval(p,(x));
    subplot(4,3,i);plot(x,y,'o');hold on;plot(x,f,'r','LineWidth',1);xlabel('{\itave}\_Δ{\itQ_{seg}}');ylabel('SOH');
    ax=gca;ax.FontSize=fs;ax.FontName='Arial';temp=num2str(roundn(corr_mean(i,2),-3));%title(['{\itρ}=' temp],'Color',[8,46,84]/255);
    xlim(roundn([min(x) max(x)],-2));
    ylim(roundn([min(y) max(y)],-2));
    xticks(roundn(linspace(min(x),max(x),4),-2));
    yticks(roundn(linspace(min(y),max(y),3),-2));
    set(gca, 'LineWidth',1);
    set(gca,'XDir','reverse');
    ax=gca;ax.FontSize=fs;ax.FontName='Arial';
    text(max(x)-0.01,min(y)+0.05,['{\itρ}=' temp],'Color',[8,46,84]/255,'linewidth',1,'FontSize',fs,'FontName','Arial');    
    RemoveSubplotWhiteArea(gca, 4, 3, floor((i-1)/3)+1, mod(i-1,3)+1); % 去除空白部分
end
figure;set(gcf,'Position',[100,100,900,800], 'color','w');Lw=1;fs=14;
for i=1:size(Sample_seg,2)
    x=deltaQ_std(:,i);y=soh;p=polyfit(x,y,1);f=polyval(p,(x));
    subplot(4,3,i);plot(x,y,'o');hold on;plot(x,f,'r','LineWidth',1);xlabel('{\itstd}\_Δ{\itQ_{seg}}');ylabel('SOH');
    ax=gca;ax.FontSize=fs;ax.FontName='Arial';temp=num2str(roundn(corr_std(i,2),-3));%title(['{\itρ}=' temp],'Color',[8,46,84]/255);
    xlim(roundn([min(x) max(x)],-2));
    ylim(roundn([min(y) max(y)],-2));
    xticks(roundn(linspace(min(x),max(x),4),-2));
    yticks(roundn(linspace(min(y),max(y),3),-2));
    set(gca, 'LineWidth',1);
    set(gca,'XDir','reverse');
    ax=gca;ax.FontSize=fs;ax.FontName='Arial';
    text(max(x)-0.005,min(y)+0.05,['{\itρ}=' temp],'Color',[8,46,84]/255,'linewidth',1,'FontSize',fs,'FontName','Arial');    
    RemoveSubplotWhiteArea(gca, 4, 3, floor((i-1)/3)+1, mod(i-1,3)+1); % 去除空白部分
end
function [Sample_seg,Label_seg]=feature(Cyc,para)
Qmax=para.Qmax;% rated capacity
V_star=para.V_star;
V_end=para.V_end;
Stride=para.Stride;%
Cha_interval=para.Cha_interval;%
Vmin=V_star;
Vmax=V_end;
seg_length=para.seg_length;
for i=1:length(Cyc)
    Seq_min=max(find(Cyc(i).V<=V_star));
    Seq_max=min(find(Cyc(i).V>=V_end));
    if ((isempty(Seq_min))||(isempty(Seq_max))||(Seq_max<=Seq_min))
        continue;
    end
    Seq=[Seq_min:Seq_max];    
    for k=Seq_min:Seq_max-1
        if(Cyc(i).V(k)>=Cyc(i).V(k+1))
            Cyc(i).V(k+1)=Cyc(i).V(k)*(1+1e-4);
        end%
    end
    V=[V_star:Cha_interval:V_end];
    Q=interp1(Cyc(i).V(Seq),Cyc(i).Q(Seq),V,'linear','extrap');
    if(length(Q)<seg_length)
        continue;
    end
    cnt1=1;
    for j=seg_length:Stride:length(V)
        Qtemp=Q(j-seg_length+1:j)-Q(j-seg_length+1);
        Sample_seg(:,cnt1,i)=[mean(V(j-seg_length+1:j));mean(Qtemp);std(Qtemp)]'; 
        Label_seg(cnt1,i)=Cyc(i).Ca/Qmax;  
        cnt1=cnt1+1;
    end
end
end