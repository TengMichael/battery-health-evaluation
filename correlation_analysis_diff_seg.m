clear;close all;clc;
run Init_HNEI.m
% run Init_NCA.m
% run Init_NMC.m 
% run Init_LFP.m
seg_max=60;
inter_num=59;
%%
corr_mean=zeros(inter_num,inter_num);
corr_std=zeros(inter_num,inter_num);
Cyc=Cell(Train_cell).Cyc;
for seg=1:inter_num
    para.Num_seq=seg_max-seg+1;
    [Sample_seg,Label_seg]=feature(Cyc,para);
    cnt=1;
    for i=1:size(Sample_seg,3)
        if (Label_seg(1,i)~=0)
            deltaQ_mean(cnt,:)=Sample_seg(2,:,i);% mean
            deltaQ_std(cnt,:)=Sample_seg(3,:,i);% std
            soh(cnt,1)=Label_seg(1,i);
            cnt=cnt+1;
        end
    end
    tem=corrcoef([soh,deltaQ_mean]);
    corr_mean(inter_num:-1:inter_num-seg+1,seg)=roundn(tem(2:end,1),-4);
    corr_mean_ave(seg)=mean(abs(tem(2:end,1)));
    tem=corrcoef([soh,deltaQ_std]);
    corr_std(inter_num:-1:inter_num-seg+1,seg)=roundn(tem(2:end,1),-4);
    corr_std_ave(seg)=mean(abs(tem(2:end,1)));
    clear deltaQ_mean deltaQ_std soh
end
corr_mean(corr_mean==0)=nan;
corr_std(corr_std==0)=nan;
%%
L=10;yellow=[243 219 32];green=[39 182 169];blue=[71 73 255];
for i=1:3
    delta=(yellow(i)-green(i))/L;
    colordata1(:,i)=[yellow(i):-delta:green(i)]/255;
end
for i=1:3
    delta=(green(i)-blue(i))/L;
    colordata2(:,i)=[green(i):-delta:blue(i)]/255;
end
colordata=[colordata1;colordata2];
colornan=[65, 72, 74]/255;
%%
close all;
tmp=[59:-1:1];
for i=1:length(tmp)
    tmpx(i)=nan;
    if (i==1||i==59||mod(i,10)==0)continue;end
    tmp(i)=nan;
end
colornan=[223, 230, 233]/255;
figure;set(gcf,'Position',[0,0,1200,600], 'color','w');Lw=1;fs=10;
h=heatmap(abs(corr_mean),'MissingDataColor','w','FontName','Arial','GridVisible','off');
colormap((colordata));h.YDisplayLabels=tmp;h.XDisplayLabels=tmpx;h.FontSize=18;
figure;set(gcf,'Position',[0,0,1200,600], 'color','w');Lw=1;fs=10;
h=heatmap(abs(corr_std),'MissingDataColor','w','FontName','Arial','GridVisible','off');
colormap((colordata));h.YDisplayLabels=tmp;h.XDisplayLabels=tmpx;h.FontSize=18;
%%
function [Sample_seg,Label_seg]=feature(Cyc,para)
Qmax=para.Qmax;% rated capacity
V_star=para.V_star;
V_end=para.V_end;
Stride=para.Stride;%
Cha_interval=para.Cha_interval;%
Vmin=V_star;
Vmax=V_end;
Num_seq=para.Num_seq;
cnt=1;
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
        end
    end
    V=[V_star:Cha_interval:V_end];
    Q=interp1(Cyc(i).V(Seq),Cyc(i).Q(Seq),V,'linear','extrap');
    if(length(Q)<Num_seq)
        continue;
    end
    cnt1=1;
    for j=Num_seq:Stride:length(V)
        Qtemp=Q(j-Num_seq+1:j)-Q(j-Num_seq+1);
        Sample_seg(:,cnt1,i)=[mean(V(j-Num_seq+1:j));mean(Qtemp);std(Qtemp)]';
        Label_seg(cnt1,i)=Cyc(i).Ca/Qmax;
        cnt1=cnt1+1;
    end
end
end