clear;close all;clc;
% input parameters  
run Init_HNEI.m
% run Init_NCA.m
% run Init_NMC.m
% run Init_LFP.m
%% use one cell for training and one cell for test
% produce training samples
Cyc=Cell(Train_cell).Cyc;
[Sample_seg,Label_seg]=feature(Cyc,para);
index=find(Label_seg(1,:)~=0);% exclude zero
Label_seg=Label_seg(:,index);
Sample_seg=Sample_seg(:,:,index);
clear Sample_seg_filter Label_seg_filter 
for i=1:size(Sample_seg,2)
    temp=reshape(Sample_seg(2,i,:),[1,size(Sample_seg(2,i,:),3)]);
    Sample_seg_filter(1,i,:)=myfilter(temp);% filtering 
    temp=Sample_seg(3,i,:);
    temp=reshape(Sample_seg(3,i,:),[1,size(Sample_seg(3,i,:),3)]);
    Sample_seg_filter(2,i,:)=myfilter(temp);% filtering 
    Label_seg_filter(i,:)=myfilter(Label_seg(i,:));% filtering
end
xt=reshape(Sample_seg_filter,[2,size(Sample_seg_filter,2)*size(Sample_seg_filter,3)]);xt=xt';
yt=reshape(Label_seg_filter,1,[]);yt=yt';
% produce test samples
Cyc=Cell(Test_cell).Cyc;
[Sample_seg,Label_seg]=feature(Cyc,para);
index=find(Label_seg(1,:)~=0);% exclude zero
Label_seg=Label_seg(:,index);
Sample_seg=Sample_seg(:,:,index);
clear Sample_seg_filter Label_seg_filter 
for i=1:size(Sample_seg,2)
    temp=reshape(Sample_seg(2,i,:),[1,size(Sample_seg(2,i,:),3)]);
    Sample_seg_filter(1,i,:)=myfilter(temp);% filtering
    temp=Sample_seg(3,i,:);
    temp=reshape(Sample_seg(3,i,:),[1,size(Sample_seg(3,i,:),3)]);
    Sample_seg_filter(2,i,:)=myfilter(temp);% filtering
    Label_seg_filter(i,:)=myfilter(Label_seg(i,:));% filtering
end
xs=reshape(Sample_seg_filter,[2,size(Sample_seg_filter,2)*size(Sample_seg_filter,3)]);xs=xs';
ys=reshape(Label_seg_filter,1,[]);ys=ys';
%% MLR
coe=regress(yt,[ones(length(xt(:,1)),1) xt]);% regression model
yt_MLR=coe(1)+xt*coe(2:end);
ys_MLR=coe(1)+xs*coe(2:end);
my_plot(xt,yt*100,yt_MLR*100);
my_plot(xs,ys*100,ys_MLR*100);
train_err = (yt_MLR-yt);disp(mean(abs(train_err)));disp(sqrt(mse(train_err)));
test_err = (ys_MLR-ys);disp(mean(abs(test_err)));disp(sqrt(mse(test_err)));
%% MLR random application
load random_seg_num.mat
L=size(Sample_seg_filter,3);
r_index=r_index(1:L);
for i=1:size(Sample_seg_filter,3)
xs_rand(i,:)=Sample_seg_filter(:,r_index(i),i);
end
ys_rand=Label_seg_filter(1,:)';
ys_rand_est=coe(1)+xs_rand*coe(2:end);
test_rand_err =ys_rand_est-ys_rand;
disp({'random test MAE',mean(abs(test_rand_err))});
disp({'random test RMSE',sqrt(mse(test_rand_err))});
my_plot_rand(xs_rand,ys_rand*100,ys_rand_est*100);
function [Sample_seg,Label_seg]=feature(Cyc,para)
%% general seting
Qmax=para.Qmax;% rated capacity
V_star=para.V_star;%
V_end=para.V_end;
Stride=para.Stride;% stride size 
Cha_interval=para.Cha_interval;% voltage interval
Vmin=V_star;
Vmax=V_end;
seg_length=para.seg_length;% length of each segment
for i=1:length(Cyc)
    Seq_min=max(find(Cyc(i).V<=V_star));
    Seq_max=min(find(Cyc(i).V>=V_end));
    if ((isempty(Seq_min))||(isempty(Seq_max))||(Seq_max<=Seq_min))
        continue;
    end
    Seq=[Seq_min:Seq_max];    
    for k=Seq_min:Seq_max-1
        if(Cyc(i).V(k)>=Cyc(i).V(k+1))
            Cyc(i).V(k+1)=Cyc(i).V(k)*(1+1e-4);% to keep monotonicity
        end
    end
    V=[V_star:Cha_interval:V_end];
    Q=interp1(Cyc(i).V(Seq),Cyc(i).Q(Seq),V,'linear','extrap');
    V=(V-Vmin)/(Vmax-Vmin);
    Q=Q/Qmax;
    if(length(Q)<seg_length)
        continue;
    end
    cnt1=1;
    for j=seg_length:Stride:length(V)
        Qtemp=Q(j-seg_length+1:j)-Q(j-seg_length+1);% calculate capacity increment
        Sample_seg(:,cnt1,i)=[mean(V(j-seg_length+1:j));mean(Qtemp);std(Qtemp)]'; 
        Label_seg(cnt1,i)=Cyc(i).Ca/Qmax;  
        cnt1=cnt1+1;
    end
end
end
function y=myfilter(x)
W = fspecial('gaussian',[5,5],2);
y=imfilter(x,W,'replicate');
end
function my_plot(xt,yt,yt_MLR)
red1=[255 111 94]/255;
green1=[126 186 162]/255;
Lw=2;fs=10;makersize=2;
x=[1:length(xt)];
figure;set(gcf,'Position',[100,300,300,200], 'color','w');
plot(x,yt_MLR,'-o','color',red1,'MarkerSize',makersize,'Linewidth',Lw/2);hold on;
plot(x,yt,'b-.','MarkerSize',makersize,'Linewidth',Lw);set(gca, 'LineWidth',1);
xlim([min(x) max(x)]);xlabel('Sample');ylabel('SOH (%)');%ylim([70 100]);
ax=gca;ax.FontSize=fs;ax.FontName='Arial';legend('MLR','Real');
end
function my_plot_rand(xt,yt,yt_MLR)
red1=[255 111 94]/255;
green1=[126 186 162]/255;
Lw=2;fs=10;makersize=2;
x=[1:length(xt)];
figure;set(gcf,'Position',[100,300,300,200], 'color','w');
plot(x,yt_MLR,'-o','color',red1,'MarkerSize',makersize,'Linewidth',Lw/2);hold on;
plot(x,yt,'b-.','MarkerSize',makersize,'Linewidth',Lw);set(gca, 'LineWidth',1);
xlim([min(x) max(x)]);xlabel('Cycle');ylabel('SOH (%)');%ylim([70 100]);
ax=gca;ax.FontSize=fs;ax.FontName='Arial';legend('MLR','Real');
end