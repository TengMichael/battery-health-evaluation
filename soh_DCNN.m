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
Sample_seg=Sample_seg(:,:,:,index);
Sample_seg(:,2,:,:)=myfilter_DCNN(Sample_seg(:,2,:,:));% filtering 
Label_seg_filter=myfilter(Label_seg);% filtering 
xt=reshape(Sample_seg,input_size(1),input_size(2),input_size(3),[]);
yt=reshape(Label_seg_filter,1,[]);yt=yt';
% produce test samples
Cyc=Cell(Test_cell).Cyc;
[Sample_seg,Label_seg]=feature(Cyc,para);
index=find(Label_seg(1,:)~=0);% exclude zero
Label_seg=Label_seg(:,index);
Sample_seg=Sample_seg(:,:,:,index);
Sample_seg(:,2,:,:)=myfilter_DCNN(Sample_seg(:,2,:,:));% filtering 
Label_seg_filter=myfilter(Label_seg);% filtering 
xs=reshape(Sample_seg,input_size(1),input_size(2),input_size(3),[]);
ys=reshape(Label_seg_filter,1,[]);ys=ys';
%% DCNN
gpuDevice(1);% GPU
layers = [...
    imageInputLayer(input_size); % input
    convolution2dLayer([floor(input_size(1)/2) 2],10);
    batchNormalizationLayer(); 
    reluLayer()% Relu
    convolution2dLayer([floor(input_size(1)/2)-2 1],10);
    batchNormalizationLayer(); 
    reluLayer()% Relu
    maxPooling2dLayer([4 1]);
    fullyConnectedLayer(1);
    batchNormalizationLayer();
    regressionLayer,...
    ];
options = trainingOptions('sgdm',...
    'MaxEpochs',30,...
    'Plots','training-progress');
net_cnn = trainNetwork(xt,yt,layers,options);
%% DCNN test
close all;
yt_DCNN = predict(net_cnn,xt);
train_err =yt_DCNN-yt;
ys_DCNN = predict(net_cnn,xs);
test_err =ys_DCNN-ys;
my_plot(xt,yt*100,yt_DCNN*100);
my_plot(xs,ys*100,ys_DCNN*100);
disp(mean(abs(train_err)));
disp(sqrt(mse(train_err)));
disp(mean(abs(test_err)));
disp(sqrt(mse(test_err)));
%% DCNN random application
load random_seg_num.mat
L=size(Sample_seg,4);
r_index=r_index(1:L);
for i=1:size(Sample_seg,4)
xs_rand(:,:,1,i)=Sample_seg(:,:,r_index(i),i);
end
ys_rand=Label_seg_filter(1,:)';
ys_rand_est = predict(net_cnn,xs_rand);
test_rand_err =ys_rand_est-ys_rand;
disp({'random test MAE',mean(abs(test_rand_err))});
disp({'random test RMSE',sqrt(mse(test_rand_err))});
my_plot_rand(xs_rand,ys_rand*100,ys_rand_est*100);
%%
function [Sample_seg,Label_seg]=feature(Cyc,para)
Qmax=para.Qmax;% rated capacity
V_star=para.V_star;
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
        Sample_seg(:,:,cnt1,i)=[V(j-seg_length+1:j);Qtemp]';
        Label_seg(cnt1,i)=Cyc(i).Ca/Qmax;
        cnt1=cnt1+1;
    end
end
end
%%
function y=myfilter(x)
W = fspecial('gaussian',[5,5],2);
y=imfilter(x,W,'replicate');
end
function y=myfilter_DCNN(x)
seq_Num=size(x,3);
seq_length=size(x,1);
for i=1:size(x,4)
   tmp(1:seq_length,i)=reshape(x(:,1,1,i),[],1);
   for j=seq_length+1:seq_Num+seq_length-1
       tmp(j,i)=x(end,1,j-seq_length+1,i); 
   end
end
for i=1:size(tmp,1)
   tmp_filter(i,:)=myfilter(tmp(i,:));
end
for i=1:size(tmp,2)
   for j=1:seq_Num
       tmp_out(:,1,j,i)=tmp_filter(j:end-seq_Num+j,i);
   end
end
y=tmp_out;
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
ax=gca;ax.FontSize=fs;ax.FontName='Arial';legend('DCNN','Real');
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
ax=gca;ax.FontSize=fs;ax.FontName='Arial';legend('DCNN','Real');
end