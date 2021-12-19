%% Initial parameters for HNEI dataset
load('SNL_NMC.mat');
para.Qmax=3.0;% rated capacity
para.V_star=3.6;
para.V_end=4.19;
para.Stride=1;% stride size
para.Cha_interval=0.01;% voltage interval
Train_cell=7;
Test_cell=9;
seq_Num=12;% number of segments
para.seg_length=round((para.V_end-para.V_star)/para.Cha_interval)+1-seq_Num+1;
input_size=[para.seg_length,2,1];