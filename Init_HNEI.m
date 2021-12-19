%% Initial parameters for HNEI dataset
load('HNEI_cell.mat');
para.Qmax=2.8;% rated capacity
para.V_star=3.7;
para.V_end=4.29;
para.Stride=1;% stride size
para.Cha_interval=0.01;% voltage interval
Train_cell=1;
Test_cell=8;
seq_Num=12;% number of segments
para.seg_length=round((para.V_end-para.V_star)/para.Cha_interval)+1-seq_Num+1;
input_size=[para.seg_length,2,1];