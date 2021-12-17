clear;close all;clc;
path='datasets\';
File =dir(fullfile(path,'*.csv'));
FileNames = {File.name}';
%% extract HNEI datasets
filename=FileNames(1:15);% select files of HNEI datasets
for i=1:length(filename)
    A=xlsread([path filename{i}]);
    A(:,1)=[];
    data.time=A(:,1);
    data.cycle=A(:,2);
    data.current=A(:,3);
    data.voltage=A(:,4);
    data.chaQ=A(:,5);
    data.disQ=A(:,6);
    HNEI_cell(i).data=data;
    HNEI_cell(i).Cellname=filename{i};
    i
end
HNEI_cell(10)=[];% exclude empty data
%% extract SNL datasets
filename=FileNames(16:end);% select files of SNL datasets
for i=1:length(filename)
    A=xlsread([path filename{i}]);
    A(:,1)=[];
    data.time=A(:,1);
    data.cycle=A(:,2);
    data.current=A(:,3);
    data.voltage=A(:,4);
    data.chaQ=A(:,5);
    data.disQ=A(:,6);
    SNL_cell(i).data=data;
    SNL_cell(i).Cellname=filename{i};
    i
end
%% extract charging data and capacity of each cycle for HNEI cells
Qmax=2.8;% Ah nominal capacity of HNEI (NMC-LCO) cells
for i=1:length(HNEI_cell)
    data=HNEI_cell(i).data;
    CycMax=max(data.cycle);
    cnt=1;clear Cyc Label;
    for j=1:CycMax
        Num=find(data.cycle==j);
        if(isempty(Num))continue;end
        tmp.current=data.current(Num);
        tmp.voltage=data.voltage(Num);
        tmp.chaQ=data.chaQ(Num);
        tmp.time=data.time(Num);
        tmp.Ca=max(data.disQ(Num));
        if(((tmp.Ca/Qmax)>1.2)||((tmp.Ca/Qmax)<0.7))continue;end
        FindNum=find(tmp.current>0);
        if(isempty(FindNum))continue;end
        During_t=tmp.time(FindNum(end))-tmp.time(FindNum(1));
        if(length(FindNum)<100||During_t>1.3e4)continue;end
        Cyc(cnt).I=tmp.current(FindNum);
        Cyc(cnt).V=tmp.voltage(FindNum);
        Cyc(cnt).Q=tmp.chaQ(FindNum);
        Cyc(cnt).t=tmp.time(FindNum);Cyc(cnt).t=Cyc(cnt).t-Cyc(cnt).t(1);
        Cyc(cnt).Ca=tmp.Ca;
        Label(cnt)=Cyc(cnt).Ca/Qmax;
        cnt=cnt+1;
    end
    Cell(i).Cyc=Cyc;Cell(i).Label=Label;Cell(i).Cellname=HNEI_cell(i).Cellname;
end
save HNEI_cell.mat Cell
%% extract charging data and capacity of each cycle for SNL cells
Qmax=1.1;% Ah nominal capacity of LFP cells
MultiCell=SNL_cell([1:15,25:30]);% LFP 100% DOD
Cell=SNL_data_extract(MultiCell,Qmax);
save LFP_cell.mat Cell
Qmax=3.2;% Ah nominal capacity of NCA cells
MultiCell=SNL_cell([31:42,49:54]);% NCA 100% DOD
Cell=SNL_data_extract(MultiCell,Qmax);
save NCA_cell.mat Cell
Qmax=3;% Ah nominal capacity of NMC cells
MultiCell=SNL_cell([55:70,81:86]);% NCM 100% DOD
Cell=SNL_data_extract(MultiCell,Qmax);
save NMC_cell.mat Cell
%%
function Cell=SNL_data_extract(MultiCell,Qmax)
for i=1:length(MultiCell)
    data=MultiCell(i).data;
    CycMax=max(data.cycle);
    cnt=1;clear Cyc Label;
    for j=1:CycMax
        Num=find(data.cycle==j);
        if(isempty(Num))continue;end
        tmp.current=data.current(Num);
        tmp.voltage=data.voltage(Num);
        tmp.chaQ=data.chaQ(Num);
        tmp.time=data.time(Num);
        tmp.Ca=max(data.disQ(Num));
        if(((tmp.Ca/Qmax)>1.2)||((tmp.Ca/Qmax)<0.5))continue;end
        FindNum=find(tmp.current>0);
        if(isempty(FindNum))continue;end
        During_t=tmp.time(FindNum(end))-tmp.time(FindNum(1));
        if(length(FindNum)<25||During_t<5e3)continue;end
        Cyc(cnt).I=tmp.current(FindNum);
        Cyc(cnt).V=tmp.voltage(FindNum);
        Cyc(cnt).Q=tmp.chaQ(FindNum);
        Cyc(cnt).t=tmp.time(FindNum);Cyc(cnt).t=Cyc(cnt).t-Cyc(cnt).t(1);
        Cyc(cnt).Ca=tmp.Ca;
        Label(cnt)=Cyc(cnt).Ca/Qmax;
        cnt=cnt+1;
    end
    Cell(i).Cyc=Cyc;Cell(i).Label=Label;Cell(i).Cellname=MultiCell(i).Cellname;
end
end