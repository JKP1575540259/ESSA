% The script conducts the single simulation experiment to intuitively show the
% applicability of our approach on analyzing the incomplete time series with 
% different percentage of missing data (10% - 60%).

% DATA_1C - DATA_6C: the incomplete time series with missing percentage of
% 10% - 60%. The missing values are marked as NaN.

% FORM_1C - FORM_6C: the associated formal errors.

L = 365; % window size
k = 5;   % reconstruction order
path = 'D:\SuppMaterial\'; % The users can change the folder
DATA_1C = cell2mat(struct2cell(load(fullfile(path,'data','Single','DATA_1C.mat'))));
DATA_2C = cell2mat(struct2cell(load(fullfile(path,'data','Single','DATA_2C.mat'))));
DATA_3C = cell2mat(struct2cell(load(fullfile(path,'data','Single','DATA_3C.mat'))));
DATA_4C = cell2mat(struct2cell(load(fullfile(path,'data','Single','DATA_4C.mat'))));
DATA_5C = cell2mat(struct2cell(load(fullfile(path,'data','Single','DATA_5C.mat'))));
DATA_6C = cell2mat(struct2cell(load(fullfile(path,'data','Single','DATA_6C.mat'))));
FORM_1C = cell2mat(struct2cell(load(fullfile(path,'data','Single','FORM_1C.mat'))));
FORM_2C = cell2mat(struct2cell(load(fullfile(path,'data','Single','FORM_2C.mat'))));
FORM_3C = cell2mat(struct2cell(load(fullfile(path,'data','Single','FORM_3C.mat'))));
FORM_4C = cell2mat(struct2cell(load(fullfile(path,'data','Single','FORM_4C.mat'))));
FORM_5C = cell2mat(struct2cell(load(fullfile(path,'data','Single','FORM_5C.mat'))));
FORM_6C = cell2mat(struct2cell(load(fullfile(path,'data','Single','FORM_6C.mat'))));

signal_true = cell2mat(struct2cell(load(fullfile(path,'data','signal_true.mat')))); % simulated true signals
t = cell2mat(struct2cell(load(fullfile(path,'data','t.mat')))); % time series in decimal year

DATA_T = {DATA_1C,DATA_2C,DATA_3C,DATA_4C,DATA_5C,DATA_6C};
FORM_T = {FORM_1C,FORM_2C,FORM_3C,FORM_4C,FORM_5C,FORM_6C};
SIGNAL_F = zeros(7022,6);
SIGNAL_I = zeros(7022,6);
for i = 1:6
    data = DATA_T{i};
    form = FORM_T{i};
    index_miss = find(isnan(data)==1);
    index_avai = find(isnan(data)==0);
    form = form(index_avai);
    data = data(index_avai);
    signal_f = ESSA_F(form,data,index_miss,L,k );
    signal_i = ESSA_I(data,index_miss,L,k );
    SIGNAL_F(:,i) = signal_f;
    SIGNAL_I(:,i) = signal_i;
end

% Notes:
% (1) plot the estimated signals by ESSA (Identity) and ESSA (Weight) as well as the true signals
% (2) the first sub-figure is the true signals
% (3) the following six sub-figures are the estimated signals under different missing percentages (10%-60%, order: left->right,top->bottom) 
figure
subplot(4,2,1)
plot(t,signal_true,'k')
legend('True signals')
set(gca,'xtick',[]);
for i = 1:6
    subplot(4,2,i+2)
    plot(t,SIGNAL_I(:,i),'g')
    hold on;
    plot(t,SIGNAL_F(:,i),'r')
    if i==1
        legend('ESSA(Identity)','ESSA(Weight)');
    end
end

% plot the difference between the estimated signals and true signals
figure
for i = 1:6
    subplot(3,2,i)
    plot(t,signal_true - SIGNAL_I(:,i),'g')
    hold on;
    plot(t,signal_true - SIGNAL_F(:,i),'r')
end






