%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    co_coincis.m                                                             %
%                                                                             %
%    该文件用于对采集到的数据进行 时间符合 运算，获得时间窗内的符合数据。     %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    可配置参数                                                               %
%    FP_DIR   --  事件数据保存路径                                            %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FP_DIR='D:\projects\Data\2015-6-15-19-42\';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    将所有的符合数据汇总保存成单个文件                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_list = dir(fullfile([FP_DIR, '*.bcoi']));
data    = load( [FP_DIR, fp_list(1).name], '-mat');
coincis = data.data;
clear data;

for fpi = 2 : length( fp_list )
    data    = load( [FP_DIR, fp_list(fpi).name], '-mat');
    coincis = [coincis; data.data];
    
    clear data;
end

sect    = length(coincis(1, :))/2;
coincis = [coincis(:, 1:sect), ...
           zeros(length(coincis(:, 1)), 1), ...
           coincis(:, sect+1:end), ...
           zeros(length(coincis(:, 1)), 1)];

fp = fopen([FP_DIR, 'coincidencs.cocis'], 'w');
fwrite(fp, coincis', 'double', 'ieee-be');
fclose(fp);

for fpi = 1 : length(fp_list)
    delete([FP_DIR, fp_list(fpi).name]); 
end  %%% for fpi =

clear all; exit;
