%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    DoNECR.m                                                         %
%                                                                     %
%    该文件用于PET系统NECR测试。                                      %
%                                                                     %
%    在进行系统散射因子测试时，要满足系统如下两个条件：               %
%     1) 随机符合率和真事件符合率的比值<1.0%                          %
%     2) 每次采集的符合计数>500,000                                   %
%                                                                     %
%    实际测量中，将满足如下条件的采集数据作为计算散射因子的           %
%    有效数据：                                                       %
%     1) 药物浓度在 1/4 Aref 和 Aref 之间， Aref为系统本底            %
%        等效活度；                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SetSysEnv();%%%%%%%

FPDIR = 'F:\projects\Digital-PET\实验记录\2015-11-13_系统NECR测试\';
isFigure     = 1;
FPEXT        = '.bin';
SinoGramFpExt= '.sin';
BkgFPExt     = '.bkg';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     以下设置需要根据每次数据进行重新配置        %%%%%
Acal    = 279 .* 3.7e4 ./ 1000;   %%% 放射源初始活度(kcps)
Tcal    = '2016-04-11 10:32';     %%% Acal测试时间10:45 
BkgAref = 2875.0557;              %%% 系统背景等效活度（kcps）
%%%%%% --END--  以下设置需要根据每次数据进行重新配置   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     以下为设备探测器结构相关参数                %%%%%
SliceCnt = 11;  AlphaCnt = 529;  RingCnt  = 529;
SinoGramDataStruct = 'float';
SinoGramCenterIndex= 265;
SinoGram04CMIndex  = [128, 138];   %%% 02cm以外的区域  
                                   %%% FOV直径为740mm
SinoGram12CMIndex  = [209, 321];   %%% 12cm以外的区域  
                                   %%% FOV直径为740mm
%%%%%% --END--  以下为设备探测器结构相关参数           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%     实际数据处理流程见下    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  1. 对背景数据进行处理                           %%%%%
fps   = dir(fullfile([FPDIR, '*', BkgFPExt]));
if isempty(fps)
    disp('There is no Background Listmode file!!!');
else
    %%% 1) 将背景Listmode数据转化成Sinogram数据
    Singles   = GetSinglesFromFile([FPDIR, fps(1).name], 20);
    Singles(:, 1) = Singles(:, 1) - 256; 
    CoSingles = reshape(Singles', ...
                        length(Singles(1, :))*2, ...
                        length(Singles(:, 1))/2)';
    clear Singles;

    SinoGram = SinoGramFn(CoSingles);
     
    fpd = fopen([FPDIR, strrep(fps(1).name, ...
                 BkgFPExt, SinoGramFpExt)], 'w');
    fwrite(fpd, SinoGram, SinoGramDataStruct);
    fclose(fpd);
    
    %%% 2) 处理Background的Sinogram数据
    fpd = fopen([FPDIR, strrep(fps(1).name, ...
                 BkgFPExt, SinoGramFpExt)], 'r');
    dat = fread(fpd, SinoGramDataStruct);
    fclose(fpd);
    
    BkgCounts = zeros(1, SliceCnt);
    SinoGrams = reshape(dat, SliceCnt, AlphaCnt, RingCnt); clear dat;
    for slice = 1 : SliceCnt
                
        %%% a)  获取单层Sinogram数据
        SinoGram = reshape(SinoGrams(slice, :, :),  AlphaCnt, RingCnt);
        if 1 == isFigure
            figure; imshow(SinoGram, [], 'InitialMagnification','fit');  
        end
        
        %%% b) 将SinoGram中心12cm以外的区域置为0
        SinoGram(:, 1:SinoGram12CMIndex(1))   = 0;
        SinoGram(:, SinoGram12CMIndex(2):end) = 0;
        if 1 == isFigure
            figure; imshow(SinoGram, [], 'InitialMagnification','fit');   
        end
        
        BkgCounts(1, slice) = sum(SinoGram(:));
    end  %%% for slice =
    
    %%% 4) 获取Background的采集开始时间T0和采集时间Tacq
    date    = sscanf(strrep(fps(1).name, '-', ''), '%s');
    Start   = datetime(fps(1).date);
    Stop    = datetime(date(1:14), 'InputFormat', 'yyyyMMddHHmmss');
    BkgT0   = fix(minutes(abs(Start - datetime(Tcal))));
    BkgTacq = fix(minutes(abs(Stop - Start)));
    
    %%% 4）背景计数率（kcps）
    Rbkg    =  BkgCounts ./  (ones(1, SliceCnt) * BkgTacq) ./ 60 ./ 1000;
    
    %%% 5) 删除Backgraound的Sinogram数据
    delete([FPDIR, strrep(fps(1).name, BkgFPExt, SinoGramFpExt)]);
end %%% if--


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  2. 获得所有待处理Listmode数据文件列表           %%%%%
fps   = dir(fullfile([FPDIR, '*', FPEXT]));
if isempty(fps)
    disp('There is no valid files!!!');
    finish;
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  3. 将本次实验的Listmode数据转换成Sinogram数据   %%%%%                              
disp(['There are ', num2str(length(fps)), ...
      ' files to be processed.']);
for fpi = 1 : length(fps)
    disp(['No.', num2str(fpi), ' file is processing...']);
    
    Singles   = GetSinglesFromFile([FPDIR, fps(fpi).name], 20);
    Singles(:, 1) = Singles(:, 1) - 256; 
    CoSingles = reshape(Singles', ...
                        length(Singles(1, :))*2, ...
                        length(Singles(:, 1))/2)';
    clear Singles;

    SinoGrams = SinoGramFn(CoSingles);
     
    fpd = fopen([FPDIR, strrep(fps(fpi).name, ...
                 FPEXT, '.sin')], 'w');
    fwrite(fpd, SinoGrams, 'float');
    fclose(fpd);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    收集本次实验数据集的采集时间                  %%%%%
%%%%%  所有文件名按照 YYYYMMDD-HHMMSS-*.bin格式命名    %%%%%
T0   = zeros(length(fps)-1, 1);  %%% 记录非背景数据的采集    
Tacq = zeros(length(fps)-1, 1);  %%% 开始时间T0和采集时长Tacq 
for fpi = 1 : length(fps)
    date  = sscanf(strrep(fps(fpi).name, '-', ''), '%s');
    Start = datetime(fps(fpi).date);
    Stop  = datetime(date(1:14), 'InputFormat', 'yyyyMMddHHmmss');
    
    T0(fpi)   = fix(minutes(abs(Start - datetime(Tcal))));
    Tacq(fpi) = fix(minutes(abs(Stop - Start)));
end
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%      获得本次实验数据集的符合计数SingleCounts    %%%%%
%%%%%         符合计数率SingleRate                     %%%%%
%%%%%         放射源活度A0 和 Aave                     %%%%%
fps   = dir(fullfile([FPDIR, '*', FPEXT]));
if isempty(fps)
    disp('There is no valid Listmode files!!!');
else
    SingleCounts = zeros(length(fps), 1);      
    for fpi = 1 : length(fps)
        SingleCounts(fpi) = fps(fpi).bytes / 20 / 2;
    end
    
    SingleRate  = SingleCounts ./ Tacq ./ 1000;  %%% kcps
    [A0, Aave]  = GetActivityFn(Acal, 0, T0, Tacq);
end
%%%%% --END-- 获得本次实验数据集的符合计数及符合计数率 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           对满足条件的Sinogram数据进行处理       %%%%%
fps   = dir(fullfile([FPDIR, '*', SinoGramFpExt]));
if isempty(fps)
    disp('There is no valid Sinogram files!!!');
else
    %%%%  I) 先对每次采集的Sinogram数据按照NEMA标准进行处理，
    %%%%     并分别获得以下数据
    Cij_RandScat    = zeros(length(fps), SliceCnt);  %%% 随机&散射事件计数
    Cij_Total       = zeros(length(fps), SliceCnt);  %%% 总事件计数
    
    SingleCounts    = zeros(length(fps), 2);         %%% 原始总计数
                                                     %%% 处理后总计数
    for fpi = 1 : length(fps)
        fpd = fopen([FPDIR, fps(fpi).name], 'r');
        dat = fread(fpd, SinoGramDataStruct);
        fclose(fpd);
        
        SinoGrams = reshape(dat, SliceCnt, AlphaCnt, RingCnt);  clear dat;
        
        SingleCounts(fpi, 1) = sum(SinoGrams(:));  %%% For debugging
        
        for slice = 1 : SliceCnt
            
            %%% 1)  获取单层Sinogram数据
            SinoGram = reshape(SinoGrams(slice, :, :),  AlphaCnt, RingCnt);
            if 1 == isFigure
                figure; imshow(SinoGram, [], 'InitialMagnification','fit');  
            end
            
            %%%  2) 将SinoGram中心12cm以外的区域置为0
            SinoGram(:, 1:SinoGram12CMIndex(1))   = 0;
            SinoGram(:, SinoGram12CMIndex(2):end) = 0;
            if 1 == isFigure
                figure; imshow(SinoGram, [], 'InitialMagnification','fit');   
            end
            SingleCounts(fpi, 2) = SingleCounts(fpi, 2) + sum(SinoGram(:));
            
            %%%  3) 将线源平移至SinoGram中心位置
            [M, I] = max(SinoGram, [], 2);  %%%  每行最大值的位置
            for row = 1 : AlphaCnt
                SinoGram(row, :) = circshift(SinoGram(row, :), ...
                             [0, SinoGramCenterIndex - I(row)]);
            end
            if 1 == isFigure
                figure; imshow(SinoGram, [], 'InitialMagnification','fit');    
            end
            
            %%%  4) 将SinoGram中各角度下的投影值相加
            %%%     由于Alpha在 n 和 (n+1) 行的每个投影值的间隔
            %%%     相差为 1，故需要跳过这些间隔
            SinoGram = sum(SinoGram, 1);   %%% 将每列相加
            if (mod(SinoGramCenterIndex, 2) == 0)
                SinoGram = SinoGram(2:2:end);
            else
                SinoGram = SinoGram(1:2:end);
            end
            if 1 == isFigure
                figure; plot(SinoGram, 'b.');      
            end
            
            %%%  5) 计算散射&随机计数、真事件计数
            C_L   = (SinoGram(SinoGram04CMIndex(1)) + ...
                     SinoGram(SinoGram04CMIndex(1) + 1)) / 2;
            C_R   = (SinoGram(SinoGram04CMIndex(2)) + ...
                     SinoGram(SinoGram04CMIndex(2) - 1)) / 2;
                 
            if 1 == isFigure
                hold on; plot(SinoGram04CMIndex(1),   SinoGram(SinoGram04CMIndex(1)), 'r*');
                hold on; plot(SinoGram04CMIndex(1)+1, SinoGram(SinoGram04CMIndex(1)+1), 'r*');  
                hold on; plot(SinoGram04CMIndex(2),   SinoGram(SinoGram04CMIndex(2)), 'r*');
                hold on; plot(SinoGram04CMIndex(2)-1, SinoGram(SinoGram04CMIndex(2)-1), 'r*');  
            end
                     
            Cij_RandScat(fpi, slice)  = (C_L + C_R) ./ 2 .* ...
                (SinoGram04CMIndex(2) - SinoGram04CMIndex(1) - 1) + ...
                sum(SinoGram(1 : SinoGram04CMIndex(1))) + ...
                sum(SinoGram(SinoGram04CMIndex(2) : end));         
            Cij_Total(fpi, slice) = sum(SinoGram);
        end  %%% for slice =
        
%         delete([FPDIR, fps(fpi).name]);
    end  %%% for fpi =
    
% % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %     %%%%%%   以下为NEMA 2007标准
% % %     %%%% II) 使用FDG活度在 1/4 Aref 和 Aref 间的数据集计算散射因子
% % %     SFi = sum(Cij_RandScat(Aave > 0.25*BkgAref & Aave < BkgAref, :), 1) ./ ...
% % %           sum(Cij_Total(   Aave > 0.25*BkgAref & Aave < BkgAref, :), 1);
% % %     SF  = sum(sum(Cij_RandScat(Aave > 0.25*BkgAref & Aave < BkgAref, :))) ./ ...
% % %           sum(sum(Cij_Total(   Aave > 0.25*BkgAref & Aave < BkgAref, :)));
% % %     
% % %     %%%% III) 根据计算的散射因子计算如下各计数率(cps)
% % %     Rij_Total       = Cij_Total ./ ...
% % %                       (Tacq * ones(1, SliceCnt));    %%% 总事件计数率
% % %     Rij_True        = (Cij_Total - Cij_RandScat) ./ ...
% % %                       (Tacq * ones(1, SliceCnt));    %%% 真事件计数率
% % %     Rij_Random      = Rij_Total - Rij_True ./ ...
% % %                       (ones(length(fps), 1) * (ones(1, SliceCnt) - SFi)); 
% % %                                                      %%% 随机事件计数率
% % %     Rij_Scatter     = ones(length(fps), 1) * ...
% % %                       (SFi ./ (ones(1, SliceCnt) - SFi)) .* Rij_True;  
% % %                                                      %%% 散射事件计数率
% % %     Rij_NECR        = Rij_True .* Rij_True ./ ...
% % %                       (Rij_Total + Rij_Random);      %%% 噪声等效计数率
% % %     %%%%%% --END--   以下为NEMA 2007标准
% % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  对NEMA标准进行校正后的方法
    %%%%    Rscatter = Rtotal - Rtrue - Rbackground
    %%%%    SF = Rscatter/(Rscatter + Rtrue)
    %%%%    NECR = Rtrue * Rtrue / Rtotal
    %%%%
    %%%% II) 使用FDG活度在 1/4 Aref 和 Aref 间的数据集计算散射因子
    AcqCnt   = nnz(Aave > 0.25*BkgAref & Aave < BkgAref);
    
    R_tot  = Cij_Total(Aave > 0.25*BkgAref & Aave < BkgAref, :) ./ ...
            (Tacq(     Aave > 0.25*BkgAref & Aave < BkgAref) * ...
             ones(1, SliceCnt)) ./ 60 ./ 1000;    %%% 总事件计数率（kcps）
    R_true = (Cij_Total(   Aave > 0.25*BkgAref & Aave < BkgAref, :) - ...
              Cij_RandScat(Aave > 0.25*BkgAref & Aave < BkgAref, :)) ./ ...
             (Tacq(Aave > 0.25*BkgAref & Aave < BkgAref) * ...
              ones(1, SliceCnt)) ./ 60 ./ 1000;    %%% 真事件计数率（kcps）
                
    R_scat = R_tot - R_true - ones(AcqCnt, 1) * Rbkg; 
    
    SFi    = sum(R_scat, 1) ./ sum(R_true + R_scat, 1);
    SF     = sum(R_scat(:)) ./ sum(sum(R_true + R_scat));

    %%%% III) 根据计算的散射因子计算如下各计数率(cps)
    Rij_Total       = Cij_Total ./ (Tacq * ones(1, SliceCnt));
    Rij_True        = (Cij_Total - Cij_RandScat) ./  ...
                      (Tacq * ones(1, SliceCnt));    %%% 真符合事件计数率
    Rij_Random      = Rij_Total - Rij_True ./ ...
                      (ones(length(fps), 1) * (ones(1, SliceCnt) - SFi)); 
                                                     %%% 随机事件计数率
    Rij_Scatter     = ones(length(fps), 1) * ...
                      (SFi ./ (ones(1, SliceCnt) - SFi)) .* Rij_True;  
                                                     %%% 散射事件计数率
    Rij_NECR        = Rij_True .* Rij_True ./ ...
                      (Rij_Total + Rij_Random);      %%% 噪声等效计数率
    %%%% --END--  对NEMA标准进行校正后的方法
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Rj_Total    = sum(Rij_Total,   2);   %%% cps
    Rj_True     = sum(Rij_True,    2);   %%% cps
    Rj_Random   = sum(Rij_Random,  2);   %%% cps
    Rj_Scatter  = sum(Rij_Scatter, 2);   %%% cps
    Rj_NECR     = sum(Rij_NECR,    2);   %%% cps
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  绘制图形                                        %%%%%
figure;   line1 = plot(Aave, Rj_Total,    'bs-'); 
hold on;  line2 = plot(Aave, Rj_True,     'bo-');
hold on;  line3 =plot(Aave, Rj_Random,   'bd-');
hold on;  line4 =plot(Aave, Rj_Scatter,  'bx-');
hold on;  line5 =plot(Aave, Rj_NECR,     'b^-');
legend([line1, line2, line3, line4, line5], ...
       'TOTAL', 'TRUE', 'RANDOM', 'SCATTER', 'NECR');
hold on;  xlabel('Acivity (kcps)'); ylabel('Counts Rate (cps)');
