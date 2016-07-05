%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    DoNECR.m                                                         %
%                                                                     %
%    ���ļ�����PETϵͳNECR���ԡ�                                      %
%                                                                     %
%    �ڽ���ϵͳɢ�����Ӳ���ʱ��Ҫ����ϵͳ��������������               %
%     1) ��������ʺ����¼������ʵı�ֵ<1.0%                          %
%     2) ÿ�βɼ��ķ��ϼ���>500,000                                   %
%                                                                     %
%    ʵ�ʲ����У����������������Ĳɼ�������Ϊ����ɢ�����ӵ�           %
%    ��Ч���ݣ�                                                       %
%     1) ҩ��Ũ���� 1/4 Aref �� Aref ֮�䣬 ArefΪϵͳ����            %
%        ��Ч��ȣ�                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SetSysEnv();%%%%%%%

FPDIR = 'F:\projects\Digital-PET\ʵ���¼\2015-11-13_ϵͳNECR����\';
isFigure     = 1;
FPEXT        = '.bin';
SinoGramFpExt= '.sin';
BkgFPExt     = '.bkg';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     ����������Ҫ����ÿ�����ݽ�����������        %%%%%
Acal    = 279 .* 3.7e4 ./ 1000;   %%% ����Դ��ʼ���(kcps)
Tcal    = '2016-04-11 10:32';     %%% Acal����ʱ��10:45 
BkgAref = 2875.0557;              %%% ϵͳ������Ч��ȣ�kcps��
%%%%%% --END--  ����������Ҫ����ÿ�����ݽ�����������   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     ����Ϊ�豸̽�����ṹ��ز���                %%%%%
SliceCnt = 11;  AlphaCnt = 529;  RingCnt  = 529;
SinoGramDataStruct = 'float';
SinoGramCenterIndex= 265;
SinoGram04CMIndex  = [128, 138];   %%% 02cm���������  
                                   %%% FOVֱ��Ϊ740mm
SinoGram12CMIndex  = [209, 321];   %%% 12cm���������  
                                   %%% FOVֱ��Ϊ740mm
%%%%%% --END--  ����Ϊ�豸̽�����ṹ��ز���           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%     ʵ�����ݴ������̼���    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  1. �Ա������ݽ��д���                           %%%%%
fps   = dir(fullfile([FPDIR, '*', BkgFPExt]));
if isempty(fps)
    disp('There is no Background Listmode file!!!');
else
    %%% 1) ������Listmode����ת����Sinogram����
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
    
    %%% 2) ����Background��Sinogram����
    fpd = fopen([FPDIR, strrep(fps(1).name, ...
                 BkgFPExt, SinoGramFpExt)], 'r');
    dat = fread(fpd, SinoGramDataStruct);
    fclose(fpd);
    
    BkgCounts = zeros(1, SliceCnt);
    SinoGrams = reshape(dat, SliceCnt, AlphaCnt, RingCnt); clear dat;
    for slice = 1 : SliceCnt
                
        %%% a)  ��ȡ����Sinogram����
        SinoGram = reshape(SinoGrams(slice, :, :),  AlphaCnt, RingCnt);
        if 1 == isFigure
            figure; imshow(SinoGram, [], 'InitialMagnification','fit');  
        end
        
        %%% b) ��SinoGram����12cm�����������Ϊ0
        SinoGram(:, 1:SinoGram12CMIndex(1))   = 0;
        SinoGram(:, SinoGram12CMIndex(2):end) = 0;
        if 1 == isFigure
            figure; imshow(SinoGram, [], 'InitialMagnification','fit');   
        end
        
        BkgCounts(1, slice) = sum(SinoGram(:));
    end  %%% for slice =
    
    %%% 4) ��ȡBackground�Ĳɼ���ʼʱ��T0�Ͳɼ�ʱ��Tacq
    date    = sscanf(strrep(fps(1).name, '-', ''), '%s');
    Start   = datetime(fps(1).date);
    Stop    = datetime(date(1:14), 'InputFormat', 'yyyyMMddHHmmss');
    BkgT0   = fix(minutes(abs(Start - datetime(Tcal))));
    BkgTacq = fix(minutes(abs(Stop - Start)));
    
    %%% 4�����������ʣ�kcps��
    Rbkg    =  BkgCounts ./  (ones(1, SliceCnt) * BkgTacq) ./ 60 ./ 1000;
    
    %%% 5) ɾ��Backgraound��Sinogram����
    delete([FPDIR, strrep(fps(1).name, BkgFPExt, SinoGramFpExt)]);
end %%% if--


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  2. ������д�����Listmode�����ļ��б�           %%%%%
fps   = dir(fullfile([FPDIR, '*', FPEXT]));
if isempty(fps)
    disp('There is no valid files!!!');
    finish;
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  3. ������ʵ���Listmode����ת����Sinogram����   %%%%%                              
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
%%%%%    �ռ�����ʵ�����ݼ��Ĳɼ�ʱ��                  %%%%%
%%%%%  �����ļ������� YYYYMMDD-HHMMSS-*.bin��ʽ����    %%%%%
T0   = zeros(length(fps)-1, 1);  %%% ��¼�Ǳ������ݵĲɼ�    
Tacq = zeros(length(fps)-1, 1);  %%% ��ʼʱ��T0�Ͳɼ�ʱ��Tacq 
for fpi = 1 : length(fps)
    date  = sscanf(strrep(fps(fpi).name, '-', ''), '%s');
    Start = datetime(fps(fpi).date);
    Stop  = datetime(date(1:14), 'InputFormat', 'yyyyMMddHHmmss');
    
    T0(fpi)   = fix(minutes(abs(Start - datetime(Tcal))));
    Tacq(fpi) = fix(minutes(abs(Stop - Start)));
end
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%      ��ñ���ʵ�����ݼ��ķ��ϼ���SingleCounts    %%%%%
%%%%%         ���ϼ�����SingleRate                     %%%%%
%%%%%         ����Դ���A0 �� Aave                     %%%%%
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
%%%%% --END-- ��ñ���ʵ�����ݼ��ķ��ϼ��������ϼ����� %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           ������������Sinogram���ݽ��д���       %%%%%
fps   = dir(fullfile([FPDIR, '*', SinoGramFpExt]));
if isempty(fps)
    disp('There is no valid Sinogram files!!!');
else
    %%%%  I) �ȶ�ÿ�βɼ���Sinogram���ݰ���NEMA��׼���д���
    %%%%     ���ֱ�����������
    Cij_RandScat    = zeros(length(fps), SliceCnt);  %%% ���&ɢ���¼�����
    Cij_Total       = zeros(length(fps), SliceCnt);  %%% ���¼�����
    
    SingleCounts    = zeros(length(fps), 2);         %%% ԭʼ�ܼ���
                                                     %%% ������ܼ���
    for fpi = 1 : length(fps)
        fpd = fopen([FPDIR, fps(fpi).name], 'r');
        dat = fread(fpd, SinoGramDataStruct);
        fclose(fpd);
        
        SinoGrams = reshape(dat, SliceCnt, AlphaCnt, RingCnt);  clear dat;
        
        SingleCounts(fpi, 1) = sum(SinoGrams(:));  %%% For debugging
        
        for slice = 1 : SliceCnt
            
            %%% 1)  ��ȡ����Sinogram����
            SinoGram = reshape(SinoGrams(slice, :, :),  AlphaCnt, RingCnt);
            if 1 == isFigure
                figure; imshow(SinoGram, [], 'InitialMagnification','fit');  
            end
            
            %%%  2) ��SinoGram����12cm�����������Ϊ0
            SinoGram(:, 1:SinoGram12CMIndex(1))   = 0;
            SinoGram(:, SinoGram12CMIndex(2):end) = 0;
            if 1 == isFigure
                figure; imshow(SinoGram, [], 'InitialMagnification','fit');   
            end
            SingleCounts(fpi, 2) = SingleCounts(fpi, 2) + sum(SinoGram(:));
            
            %%%  3) ����Դƽ����SinoGram����λ��
            [M, I] = max(SinoGram, [], 2);  %%%  ÿ�����ֵ��λ��
            for row = 1 : AlphaCnt
                SinoGram(row, :) = circshift(SinoGram(row, :), ...
                             [0, SinoGramCenterIndex - I(row)]);
            end
            if 1 == isFigure
                figure; imshow(SinoGram, [], 'InitialMagnification','fit');    
            end
            
            %%%  4) ��SinoGram�и��Ƕ��µ�ͶӰֵ���
            %%%     ����Alpha�� n �� (n+1) �е�ÿ��ͶӰֵ�ļ��
            %%%     ���Ϊ 1������Ҫ������Щ���
            SinoGram = sum(SinoGram, 1);   %%% ��ÿ�����
            if (mod(SinoGramCenterIndex, 2) == 0)
                SinoGram = SinoGram(2:2:end);
            else
                SinoGram = SinoGram(1:2:end);
            end
            if 1 == isFigure
                figure; plot(SinoGram, 'b.');      
            end
            
            %%%  5) ����ɢ��&������������¼�����
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
% % %     %%%%%%   ����ΪNEMA 2007��׼
% % %     %%%% II) ʹ��FDG����� 1/4 Aref �� Aref ������ݼ�����ɢ������
% % %     SFi = sum(Cij_RandScat(Aave > 0.25*BkgAref & Aave < BkgAref, :), 1) ./ ...
% % %           sum(Cij_Total(   Aave > 0.25*BkgAref & Aave < BkgAref, :), 1);
% % %     SF  = sum(sum(Cij_RandScat(Aave > 0.25*BkgAref & Aave < BkgAref, :))) ./ ...
% % %           sum(sum(Cij_Total(   Aave > 0.25*BkgAref & Aave < BkgAref, :)));
% % %     
% % %     %%%% III) ���ݼ����ɢ�����Ӽ������¸�������(cps)
% % %     Rij_Total       = Cij_Total ./ ...
% % %                       (Tacq * ones(1, SliceCnt));    %%% ���¼�������
% % %     Rij_True        = (Cij_Total - Cij_RandScat) ./ ...
% % %                       (Tacq * ones(1, SliceCnt));    %%% ���¼�������
% % %     Rij_Random      = Rij_Total - Rij_True ./ ...
% % %                       (ones(length(fps), 1) * (ones(1, SliceCnt) - SFi)); 
% % %                                                      %%% ����¼�������
% % %     Rij_Scatter     = ones(length(fps), 1) * ...
% % %                       (SFi ./ (ones(1, SliceCnt) - SFi)) .* Rij_True;  
% % %                                                      %%% ɢ���¼�������
% % %     Rij_NECR        = Rij_True .* Rij_True ./ ...
% % %                       (Rij_Total + Rij_Random);      %%% ������Ч������
% % %     %%%%%% --END--   ����ΪNEMA 2007��׼
% % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  ��NEMA��׼����У����ķ���
    %%%%    Rscatter = Rtotal - Rtrue - Rbackground
    %%%%    SF = Rscatter/(Rscatter + Rtrue)
    %%%%    NECR = Rtrue * Rtrue / Rtotal
    %%%%
    %%%% II) ʹ��FDG����� 1/4 Aref �� Aref ������ݼ�����ɢ������
    AcqCnt   = nnz(Aave > 0.25*BkgAref & Aave < BkgAref);
    
    R_tot  = Cij_Total(Aave > 0.25*BkgAref & Aave < BkgAref, :) ./ ...
            (Tacq(     Aave > 0.25*BkgAref & Aave < BkgAref) * ...
             ones(1, SliceCnt)) ./ 60 ./ 1000;    %%% ���¼������ʣ�kcps��
    R_true = (Cij_Total(   Aave > 0.25*BkgAref & Aave < BkgAref, :) - ...
              Cij_RandScat(Aave > 0.25*BkgAref & Aave < BkgAref, :)) ./ ...
             (Tacq(Aave > 0.25*BkgAref & Aave < BkgAref) * ...
              ones(1, SliceCnt)) ./ 60 ./ 1000;    %%% ���¼������ʣ�kcps��
                
    R_scat = R_tot - R_true - ones(AcqCnt, 1) * Rbkg; 
    
    SFi    = sum(R_scat, 1) ./ sum(R_true + R_scat, 1);
    SF     = sum(R_scat(:)) ./ sum(sum(R_true + R_scat));

    %%%% III) ���ݼ����ɢ�����Ӽ������¸�������(cps)
    Rij_Total       = Cij_Total ./ (Tacq * ones(1, SliceCnt));
    Rij_True        = (Cij_Total - Cij_RandScat) ./  ...
                      (Tacq * ones(1, SliceCnt));    %%% ������¼�������
    Rij_Random      = Rij_Total - Rij_True ./ ...
                      (ones(length(fps), 1) * (ones(1, SliceCnt) - SFi)); 
                                                     %%% ����¼�������
    Rij_Scatter     = ones(length(fps), 1) * ...
                      (SFi ./ (ones(1, SliceCnt) - SFi)) .* Rij_True;  
                                                     %%% ɢ���¼�������
    Rij_NECR        = Rij_True .* Rij_True ./ ...
                      (Rij_Total + Rij_Random);      %%% ������Ч������
    %%%% --END--  ��NEMA��׼����У����ķ���
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Rj_Total    = sum(Rij_Total,   2);   %%% cps
    Rj_True     = sum(Rij_True,    2);   %%% cps
    Rj_Random   = sum(Rij_Random,  2);   %%% cps
    Rj_Scatter  = sum(Rij_Scatter, 2);   %%% cps
    Rj_NECR     = sum(Rij_NECR,    2);   %%% cps
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ����ͼ��                                        %%%%%
figure;   line1 = plot(Aave, Rj_Total,    'bs-'); 
hold on;  line2 = plot(Aave, Rj_True,     'bo-');
hold on;  line3 =plot(Aave, Rj_Random,   'bd-');
hold on;  line4 =plot(Aave, Rj_Scatter,  'bx-');
hold on;  line5 =plot(Aave, Rj_NECR,     'b^-');
legend([line1, line2, line3, line4, line5], ...
       'TOTAL', 'TRUE', 'RANDOM', 'SCATTER', 'NECR');
hold on;  xlabel('Acivity (kcps)'); ylabel('Counts Rate (cps)');
