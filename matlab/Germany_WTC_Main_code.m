clear; clc; close all;
inputFolder = 'D:\ce23resch11018\germany_timeseries\';
baseOut     = 'D:\ce23resch11018\germany\';

plotOut   = fullfile(baseOut,'WTC_PLOTS');
cohOut    = fullfile(baseOut,'WTC_COHERENCE_VALUES');
metricOut = fullfile(baseOut,'RMI_LAG');
featOut   = fullfile(baseOut,'features');
phaseOut  = fullfile(baseOut,'phase_deg_all');
lagOut    = fullfile(baseOut,'lag_days_all');

dirs = {plotOut,cohOut,metricOut,featOut,phaseOut,lagOut};
for i = 1:numel(dirs)
    if ~exist(dirs{i},'dir'), mkdir(dirs{i}); end
end

csvFiles = dir(fullfile(inputFolder,'*.csv'));
bands = [
    2     7;
    7     14;
    14    31;
    31    183;
    183   365;
    365   730;
    730   1460;
    1460  10000
];
numBands = size(bands,1);

HF_band = [2 31];
LF_band = [183 10000];
featureNames = {};
for b = 1:numBands
    featureNames = [featureNames, ...
        sprintf("MeanWCOH_Band%d",b), ...
        sprintf("PercSig_Band%d",b), ...
        sprintf("MeanPhaseRad_Band%d",b), ...
        sprintf("MeanPhaseDeg_Band%d",b), ...
        sprintf("PhaseConc_Band%d",b)];
end
featureNames = [featureNames,"DominantPeriod"];

for k = 1:length(csvFiles)

    try
        filePath = fullfile(inputFolder,csvFiles(k).name);
        station  = erase(csvFiles(k).name,'.csv');

        fprintf('\nProcessing catchment: %s\n',station);
        data = readtable(filePath);

        X = data.discharge_mm(:);
        Y = data.precipitation_mm(:);
        dates = datetime(data.date,'InputFormat','dd-MM-yyyy');

        valid = ~isnan(X) & ~isnan(Y);
        X = X(valid); Y = Y(valid); dates = dates(valid);

        if length(X) < 200 || std(X)==0 || std(Y)==0
            fprintf('Skipped (insufficient data)\n');
            continue;
        end

        X = (X - mean(X,'omitnan')) ./ std(X,'omitnan');
        Y = (Y - mean(Y,'omitnan')) ./ std(Y,'omitnan');

        [coh, period, ~, ~, sig] = wtc(X,Y);
        coh(coh>1) = NaN;
        sig(sig>1) = 1;

        cohTbl = array2table(coh','VariableNames', ...
            cellstr(num2str(period','P_%0.2f')));
        cohTbl.Date = dates;
        cohTbl = movevars(cohTbl,'Date','before',1);
        writetable(cohTbl,fullfile(cohOut,station+"_coherence.csv"));

        figure('Color','w','Position',[100 100 1200 800]);
        wtc(X,Y);
        ax = gca;

        years = year(dates);
        uniqYears = unique(years);
        nYears = numel(uniqYears);

        tickPos = zeros(nYears,1);
        for y = 1:nYears
            tickPos(y) = find(years==uniqYears(y),1,'first');
        end

        ax.XTick = tickPos;
        ax.XTickLabel = string(1:nYears);
        ax.XTickLabelRotation = 0;

        xlabel('Year of record');
        ylabel('Period (days)');
        title(['Wavelet Coherence: ',station],'Interpreter','none');
        set(gca,'FontSize',12)

        saveas(gcf,fullfile(plotOut,station+"_WTC.jpg"));
        close
        WXY = xwt(X,Y);
        phase = angle(WXY);
        phaseDeg = phase * 180/pi;

        numScales = length(period);
        lagDays = NaN(size(phase));

        for s = 1:numScales
            lagDays(s,:) = (phase(s,:) / (2*pi)) * period(s);
        end

        phaseTbl = array2table(phaseDeg);
        phaseTbl = addvars(phaseTbl,period(:), ...
            'Before',1,'NewVariableNames','PeriodDays');
        writetable(phaseTbl,fullfile(phaseOut,station+"_phase_deg.csv"));
        lagTbl = array2table(lagDays);
        lagTbl = addvars(lagTbl,period(:), ...
            'Before',1,'NewVariableNames','PeriodDays');
        writetable(lagTbl,fullfile(lagOut,station+"_lag_days.csv"));

        meanPhase = NaN(numScales,1);
        phaseConc = NaN(numScales,1);

        for s = 1:numScales
            ph = phase(s,:);
            ph = ph(~isnan(ph));
            if ~isempty(ph)
                c = exp(1i*ph);
                meanPhase(s) = angle(mean(c));
                phaseConc(s) = abs(mean(c));
            end
        end
        featureRow = [];

        for b = 1:numBands
            idx = period>=bands(b,1) & period<bands(b,2);

            meanCoh = mean(coh(idx,:),'all','omitnan');
            percSig = 100 * sum(sig(idx,:),'all')/numel(sig(idx,:));

            bandPhRad = mean(meanPhase(idx),'omitnan');
            bandPhDeg = bandPhRad*180/pi;
            bandConc  = mean(phaseConc(idx),'omitnan');

            featureRow = [featureRow ...
                meanCoh percSig bandPhRad bandPhDeg bandConc];
        end

        [~,domIdx] = max(mean(coh,2,'omitnan'));
        dominantPeriod = period(domIdx);

        featureRow = [featureRow dominantPeriod];

        FT = array2table(featureRow,'VariableNames',featureNames);
        FT.Catchment = string(station);
        FT = movevars(FT,'Catchment','before',1);
        writetable(FT,fullfile(featOut,station+"_features.csv"));
        HF_idx = period>=HF_band(1) & period<HF_band(2);
        LF_idx = period>=LF_band(1) & period<LF_band(2);

        HF_mean = mean(coh(HF_idx,:),'all','omitnan');
        LF_mean = mean(coh(LF_idx,:),'all','omitnan');
        RMI     = HF_mean / LF_mean;

        meanLag = mean(lagDays(:),'omitnan');
        meanPhR = mean(phase(:),'omitnan');
        meanPhD = meanPhR*180/pi;

        M = table(string(station),HF_mean,LF_mean,RMI, ...
                  meanLag,meanPhR,meanPhD, ...
            'VariableNames',{'Catchment','HF_Mean','LF_Mean','RMI', ...
                             'MeanLag_days','MeanPhase_rad','MeanPhase_deg'});
        writetable(M,fullfile(metricOut,station+"_RMI_LAG.csv"));

        fprintf(' Completed %s\n',station);

    catch ME
        fprintf(' ERROR %s : %s\n',station,ME.message);
    end
end

fprintf('\nALL GERMANY CATCHMENTS PROCESSED SUCCESSFULLY\n');


