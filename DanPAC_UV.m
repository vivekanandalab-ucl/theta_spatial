%% Script to assess phase-amplitude coupling in the intracranial data from
%  the spatial memory tas. Daniel Bush, UCL (2020) drdanielbush@gmail.com

%  Provide some settings
folder          = '/Users/umie/Documents/Matlab/New Cue Period Data/';
folder2         = '/Users/umie/Documents/Matlab/old_cue/';
patientNos      = [2 5 6 8 9 10 13 14 15 16];
region          = 'RealHPC';
keepTrials      = 'no';         % Average MI values across individual trials ('yes'), or concatenate trials before the MI calculation ('no')?
startTime       = -5;           % Start time of the data epochs
baseWin         = [-3 0];       % Baseline time window (s)
cueWin          = [0 3];        % Cue time window (s)
phaseBand       = [2 12];       % Wide phase frequency band (Hz)
lowPhase        = [2 5];
highPhase       = [6 11];
ampBand         = [20 120];     % Wide amplitude frequency band (Hz)
lowAmp          = [35 50];
highAmp         = [60 110];
nFreqPoints     = 50;           % Number of frequency bands
freqRange       = [min(phaseBand) max(ampBand)];
freqBands       = linspace(log(freqRange(1)),log(freqRange(2)),nFreqPoints); clear nFreqPoints freqRange
freqBands       = exp(freqBands);

%  Cycle through each patient
for p           = 1 : length(patientNos)
    
    % Load the raw data and extract trial information
    filename     = ['Patient ' int2str(patientNos(p)) '_Cue_' region '_-2to4s_UmeshERPData.mat'];
    load([folder2 filename],'Data'); clear filename
    errors       = Data.Error; clear Data
    filename     = ['Patient ' int2str(patientNos(p)) '_Cue_' region '_-5to5s_UmeshERPData.mat'];
    load([folder filename]); clear filename
    Fs           = (unique(Data.DataPoints)-1) / range(Data.Timebase);
    trials   	 = [Data.Time(:,1:2)*Fs startTime*Fs*ones(size(Data.Time,1),1)];
    channels     = Data.Contacts;
    firstChan  	= true;
    
    % Loop through each channel...
    for c    	= 1 : length(channels)
        
        % ...and, if there are enough artefact free trials on that
        % channel...        
        if sum(Data.Artifacts(:,c)==0)>1
            
            % Pre-process (artefact free trials only)
            cfg.dataset     = [folder Data.Patient];
            cfg.channel     = channels(c);
            cfg.trl         = trials(Data.Artifacts(:,c)==0,:);
            cfg             = ft_definetrial(cfg);
            epoch_data      = ft_preprocessing(cfg); clear cfg
            clc
            
            % Downsample
            cfg.resamplefs  = Settings.downFs;
            cfg.detrend     = 'yes';
            cfg.demean      = 'yes';
            epoch_data      = ft_resampledata(cfg,epoch_data); clear cfg
            clc
            
            % Generate the time-frequency representation
            cfg.method      = 'wavelet';
            cfg.width       = 5;
            cfg.trials      = 'all';
            cfg.output      = 'fourier';
            cfg.foi         = freqBands;
            cfg.keeptrials  = 'yes';
            cfg.toi         = epoch_data.time{1}(epoch_data.time{1}>=cueWin(1) & epoch_data.time{1}<=cueWin(2));
            cuePower        = ft_freqanalysis(cfg, epoch_data);
            cfg.toi         = epoch_data.time{1}(epoch_data.time{1}>=baseWin(1) & epoch_data.time{1}<=baseWin(2));
            basePower       = ft_freqanalysis(cfg, epoch_data); clear cfg epoch_data
            clc
            
            % Compute modulation index
            cfg.method      = 'mi';
            cfg.keeptrials  = keepTrials;
            cfg.freqlow     = phaseBand;
            cfg.freqhigh    = ampBand;
            CFC             = ft_crossfrequencyanalysis2(cfg,cuePower,cuePower);
            baseCFC         = ft_crossfrequencyanalysis2(cfg,basePower,basePower); clear basePower
            
            %now z score
            zCFC                = mean(cat(1,CFC.crsspctrm,baseCFC.crsspctrm),1);
            zStdCFC             = std(cat(1,CFC.crsspctrm,baseCFC.crsspctrm),[],1);
            CFC.crsspctrm       = (CFC.crsspctrm - repmat(zCFC,[size(CFC.crsspctrm,1) 1 1 1])) ./ repmat(zStdCFC,[size(CFC.crsspctrm,1) 1 1 1]);
            baseCFC.crsspctrm   = (baseCFC.crsspctrm - repmat(zCFC,[size(baseCFC.crsspctrm,1) 1 1 1])) ./ repmat(zStdCFC,[size(baseCFC.crsspctrm,1) 1 1 1]);
    
            clc                            
            
            % Store the output for that participant
            if firstChan
                if p == 1
                    out.phsFreq     = CFC.freqlow;
                    out.ampFreq     = CFC.freqhigh;
                    out.cuePAC      = nan(length(out.phsFreq),length(out.ampFreq),length(patientNos));
                    out.basePAC     = nan(length(out.phsFreq),length(out.ampFreq),length(patientNos));
                    out.goodPAC     = nan(length(out.phsFreq),length(out.ampFreq),length(patientNos));
                    out.badPAC      = nan(length(out.phsFreq),length(out.ampFreq),length(patientNos));
                end
                cuePAC              = nan(length(out.phsFreq),length(out.ampFreq),length(channels));
                basePAC             = nan(length(out.phsFreq),length(out.ampFreq),length(channels));
                goodPAC             = nan(length(out.phsFreq),length(out.ampFreq),length(channels));
                badPAC              = nan(length(out.phsFreq),length(out.ampFreq),length(channels));
                firstChan           = false;
            end
            switch keepTrials
                case 'yes'
                    cuePAC(:,:,c)	= permute(mean(CFC.crsspctrm),[3 4 1 2]);
                    basePAC(:,:,c)  = permute(mean(baseCFC.crsspctrm),[3 4 1 2]);
                case 'no'
                    cuePAC(:,:,c)	= permute(CFC.crsspctrm,[2 3 1]);
                    basePAC(:,:,c)  = permute(baseCFC.crsspctrm,[2 3 1]); 
            end
            clear cfg CFC baseCFC
            
            % Identify good and bad trials
            errors_trial    = errors(Data.Artifacts(:,c)==0);
            %goodTrials      = errors_trial <= median(errors_trial);
            %badTrials       = errors_trial >  median(errors_trial); clear errors_trial
            goodTrials      = errors_trial <= prctile(errors_trial,40);
            badTrials       = errors_trial > prctile(errors_trial,60); clear errors_trial
            
            
            
            % ...if there are a sufficient number
            if sum(goodTrials)>1 && sum(badTrials)>1
                
                % Extract the power data for good and bad trials separately
                goodPower   = cuePower;
                goodPower.fourierspctrm     = goodPower.fourierspctrm(goodTrials,:,:,:);
                goodPower.cumtapcnt         = goodPower.cumtapcnt(goodTrials,:);
                badPower    = cuePower;
                badPower.fourierspctrm      = badPower.fourierspctrm(badTrials,:,:,:);
                badPower.cumtapcnt          = badPower.cumtapcnt(badTrials,:);
                
                % Compute modulation index for each set of trials
                cfg.method      = 'mi';
                cfg.keeptrials  = keepTrials;
                cfg.freqlow     = phaseBand;
                cfg.freqhigh    = ampBand;
                goodCFC         = ft_crossfrequencyanalysis2(cfg,goodPower,goodPower); clear goodPower
                badCFC          = ft_crossfrequencyanalysis2(cfg,badPower,badPower); clear badPower
                
                %now z score
                goodCFC.crsspctrm  = (goodCFC.crsspctrm - repmat(zCFC,[size(goodCFC.crsspctrm,1) 1 1 1])) ./ repmat(zStdCFC,[size(goodCFC.crsspctrm,1) 1 1 1]);
                badCFC.crsspctrm   = (badCFC.crsspctrm - repmat(zCFC,[size(badCFC.crsspctrm,1) 1 1 1])) ./ repmat(zStdCFC,[size(badCFC.crsspctrm,1) 1 1 1]); clear zCFC zSTDCFC
    
                switch keepTrials
                    case 'yes'
                        goodPAC(:,:,c)	= permute(mean(goodCFC.crsspctrm),[3 4 1 2]);
                        badPAC(:,:,c)   = permute(mean(badCFC.crsspctrm),[3 4 1 2]);
                    case 'no'
                        goodPAC(:,:,c)	= permute(goodCFC.crsspctrm,[2 3 1]);
                        badPAC(:,:,c)   = permute(badCFC.crsspctrm,[2 3 1]);
                end                
                clear goodCFC badCFC cfg
                clc
            end
            clear goodTrials badTrials cuePower                    
                    
        end
        
    end
    clear c firstChan Data Settings channels trials errors Fs
    
    % Store the average data for that participant
    out.cuePAC(:,:,p)       = nanmean(cuePAC,3); clear cuePAC
    out.basePAC(:,:,p)      = nanmean(basePAC,3); clear basePAC
    out.goodPAC(:,:,p)      = nanmean(goodPAC,3); clear goodPAC
    out.badPAC(:,:,p)       = nanmean(badPAC,3); clear badPAC
        
end
clear p ampBand phaseBand freqBands startTime dataFold

%now get mean of low/high theta phase and low/high gamma amplitude
flpInds       = [find(out.phsFreq>=lowPhase(1),1,'first') find(out.phsFreq<=lowPhase(2),1,'last')];
fhpInds       = [find(out.phsFreq>=highPhase(1),1,'first') find(out.phsFreq<=highPhase(2),1,'last')];
flaInds       = [find(out.ampFreq>=lowAmp(1),1,'first') find(out.ampFreq<=lowAmp(2),1,'last')];
fhaInds       = [find(out.ampFreq>=highAmp(1),1,'first') find(out.ampFreq<=highAmp(2),1,'last')];

low_lowcuePAC  = squeeze(mean(mean(out.cuePAC(flpInds(1):flpInds(2),flaInds(1):flaInds(2),:)))); 
low_highcuePAC = squeeze(mean(mean(out.cuePAC(flpInds(1):flpInds(2),fhaInds(1):fhaInds(2),:)))); 
high_lowcuePAC = squeeze(mean(mean(out.cuePAC(fhpInds(1):fhpInds(2),flaInds(1):flaInds(2),:)))); 
high_highcuePAC = squeeze(mean(mean(out.cuePAC(fhpInds(1):fhpInds(2),fhaInds(1):fhaInds(2),:)))); 

low_lowbasePAC  = squeeze(mean(mean(out.basePAC(flpInds(1):flpInds(2),flaInds(1):flaInds(2),:)))); 
low_highbasePAC = squeeze(mean(mean(out.basePAC(flpInds(1):flpInds(2),fhaInds(1):fhaInds(2),:)))); 
high_lowbasePAC = squeeze(mean(mean(out.basePAC(fhpInds(1):fhpInds(2),flaInds(1):flaInds(2),:)))); 
high_highbasePAC = squeeze(mean(mean(out.basePAC(fhpInds(1):fhpInds(2),fhaInds(1):fhaInds(2),:)))); 

low_lowgoodPAC  = squeeze(mean(mean(out.goodPAC(flpInds(1):flpInds(2),flaInds(1):flaInds(2),:)))); 
low_highgoodPAC = squeeze(mean(mean(out.goodPAC(flpInds(1):flpInds(2),fhaInds(1):fhaInds(2),:)))); 
high_lowgoodPAC = squeeze(mean(mean(out.goodPAC(fhpInds(1):fhpInds(2),flaInds(1):flaInds(2),:)))); 
high_highgoodPAC = squeeze(mean(mean(out.goodPAC(fhpInds(1):fhpInds(2),fhaInds(1):fhaInds(2),:)))); 

low_lowbadPAC  = squeeze(mean(mean(out.badPAC(flpInds(1):flpInds(2),flaInds(1):flaInds(2),:)))); 
low_highbadPAC = squeeze(mean(mean(out.badPAC(flpInds(1):flpInds(2),fhaInds(1):fhaInds(2),:)))); 
high_lowbadPAC = squeeze(mean(mean(out.badPAC(fhpInds(1):fhpInds(2),flaInds(1):flaInds(2),:)))); 
high_highbadPAC = squeeze(mean(mean(out.badPAC(fhpInds(1):fhpInds(2),fhaInds(1):fhaInds(2),:)))); 

%save([folder2 '\' region '_stat_mi.mat'],'low_lowcuePAC','low_highcuePAC','high_lowcuePAC','high_highcuePAC','low_lowbasePAC','low_highbasePAC','high_lowbasePAC','high_highbasePAC','low_lowgoodPAC','low_highgoodPAC','high_lowgoodPAC', 'high_highgoodPAC','low_lowbadPAC','low_highbadPAC','high_lowbadPAC', 'high_highbadPAC');
    
%  Plot a summary of the output
figure
datRange        = ceil(max(max(abs(nanmean(out.cuePAC-out.basePAC,3).*10000))))./10000;
datRange        = [-0.3 datRange];
subplot(1,2,2)
contourf(out.phsFreq,out.ampFreq,nanmean(out.cuePAC-out.basePAC,3)',linspace(datRange(1),datRange(2),50),'linestyle','none')
colormap jet, axis square, set(gca,'FontSize',24), xlabel('Phase Frequency (Hz)','FontSize',36), ylabel('Amplitude Frequency (Hz)','FontSize',36)
colorbar, title('good - bad Cue PAC','FontSize',36), set(gca,'clim',datRange)
datRange        = ceil(max(max(abs(nanmean(out.goodPAC-out.badPAC,3).*10000))))./10000;
datRange        = [-datRange datRange];
subplot(1,2,1)
contourf(out.phsFreq,out.ampFreq,nanmean(out.goodPAC-out.badPAC,3)',linspace(datRange(1),datRange(2),50),'linestyle','none')
colormap jet, axis square, set(gca,'FontSize',24), xlabel('Phase Frequency (Hz)','FontSize',36), ylabel('Amplitude Frequency (Hz)','FontSize',36)
colorbar, title('Cue - Baseline PAC','FontSize',36), set(gca,'clim',datRange)
clear datRange

figure
datRange        = ceil(max(max(abs(nanmean(out.cuePAC,3).*10000))))./10000;
datRange        = [-datRange datRange];
subplot(2,2,1)
contourf(out.phsFreq,out.ampFreq,nanmean(out.cuePAC,3)',linspace(datRange(1),datRange(2),50),'linestyle','none')
colormap jet, axis square, set(gca,'FontSize',24), xlabel('Phase Frequency (Hz)','FontSize',36), ylabel('Amplitude Frequency (Hz)','FontSize',36)
colorbar, title('Cue PAC','FontSize',36), set(gca,'clim',datRange)
datRange        = ceil(max(max(abs(nanmean(out.basePAC,3).*10000))))./10000;
datRange        = [-datRange datRange];
subplot(2,2,2)
contourf(out.phsFreq,out.ampFreq,nanmean(out.basePAC,3)',linspace(datRange(1),datRange(2),50),'linestyle','none')
colormap jet, axis square, set(gca,'FontSize',24), xlabel('Phase Frequency (Hz)','FontSize',36), ylabel('Amplitude Frequency (Hz)','FontSize',36)
colorbar, title('Baseline PAC','FontSize',36), set(gca,'clim',datRange)
datRange        = ceil(max(max(abs(nanmean(out.goodPAC,3).*10000))))./10000;
datRange        = [-datRange datRange];
subplot(2,2,3)
contourf(out.phsFreq,out.ampFreq,nanmean(out.goodPAC,3)',linspace(datRange(1),datRange(2),50),'linestyle','none')
colormap jet, axis square, set(gca,'FontSize',24), xlabel('Phase Frequency (Hz)','FontSize',36), ylabel('Amplitude Frequency (Hz)','FontSize',36)
colorbar, title('Baseline PAC','FontSize',36), set(gca,'clim',datRange)
datRange        = ceil(max(max(abs(nanmean(out.goodPAC,3).*10000))))./10000;
datRange        = [-datRange datRange];
subplot(2,2,4)
contourf(out.phsFreq,out.ampFreq,nanmean(out.badPAC,3)',linspace(datRange(1),datRange(2),50),'linestyle','none')
colormap jet, axis square, set(gca,'FontSize',24), xlabel('Phase Frequency (Hz)','FontSize',36), ylabel('Amplitude Frequency (Hz)','FontSize',36)
colorbar, title('Good - Bad Cue PAC','FontSize',36), set(gca,'clim',datRange)
clear datRange

%stats
statlowFreqs   = [6 9];
stathighFreqs  = [60 110]; 

flInds       = [find(out.phsFreq>=statlowFreqs(1),1,'first') find(out.phsFreq<=statlowFreqs(2),1,'last')];
fhInds       = [find(out.ampFreq>=stathighFreqs(1),1,'first') find(out.ampFreq<=stathighFreqs(2),1,'last')];
cue_base    = out.cuePAC-out.basePAC;
intPower    = squeeze(mean(mean(cue_base(flInds(1):flInds(2),fhInds(1):fhInds(2),:)))); 
%goodPower   = squeeze(mean(mean(allgoodCFC(flInds(1):flInds(2),fhInds(1):fhInds(2),:)))); 
%badPower    = squeeze(mean(mean(allbadCFC(flInds(1):flInds(2),fhInds(1):fhInds(2),:)))); clear fInds tInds

[h,p,~,stats] = ttest(intPower);