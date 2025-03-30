%% Script to plot mean high frequency amplitude by theta phase in human
%  iEEG data, for Umesh
%  Daniel Bush, UCL (2018) drdanielbush@gmail.com

%  Provide some parameters for the analysis
folder              = '/Users/umie/Documents/Matlab/New Cue Period Data/';
folder2             = '/Users/umie/Documents/Matlab/old_cue/';          %patientNos          = [2 4 5 6 8 9 10 12 13 14 15 16 19];
patientNos          = [2 5 6 8 9 10 13 14 15];                                  % List of patient IDs to analyse
cueWin              = [6 8];                                            % Time window for the analysis
region              = 'RealHPC';                                                % Region to be analysed
fOrder              = 2;                                                        % Filter order (use 2 as a default)
freqBands           = [30:1:50,50:2:100,100:5:120];
phaseBand           = [2 5];                                                    % Frequency band for phase (i.e. theta band)
ampBand1            = [30 50];                                                  % Frequency band for amplitude 1 (i.e. low gamma band)
ampBand2            = [60 120];                                                 % Frequency band for ampitude 2 (i.e. high gamma band)
phaseBins           = 15;                                                       % Number of phase bins (use 30 as a default)
downFs              = 512;

%  Provide settings for the wavelet transform
Settings.tWin       = cueWin;
Settings.freqBand   = [30 120];
Settings.logTrans   = 1;
Settings.baseline   = [];
Settings.zscore     = 1;

%  Then cycle through each patient
for p               = 1 : length(patientNos)
    
    % Load the raw data and extract trial information
    filename        = ['Patient ' int2str(patientNos(p)) '_Cue_' region '_-2to4s_UmeshERPData.mat'];
    load([folder2 filename],'Data'); clear filename
    errors          = Data.Error; clear Data
    filename        = ['Patient ' int2str(patientNos(p)) '_Cue_' region '_-5to5s_UmeshERPData.mat'];
    load([folder filename],'Data');
    filename        = [folder Data.Patient];
    header          = read_edf_file('header',filename);
    trials          = Data.Artifacts;
    Artifacts       = Data.Artifacts;
    Fs              = (unique(Data.DataPoints)-1) / range(Data.Timebase);
    trials          = [Data.Time(:,1:2)*Fs trials];
    channels        = Data.Contacts; clear Data
    clear header
    
    for c               = 1:length(channels) % ADD FOR C = 1 : LENGTH(CHANNELS)
                
        % Pre-process
        cfg             = [];
        cfg.dataset     = filename;
        cfg.channel     = channels(c);
        cfg.trl         = trials;
        cfg             = ft_definetrial(cfg);
        epoch_data      = ft_preprocessing(cfg);
        
        % Reduce to artefact free trials
        trial           = find(Artifacts(:,c)==0);
        cfg             = [];
        errors_trial    = errors(trial);
        goodTrials      = errors_trial <= median(errors_trial);
        badTrials       = errors_trial > median(errors_trial); clear errors_trial
        cfg.trials      = trial;
        Sel_epoch_data	= ft_preprocessing(cfg,epoch_data); clear epoch_data
        
        %downsample
        cfg = [];
        cfg.resamplefs  = downFs;
        cfg.detrend     = 'yes';
        cfg.demean      = 'yes';
        Sel_epoch_data  = ft_resampledata(cfg,Sel_epoch_data); clear cfg
        
        % Compute the time-frequency representations
        cfg             = [];
        %cfg.channel     = channels(c); % REPLACE 'CHANNELS' WITH 'CHANNELS(C)'
        cfg.method      = 'wavelet';    
        cfg.width       = 5;
        cfg.trials      = 'all'; 
        cfg.keeptrials  = 'yes';
        cfg.output      = 'pow';
        cfg.foi         = freqBands;
        cfg.t_ftimwin(1:length(cfg.foi)) = 0.2;
        cfg.tapsmofrq(1:length(cfg.foi)) = 10;
        cfg.toi         = 5:1/512:8;    
        Power           = ft_freqanalysis(cfg, Sel_epoch_data); clear cfg 
        Power.powspctrm = log(Power.powspctrm);
        
        
        % Filter the data in the low frequency band
        [b,a]           = butter(fOrder,phaseBand.*2/downFs,'bandpass');
        timeWin         = [find(Sel_epoch_data.time{1}>=cueWin(1),1,'first') find(Sel_epoch_data.time{1}<=cueWin(2),1,'last')];
        for t           = 1 : length(Sel_epoch_data.trial)
            Sel_epoch_data.phaseSig{t}      = filtfilt(b,a,Sel_epoch_data.trial{t}')';
            Sel_epoch_data.phaseBin{t}      = angle(hilbert(Sel_epoch_data.phaseSig{t}'))';
            [~,Sel_epoch_data.phaseBin{t}]  = histc(Sel_epoch_data.phaseBin{t},linspace(-pi,pi,phaseBins+1));
            Sel_epoch_data.phaseBin{t}      = Sel_epoch_data.phaseBin{t}(:,timeWin(1):timeWin(2));
            Sel_epoch_data.phaseSig{t}      = Sel_epoch_data.phaseSig{t}(:,timeWin(1):timeWin(2));            
        end
        clear b a t
        
        % Filter the data in each gamma frequency band
        % Filter the data in the low frequency band
        [b,a]           = butter(fOrder,ampBand1.*2/downFs,'bandpass');        
        for t           = 1 : length(Sel_epoch_data.trial)
            Sel_epoch_data.amp1{t}          = filtfilt(b,a,Sel_epoch_data.trial{t}')';
            Sel_epoch_data.amp1{t}          = abs(hilbert(Sel_epoch_data.amp1{t}'))';            
            Sel_epoch_data.amp1{t}          = Sel_epoch_data.amp1{t}(:,timeWin(1):timeWin(2));            
            if Settings.logTrans
                Sel_epoch_data.amp1{t}   	= log(Sel_epoch_data.amp1{t});            
            end
        end
        clear b a t
        
        % Filter the data in the low frequency band
        [b,a]           = butter(fOrder,ampBand2.*2/downFs,'bandpass');        
        for t           = 1 : length(Sel_epoch_data.trial)
            Sel_epoch_data.amp2{t}          = filtfilt(b,a,Sel_epoch_data.trial{t}')';
            Sel_epoch_data.amp2{t}          = abs(hilbert(Sel_epoch_data.amp2{t}'))';            
            Sel_epoch_data.amp2{t}          = Sel_epoch_data.amp2{t}(:,timeWin(1):timeWin(2));     
            if Settings.logTrans
                Sel_epoch_data.amp2{t}   	= log(Sel_epoch_data.amp2{t});            
            end
        end
        clear b a t timeWin
                
        % Extract the time-frequency representation for each trial / channel,
        % and then average the amplitude in each frequency across each phase
        % bin        
        for t                   = 1 : length(Sel_epoch_data.trial)
            
            
            %[Power,Freq,~]      = DanWavelet(Sel_epoch_data.trial{t},Sel_epoch_data.time{t},downFs,Settings);
                              
            % Assign memory for all patient data
            if p==1
                fAxis           = Power.freq;
                pwrByPhase      = zeros(length(Power.freq),phaseBins,length(patientNos));
                goodpwrByPhase  = zeros(length(Power.freq),phaseBins,length(patientNos));
                badpwrByPhase   = zeros(length(Power.freq),phaseBins,length(patientNos));
                thetaByPhase    = zeros(phaseBins,length(patientNos));
                gamma1ByPhase   = zeros(phaseBins,length(patientNos));
                gamma2ByPhase   = zeros(phaseBins,length(patientNos));                
                gammaCorr       = zeros(length(patientNos),1);
                trialNos        = zeros(length(patientNos),1);
            end
            
            % Assign memory for this patient's data
            if c == 1
                pPwrByPhase     = zeros(length(Power.freq),phaseBins,size(Sel_epoch_data.trial{t},1));
                goodpPwrByPhase = zeros(length(Power.freq),phaseBins,size(Sel_epoch_data.trial{t},1));
                badpPwrByPhase  = zeros(length(Power.freq),phaseBins,size(Sel_epoch_data.trial{t},1));
                pThetaByPhase   = zeros(phaseBins,length(channels));
                pGamma1ByPhase  = zeros(phaseBins,length(channels));
                pGamma2ByPhase  = zeros(phaseBins,length(channels));
                pGammaCorr      = zeros(length(channels),1);                
            end
            
            % Assign memory for this channel's data
            if t	== 1
                trialPowByPhase	= nan(length(Power.freq),phaseBins,length(Sel_epoch_data.trial));
                tThetaByPhase   = zeros(phaseBins,length(Sel_epoch_data.trial));
                tGamma1ByPhase  = zeros(phaseBins,length(Sel_epoch_data.trial));
                tGamma2ByPhase  = zeros(phaseBins,length(Sel_epoch_data.trial));
                tGammaCorr      = zeros(length(Sel_epoch_data.trial),1);
                if c == 1
                    pTrialNos	= zeros(length(channels),1);
                end
                pTrialNos(c,1)  = length(Sel_epoch_data.trial);
            end
            
            % Split the trial into theta cycles
            [~,cycLocs]         = findpeaks(Sel_epoch_data.phaseSig{t});
            cycInds             = zeros(length(Sel_epoch_data.phaseSig{t}),1);
            cycInds([1 cycLocs+1]) = 1;  clear cycLocs
            cycInds             = cumsum(cycInds);
            
            % Compute low and high gamma power for each cycle
       
            lowGammaP           = squeeze(mean(Power.powspctrm(t,:,Power.freq>=ampBand1(1) & Power.freq<=ampBand1(2),(1:length(cycInds)))));
            lowGammaP           = accumarray(cycInds,lowGammaP,[],@mean);
            highGammaP          = squeeze(mean(Power.powspctrm(t,:,Power.freq>=ampBand2(1) & Power.freq<=ampBand2(2),(1:length(cycInds)))));
            highGammaP          = accumarray(cycInds,highGammaP,[],@mean); 
            
            % Compute the correlation between low and high gamma power on a
            % cycle by cycle basis for that trial
            b                   = glmfit(lowGammaP,highGammaP);
            tGammaCorr(t,1)     = b(2); clear b lowGammaP highGammaP
            
            % Cycle through each frequency, compute average power in each
            % phase bin
            Pow = squeeze(mean(Power.powspctrm(:,:,:,(1:length(cycInds)),1))); clear cycInds
            for f   = 1 : length(fAxis)
                trialPowByPhase(f,:,t)	= accumarray(Sel_epoch_data.phaseBin{t}',Pow(f,:)',[phaseBins 1],@mean)';
            end
            clear f Pow Freq 
            
            % Also compute mean theta signal and gamma amplitude in each phase bin
            tThetaByPhase(:,t)  = accumarray(Sel_epoch_data.phaseBin{t}',Sel_epoch_data.phaseSig{t}',[phaseBins 1],@mean)';
            tGamma1ByPhase(:,t) = accumarray(Sel_epoch_data.phaseBin{t}',Sel_epoch_data.amp1{t}',[phaseBins 1],@mean)';
            tGamma2ByPhase(:,t) = accumarray(Sel_epoch_data.phaseBin{t}',Sel_epoch_data.amp2{t}',[phaseBins 1],@mean)';
            
        end
        
        % Compute the mean z-scored power in each frequency and each phase bin
        if Settings.zscore                  
            trialPowByPhase     = (trialPowByPhase - repmat(mean(trialPowByPhase,2),[1 size(trialPowByPhase,2)])) ./ repmat(std(trialPowByPhase,[],2),[1 size(trialPowByPhase,2)]);                                    
            tGamma1ByPhase      = (tGamma1ByPhase - repmat(mean(tGamma1ByPhase),size(tGamma1ByPhase,1),1)) ./ repmat(std(tGamma1ByPhase),size(tGamma1ByPhase,1),1);
            tGamma2ByPhase      = (tGamma2ByPhase - repmat(mean(tGamma2ByPhase),size(tGamma2ByPhase,1),1)) ./ repmat(std(tGamma2ByPhase),size(tGamma2ByPhase,1),1);
        end
            
        % Average over all artefact free trials for that channel      
        pPwrByPhase(:,:,c)   	= mean(trialPowByPhase,3);        
        goodpPwrByPhase(:,:,c)	= mean(trialPowByPhase(:,:,goodTrials),3);
        badpPwrByPhase(:,:,c)   = mean(trialPowByPhase(:,:,badTrials),3); clear trialPowByPhase Power
        pThetaByPhase(:,c)      = mean(tThetaByPhase,2); clear tThetaByPhase
        pGamma1ByPhase(:,c)     = mean(tGamma1ByPhase,2); clear tGamma1ByPhase
        pGamma2ByPhase(:,c)     = mean(tGamma2ByPhase,2); clear tGamma2ByPhase
        pGammaCorr(c,1)         = mean(tGammaCorr); clear tGammaCorr
        
    end
    clear trial goodTrials badTrials filename channels c Artifacts errors trials Sel_epoch_data
    
    %Average data over channels for that participant
    pwrByPhase(:,:,p)           = mean(pPwrByPhase,3); clear t c sbjFirst %pPwrByPhase epoch_data trials
    goodpwrByPhase(:,:,p)       = mean(goodpPwrByPhase,3);
    badpwrByPhase(:,:,p)        = mean(badpPwrByPhase,3);
    thetaByPhase(:,p)           = mean(pThetaByPhase,2);
    gamma1ByPhase(:,p)          = mean(pGamma1ByPhase,2);
    gamma2ByPhase(:,p)          = mean(pGamma2ByPhase,2);
    gammaCorr(p,1)              = mean(pGammaCorr); 
    clear pPwrByPhase goodpPwrByPhase badpPwrByPhase pGamma1ByPhase pGamma2ByPhase pThetaByPhase pGammaCorr
    clc
end


%% Establish whether low and high gamma power are negatively correlated on a trial by trial basis
[~,p,~,stats]       = ttest(gammaCorr); % Note that we want a negative t-statistic!


%% Look at distribution of power by phase across participants
pAxis               = linspace(-pi+pi/phaseBins,pi-pi/phaseBins,phaseBins); clear allFirst phaseBins

gamma1ByPhaseB      = squeeze(mean(pwrByPhase(fAxis>=ampBand1(1) & fAxis<=ampBand1(2),:,:)));
gamma2ByPhaseB      = squeeze(mean(pwrByPhase(fAxis>=ampBand2(1) & fAxis<=ampBand2(2),:,:)));
 
allPrefPhase1       = circ_mean(repmat(pAxis',1,length(patientNos)),gamma1ByPhase);
allPrefPhase2       = circ_mean(repmat(pAxis',1,length(patientNos)),gamma2ByPhase);
allPrefPhase1B      = circ_mean(repmat(pAxis',1,length(patientNos)),gamma1ByPhaseB);
allPrefPhase2B      = circ_mean(repmat(pAxis',1,length(patientNos)),gamma2ByPhaseB);

[allP1,allZ1]     	= circ_rtest(allPrefPhase1);
[allP2,allZ2]     	= circ_rtest(allPrefPhase2);

[allP1B,allZ1B]     = circ_rtest(allPrefPhase1B);
[allP2B,allZ2B]     = circ_rtest(allPrefPhase2B);

[pval, table]       = circ_wwtest([allPrefPhase1' ; allPrefPhase2'], [ones(length(patientNos),1) ; 2*ones(length(patientNos),1)]);
[pvalB, tableB]     = circ_wwtest([allPrefPhase1B' ; allPrefPhase2B'], [ones(length(patientNos),1) ; 2*ones(length(patientNos),1)]);

%Plot the output (for participants with > 10 trials only!)
figure
contourf(pAxis,fAxis,mean(pwrByPhase,3),100,'linestyle','none')
%axis square
%set(gca,'clim',[-0.1 0.1])
set(gca,'linewidth',1.5)
%set(gca,'Xticklabel',[], 'Yticklabel', []) 
%set(gca,'FontSize',24)
%xlabel('Phase (rad)','FontSize',24)
%ylabel('Frequency (Hz)','FontSize',24)
c = colorbar;
c = colorbar;
c.LineWidth = 1.5;
set(c,'Yticklabel', [])
width=500;
height=300;
set(gcf,'position',[0,0,width,height])
%title('All trials','FontSize',24)
%ylabel(c,'Power','FontSize',24); clear c
set(gca,'YScale','log')
colormap jet

%figure
figure
shadedErrorBar(pAxis,mean(thetaByPhase,2),std(thetaByPhase,[],2)./sqrt(length(patientNos)))
%axis square
xlim ([-3 3])
set(gca,'linewidth',1.5)
set(gca,'Xticklabel',[], 'Yticklabel', []) 
%set(gca,'FontSize',18)
width=500;
height=150;
set(gcf,'position',[0,0,width,height])
%title('Average Theta Signal','FontSize',24)
%xlabel('Phase (rad)','FontSize',24)
%ylabel('Mean Amplitude (uV)','FontSize',24)

figure
subplot(1,2,1)
contourf(pAxis,fAxis,mean(goodpwrByPhase,3),50,'linestyle','none')
%axis square
set(gca,'clim',[-0.2 0.2])
set(gca,'FontSize',24)
xlabel('Phase (rad)','FontSize',24)
ylabel('Frequency (Hz)','FontSize',24)
c = colorbar;
title('Good trials','FontSize',24)
ylabel(c,'Amplitude (z-scored)','FontSize',24); clear c
colormap jet

subplot(1,2,2)
contourf(pAxis,fAxis,mean(badpwrByPhase,3),50,'linestyle','none')
axis square
set(gca,'clim',[-0.2 0.2])
set(gca,'FontSize',24)
xlabel('Phase (rad)','FontSize',24)
ylabel('Frequency (Hz)','FontSize',24)
c = colorbar;
title('Bad trials','FontSize',24)
ylabel(c,'Amplitude (z-scored)','FontSize',24); clear c
colormap jet

figure
subplot(2,2,1)
shadedErrorBar(pAxis,mean(gamma1ByPhase,2),std(gamma1ByPhase,[],2)./sqrt(length(patientNos)))
axis square
set(gca,'FontSize',18)
title('Low Gamma','FontSize',24)
ylabel('Mean Amplitude','FontSize',24)

subplot(2,2,2)
shadedErrorBar(pAxis,mean(gamma2ByPhase,2),std(gamma2ByPhase,[],2)./sqrt(length(patientNos)))
axis square
set(gca,'FontSize',18)
title('High Gamma','FontSize',24)

subplot(2,2,3)
shadedErrorBar(pAxis,mean(gamma1ByPhaseB,2),std(gamma1ByPhaseB,[],2)./sqrt(length(patientNos)))
axis square
set(gca,'FontSize',18)
xlabel('Phase (rad)','FontSize',24)
ylabel('Mean Amplitude','FontSize',24)

subplot(2,2,4)
shadedErrorBar(pAxis,mean(gamma2ByPhaseB,2),std(gamma2ByPhaseB,[],2)./sqrt(length(patientNos)))
axis square
set(gca,'FontSize',18)
xlabel('Phase (rad)','FontSize',24)