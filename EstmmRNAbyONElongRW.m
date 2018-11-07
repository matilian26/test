%% Estimating expected mRNA based on FractON and elongation rate estimate
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'Kr2xProxEmpty';'KrBoth';'KrBothEmpty';'KrDist32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'HbEmpty'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgTimes=[];
    Avgnc14Times=[];
    FluoTimes=[];
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        filename=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(filename);
        AvgTimes=[AvgTimes,ElapsedTime(end)];
        Avgnc14Times=[Avgnc14Times, (ElapsedTime(end)-ElapsedTime(nc14))];
        nc_number=[CompiledParticles.nc];
        CompiledParticles_14=CompiledParticles(nc_number==14);
        for n=1:length(CompiledParticles_14)
            FluoTimes=[FluoTimes, (ElapsedTime(CompiledParticles_14(n).Frame(end))-ElapsedTime(CompiledParticles_14(n).Frame(1)))];
        end
    end
    AverageTime(cc).AvgAllTimes=nanmean([AvgTimes]); 
    AvgTime(cc).nc14Avg=nanmean(Avgnc14Times);
    AvgTime(cc).AvgFluoLength=nanmean(FluoTimes);
end

%% Estimate mRNA production using fractON and RNAP elongation rate estimate
FractFilename=[DropboxFolder filesep 'Constructs' filesep 'FractON.mat'];
load(FractFilename);
ElongationRate=2000; %range of 1.5-6kb/minute
PolySpacing=(1/150);
for cc=1:length(ConstructList)
    MinAct=(AvgFractONAllAP(cc).AvgFractON .*50);%(AvgTime(cc).nc14Avg));
    %TranscriptsMade=(1.5/5.6).*(AvgTime(cc).AvgFluoLength);
    TransPerMinute=(ElongationRate * PolySpacing)
    TranscriptsMade=(MinAct .* TransPerMinute);
    %TranscriptsMade=(1.5/5.6)*MinAct; % 1.5kb/min and 5.6kb length gene
    %KbTranscribed=MinAct .* 3.5;  % estm of 1.5kb/min elongation
    Transcripts(cc).NumTrans=TranscriptsMade;
    %Transcripts(cc).NumTrans=TranscriptsMade .*20 %56 is estimate of # polymerases could be on gene body saying each takes up 100nt
    %Transcripts(cc).NumTrans=(KbTranscribed ./ 5.6).*20; %20 is estm based
    %on Convoy RNA pol paper spacing/rate data
end
 TotmRNACalc=[DropboxFolder filesep 'Constructs' filesep 'AllTotalmRNAProd.mat'];
 load(TotmRNACalc);
 AmpCalc=[DropboxFolder filesep 'Constructs' filesep 'BurstAmplitude'];
 load(AmpCalc);
 FRNAP=input('Which FRNAP value to use?');
 ConvVal=FRNAP*15; %15 being estm of nc13 length;

for cc=1:length(ConstructList)
    figure
    %plot(([AvgAmpAllAP(cc).AvgAmp]./377));
    plot(([AvgProdAllAP(cc).AvgProd]./ConvVal)) %3246 is based on max slope AU values divided by loading rates in (Garcia 2013)- max they saw was ~20RNAP/min --->10RNAP/timept 
    hold on 
    plot(Transcripts(cc).NumTrans,'--');
    title(ConstructList{cc});
end