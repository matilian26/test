%% just take derivatives of smoothed traces 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrBothSep32C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

% go through each embryo of each construct
for cc=1:length(ConstructList)
     Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    ConRawSlopes=[];
    ConPosSlopes=[];
    ConNegSlopes=[];
    ConFluoPoints=[];
    for ee=1:Datalength(cc)
        PrefixName=Data(ee).Prefix;
        Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        Filename2=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat']
        load(Filename);
        load(Filename2);
        EmbRawSlopes=[];
        EmbPosSlopes=[];
        EmbNegSlopes=[];
        EmbFluoPoints=[];
%specifically look at bursts in nc14
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);

for n=1:length(SpotDiff)
    if length(SpotDiff(n).Frame1) >=3 & length(SpotDiff(n).SpotOne)==length(ElapsedTime)
    SmoothParticles=smooth(ElapsedTime,SpotDiff(n).SpotOne,0.1,'lowess');
    Derivs(n).RawSlopes=diff(SmoothParticles)./[diff(ElapsedTime)]';
    Derivs(n).PositiveSlopes=[Derivs(n).RawSlopes(Derivs(n).RawSlopes >=0)];
    Derivs(n).NegSlopes=[Derivs(n).RawSlopes(Derivs(n).RawSlopes < 0)];
    Fluos(n).FluoPoints=[SpotDiff(n).SpotOne];
    EmbRawSlopes=[EmbRawSlopes; [Derivs(n).RawSlopes]];
    EmbPosSlopes=[EmbPosSlopes; [Derivs(n).PositiveSlopes]];
    EmbNegSlopes=[EmbNegSlopes; [Derivs(n).NegSlopes]];
    EmbFluoPoints=[EmbFluoPoints; [Fluos(n).FluoPoints]];
    end
end
AllSlopes(cc).Embryo(ee).RawSlopes=[EmbRawSlopes];
AllSlopes(cc).Embryo(ee).PosSlopes=[EmbPosSlopes];
AllSlopes(cc).Embryo(ee).NegSlopes=[EmbNegSlopes];
ConRawSlopes=[ConRawSlopes; EmbRawSlopes];
ConPosSlopes=[ConPosSlopes; EmbPosSlopes];
ConNegSlopes=[ConNegSlopes; EmbNegSlopes];
ConFluoPoints=[ConFluoPoints; EmbFluoPoints];
    end
    AllSlopes(cc).AllRawSlopes=[ConRawSlopes];
    AllSlopes(cc).AllPosSlopes=ConPosSlopes;
    AllSlopes(cc).AllNegSlopes=ConNegSlopes;
    AllSlopes(cc).AllFluoPoints=ConFluoPoints;
end

clear SmoothParticles

%% Try plotting traces with different cutoff values 
for n=1:length(SpotDiff)
    if length(SpotDiff(n).Frame1) >=3 & length(SpotDiff(n).SpotOne)==length(ElapsedTime)
       %try limiting frames to nc14 to reduce extreme smoothing
       SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:end),SpotDiff(n).SpotOne(nc14:end),0.1,'lowess');
    %SmoothParticles(n).Smoothed=smooth(ElapsedTime,SpotDiff(n).SpotOne,0.1,'lowess');
    SmoothParticles(n).Smoothed(SmoothParticles(n).Smoothed<0)=0;
    end
end
%Q=input('Choose spot to look at');
P=input('Choose ON threshold');
Q=input('Choose OFF threshold');
for ss=1:15:length(SmoothParticles)
    if ~isempty(SmoothParticles(ss).Smoothed);
figure 
plot(SmoothParticles(ss).Smoothed);
hold on 
dx=diff(SmoothParticles(ss).Smoothed)./([diff(ElapsedTime(nc14:end))]');
%dx=diff(SmoothParticles(ss).Smoothed)./([diff(ElapsedTime)]');
ONThresholdSlope=dx>=P;
OFFThresholdSlope=dx<=Q;
AboveLine=[SmoothParticles(ss).Smoothed];
AboveLine(~ONThresholdSlope)=nan;
AboveLine(OFFThresholdSlope)=nan;
BelowLine=[SmoothParticles(ss).Smoothed];
BelowLine(ONThresholdSlope)=nan;
BelowLine(~OFFThresholdSlope)=nan;
plot(AboveLine,'r','Linewidth',2);
plot(BelowLine,'g','Linewidth',2);
title(['spot',num2str(ss),' ',num2str(P),' ','ON threshold',num2str(Q),'OFF threshold']);
legend('trace','ON','OFF');
    end
end

