%% Burst calling based on slope threshold 

%Look at nc14 only
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);

for n=1:length(SpotDiff)
    if length(SpotDiff(n).Frame1) >=3 & length(SpotDiff(n).SpotOne)==length(ElapsedTime)
       %try limiting frames to nc14 to reduce extreme smoothing
       SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:end),SpotDiff(n).SpotOne(nc14:end),0.1,'lowess');
    %fluorescence should never be negative but sometimes with lots of 0s
    %looks like the smoothing f(x) makes it so
       SmoothParticles(n).Smoothed(SmoothParticles(n).Smoothed<0)=0;
       SmoothParticles(n).RawSlope=diff(SmoothParticles(n).Smoothed)./([diff(ElapsedTime(nc14:end))]');
%Slope thresholds
ONThreshold=1500; %2000;  %6/20/18 thinking ~288AU=1 polymerase so 1500 ~ 5ish polymerases
OFFThreshold=-1500; %-2000;
ONFluoLimit=900;
OFFFluoLimit=250;
       %OFFThresholdSlope=dx<=Q;
% AboveLine=[SmoothParticles(ss).Smoothed];
% AboveLine(~ONThresholdSlope)=nan;
% AboveLine(OFFThresholdSlope)=nan;
% BelowLine=[SmoothParticles(ss).Smoothed];
% BelowLine(ONThresholdSlope)=nan;
% BelowLine(~OFFThresholdSlope)=nan;
StartFramePool=[];
StartSlopePool=[];
EndSlopePool=[];
Starts=[find(SmoothParticles(n).RawSlope>=ONThreshold)];
Starts=min(Starts);
StartFramePool=[StartFramePool,Starts];
EndPts=[find(SmoothParticles(n).RawSlope<=OFFThreshold)];
EndPts=min(EndPts);
EndFramePool=[1000,EndPts];
if ~isempty(EndPts) & (~isempty(Starts))
while EndFramePool(end) ~= EndFramePool(end-1)
    for ss=EndFramePool(end):(length(SmoothParticles(n).RawSlope)-3)
        if SmoothParticles(n).RawSlope(ss) >=ONThreshold | (SmoothParticles(n).Smoothed(ss:(ss+3)) >= ONFluoLimit &(SmoothParticles(n).RawSlope(ss) >OFFThreshold))
            StartFramePool=[StartFramePool,ss];
            StartSlopePool=[StartSlopePool, SmoothParticles(n).RawSlope(ss)];
            break
%         elseif ss <= (length(SmoothParticles(n).RawSlope)-3)
%             if SmoothParticles(n).Smoothed(ss:(ss+3)) >= ONFluoLimit & (SmoothParticles(n).RawSlope(ss) > OFFThreshold)
%                 StartFramePool=[StartFramePool,ss];
%                 StartSlopePool=[StartSlopePool, SmoothParticles(n).RawSlope(ss)];
%                 break
%             end
        end
    end
    for ss=StartFramePool(end):(length(SmoothParticles(n).RawSlope)-3)
        if SmoothParticles(n).RawSlope(ss) <=OFFThreshold | (SmoothParticles(n).Smoothed(ss:(ss+3)) <= OFFFluoLimit)
            EndFramePool=[EndFramePool, ss];
            EndSlopePool=[EndSlopePool,SmoothParticles(n).RawSlope(ss)];
            break
%         elseif ss <= (length(SmoothParticles(n).RawSlope)-3)
%             if SmoothParticles(n).Smoothed(ss:(ss+3))<= OFFFluoLimit
%                 EndFramePool=[EndFramePool, ss];
%                 EndSlopePool=[EndSlopePool, SmoothParticles(n).RawSlope(ss)];
%                 break
%             end
        end
    end
    if sum(SmoothParticles(n).RawSlope(EndFramePool(end):end)>=ONThreshold)==0
%         EndFramePool(end+1)=EndFramePool(end);
%         StartFramePool(end+1)=StartFramePool(end);
        PotentialPeaks(n).Starts=[StartFramePool];
        PotentialPeaks(n).Ends=[EndFramePool(2:end)];
        break
    elseif sum(SmoothParticles(n).RawSlope(StartFramePool(end):end)<=OFFThreshold)==0
%         EndFramePool(end+1)=EndFramePool(end);
%         StartFramePool(end+1)=StartFramePool(end);
        PotentialPeaks(n).Starts=[StartFramePool];
        
        PotentialPeaks(n).Ends=[EndFramePool(2:end)];
        break
    end
end
if length(PotentialPeaks(n).Ends)>1 & (PotentialPeaks(n).Ends(end)==PotentialPeaks(n).Ends(end-1))
    PotentialPeaks(n).Ends=[PotentialPeaks(n).Ends(1:end-1)];
end
if length(PotentialPeaks(n).Starts)>1 & (PotentialPeaks(n).Starts(end)==PotentialPeaks(n).Starts(end-1))
    PotentialPeaks(n).Starts=[PotentialPeaks(n).Starts(1:end-1)];
end
BurstProperties(n).SmoothTrace=SmoothParticles(n).Smoothed;
% not sure if should say frame is equal to slope frame or +1, for rn
% setting as equal 
%If start frames > end frames, ignore the last start since we can't say how
%long the burst lasts
if length(PotentialPeaks(n).Starts) > (length(PotentialPeaks(n).Ends))
    BurstProperties(n).Duration=ElapsedTime(PotentialPeaks(n).Ends)-ElapsedTime(PotentialPeaks(n).Starts(1:end-1));
    BurstProperties(n).Duration=BurstProperties(n).Duration(BurstProperties(n).Duration>=0);
    BurstProperties(n).ONFrames=[PotentialPeaks(n).Starts(1:end-1)];
else
BurstProperties(n).Duration=ElapsedTime(PotentialPeaks(n).Ends)-ElapsedTime(PotentialPeaks(n).Starts);
BurstProperties(n).Duration=BurstProperties(n).Duration(BurstProperties(n).Duration>=0);
BurstProperties(n).ONFrames=[PotentialPeaks(n).Starts];
end

BurstProperties(n).OFFFrames=[PotentialPeaks(n).Ends];
if PotentialPeaks(n).Ends(1)==1
    BurstProperties(n).BurstAmplitude=[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends(2:end)]-1)]';
else
BurstProperties(n).BurstAmplitude=[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends]-1)]';
end
BurstProperties(n).FirstTimeOn=ElapsedTime(PotentialPeaks(n).Starts(1));
BurstProperties(n).NBursts=length(BurstProperties(n).BurstAmplitude);
BurstProperties(n).Frequency=(BurstProperties(n).NBursts)/(ElapsedTime(end)-ElapsedTime(PotentialPeaks(n).Starts(1)));
BurstProperties(n).FrequencyAct=(BurstProperties(n).NBursts)/(ElapsedTime(PotentialPeaks(n).Ends(end))-ElapsedTime(PotentialPeaks(n).Starts(1)));
BurstProperties(n).APBin=SpotDiff(n).APBin;
BurstProperties(n).Nucleus=SpotDiff(n).Nucleus;
BurstProperties(n).TotalmRNA=SpotDiff(n).TotalmRNAOne;
BurstProperties(n).FractON=sum(BurstProperties(n).Duration)/(length(ElapsedTime(nc14:end))); 


end
    end
end
