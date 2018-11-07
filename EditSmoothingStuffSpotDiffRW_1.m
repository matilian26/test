%% Smooth fluorescence trace and find potential peaks/troughs in nc14

%specifically look at bursts in nc14
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);
Nucinfo=[CompiledParticles_14.Nucleus];
for n=1:length(SpotDiff)
    CompParEquiv=CompiledParticles_14(Nucinfo==SpotDiff(n).Nucleus);
    if length(SpotDiff(n).Frame1) >=3 & (length(SpotDiff(n).SpotOne)==length(ElapsedTime))       
        PotentialPeak(n).Frame=[];
        PotentialPeak(n).Fluo=[];
        PotentialTrough(n).Frame=[];
        PotentialTrough(n).Fluo=[];
        SmoothParticles=smoothdata([SpotDiff(n).SpotOne],'lowess','omitnan','SmoothingFactor',0.1);
        %SmoothParticles=smooth(ElapsedTime,SpotDiff(n).SpotOne,0.1,'lowess');
        Baseline=min(SmoothParticles);
        Tolerance=0.5*10^4;
        for ii=1:length(SpotDiff(n).SpotOne)
            %dealing with endpoints
            if (ii==1) & (SmoothParticles(ii+1) >= SmoothParticles(ii)) & (SmoothParticles(ii)-Baseline <= Tolerance)
                PotentialTrough(n).Frame=[PotentialTrough(n).Frame, ii];
                PotentialTrough(n).Fluo= [PotentialTrough(n).Fluo,SmoothParticles(ii)];
            end
            if ii==length(SpotDiff(n).SpotOne) & SmoothParticles(ii-1) >= SmoothParticles(ii)
                PotentialTrough(n).Frame=[PotentialTrough(n).Frame, ii];
                PotentialTrough(n).Fluo=[PotentialTrough(n).Fluo, SmoothParticles(ii)];
            end

            if (ii~=1) & (ii~=length(SpotDiff(n).SpotOne))
                if (SmoothParticles(ii) >= SmoothParticles(ii-1)) & (SmoothParticles(ii)>= SmoothParticles(ii+1))
                    PotentialPeak(n).Frame=[PotentialPeak(n).Frame,ii];
                    PotentialPeak(n).Fluo=[PotentialPeak(n).Fluo,SmoothParticles(ii)];
                elseif (SmoothParticles(ii) <= SmoothParticles(ii-1)) & (SmoothParticles(ii) <= SmoothParticles(ii+1))
                    PotentialTrough(n).Frame=[PotentialTrough(n).Frame,ii];
                    PotentialTrough(n).Fluo=[PotentialTrough(n).Fluo,SmoothParticles(ii)];
%                 elseif (SmoothParticles(ii) <= SmoothParticles(ii-1)) & (SmoothParticles(ii)<= (Tolerance-Baseline)) 
%                     PotentialTrough(n).Frame=[PotentialTrough(n).Frame,ii];
%                     PotentialTrough(n).Fluo=[PotentialTrough(n).Fluo,SmoothParticles(ii)];
%                 elseif (SmoothParticles(ii) <= SmoothParticles(ii+1)) & (SmoothParticles(ii) <= (Tolerance-Baseline))
%                     PotentialTrough(n).Frame=[PotentialTrough(n).Frame,ii];
%                     PotentialTrough(n).Fluo=[PotentialTrough(n).Fluo, SmoothParticles(ii)];
                end
            end
        end

        %
    %     plot(ElapsedTime(PotentialPeak(n).Frame),PotentialPeak(n).Fluo,'og')
    %     plot(ElapsedTime(PotentialTrough(n).Frame),PotentialTrough(n).Fluo,'oc')

        % Determine and keep real troughs 
        numTrough=[1500,1000];
        while(numTrough(end-1) ~= numTrough(end))
            if exist('tempTroughFluo','var')
                TroughPoolFluo=tempTroughFluo;
                TroughPoolFrame=tempTroughFrame;
            else 
                TroughPoolFluo=PotentialTrough(n).Fluo;
                TroughPoolFrame=PotentialTrough(n).Frame;
            end
            tempTroughFluo = [];
            tempTroughFrame = [];
            PeakStatus=[];
            for ii=1:(length(TroughPoolFluo)-1)
                MaxTrough=max(TroughPoolFluo([ii:ii+1]));
                TroughHeight=MaxTrough-Baseline;
                PeakHeight=PotentialPeak(n).Fluo(ii);
                if (PeakHeight-Baseline) >= 2*TroughHeight
                    PeakStatus=[PeakStatus, 1];    
                    ExistFrame=find(tempTroughFrame==TroughPoolFrame(ii));
                    if isempty(ExistFrame)
                        tempTroughFluo=[tempTroughFluo,TroughPoolFluo(ii:ii+1)];
                        tempTroughFrame = [tempTroughFrame,TroughPoolFrame(ii:ii+1)];
                    else
                        tempTroughFluo=[tempTroughFluo,TroughPoolFluo(ii+1)];
                        tempTroughFrame=[tempTroughFrame,TroughPoolFrame(ii+1)];
                    end
                else
                    PeakStatus=[PeakStatus,0];
                    ExistFrame=find(tempTroughFrame==TroughPoolFrame(ii));
                    if isempty(ExistFrame)
                        tempTroughFluo=[tempTroughFluo,TroughPoolFluo(ii:ii+1)];
                        tempTroughFrame = [tempTroughFrame,TroughPoolFrame(ii:ii+1)];
                    else
                        tempTroughFluo=[tempTroughFluo,TroughPoolFluo(ii+1)];
                        tempTroughFrame=[tempTroughFrame,TroughPoolFrame(ii+1)];
                    end
                end
                if ii >= 2
                    if (PeakStatus(ii-1)==0 & PeakStatus(ii)==0) 
                        tempTroughFluo(ii)=[];
                        tempTroughFrame(ii)=[];
                    end
                end


            end
            numTrough=[numTrough, length(tempTroughFluo)];
        end
        
        temp2TroughFluo = [];
        temp2TroughFrame = [];
        PeakStatus=[];
        BurstAmplitude=[];
        for ii=1:(length(tempTroughFluo)-1)
             MaxTrough=max(tempTroughFluo([ii:ii+1]));
            TroughHeight=MaxTrough-Baseline;
            %finding frame of peak inbetween this trough and next one 
            Frameindex=find([PotentialPeak(n).Frame >= tempTroughFrame(ii) & PotentialPeak(n).Frame <= tempTroughFrame(ii+1)]);

            PeakHeight=max(PotentialPeak(n).Fluo(Frameindex));
            if (PeakHeight-TroughHeight) >= 2*TroughHeight
                PeakStatus=[PeakStatus, 1];  
                BurstAmplitude=[BurstAmplitude,PeakHeight];
                ExistFrame=find(temp2TroughFrame==tempTroughFrame(ii));
                if isempty(ExistFrame)
                    temp2TroughFluo=[temp2TroughFluo,tempTroughFluo(ii:ii+1)];
                    temp2TroughFrame = [temp2TroughFrame,tempTroughFrame(ii:ii+1)];
                else
                    temp2TroughFluo=[temp2TroughFluo,tempTroughFluo(ii+1)];
                    temp2TroughFrame=[temp2TroughFrame,tempTroughFrame(ii+1)];
                end
            else
                PeakStatus=[PeakStatus,0];
            end

            if ii >= 2
                if ii==2 || ii==(length(tempTroughFluo)-1)
                    if (PeakStatus(ii-1)==1) & (PeakStatus(ii)==0)
                        maxtroughindex=find([tempTroughFluo==max(tempTroughFluo([ii,ii+1]))]);
                    if maxtroughindex==ii
                        temp2TroughFluo(end)=tempTroughFluo(ii+1);
                        temp2TroughFrame(end)=tempTroughFrame(ii+1);
                    end
                    end
                    if (PeakStatus(ii-1)==0) & (PeakStatus(ii)==1)
                        maxtroughindex=find([tempTroughFluo==max(tempTroughFluo([ii-1,ii]))]);
                        if maxtroughindex==ii
                            temp2TroughFluo(end-1)=tempTroughFluo(ii-1);
                            temp2TroughFrame(end-1)=tempTroughFrame(ii-1);
                        end
                    end
                end


                if (PeakStatus(ii-1)==0 && PeakStatus(ii)==1) 
                    maxtroughindex=find([tempTroughFluo==max(tempTroughFluo([ii-1,ii]))]);
                    if maxtroughindex==ii
                        temp2TroughFluo(end-1)=tempTroughFluo(ii-1);
                        temp2TroughFrame(end-1)=tempTroughFrame(ii-1);
                    end
                end
            end


        end
        %record the frames of each burst peak
        %RealBursts=PotentialPeak(n).Fluo(PotentialPeak(n).Fluo~=0);
        BurstyFrame=find(PotentialPeak(n).Fluo(ismember(PotentialPeak(n).Fluo,BurstAmplitude)));
        BurstProperties(n).BurstFrames=PotentialPeak(n).Frame(BurstyFrame);
        BurstProperties(n).BtwnBurstTime=(ElapsedTime(BurstProperties(n).BurstFrames(2:end))-ElapsedTime(BurstProperties(n).BurstFrames(1:end-1)));
        
            %Record other properites 
        BurstProperties(n).MeanAP=SpotDiff(n).MeanAP;
        BurstProperties(n).TotalmRNA=SpotDiff(n).TotalmRNAOne;
        BurstProperties(n).TotalmRNAError=CompiledParticles_14(n).TotalmRNAError;
        BurstProperties(n).BurstAmplitude=BurstAmplitude;
        BurstProperties(n).Troughs=temp2TroughFluo;
        BurstProperties(n).TroughFrame=temp2TroughFrame;
        BurstProperties(n).Duration=(ElapsedTime(temp2TroughFrame(2:end))-ElapsedTime(temp2TroughFrame(1:end-1)));
        BurstProperties(n).Baseline=Baseline;
        BurstProperties(n).NBursts=length(BurstAmplitude);
        BurstProperties(n).Nucleus=SpotDiff(n).Nucleus;
        BurstProperties(n).SpotLife=(ElapsedTime(CompiledParticles_14(n).Frame(end))-(ElapsedTime(CompiledParticles_14(n).Frame(1))));
        BurstProperties(n).FirstTimeOn=(ElapsedTime(CompiledParticles_14(n).Frame(1))- ElapsedTime(nc14));
        BurstProperties(n).FluoFrames=CompiledParticles_14(1).Frame;
        BurstProperties(n).Fluo=CompiledParticles_14(n).Fluo;
        BurstProperties(n).FullFluoOne=SpotDiff(n).SpotOne;
        BurstProperties(n).SmoothFluo=SmoothParticles;
        
        
%         for bb=1:length(BurstProperties(n).BurstFrames)
%             BurstProperties(n).Decay(bb)=(ElapsedTime(BurstProperties(n).BurstFrames(bb))-ElapsedTime(BurstProperties(n).TroughFrame(bb+1)));
%             BurstProperties(n).Onset(bb)=(ElapsedTime(BurstProperties(n).BurstFrames(bb))-ElapsedTime(BurstProperties(n).TroughFrame(bb)));
%         end
        if length(BurstAmplitude) ~=0
            BurstProperties(n).Frequency=length(BurstAmplitude)/(ElapsedTime(temp2TroughFrame(end))-ElapsedTime(temp2TroughFrame(1)));
            BurstProperties(n).Frequency2=length(BurstAmplitude)/(ElapsedTime(end)-ElapsedTime(nc14)); %number of burst per time of nc14
        else 
            BurstProperties(n).Frequency=0;
            BurstProperties(n).Frequency2=0;
        end
    end
    
clear 'tempTroughFluo'
end

%get rid of empty nucleus values
for ii=1:length(BurstProperties)
    if isempty(BurstProperties(ii).Nucleus)
    BurstProperties(ii).Nucleus=nan;
    end
end
%add up total mRNA from both spots in a nucleus
BurstNuclei=unique([BurstProperties.Nucleus]);
for nn=1:length(BurstNuclei)
    tempnuclei=find([BurstProperties.Nucleus]==BurstNuclei(nn));
    if length(tempnuclei) > 1
        BurstProperties(tempnuclei(1)).BothAllelesProd=((BurstProperties(tempnuclei(1)).TotalmRNA) + (BurstProperties(tempnuclei(2)).TotalmRNA));
        BurstProperties(tempnuclei(2)).BothAllelesProd=((BurstProperties(tempnuclei(1)).TotalmRNA) + (BurstProperties(tempnuclei(2)).TotalmRNA));   
        BurstProperties(tempnuclei(1)).BothAllelesError=((BurstProperties(tempnuclei(1)).TotalmRNAError) + (BurstProperties(tempnuclei(2)).TotalmRNAError));
        BurstProperties(tempnuclei(2)).BothAllelesError=((BurstProperties(tempnuclei(1)).TotalmRNAError) + (BurstProperties(tempnuclei(2)).TotalmRNAError));
    elseif length(tempnuclei)==1
        BurstProperties(tempnuclei(1)).BothAllelesProd=BurstProperties(tempnuclei(1)).TotalmRNA;
        BurstProperties(tempnuclei(1)).BothAllelesError=BurstProperties(tempnuclei(1)).TotalmRNAError;
    end
end

for tt=1:length(BurstProperties)
    if ~isempty(BurstProperties(tt).Duration)
    APEstm(tt)=round(BurstProperties(tt).MeanAP,3);
    for jj=1:length(APbinID)
        if APEstm(tt) < APbinID(jj)
            BurstProperties(tt).APBin=APbinID(jj);
            break;
        end
    end
    end
end

for tt=1:length(BurstProperties)
    if isempty([BurstProperties(tt).APBin])
        BurstProperties(tt).APBin=nan;
    end
end

%save bursting info 
%save('BurstPropertiesnc14','BurstProperties')