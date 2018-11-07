%% Smooth fluorescence trace and find potential peaks/troughs in nc14

%specifically look at bursts in nc14
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);

for n=1:length(CompiledParticles_14)

    if length(CompiledParticles_14(n).Fluo) >=3        
        PotentialPeak(n).Frame=[];
        PotentialPeak(n).Fluo=[];
        PotentialTrough(n).Frame=[];
        PotentialTrough(n).Fluo=[];
        SmoothParticles=smooth(ElapsedTime(CompiledParticles_14(n).Frame),CompiledParticles_14(n).Fluo,0.1,'lowess');
        Baseline=min(SmoothParticles);
        Tolerance=0.5*10^4;
        for ii=1:length(CompiledParticles_14(n).Frame)
            %dealing with endpoints
            if (ii==1) & (SmoothParticles(ii+1) >= SmoothParticles(ii)) & (SmoothParticles(ii)-Baseline <= Tolerance)
                PotentialTrough(n).Frame=[PotentialTrough(n).Frame, CompiledParticles_14(n).Frame(ii)];
                PotentialTrough(n).Fluo= [PotentialTrough(n).Fluo,SmoothParticles(ii)];
            end
            if ii==length(CompiledParticles_14(n).Frame) & SmoothParticles(ii-1) >= SmoothParticles(ii)
                PotentialTrough(n).Frame=[PotentialTrough(n).Frame, CompiledParticles_14(n).Frame(ii)];
                PotentialTrough(n).Fluo=[PotentialTrough(n).Fluo, SmoothParticles(ii)];
            end

            if (ii~=1) & (ii~=length(CompiledParticles_14(n).Frame))
                if (SmoothParticles(ii) >= SmoothParticles(ii-1)) & (SmoothParticles(ii)>= SmoothParticles(ii+1))
                    PotentialPeak(n).Frame=[PotentialPeak(n).Frame,CompiledParticles_14(n).Frame(ii)];
                    PotentialPeak(n).Fluo=[PotentialPeak(n).Fluo,SmoothParticles(ii)];
                elseif (SmoothParticles(ii) <= SmoothParticles(ii-1)) & (SmoothParticles(ii) <= SmoothParticles(ii+1))
                    PotentialTrough(n).Frame=[PotentialTrough(n).Frame,CompiledParticles_14(n).Frame(ii)];
                    PotentialTrough(n).Fluo=[PotentialTrough(n).Fluo,SmoothParticles(ii)];
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
            if (PeakHeight-Baseline) >= 2*TroughHeight
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
%         if length(BurstAmplitude)==0 & (~isempty(CompiledParticles_14(n).TotalmRNA)) %Deal with spots where made mRNA but high end point makes not recognize as a burst
%             temp2TroughFrame=[CompiledParticles_14(n).Frame(1), CompiledParticles_14(n).Frame(end)];
%             temp2TroughFluo=[CompiledParticles_14(n).Fluo(1), CompiledParticles_14(n).Fluo(end)];
%             BurstAmplitude=max([CompiledParticles_14(n).Fluo]);
%         end
        if length(BurstAmplitude)==0 & (~isempty(CompiledParticles_14(n).TotalmRNA))
            temp2TroughFluo=[(min([CompiledParticles_14(n).Fluo])), CompiledParticles_14(n).Fluo(end)];
            minFluoFrame=find([CompiledParticles_14(n).Fluo]==min([CompiledParticles_14(n).Fluo]));
            temp2TroughFrame=[CompiledParticles_14(n).Frame(minFluoFrame), CompiledParticles_14(n).Frame(end)];
            BurstAmplitude=max([CompiledParticles_14(n).Fluo(minFluoFrame:end)]);
        end
        
        BurstProperties(n).MeanAP=CompiledParticles_14(n).MeanAP;
        BurstProperties(n).TotalmRNA=CompiledParticles_14(n).TotalmRNA;
        BurstProperties(n).TotalmRNAError=CompiledParticles_14(n).TotalmRNAError;
        BurstProperties(n).BurstAmplitude=BurstAmplitude;
        BurstProperties(n).Troughs=temp2TroughFluo;
        BurstProperties(n).TroughFrame=temp2TroughFrame;
        BurstProperties(n).Duration=(ElapsedTime(temp2TroughFrame(2:end))-ElapsedTime(temp2TroughFrame(1:end-1)));
        BurstProperties(n).Baseline=Baseline;
        BurstProperties(n).NBursts=length(BurstAmplitude);
        BurstProperties(n).Nucleus=CompiledParticles_14(n).Nucleus;
        BurstProperties(n).SpotLife=(ElapsedTime(CompiledParticles_14(n).Frame(end))-(ElapsedTime(CompiledParticles_14(n).Frame(1))));
        BurstProperties(n).FirstTimeOn=(ElapsedTime(CompiledParticles_14(n).Frame(1))- ElapsedTime(nc14));
        BurstProperties(n).FluoFrames=CompiledParticles_14(n).Frame;
        BurstProperties(n).Fluo=CompiledParticles_14(n).Fluo;
        
        if length(BurstAmplitude) ~=0
            BurstProperties(n).Frequency=length(BurstAmplitude)/(ElapsedTime(temp2TroughFrame(end))-ElapsedTime(temp2TroughFrame(1)));
            %BurstProperties(n).Frequency2=length(BurstAmplitude)/(ElapsedTime(end)-ElapsedTime(nc14)); %number of burst per time of nc14
            BurstProperties(n).Frequency3=length(BurstAmplitude)/(ElapsedTime(CompiledParticles_14(n).Frame(end))-ElapsedTime(CompiledParticles_14(n).Frame(1))); %over time of first recorded fluo to last 
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
