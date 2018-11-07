n = 350;
figure;
plot(ElapsedTime(CompiledParticles(n).Frame), CompiledParticles(n).Fluo)
hold on
plot(ElapsedTime(CompiledParticles(n).Frame), smooth(ElapsedTime(CompiledParticles(n).Frame),CompiledParticles(n).Fluo,0.1,'lowess'),'*')

%% Smooth fluorescence trace and find potential peaks/troughs 
for n=1:length(CompiledParticles)

    if length(CompiledParticles(n).Fluo) >=3        
        PotentialPeak(n).Frame=[];
        PotentialPeak(n).Fluo=[];
        PotentialTrough(n).Frame=[];
        PotentialTrough(n).Fluo=[];
        SmoothParticles=smooth(ElapsedTime(CompiledParticles(n).Frame),CompiledParticles(n).Fluo,0.1,'lowess');
        Baseline=min(SmoothParticles);
        Tolerance=0.5*10^4;
        for ii=1:length(CompiledParticles(n).Frame)
            if (ii==1) & (SmoothParticles(ii+1) >= SmoothParticles(ii)) & (SmoothParticles(ii)-Baseline <= Tolerance)
                PotentialTrough(n).Frame=[PotentialTrough(n).Frame, CompiledParticles(n).Frame(ii)];
                PotentialTrough(n).Fluo= [PotentialTrough(n).Fluo,SmoothParticles(ii)];
            end
            if ii==length(CompiledParticles(n).Frame) & SmoothParticles(ii-1) >= SmoothParticles(ii)
                PotentialTrough(n).Frame=[PotentialTrough(n).Frame, CompiledParticles(n).Frame(ii)];
                PotentialTrough(n).Fluo=[PotentialTrough(n).Fluo, SmoothParticles(ii)];
            end

            if (ii~=1) & (ii~=length(CompiledParticles(n).Frame))
                if (SmoothParticles(ii) >= SmoothParticles(ii-1)) & (SmoothParticles(ii)>= SmoothParticles(ii+1))
                    PotentialPeak(n).Frame=[PotentialPeak(n).Frame,CompiledParticles(n).Frame(ii)];
                    PotentialPeak(n).Fluo=[PotentialPeak(n).Fluo,SmoothParticles(ii)];
                elseif (SmoothParticles(ii) <= SmoothParticles(ii-1)) & (SmoothParticles(ii) <= SmoothParticles(ii+1))
                    PotentialTrough(n).Frame=[PotentialTrough(n).Frame,CompiledParticles(n).Frame(ii)];
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
                if PeakHeight >= 2*TroughHeight
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
            Frameindex=find([PotentialPeak(n).Frame >= tempTroughFrame(ii) & PotentialPeak(n).Frame <= tempTroughFrame(ii+1)]);

            PeakHeight=max(PotentialPeak(n).Fluo(Frameindex));
            if PeakHeight >= 2*TroughHeight
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
        BurstProperties(n).BurstAmplitude=BurstAmplitude;
        BurstProperties(n).Troughs=temp2TroughFluo;
        BurstProperties(n).TroughFrame=temp2TroughFrame;
        BurstProperties(n).Duration=(ElapsedTime(temp2TroughFrame(2:end))-ElapsedTime(temp2TroughFrame(1:end-1)));
        if length(BurstAmplitude) ~=0
            BurstProperties(n).Frequency=length(BurstAmplitude)/(ElapsedTime(temp2TroughFrame(end))-ElapsedTime(temp2TroughFrame(1)));
        else 
            BurstProperties(n).Frequency=0;
        end
    end
    
clear 'tempTroughFluo'
end

save('BurstProperties','BurstProperties')




