%% plot mean fluorescence of a construct across space at it's peak exp time in NC13 vs NC14
KrConstructs

for ii=1:length(KrData)
    figure
    [M,I]=max(KrData(ii).MeanNC13);  %find which time point has the highest exp in NC13
    I=I(I~=1);
    if ~isempty(I)
    Time13Use=round(mean(I));
    else
        Time13Use=20;
    end
    [M,II]=max(KrData(ii).MeanNC14);
    II=II(II~=1);
    Time14Use=round(mean(II));
    plot(KrData(ii).MeanNC13(Time13Use,:))
    hold on 
    plot(KrData(ii).MeanNC14(Time14Use,:));
    legend('NC13','NC14','Location','best');
    title([ConstructList{ii},'NC13',' ', num2str(Time13Use),' ','NC14',' ',num2str(Time14Use)]);
    ylim([0 40000]);
    xlim([15 30])
    xlabel('AP bin');
    ylabel('Mean Fluorescence');
end
    