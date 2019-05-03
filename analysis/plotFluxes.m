% workerID=4;

text_leg=cell(4,1);

for workerID=1:4;
    
    time=Coupled(workerID).OUT(5).timestamp;
    FluxRight=Coupled(workerID).OUT(5).lateral.water_fluxes(workerID,workerID+1,:);
    FluxRight=permute(FluxRight,[3 2 1]);
    fprintf('Sum gained by w%1.0f : %3.2e m\n',workerID,nansum(FluxRight))
    dat=FluxRight~=0;
    FluxRight=FluxRight(dat);
    timeRight=time(dat);
    
    plot(timeRight, FluxRight)
    hold on
    text_leg{workerID,1}=sprintf('Line %1.0f',workerID);
    
end
legend(text_leg)
datetick
hold off