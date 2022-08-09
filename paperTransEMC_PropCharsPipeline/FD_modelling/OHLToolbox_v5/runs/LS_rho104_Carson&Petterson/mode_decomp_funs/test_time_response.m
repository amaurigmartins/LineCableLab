SampleTime=1e-5;
TotalSampleNumber=10^6;
InputTime=double((1:TotalSampleNumber)')*SampleTime;

%InputSignal=ones(1,TotalSampleNumber);

InputSignal=zeros(1,TotalSampleNumber);
InputSignal(1)=10^20;

[y,t]=timeresp(fit(1,4),InputSignal,SampleTime);

plot(t,y)