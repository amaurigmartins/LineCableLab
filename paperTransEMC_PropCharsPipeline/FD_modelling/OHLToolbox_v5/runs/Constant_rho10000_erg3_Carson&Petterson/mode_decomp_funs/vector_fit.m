function [Ti_vf,fit,check]=vector_fit(Ti,ord,num_files,freq)

Ti_temp=zeros(num_files,ord);
Ti_vf_temp=zeros(num_files,ord);
Ti_vf=zeros(num_files,ord^2);

check=zeros(ord,ord);

fprintf('Rational fit\n');
for k=1:ord
    fprintf('mode #%i\n',k);
    for o=1:ord
        Ti_temp(:,o)=Ti(:,k+(o-1)*ord);
    end
        h=rationalfit(freq,Ti_temp,-100,[],[],false,10,24,true);
        fit(k,1:ord)=h(1:ord); % Προσοχή! Δεν έχει γίνει preallocation!
    for o=1:ord
        [Ti_vf_temp(:,o),freq]=freqresp(h(o),freq);
    end
    for o=1:ord
        Ti_vf(:,k+(o-1)*ord)=Ti_vf_temp(:,o);
    end
end
    
for k=1:ord
    for o=1:ord
        check(k,o)=ispassive(fit(k,o));
    end
end