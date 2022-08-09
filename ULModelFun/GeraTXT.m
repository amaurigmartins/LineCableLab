% GeraTXT: gera arquivo .txt com informações dos parâmetros da linha e ajuste vetorial.

filename='C:\ATP\work\fitULM000001.txt';

fid = fopen(filename,'wt');
fprintf(fid,'%d\n',NFase); %numero de fases
fprintf(fid,'%d\n',NFase); %numero de modos
fprintf(fid,'%d\n',size(polesTr,1)); %numero de polos Yc
fprintf(fid,'%d\n',size(pPj,1)); %numero de polos A

%polYcvet
for ii = 1:(length(polesTr))
    fprintf(fid,'%.16e\t%.16e\n',real(polesTr(ii)),imag(polesTr(ii)));
end
%resYc
for kk = 1:(size(rYc,3))
    for ii = 1:(size(rYc,1))
        for jj = 1:(size(rYc,2))
            fprintf(fid,'%.16e\t%.16e\n',real(rYc(ii,jj,kk)),imag(rYc(ii,jj,kk)));
        end
    end
end
%polAvet
for ii = 1:(length(pPj))
    fprintf(fid,'%.16e\t%.16e\n',real(pPj(ii)),imag(pPj(ii)));
end
%resA
for mm = 1:(size(rA,4))
    for kk = 1:(size(rA,3))
        for ii = 1:(size(rA,1))
            for jj = 1:(size(rA,2))
                fprintf(fid,'%.16e\t%.16e\n',real(rA(ii,jj,kk,mm)),imag(rA(ii,jj,kk,mm)));
            end
        end
    end
end
% taus otimizados
for jj = 1:(length(tauOtim))
    fprintf(fid,'%.16e\n',tauOtim(jj));
end
% resid0Yc
for ii = 1:(size(resid0Yc,1))
    for jj = 1:(size(resid0Yc,2))
        fprintf(fid,'%.16e\t%.16e\n',real(resid0Yc(ii,jj)),imag(resid0Yc(ii,jj)));
    end
end
fclose(fid);