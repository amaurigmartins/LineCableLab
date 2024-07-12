function [] = makeATPXML(f,Z,Y,jobid,currPath)

Zprint=Z;
Yprint=Y;
freq=f;
length=1; %1 meter because TOOLBOX output: Z [Ohm/m] | Y [Mho/m]

%Hardcoded Y Printing option
%ATP='C';
ATP='G+Bi';

%XML Printing
filename = fullfile(currPath,[jobid '_ZYATP.xml']);
fid = fopen(filename,'wt');

%Header
if ATP == 'C'
    fprintf(fid,'<ZY NumPhases="%d" Length="%.4f" ZFmt="R+Xi" YFmt="C">\n',size(Zprint,1),length);
else
    fprintf(fid,'<ZY NumPhases="%d" Length="%.4f" ZFmt="R+Xi" YFmt="G+Bi">\n',size(Zprint,1),length);
end

%Z Printhing
for kk = 1:size(freq,1);    
    fprintf(fid,'  <Z Freq="%.16E">\n',freq(kk));
        for ii = 1:size(Zprint,1);
            for jj = 1:size(Zprint,2);
                if jj ~= size(Zprint,2)
                    fprintf(fid,'%s,',num2str(Zprint(ii,jj,kk),'%.16E'));
                else
                    fprintf(fid,'%s',num2str(Zprint(ii,jj,kk),'%.16E'));
                end
            end
            fprintf(fid,'\n');
        end
fprintf(fid,'  </Z>\n');
end

%Y Printhing
if ATP == 'C'
    fprintf(fid,'  <Y Freq="%.16E">\n',freq(1));
    for ii = 1:size(Yprint,1);
        for jj = 1:size(Yprint,2);
            if jj ~= size(Yprint,2)
                fprintf(fid,'%s,',num2str(imag(Yprint(ii,jj,1)/(2*pi*freq(1))),'%.16E'));
            else
                fprintf(fid,'%s',num2str(imag(Yprint(ii,jj,1)/(2*pi*freq(1))),'%.16E'));
            end
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'  </Y>\n');    
else
    for kk = 1:size(freq,1);    
        fprintf(fid,'  <Y Freq="%.16E">\n',freq(kk));
            for ii = 1:size(Yprint,1);
                for jj = 1:size(Yprint,2);
                    if jj ~= size(Yprint,2)
                        fprintf(fid,'%s,',num2str(Yprint(ii,jj,kk),'%.16E'));
                    else
                        fprintf(fid,'%s',num2str(Yprint(ii,jj,kk),'%.16E'));
                    end
                end
                fprintf(fid,'\n');
            end
    fprintf(fid,'  </Y>\n');
    end
end

fprintf(fid,'</ZY>');
fclose(fid);
end