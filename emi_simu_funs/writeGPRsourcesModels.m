function [outM] = writeGPRsourcesModels(G,NodeT,NodeP,N_LT,N_DT)
%% Criação do cartão ATP models

    %Teste para compatibilidade de tamanho dos vetores G, numT e NumP    
    [numT,numP] = size(G);

%     if N_LT ~= numP || (N_DT+1) ~= numT
%         error('Wrong G dimensions.');
%         return;
%     end

    %inicialização de variáveis de Saída
    outM = [];
    outT = [];
    
%     %Início da escrita do cartão Models
%     outM = strcat(outM,'/MODELS\n');
%     outM = strcat(outM,'MODELS\n');
    outM = strcat(outM,'INPUT\n');

    %Escrita do trecho INPUT, da Models
    for k =1:size(NodeT,2)
        outM = strcat(outM,sprintf(' MM%04d {i(%s)}',k,NodeT{k}),'\n');
    end

%    Escrita do Trecho OUTPUT da Models
    outM = strcat(outM,'OUTPUT\n');

    for k = 1:size(NodeP,2)
        outM = strcat(outM,sprintf(' G%04d',k),'\n');
    end

    %Definição do nome da Models
    outM = strcat(outM,'MODEL GROUNDCALC\nINPUT');

    %Identificação das entradas da Models
    for k = 1:(size(NodeT,2))
        if k == 1
            outM = strcat(outM,sprintf(' F%d,',k));
        elseif k == size(NodeT,2)
            outM = strcat(outM,sprintf('F%d',k));
        elseif mod(k,5) == 0
            outM = strcat(outM,sprintf('F%d,',k),'\n');    
        else
            outM = strcat(outM,sprintf('F%d,',k));
        end
    end

    %Identificação das saídas do bloco Models
    outM = strcat(outM,'\nOUTPUT\x20');

    cont = 1;
    for k = 1:size(NodeP,2)
        if k == (9*cont - 8)  && k ~= size(NodeP,2)          %Imprime espaço após a quebra de linha (Para múltiplos de 9)
            outM = strcat(outM,sprintf('OUT%04d,',k));
            cont = cont+1;
        elseif mod(k,9) == 0                            %Imprime a quebra de linha em k igual a multiplos de 9
            outM = strcat(outM,sprintf('OUT%04d,',k),'\n\x20\x20\x20\x20\x20\x20\x20'); %,sprintf(' ')
        elseif k == size(NodeP,2) 
            outM = strcat(outM,sprintf('OUT%04d',k),'\n\x20\x20\x20\x20\x20\x20\x20');
        else
            outM = strcat(outM,sprintf('OUT%04d,',k));
        end
    end

    %Identificação das Variáveis do Bloco
    for k = 1:size(NodeP,2)
        if k == 1
            outM = strcat(outM,sprintf('\nVAR OUT%04d',k),'\n');
        else 
            outM = strcat(outM,sprintf('    OUT%04d',k),'\n');
        end
    end

    %Inicialização das variáveis do bloco
    outM = strcat(outM,sprintf('    G[1..%d]\nINIT\n',numT*numP),'\n');

    for k = 1:size(NodeP,2)
       outM = strcat(outM,sprintf('   OUT%04d:= 0.0',k),'\n');
    end

    %Conversão dos valores de uma matriz quadrada para uma de linha única
    %Devido a restrições do ATP
    for i = 1:numT
        for j = 1:numP
            G_models((i-1)*numP+j) = G(i,j);
        end
    end
    
    %Inicialização e inserção dos valores da função de Green no bloco
    for k = 1:numT*numP
        if k == numT*numP
            outM = strcat(outM,sprintf('   G[%d]:= %d',k,G_models(k)),'\n','ENDINIT\n');
        else
            outM = strcat(outM,sprintf('   G[%d]:= %d',k,G_models(k)),'\n');
        end
    end

    %Escrita da execução feita pelo Bloco Models
    outM = strcat(outM,'EXEC\n');

    
    %Escrita da Soma com inclusão dos valores da Função de Green
    cont_g = 1;
    
    for k = 1:(size(NodeP,2))
        cont2 = 2;
        outM = strcat(outM,sprintf('  OUT%04d:=',k));
        for j = 1:(size(NodeT,2))
            if j == 1
                outM = strcat(outM,sprintf(' G[%d]*F%d+',cont_g,j));
                cont_g = cont_g + 1;
            elseif j == ((size(NodeT,2)))
                if j == (4*cont2 - 3)
                    tt='\t\x20\x20\x20\x20';
                else
                    tt='';
                end
                outM = strcat(outM,tt,sprintf('G[%d]*F%d',cont_g,j),'\n');
                cont_g = cont_g + 1;
            elseif mod(j,4) == 0
                outM = strcat(outM,sprintf('G[%d]*F%d+',cont_g,j),'\n');
                cont_g = cont_g + 1;
            elseif j == (4*cont2 - 3)
                outM = strcat(outM,'\t\x20\x20\x20',sprintf('\x20G[%d]*F%d+',cont_g,j));
                cont_g = cont_g + 1;
                cont2 = cont2 + 1;
            else
                outM = strcat(outM,sprintf('G[%d]*F%d+',cont_g,j));
                cont_g = cont_g + 1;
            end
        end    

    end

    %Finalização da execução e do Bloco Models
    outM = strcat(outM, 'ENDEXEC\n', 'ENDMODEL\n','USE GROUNDCALC AS DEFAULT\n','INPUT\n');

    %COnexão entre a saída da Models e o circuito, por meio dos nomes dos
    %nós
    
    for i = 1:(size(NodeT,2))
        outM = strcat(outM, sprintf(' F%d:= MM%04d',i,i),'\n');
    end
    
    outM = strcat(outM, 'OUTPUT\n');
    
    for i = 1:(size(NodeP,2))
        outM = strcat(outM, sprintf(' G%04d:= OUT%04d',i,i),'\n');
    end
    
    outM = strcat(outM, 'ENDUSE\n');%,'ENDMODELS\n');
    
    
outM=sprintf(outM);

end