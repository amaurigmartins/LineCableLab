function [] = exportCrossSection2XML(app,path)
            %output file name
            jobid=app.JobIDEditField.Value;
            fname=fullfile(path,['ATPDraw_' jobid '.xml']);

            % DECLARE SOME USEFUL CONSTANTS AND VARIABLES
            DEFTSTART=-1;
            DEFTSTOP=1000;
            freq=60;
            Leq = app.LinelengthmEditField.Value;
            rho=app.ResistivitymEditField.Value;

            linecable_data = cell2matButNotStupid(app,app.UITable.Data);

            numphases=size(linecable_data,1);
            TopLeftX=200; % 200 + 240 LCCs per row give the optimal layout
            TopLeftY=0; %will need this later to place the LCCs

            % writes ATPdraw header
            outstruct.ApplicationAttribute='ATPDraw';
            outstruct.VersionAttribute=7;
            outstruct.VersionXMLAttribute=1;
            outstruct.header.TimestepAttribute=1e-6;
            outstruct.header.TmaxAttribute=.1;
            outstruct.header.XOPTAttribute=0;
            outstruct.header.COPTAttribute=0;
            outstruct.header.SysFreqAttribute=freq;
            outstruct.header.TopLeftXAttribute=TopLeftX;
            outstruct.header.TopLeftYAttribute=TopLeftY;

            % now the damn components section
            offsetX=80;
            offsetY=360;
            tower=1;
            LCC_per_row=20; %240;
            rowcount=1; %counts number of times that LCC_per_row is reached
            colcount=1; % will be reset to 1 every LCC_per_row
            conn_idx=1;

            numblocks=1;

            for i=1:numblocks
                outstruct.objects.comp(i).NameAttribute='LCC';
                outstruct.objects.comp(i).IdAttribute=string([jobid '_' num2str(i)]);
                outstruct.objects.comp(i).CapanglAttribute=90;
                outstruct.objects.comp(i).CapPosXAttribute=-10;
                outstruct.objects.comp(i).CapPosYAttribute=-25;
                outstruct.objects.comp(i).comp_content.PosXAttribute=TopLeftX+(colcount*offsetX);
                outstruct.objects.comp(i).comp_content.PosYAttribute=TopLeftY+(rowcount*offsetY);

                colcount=colcount+1;
                thiscaption='';
                outstruct.objects.comp(i).CaptionAttribute=thiscaption;

                if mod(i,LCC_per_row)==0
                    colcount=1;
                    rowcount=rowcount+1;
                end

                outstruct.objects.comp(i).comp_content.NumPhasesAttribute=numphases;
                outstruct.objects.comp(i).comp_content.IconAttribute='default';
                outstruct.objects.comp(i).comp_content.SinglePhaseIconAttribute='true';

                %EACH ONE OF THESE BIG BOYS NEEDS Nph NODES
                y0=-20;
                for k=1:numphases %input nodes
                    y0=y0+10;
                    outstruct.objects.comp(i).comp_content.node(k).NameAttribute=string(sprintf('IN%d',k));
                    outstruct.objects.comp(i).comp_content.node(k).ValueAttribute=string(sprintf('C%dSND',k));
                    % outstruct.objects.comp(i).comp_content.node(k).ValueAttribute=string(sprintf('C%d%04d',k,i));
                    outstruct.objects.comp(i).comp_content.node(k).UserNamedAttribute='true';
                    %         outstruct.objects.comp(i).comp_content.node(k).NumPhasesAttribute=1;
                    outstruct.objects.comp(i).comp_content.node(k).KindAttribute=k;
                    outstruct.objects.comp(i).comp_content.node(k).PosXAttribute=-20;
                    outstruct.objects.comp(i).comp_content.node(k).PosYAttribute=y0;
                    outstruct.objects.comp(i).comp_content.node(k).NamePosXAttribute=0;
                    outstruct.objects.comp(i).comp_content.node(k).NamePosYAttribute=0;
                end

                y0=-20;
                for k=numphases+1:2*numphases %output nodes
                    w=k-numphases;
                    y0=y0+10;
                    outstruct.objects.comp(i).comp_content.node(k).NameAttribute=string(sprintf('OUT%d',w));
                    outstruct.objects.comp(i).comp_content.node(k).ValueAttribute=string(sprintf('C%dRCV',w));
                    % outstruct.objects.comp(i).comp_content.node(k).ValueAttribute=string(sprintf('C%d%04d',w,i+1));
                    outstruct.objects.comp(i).comp_content.node(k).UserNamedAttribute='true';
                    %         outstruct.objects.comp(i).comp_content.node(k).NumPhasesAttribute=1;
                    outstruct.objects.comp(i).comp_content.node(k).KindAttribute=w;
                    outstruct.objects.comp(i).comp_content.node(k).PosXAttribute=20;
                    outstruct.objects.comp(i).comp_content.node(k).PosYAttribute=y0;
                    outstruct.objects.comp(i).comp_content.node(k).NamePosXAttribute=0;
                    outstruct.objects.comp(i).comp_content.node(k).NamePosYAttribute=0;
                end

                % yeah this fucking sucks and probably can be done better
                outstruct.objects.comp(i).comp_content.data(1).NameAttribute='Length';
                outstruct.objects.comp(i).comp_content.data(1).ValueAttribute=Leq; %yep, it's here
                outstruct.objects.comp(i).comp_content.data(2).NameAttribute='Freq';
                outstruct.objects.comp(i).comp_content.data(2).ValueAttribute=freq; %yep, it's here
                outstruct.objects.comp(i).comp_content.data(3).NameAttribute='Grnd resis';
                outstruct.objects.comp(i).comp_content.data(3).ValueAttribute=rho; %yep, it's here
                %actual LCC data
                outstruct.objects.comp(i).LCC.NumPhasesAttribute=numphases;
                outstruct.objects.comp(i).LCC.IconLengthAttribute='true';
                outstruct.objects.comp(i).LCC.LineCablePipeAttribute=2; %probably 2=single core
                outstruct.objects.comp(i).LCC.ModelTypeAttribute=1; %0 = bergeron, 1 = pi
                outstruct.objects.comp(i).LCC.cable_header.InAirGrndAttribute=1;
                outstruct.objects.comp(i).LCC.cable_header.MatrixOutputAttribute='true';
                for k=1:numphases
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).NumCondAttribute=1;
                    if isnan(linecable_data(k,8))
                        rrout=linecable_data(k,5);
                    else
                        rrout=linecable_data(k,8);
                    end
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).RoutAttribute=rrout;
                    horzdist=linecable_data(k,2);
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).PosXAttribute=horzdist;
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).PosYAttribute=linecable_data(k,3);
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).conductor.RinAttribute=linecable_data(k,4);
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).conductor.RoutAttribute=linecable_data(k,5);
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).conductor.rhoAttribute=linecable_data(k,6);
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).conductor.muCAttribute=linecable_data(k,7);
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).conductor.muIAttribute=linecable_data(k,9);
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).conductor.epsIAttribute=linecable_data(k,10);
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).conductor.CextAttribute=0;
                    outstruct.objects.comp(i).LCC.cable_header.cable(k).conductor.GextAttribute=0;
                end

                % writes ATPDraw variables
                outstruct.variables.NumSimAttribute=1;
                outstruct.variables.IOPCVPAttribute=0;
                outstruct.variables.UseParserAttribute='false';
                % datastr={}; %does not seem to work when loading from XML
                % datastr{1,1}='RT=10.0 ';
                % datastr{end+1,1}='RC=1.0E-5 ';
                % outstruct.variables.VarStr = strjoin(datastr,char(10));

                %output to XML and fix non-conforming tags
                writestruct(outstruct, fname,'StructNodeName','project');
                if isfile(fname)
                    msg=sprintf('File %s created sucessfully!',['ATPDraw_' jobid '.xml']);
                    f = uimsgbox(app, msg, 'It works!','success');
                end
                % func_fix_text_labels(fname, fname)
                % func_fix_var_tags(fname, fname)
            end
        end