function [] = buildEMISimu(app,path)
addpath('emi_simu_funs')

% Read and pre-process geometry
jobid=app.JobIDEditField.Value;

freq=app.FrequencyHzEditField.Value;

LINE_FILE_NAME = app.SourceCSVEditField.Value;
src=readmatrix(LINE_FILE_NAME);

PIPE_FILE_NAME = app.TargetCSVEditField.Value;
tgt=readmatrix(PIPE_FILE_NAME);

rho=app.ResistivitymEditField.Value; %soil resistivity

minTargetLength=app.SubdivisionmEditField.Value; %pipeline will be subdivided every minTargetLength meters

[GeometricData, LCCData] = calcCouplingRegionsParams(src,tgt,minTargetLength); %don't try to understand what this is, it's just... WEIRD

src=GeometricData(:,1:2);
src=src(LCCData(:,4)==1,:);
tgt=GeometricData(:,3:4);

% Define ground potential coefficients
G=calcGreenFuns(src,tgt,rho);

% Cross-section
linecable_data = cell2matButNotStupid(app,app.UITable.Data);
tgt_idx=str2double(app.TargetphaseDropDown.Value);
linecable_data(tgt_idx,2)=nan;
sw_idx=str2double(app.ShieldwiresListBox.Value);

phaseNum=linecable_data(:,1);
uniquePhases = unique(phaseNum);
% Exclude indices given in tgt_idx and sw_idx
% exclude_idx = unique([tgt_idx; sw_idx]);
% activePhases = setdiff(uniquePhases, exclude_idx);
% Generate unique colors because it can be done
colors = lines(length(uniquePhases));
% Convert RGB to Hexadecimal
hexColors = rgb2hex(colors);
hexColors = strrep(hexColors,'#','');
% Convert to Windows color codes
winColors = hex2dec(hexColors);
winColors(sw_idx)=65280; %shield wires are green...
winColors(tgt_idx)=16711680; % the target is blue, this is dumb and so are you

% Initialize the cell array for bundled conductors
bundled_conductors = cell(length(uniquePhases), 1);
% Group conductors by phase number
for i = 1:length(uniquePhases)
    phase = uniquePhases(i);
    bundled_conductors{i} = find(phaseNum == phase);
end


% BUILD AN ATPDRAW FILE DESCRIBING THE EMI CIRCUIT
% <<<<<<<<<----------- don't change anything below this line

% don't bother, this is for debug purposes only
% sample= readstruct('cruz45_teste_len100_rho1000.xml', 'FileType','xml');


%output file name
fname=['EMISimu_' jobid '.xml'];

numblocks=size(LCCData,1);
numtowers=sum(LCCData(:,4));

% DECLARE SOME USEFUL CONSTANTS AND VARIABLES
nodeTGnd=sprintf('P%d',sw_idx(1));
nodePGnd=sprintf('P%d',tgt_idx);
DEFTSTART=-1;
DEFTSTOP=1000;

numphases=size(linecable_data,1);

TopLeftX=200; % 200 + 240 LCCs per row give the optimal layout
TopLeftY=0; %will need this later to place the LCCs

% writes ATPdraw header
outstruct.ApplicationAttribute='ATPDraw';
outstruct.VersionAttribute=7;
outstruct.VersionXMLAttribute=1;
outstruct.header.TimestepAttribute=1/app.SamplingrateHzEditField.Value;
outstruct.header.TmaxAttribute=app.SimulationtimesEditField.Value;
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

for i=1:numblocks
    deq=LCCData(i,1);
    Leq=LCCData(i,2);
    hastower=logical(LCCData(i,4));

    outstruct.objects.comp(i).NameAttribute='LCC';
    outstruct.objects.comp(i).IdAttribute=string(['LCCSEC' num2str(i)]);
    outstruct.objects.comp(i).CapanglAttribute=90;
    outstruct.objects.comp(i).CapPosXAttribute=-10;
    outstruct.objects.comp(i).CapPosYAttribute=-25;
    outstruct.objects.comp(i).comp_content.PosXAttribute=TopLeftX+(colcount*offsetX);
    outstruct.objects.comp(i).comp_content.PosYAttribute=TopLeftY+(rowcount*offsetY);
    colcount=colcount+1;
    thiscaption='';
    if hastower
        thiscaption=string(['T' num2str(tower)]);
        % tower resistors
        towRstruct.objects.comp(tower).NameAttribute='RESISTOR';
        towRstruct.objects.comp(tower).CaptionAttribute='';%sprintf('T%04d',tower);
        towRstruct.objects.comp(tower).IdAttribute=[];
        towRstruct.objects.comp(tower).LCC=[];
        towRstruct.objects.comp(tower).CapanglAttribute=0;
        towRstruct.objects.comp(tower).CapPosXAttribute=20;
        towRstruct.objects.comp(tower).CapPosYAttribute=5;
        towRstruct.objects.comp(tower).comp_content.AngleAttribute=270;
        towRstruct.objects.comp(tower).comp_content.PosXAttribute=TopLeftX+((colcount-1)*offsetX)+offsetX-20;
        towRstruct.objects.comp(tower).comp_content.PosYAttribute=TopLeftY+(rowcount*offsetY)+numphases*10+60;
        towRstruct.objects.comp(tower).comp_content.OutputAttribute=1; %1 for current
        towRstruct.objects.comp(tower).comp_content.IconAttribute='default';
        towRstruct.objects.comp(tower).comp_content.node(1).NameAttribute='From';
        towRstruct.objects.comp(tower).comp_content.node(1).ValueAttribute=sprintf('T%04d',tower);  %%%%%%%%%%%%%%% CHECKME FOR GENERALIZATION
        towRstruct.objects.comp(tower).comp_content.node(1).UserNamedAttribute=true;
        towRstruct.objects.comp(tower).comp_content.node(1).KindAttribute=1;
        towRstruct.objects.comp(tower).comp_content.node(1).PosXAttribute=-20;
        towRstruct.objects.comp(tower).comp_content.node(1).PosYAttribute=0;
        towRstruct.objects.comp(tower).comp_content.node(1).NamePosXAttribute=10;
        towRstruct.objects.comp(tower).comp_content.node(1).NamePosYAttribute=5;
        towRstruct.objects.comp(tower).comp_content.node(2).NameAttribute='To';
        towRstruct.objects.comp(tower).comp_content.node(2).ValueAttribute='';
        towRstruct.objects.comp(tower).comp_content.node(2).UserNamedAttribute=false;
        towRstruct.objects.comp(tower).comp_content.node(2).KindAttribute=1;
        towRstruct.objects.comp(tower).comp_content.node(2).PosXAttribute=20;
        towRstruct.objects.comp(tower).comp_content.node(2).PosYAttribute=0;
        towRstruct.objects.comp(tower).comp_content.node(2).NamePosXAttribute=10;
        towRstruct.objects.comp(tower).comp_content.node(2).NamePosYAttribute=5;
        towRstruct.objects.comp(tower).comp_content.node(2).GroundAttribute=1;
        towRstruct.objects.comp(tower).comp_content.data(1).NameAttribute='R';
        towRstruct.objects.comp(tower).comp_content.data(1).UnitAttribute='Ohm';
        towRstruct.objects.comp(tower).comp_content.data(1).ValueAttribute='RT';
        % current measuring probe
        towSWstruct.objects.comp(tower).NameAttribute='SWMEAS';
        towSWstruct.objects.comp(tower).CaptionAttribute='';%'→ TACS';
        towSWstruct.objects.comp(tower).IdAttribute=[];
        towSWstruct.objects.comp(tower).LCC=[];
        towSWstruct.objects.comp(tower).CapanglAttribute=0;
        towSWstruct.objects.comp(tower).CapPosXAttribute=10;
        towSWstruct.objects.comp(tower).CapPosYAttribute=5;
        towSWstruct.objects.comp(tower).comp_content.AngleAttribute=270;
        towSWstruct.objects.comp(tower).comp_content.PosXAttribute=TopLeftX+((colcount-1)*offsetX)+offsetX-20;
        towSWstruct.objects.comp(tower).comp_content.PosYAttribute=TopLeftY+(rowcount*offsetY)+numphases*10+20;
        towSWstruct.objects.comp(tower).comp_content.IconAttribute='default';
        towSWstruct.objects.comp(tower).comp_content.node(1).NameAttribute='SWF';
        towSWstruct.objects.comp(tower).comp_content.node(1).ValueAttribute=sprintf('%s%04d',nodeTGnd,i);  %%%%%%%%%%%%%%% CHECKME FOR GENERALIZATION
        towSWstruct.objects.comp(tower).comp_content.node(1).UserNamedAttribute=true;
        towSWstruct.objects.comp(tower).comp_content.node(1).KindAttribute=1;
        towSWstruct.objects.comp(tower).comp_content.node(1).PosXAttribute=-20;
        towSWstruct.objects.comp(tower).comp_content.node(1).PosYAttribute=0;
        towSWstruct.objects.comp(tower).comp_content.node(1).NamePosXAttribute=10;
        towSWstruct.objects.comp(tower).comp_content.node(1).NamePosYAttribute=5;
        towSWstruct.objects.comp(tower).comp_content.node(2).NameAttribute='SWT';
        towSWstruct.objects.comp(tower).comp_content.node(2).ValueAttribute=sprintf('T%04d',tower);
        towSWstruct.objects.comp(tower).comp_content.node(2).UserNamedAttribute=false;
        towSWstruct.objects.comp(tower).comp_content.node(2).KindAttribute=1;
        towSWstruct.objects.comp(tower).comp_content.node(2).PosXAttribute=20;
        towSWstruct.objects.comp(tower).comp_content.node(2).PosYAttribute=0;
        towSWstruct.objects.comp(tower).comp_content.node(2).NamePosXAttribute=10;
        towSWstruct.objects.comp(tower).comp_content.node(2).NamePosYAttribute=5;
        tower=tower+1;
    end
    outstruct.objects.comp(i).CaptionAttribute=thiscaption;


    % pipeline coating resistances
    pipRstruct.objects.comp(i).NameAttribute='RESISTOR';
    pipRstruct.objects.comp(i).CaptionAttribute='';
    pipRstruct.objects.comp(i).IdAttribute=[];
    pipRstruct.objects.comp(i).LCC=[];
    pipRstruct.objects.comp(i).CapanglAttribute=0;
    pipRstruct.objects.comp(i).CapPosXAttribute=20;
    pipRstruct.objects.comp(i).CapPosYAttribute=5;
    pipRstruct.objects.comp(i).comp_content.AngleAttribute=270;
    pipRstruct.objects.comp(i).comp_content.PosXAttribute=TopLeftX+((colcount-1)*offsetX)+20;
    pipRstruct.objects.comp(i).comp_content.PosYAttribute=TopLeftY+(rowcount*offsetY)+numphases*10+20;
    pipRstruct.objects.comp(i).comp_content.OutputAttribute=0;
    pipRstruct.objects.comp(i).comp_content.IconAttribute='default';
    pipRstruct.objects.comp(i).comp_content.node(1).NameAttribute='From';
    pipRstruct.objects.comp(i).comp_content.node(1).ValueAttribute=sprintf('%s%04d',nodePGnd,i);  %%%%%%%%%%%%%%% CHECKME FOR GENERALIZATION
    pipRstruct.objects.comp(i).comp_content.node(1).UserNamedAttribute=true;
    pipRstruct.objects.comp(i).comp_content.node(1).KindAttribute=1;
    pipRstruct.objects.comp(i).comp_content.node(1).PosXAttribute=-20;
    pipRstruct.objects.comp(i).comp_content.node(1).PosYAttribute=0;
    pipRstruct.objects.comp(i).comp_content.node(1).NamePosXAttribute=10;
    pipRstruct.objects.comp(i).comp_content.node(1).NamePosYAttribute=5;
    pipRstruct.objects.comp(i).comp_content.node(2).NameAttribute='To';
    pipRstruct.objects.comp(i).comp_content.node(2).ValueAttribute=sprintf('G%04d',i);
    pipRstruct.objects.comp(i).comp_content.node(2).UserNamedAttribute=false;
    pipRstruct.objects.comp(i).comp_content.node(2).KindAttribute=1;
    pipRstruct.objects.comp(i).comp_content.node(2).PosXAttribute=20;
    pipRstruct.objects.comp(i).comp_content.node(2).PosYAttribute=0;
    pipRstruct.objects.comp(i).comp_content.node(2).NamePosXAttribute=10;
    pipRstruct.objects.comp(i).comp_content.node(2).NamePosYAttribute=5;
    pipRstruct.objects.comp(i).comp_content.node(2).GroundAttribute=0;
    pipRstruct.objects.comp(i).comp_content.data.NameAttribute='R';
    pipRstruct.objects.comp(i).comp_content.data.UnitAttribute='Ohm';
    
    target_row = linecable_data(ismember(linecable_data(:, 1), tgt_idx),:);
    Rout=target_row(1,5);
    RhoCoat=app.CoatingresistivitymEditField.Value;
    if isnan(target_row(1,8))
        ThickCoat=0;
        RelPermCoat=1;
    else
        ThickCoat=target_row(1,8)-target_row(1,5);
        RelPermCoat=target_row(1,10);
    end
    

    ZC=shuntImpedanceSES(freq,Rout,RhoCoat,ThickCoat,RelPermCoat,Leq);
    pipRstruct.objects.comp(i).comp_content.data.ValueAttribute=real(ZC);

    % pipeline external GPR source
    pipGstruct.objects.comp(i).NameAttribute='TACSSOUR';
    pipGstruct.objects.comp(i).CaptionAttribute='';
    pipGstruct.objects.comp(i).IdAttribute=[];
    pipGstruct.objects.comp(i).LCC=[];
    pipGstruct.objects.comp(i).CapanglAttribute=0;

    pipGstruct.objects.comp(i).CapPosXAttribute=20;
    pipGstruct.objects.comp(i).CapPosYAttribute=5;
    pipGstruct.objects.comp(i).comp_content.AngleAttribute=90;
    pipGstruct.objects.comp(i).comp_content.PosXAttribute=TopLeftX+((colcount-1)*offsetX)+20;
    pipGstruct.objects.comp(i).comp_content.PosYAttribute=TopLeftY+(rowcount*offsetY)+numphases*10+60;
    pipGstruct.objects.comp(i).comp_content.IconAttribute='default';
    pipGstruct.objects.comp(i).comp_content.node(1).NameAttribute='TACS';
    pipGstruct.objects.comp(i).comp_content.node(1).ValueAttribute=sprintf('G%04d',i);  %%%%%%%%%%%%%%% CHECKME FOR GENERALIZATION
    pipGstruct.objects.comp(i).comp_content.node(1).UserNamedAttribute=true;
    pipGstruct.objects.comp(i).comp_content.node(1).KindAttribute=1;
    pipGstruct.objects.comp(i).comp_content.node(1).PosXAttribute=20;
    pipGstruct.objects.comp(i).comp_content.node(1).PosYAttribute=0;
    pipGstruct.objects.comp(i).comp_content.node(1).NamePosXAttribute=10;
    pipGstruct.objects.comp(i).comp_content.node(1).NamePosYAttribute=5;
    pipGstruct.objects.comp(i).comp_content.data(1).NameAttribute='U/I';
    pipGstruct.objects.comp(i).comp_content.data(1).UnitAttribute='';
    pipGstruct.objects.comp(i).comp_content.data(1).ValueAttribute=0;
    pipGstruct.objects.comp(i).comp_content.data(2).NameAttribute='TStart';
    pipGstruct.objects.comp(i).comp_content.data(2).UnitAttribute='s';
    pipGstruct.objects.comp(i).comp_content.data(2).ValueAttribute=-1;
    pipGstruct.objects.comp(i).comp_content.data(3).NameAttribute='TStop';
    pipGstruct.objects.comp(i).comp_content.data(3).UnitAttribute='s';
    pipGstruct.objects.comp(i).comp_content.data(3).ValueAttribute=DEFTSTOP;
    
    
    for j=1:numphases
        outstruct.objects.conn(conn_idx).conn_content.NumPhasesAttribute=1;
        outstruct.objects.conn(conn_idx).conn_content.Pos1XAttribute=TopLeftX+((colcount-1)*offsetX)+20;
        outstruct.objects.conn(conn_idx).conn_content.Pos1YAttribute=TopLeftY+(rowcount*offsetY)-20+j*10;
        outstruct.objects.conn(conn_idx).conn_content.Pos2XAttribute=TopLeftX+((colcount-1)*offsetX)+offsetX-20;
        outstruct.objects.conn(conn_idx).conn_content.Pos2YAttribute=TopLeftY+(rowcount*offsetY)-20+j*10;

        for jj=1:numel(bundled_conductors)
            % thiscolor='0'; %dont know what to do, paint it black
            if ismember(j,bundled_conductors{jj})
                thiscolor=winColors(jj);
            end
        end

        outstruct.objects.conn(conn_idx).conn_content.ColorAttribute=thiscolor;
        conn_idx=conn_idx+1;
    end

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
        outstruct.objects.comp(i).comp_content.node(k).ValueAttribute=string(sprintf('P%d%04d',k,i-1));
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
        outstruct.objects.comp(i).comp_content.node(k).ValueAttribute=string(sprintf('P%d%04d',w,i)); %string([linecable_data{w,12} sprintf('%04d',i+1)]);
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
            outstruct.objects.comp(i).LCC.cable_header.cable(k).RoutAttribute=linecable_data(k,5);
        else
            outstruct.objects.comp(i).LCC.cable_header.cable(k).RoutAttribute=linecable_data(k,8);
        end
        if isnan(linecable_data(k,2)) %modified this to autodetect what is the target conductor
            horzdist=deq;
        else
            horzdist=linecable_data(k,2);
        end
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
end

% LOCAL AND REMOTE TERMINATIONS
i=i+1;
if ~app.MatchwithcharacteristicimpedanceCheckBox.Value
    SendZVal=str2double(app.UITableTermImpedances.Data.SendZ);
    RecZVal=str2double(app.UITableTermImpedances.Data.RecZ);
else
    SendZVal=1e-5*ones(numphases,1);
    RecZVal=1e-5*ones(numphases,1);
end

for j = 1:numphases
    % Find the corresponding bundle
    for jj = 1:length(bundled_conductors)
        if ismember(j, bundled_conductors{jj})
            % Check if the current element is the smallest in the bundle
            bundle_indices=bundled_conductors{jj};
            if j == bundle_indices(1)
                outstruct.objects.comp(i).NameAttribute='RESISTOR';
                outstruct.objects.comp(i).CaptionAttribute='';
                outstruct.objects.comp(i).IdAttribute=[];
                outstruct.objects.comp(i).LCC=[];
                outstruct.objects.comp(i).CapanglAttribute=0;
                outstruct.objects.comp(i).CapPosXAttribute=20;
                outstruct.objects.comp(i).CapPosYAttribute=5;
                outstruct.objects.comp(i).comp_content.AngleAttribute=0;
                outstruct.objects.comp(i).comp_content.PosXAttribute=TopLeftX+1*offsetX-40;
                outstruct.objects.comp(i).comp_content.PosYAttribute=TopLeftY+(1*offsetY)+(j-1)*10-10;
                outstruct.objects.comp(i).comp_content.OutputAttribute=0;
                outstruct.objects.comp(i).comp_content.IconAttribute='default';
                outstruct.objects.comp(i).comp_content.node(1).NameAttribute='From';
                outstruct.objects.comp(i).comp_content.node(1).ValueAttribute=sprintf('%s','');  %%%%%%%%%%%%%%% CHECKME FOR GENERALIZATION
                outstruct.objects.comp(i).comp_content.node(1).UserNamedAttribute=true;
                outstruct.objects.comp(i).comp_content.node(1).KindAttribute=1;
                outstruct.objects.comp(i).comp_content.node(1).PosXAttribute=-20;
                outstruct.objects.comp(i).comp_content.node(1).PosYAttribute=0;
                outstruct.objects.comp(i).comp_content.node(1).NamePosXAttribute=10;
                outstruct.objects.comp(i).comp_content.node(1).NamePosYAttribute=5;
                outstruct.objects.comp(i).comp_content.node(2).NameAttribute='To';
                outstruct.objects.comp(i).comp_content.node(2).ValueAttribute=sprintf('P%d%04d',j,0);
                outstruct.objects.comp(i).comp_content.node(2).UserNamedAttribute=false;
                outstruct.objects.comp(i).comp_content.node(2).KindAttribute=1;
                outstruct.objects.comp(i).comp_content.node(2).PosXAttribute=20;
                outstruct.objects.comp(i).comp_content.node(2).PosYAttribute=0;
                outstruct.objects.comp(i).comp_content.node(2).NamePosXAttribute=10;
                outstruct.objects.comp(i).comp_content.node(2).NamePosYAttribute=5;
                outstruct.objects.comp(i).comp_content.node(2).GroundAttribute=0;
                outstruct.objects.comp(i).comp_content.data.NameAttribute='R';
                outstruct.objects.comp(i).comp_content.data.UnitAttribute='Ohm';
                outstruct.objects.comp(i).comp_content.data.ValueAttribute=SendZVal(jj);
                i=i+1;
                % else
                %     fprintf('This is the secondary conductor in a bundle\n');
            end
            break;
        end
    end
end

for j = 1:numphases
    % Find the corresponding bundle
    for jj = 1:length(bundled_conductors)
        if ismember(j, bundled_conductors{jj})
            % Check if the current element is the smallest in the bundle
            bundle_indices=bundled_conductors{jj};
            if j == bundle_indices(1)
                outstruct.objects.comp(i).NameAttribute='RESISTOR';
                outstruct.objects.comp(i).CaptionAttribute='';
                outstruct.objects.comp(i).IdAttribute=[];
                outstruct.objects.comp(i).LCC=[];
                outstruct.objects.comp(i).CapanglAttribute=0;
                outstruct.objects.comp(i).CapPosXAttribute=20;
                outstruct.objects.comp(i).CapPosYAttribute=5;
                outstruct.objects.comp(i).comp_content.AngleAttribute=0;
                outstruct.objects.comp(i).comp_content.PosXAttribute= TopLeftX+((colcount-1)*offsetX)+offsetX;
                outstruct.objects.comp(i).comp_content.PosYAttribute= TopLeftY+(rowcount*offsetY)-20+j*10;
                outstruct.objects.comp(i).comp_content.OutputAttribute=0;
                outstruct.objects.comp(i).comp_content.IconAttribute='default';
                outstruct.objects.comp(i).comp_content.node(1).NameAttribute='From';
                outstruct.objects.comp(i).comp_content.node(1).ValueAttribute=sprintf('P%d%04d',j,numblocks);
                outstruct.objects.comp(i).comp_content.node(1).UserNamedAttribute=true;
                outstruct.objects.comp(i).comp_content.node(1).KindAttribute=1;
                outstruct.objects.comp(i).comp_content.node(1).PosXAttribute=-20;
                outstruct.objects.comp(i).comp_content.node(1).PosYAttribute=0;
                outstruct.objects.comp(i).comp_content.node(1).NamePosXAttribute=10;
                outstruct.objects.comp(i).comp_content.node(1).NamePosYAttribute=5;
                outstruct.objects.comp(i).comp_content.node(2).NameAttribute='To';
                outstruct.objects.comp(i).comp_content.node(2).ValueAttribute=sprintf('%s','');  %%%%%%%%%%%%%%% CHECKME FOR GENERALIZATION
                outstruct.objects.comp(i).comp_content.node(2).UserNamedAttribute=false;
                outstruct.objects.comp(i).comp_content.node(2).KindAttribute=1;
                outstruct.objects.comp(i).comp_content.node(2).PosXAttribute=20;
                outstruct.objects.comp(i).comp_content.node(2).PosYAttribute=0;
                outstruct.objects.comp(i).comp_content.node(2).NamePosXAttribute=10;
                outstruct.objects.comp(i).comp_content.node(2).NamePosYAttribute=5;
                outstruct.objects.comp(i).comp_content.node(2).GroundAttribute=0;
                outstruct.objects.comp(i).comp_content.data.NameAttribute='R';
                outstruct.objects.comp(i).comp_content.data.UnitAttribute='Ohm';
                outstruct.objects.comp(i).comp_content.data.ValueAttribute=RecZVal(jj);
                i=i+1;
                % else
                %     fprintf('This is the secondary conductor in a bundle\n');
            end
            break;
        end
    end
end

% ADDITIONAL CARDS
% i=i+1; % already incremented above
crd=0;bspac=80;
datastr={};
outstruct.objects.comp(i).NameAttribute='ADDITIONAL';
outstruct.objects.comp(i).CaptionAttribute='TACS/MODELS SETUP';
outstruct.objects.comp(i).CapanglAttribute=0;
outstruct.objects.comp(i).CapPosXAttribute=-35;
outstruct.objects.comp(i).CapPosYAttribute=35;
outstruct.objects.comp(i).comp_content.PosXAttribute=TopLeftX+crd*2*bspac+(bspac+10);
outstruct.objects.comp(i).comp_content.PosYAttribute=TopLeftY+2*bspac;
outstruct.objects.comp(i).comp_content.IconAttribute='default';
outstruct.objects.comp(i).comp_content.data.NameAttribute='Kind';
outstruct.objects.comp(i).comp_content.data.UnitAttribute  ='';
outstruct.objects.comp(i).comp_content.data.ValueAttribute  = 1; %1 = TACS
datastr{1,1}=        'C TACS CURRENT INPUT TO MODELS';
for t=1:numtowers
    %'91XXXXX                                                            -1.      1.E3';
    datastr{end+1,1}=sprintf('91T%04d                                                            -1.      1.E3',t);
end
datastr{end+1,1}='C MODELS OUTPUT INTO TACS SOURCE';
for t=1:numtowers
    %'27XXXXX';
    datastr{end+1,1}=sprintf('27G%04d',t);
end
outstruct.objects.comp(i).comp_content.datastring = strjoin(datastr,char(10));

i=i+1;crd=crd+1;
datastr={};
outstruct.objects.comp(i).NameAttribute='ADDITIONAL';
outstruct.objects.comp(i).CaptionAttribute='BUNDLED PHASES';
outstruct.objects.comp(i).CapanglAttribute=0;
outstruct.objects.comp(i).CapPosXAttribute=-35;
outstruct.objects.comp(i).CapPosYAttribute=35;
outstruct.objects.comp(i).comp_content.PosXAttribute=TopLeftX+crd*2*bspac+(bspac+10);
outstruct.objects.comp(i).comp_content.PosYAttribute=TopLeftY+2*bspac;
outstruct.objects.comp(i).comp_content.IconAttribute='default';
outstruct.objects.comp(i).comp_content.data.NameAttribute='Kind';
outstruct.objects.comp(i).comp_content.data.UnitAttribute  ='';
outstruct.objects.comp(i).comp_content.data.ValueAttribute  = 3; %3 = BRANCH
datastr{1,1}=    'C SHORT-CIRCUIT BUNDLED PHASES';
datastr{end+1,1}='C < n1 >< n2 ><ref1><ref2>< R  >< L  >< C  >';
for k=0:1:numblocks       %'  50005A50005B              1E-6'
    for j=1:length(sw_idx)-1
        datastr{end+1,1}=sprintf('  P%d%04dP%d%04d            1.0E-5           ',sw_idx(j),k,sw_idx(j+1),k); %fucking lazy half assed solution, needs to improve for generality
    end

    for j=1:numel(bundled_conductors)
        thisbundle=bundled_conductors{j};
        if length(thisbundle)>1
            for jj=1:length(thisbundle)-1
                datastr{end+1,1}=sprintf('  P%d%04dP%d%04d            1.0E-5           ',thisbundle(jj),k,thisbundle(jj+1),k); %fucking lazy half assed solution, needs to improve for generality
            end
        end
    end
end
outstruct.objects.comp(i).comp_content.datastring = strjoin(datastr,char(10));

i=i+1;crd=crd+1;
datastr={};
outstruct.objects.comp(i).NameAttribute='ADDITIONAL';
outstruct.objects.comp(i).CaptionAttribute='REQUEST OUTPUTS';
outstruct.objects.comp(i).CapanglAttribute=0;
outstruct.objects.comp(i).CapPosXAttribute=-35;
outstruct.objects.comp(i).CapPosYAttribute=35;
outstruct.objects.comp(i).comp_content.PosXAttribute=TopLeftX+crd*2*bspac+(bspac+10);
outstruct.objects.comp(i).comp_content.PosYAttribute=TopLeftY+2*bspac;
outstruct.objects.comp(i).comp_content.IconAttribute='default';
outstruct.objects.comp(i).comp_content.data.NameAttribute='Kind';
outstruct.objects.comp(i).comp_content.data.UnitAttribute  ='';
outstruct.objects.comp(i).comp_content.data.ValueAttribute  = 8; %8 = REQUEST
datastr{1,1}=    'C OUTPUT NODE VOLTAGES';
for p=1:numblocks
    datastr{end+1,1}=sprintf('00%s%04d',nodePGnd,p);
    datastr{end+1,1}=sprintf('00G%04d',p); %dont ever again forget the trailing zeroes you fucking potato
end
outstruct.objects.comp(i).comp_content.datastring = strjoin(datastr,char(10));


% NOW IT'S TIME TO FACE THE MIGHTY BEAST - THE GPR MODELS OF THE DAMNED
nodeT=sprintfc('T%04d',[1:numtowers]);
nodeP=sprintfc('G%04d',[1:numblocks]);
MODstr=writeGPRsourcesModels(G,nodeT,nodeP,numtowers,numblocks);
i=i+1;crd=crd+1;
outstruct.objects.comp(i).NameAttribute='ADDITIONAL';
outstruct.objects.comp(i).CaptionAttribute='GROUND POTENTIAL RISE';
outstruct.objects.comp(i).CapanglAttribute=0;
outstruct.objects.comp(i).CapPosXAttribute=-35;
outstruct.objects.comp(i).CapPosYAttribute=35;
outstruct.objects.comp(i).comp_content.PosXAttribute=TopLeftX+crd*2*bspac+(bspac+10);
outstruct.objects.comp(i).comp_content.PosYAttribute=TopLeftY+2*bspac;
outstruct.objects.comp(i).comp_content.IconAttribute='default';
outstruct.objects.comp(i).comp_content.data.NameAttribute='Kind';
outstruct.objects.comp(i).comp_content.data.UnitAttribute  ='';
outstruct.objects.comp(i).comp_content.data.ValueAttribute  = 2; %2 = MODELS
outstruct.objects.comp(i).comp_content.datastring = MODstr;



% adds a sample text - TEXTS MUST COME AFTER THE COMPONENTS SECTION. You
% have been warned.

% i=1;
% outstruct.objects.text(i).text_content.PosXAttribute=TopLeftX;
% outstruct.objects.text(i).text_content.PosYAttribute=TopLeftY;
% outstruct.objects.text(i).text_content.FontNameAttribute='Microsoft Sans Serif';
% outstruct.objects.text(i).text_content.FontSizeAttribute=12;
% outstruct.objects.text(i).text_content.FontStyleAttribute=2;
% outstruct.objects.text(i).text_content.ColorAttribute=0;
% outstruct.objects.text(i).text_content.OrientationAttribute=0;
% outstruct.objects.text(i).text_content.BckColAttribute=16777215;
% outstruct.objects.text(i).text_content.FrmColAttribute=16777215;
% % outstruct.objects.text(i).text_content.Text=['Example label. Enforce a carriage return character to ↵ break line, or maybe a char(10)' char(10) 'like this.'];
% outstruct.objects.text(i).text_content.Text=['(TopLeftX,TopLeftY)'];
% i=i+1;
% outstruct.objects.text(i).text_content.PosXAttribute=offsetX;
% outstruct.objects.text(i).text_content.PosYAttribute=offsetY;
% outstruct.objects.text(i).text_content.FontNameAttribute='Microsoft Sans Serif';
% outstruct.objects.text(i).text_content.FontSizeAttribute=12;
% outstruct.objects.text(i).text_content.FontStyleAttribute=2;
% outstruct.objects.text(i).text_content.ColorAttribute=0;
% outstruct.objects.text(i).text_content.OrientationAttribute=0;
% outstruct.objects.text(i).text_content.BckColAttribute=16777215;
% outstruct.objects.text(i).text_content.FrmColAttribute=16777215;
% % outstruct.objects.text(i).text_content.Text=['Example label. Enforce a carriage return character to ↵ break line, or maybe a char(10)' char(10) 'like this.'];
% outstruct.objects.text(i).text_content.Text=['(offsetX,offsetY)'];

% writes ATPDraw variables
outstruct.variables.NumSimAttribute=1;
outstruct.variables.IOPCVPAttribute=0;
outstruct.variables.UseParserAttribute='false';
datastr={}; %does not seem to work when loading from XML
datastr{1,1}=sprintf('RT=%1.2f ',app.TowerresistanceEditField.Value);
% datastr{end+1,1}='RC=1.0E-5 ';
outstruct.variables.VarStr = strjoin(datastr,char(10));


%MERGE ALL COMPONENT STRUCTS
outstruct.objects.comp = [outstruct.objects.comp towRstruct.objects.comp towSWstruct.objects.comp ...
    pipRstruct.objects.comp pipGstruct.objects.comp];

%output to XML and fix non-conforming tags
fullfname=fullfile(path,fname);

writestruct(outstruct, fullfname,'StructNodeName','project');
func_fix_text_labels(fullfname, fullfname)
func_fix_var_tags(fullfname, fullfname)
func_replace_string(fullfname, fullfname, '1e-05', '1.0E-5')

if isfile(fullfname)
    msg=sprintf('File %s created sucessfully!',fname);
    f = uimsgbox(app, msg, 'It works!','success');
end

end