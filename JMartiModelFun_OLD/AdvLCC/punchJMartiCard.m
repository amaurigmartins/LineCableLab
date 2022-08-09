function [out] = punchJMartiCard(NP, polesZc, resZc, polesA1, resA1, Ti)
% ESCOPO:
%lembrar que cada instrução BRANCH começa com -N_modo

%polesZc = {NORD InfVal [p_1, p_2, p_3, ..., p_NORD] } modo 1
%          {NORD InfVal [p_1, p_2, p_3, ..., p_NORD] } modo 2 ...
%polesZc tem NP linhas x 3 colunas - terceira coluna é um vetor [1 x NORD]

%resZc = {NORD [r_1, r_2, r_3, ..., r_NORD] } modo 1
%        {NORD [r_1, r_2, r_3, ..., r_NORD] } modo 2 ...
%resZc tem NP linhas x 3 colunas - terceira coluna é um vetor [1 x NORD]

% repita para polesA1 e resA1

% quando printar os polos, checar se o contador é múltiplo de 3 e se for
% quebrar a linha

% Ti é uma matriz 2*NP x 2*NP
% linhas impares contém parte real, linhas pares -> parte imaginaria
% termos da matriz com 12 posições
% real_modo1 real_modo2 real_modo3 -> modos evoluem ao longo das colunas
% im_modo1 im_modo2 im_modo3 -> modos evoluem ao longo das colunas
% no maximo 6 por linha

% consultar pag321 do PDF do rule book raiz

% a principio, escrever os numeros com 25 posições

% escrever out em uma string


SKIP = 2;
PDT0 = '  0.00';
inPref = 'IN_';
inSuf = 'OUT';

seq = 'ABCDEFGHIJKLMNOPQRSTUWVXYZ';

out_params = [];
out_matrixTi = [];

Real_Ti = real(Ti);
Imag_Ti = imag(Ti);

for i = 1: NP
    
    node1 = sprintf('%s__%s',inPref,seq(i));
    node2 = sprintf('%s__%s',inSuf,seq(i));
    
    %% BRANCH Line
    
    out_branch=blanks(80);
    out_branch(1:2) = sprintf('-%d',i);
    out_branch(3:8)= node1;
    out_branch(9:14)= node2;
    out_branch(27:32) = sprintf('%5d.',SKIP);
    out_branch(33:38) = sprintf('%s',PDT0);
    out_branch(53:54) = '-2';
    out_branch(55:56) = sprintf('%2d',NP);
    
    %% CHARACTERISTIC IMPEDANCE ZC
    
    NORD_PolesZc = polesZc{i,1};
    InfVal_Zc = polesZc{i,2};
    
    out_infValZc = blanks(80);
    out_infValZc(5:8) = sprintf('%4d',NORD_PolesZc);
    out_infValZc(11:40) = modprin(30,InfVal_Zc);
    
    % ZC's POLES:
    
    out_PolesZc = [];
    
    for j = 1:NORD_PolesZc       
        out1 = modprin(23,polesZc{i,3}(1,j));
        out_PolesZc = [out_PolesZc '   ' out1];
        
        if divisible(j,3) && j ~= NORD_PolesZc
            out_PolesZc = strcat(out_PolesZc,'\n');
        end 
    end
    
    % ZC's RESIDUOS:
    
    NORD_ResZc = resZc{i,1};
    
    out_ResZc =[];
    
    for j = 1: NORD_ResZc
        out1 = modprin(23,resZc{i,2}(1,j));
        out_ResZc = [out_ResZc '   ' out1];
        
        if divisible(j,3) && j ~= NORD_ResZc
            out_ResZc = strcat(out_ResZc,'\n');
        end 
    end
    
    % TOTAL ZC's PARAMS:
    
    out_Zc = strcat(out_PolesZc,'\n',out_ResZc);
    
    %% PROPAGATION CONSTANT A1
    
    NORD_PolesA1 = polesA1{i,1};
    InfVal_A1 = polesA1{i,2};
    
    out_infValA1 = blanks(80);
    out_infValA1(5:8) = sprintf('%4d',NORD_PolesA1);
    out_infValA1(11:40) = modprin(30,InfVal_A1);
    
    % A1's POLES:
    
    out_PolesA1 = [];
    
    for j = 1:NORD_PolesA1       
        out1 = modprin(23,polesA1{i,3}(1,j));
        out_PolesA1 = [out_PolesA1 '   ' out1];
        
        if divisible(j,3) && j ~= NORD_PolesA1
            out_PolesA1 = strcat(out_PolesA1,'\n');
        end 
    end
    
    % A1's RESIDUOS:
    
    NORD_ResA1 = resA1{i,1};
    
    out_ResA1 =[];
    
    for j = 1: NORD_ResA1
        out1 = modprin(23,resA1{i,2}(1,j));
        out_ResA1 = [out_ResA1 '   ' out1];
        
        if divisible(j,3) && j ~= NORD_ResA1
            out_ResA1 = strcat(out_ResA1,'\n');
        end 
    end
    
    % TOTAL A1's PARAMS:
    
    out_A1 = strcat(out_PolesA1,'\n',out_ResA1);
    
    out_params = strcat(out_params,out_branch,'\n',out_infValZc,'\n',out_Zc,'\n',out_infValA1,'\n',out_A1,'\n');
    
    %% Ti MATRIX
    
    out_RealTi = [];
    out_ImagTi = [];
    
    for j = 1: NP
        out_RealTi = [out_RealTi '  ' modprin(10,Real_Ti(i,j))];
        out_ImagTi = [out_ImagTi '  ' modprin(10,Imag_Ti(i,j))];
        
        if divisible(j,6) && j ~= NP
            out_RealTi = strcat(out_RealTi,'\n');
            out_ImagTi = strcat(out_ImagTi,'\n');
        end
    end
    
    out_matrixTi = strcat(out_matrixTi,out_RealTi,'\n',out_ImagTi,'\n');
    
end

out = strcat(out_params,out_matrixTi);

end



%% AUXILIARY FUNCTIONS:

function [out] = modprin(len,num)

y=prin(['%' sprintf('%d',len) 'V'], num);

k = strfind(y,'e');

if ~isempty(k)
    
    y=prin(['%' sprintf('%d',len) 'v'], num);
    k = strfind(y,'e');
    substr=y(1:k-1);
    kk=strfind(substr,'.');
    
    if ~isempty(kk)
        num_after_dot=length(substr(kk+1:end));
        pow=y(k+1:end);
        newpow=num2str(str2num(pow)-num_after_dot + 1);
        substr=erase(substr,'.');
        substr = substr(1:end-1);
        abspow = abs(str2num(pow));
        absnewpow = abs(str2num(newpow));
        if (abspow < 10)
            if(absnewpow > 9)
                substr = substr(1:end-1);
                newpow=num2str(str2num(pow)-num_after_dot +2);
            end
        end
    end  
        y = [substr '.e' newpow];
end

if str2num(y) > 0 && str2num(y) < 1 && strfind(y,'.') == 1
    y = ['0' y(1:end-1)];
elseif str2num(y) < 0 && str2num(y) > -1 && strfind(y,'.') == 2
    y = ['-0' y(2:end-1)];
end

if len ~=length(y)
 y = [' ' y];
end



out=y;

end

function s = prin(fmt,varargin)
%
% prin.m:  An alternative to sprintf() & fprintf() - version 16Feb17
%          Calls Pftoa.m
%
% Calling sequence: ------------------------------
% s = prin(FormatString,OptionalArguments)
%               or
% s = prin(FileID,FormatString,OptionalArguments)
%
% FormatString '%nV' uses exactly n characters
% FormatString '%nW' uses at most n characters
% Formats v,w are similar except that decimal points are not counted
%
% For a more complete description, type "prin" at the command prompt.
%     (This will open prin.pdf)
%
% Author:  Paul Mennen (paul@mennen.org)
%          Copyright (c) 2017, Paul Mennen

s = '';
if ~nargin % open the documentation file if called without arguments
  if exist('s') open(feval('which','prin.pdf')); return; end;
end;
if iscell(fmt) & nargin==1                                  % special form used to set pretty print defaults
  if length(fmt) setappdata(0,'prin',fmt);                  % set new defaults
  else           rmappdata(0,'prin');                       % revert to defaults
  end;
  return;
end;
if isnumeric(fmt)                                           % here if the first argument is numeric?
  if nargin==1                                              % single numeric argument special case (array pretty print)
    f = '%+7W  ';  hdr = '%3d: ';  epr = 10;                % default format, header format, and number of entries per row
    v = getappdata(0,'prin');  lv = length(v);              % get f,epr,hdr from prin appdata if it exists
    if lv  f=v{1};  if lv>1 epr=v{2}; end;  if lv>2 hdr=v{3}; end;
           if isnumeric(f) f = sprintf('%%+%dW  ',f); end;
    end;
    f = sprintf('%s%d{%s}\\n',hdr,epr,f);                    % build prin format string
    v = fmt(:);  n = length(v);  pad = ceil(n/epr)*epr - n;  % number of pad elements to make a multiple of 10
    c = (n+pad)/epr;  head = epr*(1:c) - epr + 1;            % number of columns to print, header row
    v = [head; reshape([v; zeros(pad,1)],epr,c)];  v = v(:); % reshape input
    s = prin(f,v(1:end-pad));
    return;
  end;
  s = prin(varargin{:});                                    % come here if the 1st arguement a FileID or a file name
  if length(fmt)>1                                          % is the first argument a file name?
       if fmt(1)>0 p = 'at'; else p = 'wt'; fmt=-fmt; end;  % yes, come here (+/- = append/write)
       fmt = char(fmt);  nc = fmt(1)==' ';  if nc fmt(1)=''; end; % don't close if space in filename
       fmt = fopen(fmt,p); fwrite(fmt,s);
       if nc setappdata(0,'FIDp',fmt); else fclose(fmt); end;
  else v = version;                                           % no, come here if it's a FileID
       if v(1)=='6' & (fmt==1 | fmt==2)                       % fwrite to standard out doesn't 
            fprintf(fmt,strrep(strrep(s,'%','%%'),'\','\\')); % work in ver 6.1
       elseif fmt<1 f = getappdata(0,'FIDp');  fwrite(f,s);  if fmt<0 fclose(f); end;
       else fwrite(fmt,s);
       end;
  end;
  return;
end;
fmt = strrep(strrep(fmt,'\{','\173'),'\}','\175'); % replace \{ and \} with octal codes
p = find(fmt=='{');                     % search for last repeat block
while length(p)                         % keep going until all {} repeat blocks are expanded
  p = p(end);  v = p;                   % p points to left bracket, v points to repeat count
  q = find(fmt(p:end)=='}');
  if isempty(q) disp('unmatched { in format string'); return; end;
  q = q(1)+p;                           % point to character after matching end brace
  bb = fmt(p+1:q-2);                    % contents between the braces
  r = fmt(max(1,p-1))-48;               % get single digit repeat count
  if r>=0 & r<=9                        % valid repeat count.
       v = v-1;
       if p>2 r2 = fmt(p-2)-48;                                 % is this a 2 digit repeat count?
              if r2>=0 & r2<=9  v=v-1; r=r+10*r2; end;          % if yes, compute new repeat count
       end;
       bb = repmat(bb,1,r);                                     % replicate the repeat block
       e = find(bb=='!');                                       % find all the explanation points
       if length(e) bb([e e(end)+1:end]) = []; end;             % remove all ! & text after the last !
       fmt = [fmt(1:v-1) bb fmt(q:end)];                        % insert the replicated text
  else if length(find(bb=='%'))~=1 | r==46 | r==47              % no valid repeat count
            fmt = [fmt(1:p-1) '\173' bb '\175' fmt(q:end)];     % assume normal text (super/subscripts)
       else fmt = [fmt(1:p-1) 'LbRaCe' bb 'RbRaCe' fmt(q:end)]; % assume vector format (exactly one %)
       end;
  end;
  p = find(fmt=='{');                   % cycle thru all repeat blocks (back to front order)
end;  % end while length(p)
fmt = strrep(strrep(fmt,'LbRaCe','{'),'RbRaCe','}'); % change back to braces around the delimiters
p = find(fmt=='%');     % search for '%' format codes
q = find(diff(p)==1);   % search for '%%'
while length(q) p(q(1):q(1)+1) = []; q = find(diff(p)==1); end; % remove all '%%' from the list
n = length(p);          % number of format codes
nn = length(fmt);       % length of fmt argument
codes = 'vVwWscdiouxXfeEgG';
c = cell(1,n+1);        % will contain all the characters of fmt with format strings removed
f = cell(1,n);          % will contain all the format strings contain in fmt
ws = zeros(1,n);        % used to save the code index associated with f{n}
g=f; d=f; e=ws;         % f,d has format text surrounding bracketed format string f{k}. e has position of '!'
cb = 1;                 % point to beginning of next c{} string
for k=1:n               % extract the n format strings into f{1} to f{n}
  q = p(k);  w = '';    % q is the location of the % sign for this format string
  c{k} = fmt(cb:q-1);   % save format text between f{k-1} and f{k} (not including format strings)
  lb = find(c{k}=='{'); % pointer to left bracket
  if length(lb)
    lb = lb(1);
    g{k} = sprintf(c{k}(lb+1:end));  % remove first part of bracketed format from c and put it in g
    c{k} = c{k}(1:lb-1);
  end;
  for j=q+1:nn         % search for the format code (i.e. the end of the format string)
    w = find(codes==fmt(j));
    if length(w) break; end;
  end;
  if isempty(w) disp(sprintf('sprint: Unknown format code starting with %s',fmt(q:end))); return; end;
  ws(k) = w;               % save format code index
  f{k} = ['%' fmt(q+1:j)]; % save format string from the '%' to the format code
  cb = j+1;
  if length(lb)
    w = find(fmt(cb:end)=='}');
    if isempty(w) disp('unmatched { in format string'); return; end;
    w = w(1) + cb;                  % point to character after the matching right brace
    d{k} = sprintf(fmt(cb:w-2));    % the remaining portion of the bracket vecotr format
    ek = find(d{k}=='!');           % search for delimiter character
    if length(ek) d{k} = strrep(strrep(d{k},'!row','! ~, '),'!col','! ~; ');
                  ek=ek(1); e(k)=ek; d{k}(ek)=[];    % if found, record its position and remove it
    else          e(k) = length(d{k})+1;             % nonzero indicates a bracket vector format
    end;
    cb = w;
  end;
end;  % end for k=1:n
if cb<=nn c{n+1} = sprintf(fmt(cb:nn)); end;         % save format text that follows the last format string
for k=1:n c{k} = sprintf(c{k}); end;                 % sprintf conversions (e.g. \t to tab character)
s = c{1};  q = 1;                                    % q points to the next format string
na = length(varargin);  k = 0;                       % k points to the argument being processed
if na & ~n disp('prin() warning: No format codes. Variables not converted.'); na=0; end;
while k < na                                         % cycle thru all the arguments
  k = k+1;  arg = varargin{k};                       % get next argument
  while iscell(arg)
    varargin = [arg(:)' varargin(k+1:end)];          % if it's a cell array, append it to the
    k=1;  arg = varargin{1};  na = length(varargin); % front of the remaining arguments
  end;
  if e(q)>0
    arg = arg(:);  nb = length(arg);                 % bracket vector format comes here
    for j=1:nb                                       % loop thru each number in the vector
      ft = Pftoa(f{q},arg(j));                       % convert the next number in the vector to ascii
      if j<nb s = [s g{q} ft d{q}];                  % surround the number with left and right text (g,d)
      else    s = [s g{q} ft d{q}(1:e(q)-1) c{q+1}]; % for last vector element, don't include the delimiter
      end;
    end;
    q = q + 1;                                       % advance to next format string
    if q>n & k<na q=1; s=[s c{1}]; end;              % reuse the format string from the beginning
  else                                               % here for not a bracket vector format
    if ws(q)==5 s = [s sprintf(f{q},arg) c{q+1}];    % here for %s format
                q = q + 1;                           % advance to next format string
                if q>n & k<na q=1; s=[s c{1}]; end;  % reuse the format string
    else                                             % here for all other formats (except %s)
      arg = arg(:);  nb = length(arg);
      for j=1:nb                                     % cycle thru each element of next argument
        s = [s Pftoa(f{q},arg(j)) c{q+1}];           % conversion with Pftoa
        q = q + 1;                                   % advance to next format string
        if q>n & (j<nb | k<na) q=1; s=[s c{1}]; end; % reuse the format string from the beginning
      end;    % end for j=1:nb
    end;      % end if ws(q)==5
  end;        % end if e(k)>=0
end;          % end while k <= na
pr = findstr(s,' ~, '); pc = findstr(s,' ~; ');      % find the cell array row & column separators
n = length(pr) + length(pc);
if ~n return; end;                                   % return if cell seperators were not used
pp = sort([pr pc+.1 length(s)+1]);  p = round(pp);   % fractional part indicates column marker
row = 1;  col = 1;  c = s;  s = [];  p2 = -3;
for k = 1:n+1                                        % break the output string into a cell array
  p1 = p2;  p2 = p(k);                               % p1/p2 are previous and current separator locations
  s{row,col} = c(p1+4:p2-1);                         % exclude the separator characters
  if p2==pp(k) col=col+1; else col=1; row=row+1; end; % advance column or row as appropriate
end;
end
% end function prin

function s = Pftoa(fmtstr,val) % floating point to ascii convertion - called by prin.m
%
% Pftoa.m:  Alternative number conversion formatting - version 01Jan17
% Author:   Paul Mennen (paul@mennen.org)
%           Copyright (c) 2017, Paul Mennen
%
% function s = Pftoa(fmtstr,val)
% returns a string representing the number val using the format specified by fmtstr
%.
% fmtstr: format description string
% val:    the number to be converted to ascii
%
%
% fmtstr in the form '%nV' --------------------------------------------
% n: the field width
% s: the string representation of x with the maximum resolution possible
%    while using at exactly n characters.
%
% fmtstr in the form '%nv' ---------------------------------------------
% n: the field width, not counting decimal point (since it's so skinny)
% s: the string representation of x with the maximum resolution possible
%    while using at exactly n+1 characters. (If a decimal point is not
%    needed, then only n characters will be used).
%
% fmtstr in the form '%nW' ---------------------------------------------
% n: the field width
% s: the string representation of x with the maximum resolution possible
%    while using at most n characters.
%
% fmtstr in the form '%nw' ---------------------------------------------
% n: the field width, not counting decimal point (since it's so skinny)
% s: the string representation of x with the maximum resolution possible
%    while using at most n+1 characters.
%
% For any of the VWvw formats, if the field width is too small to allow
% even one significant digit, then '*' is returned.

% If the format code is not one of the four characters VWvw then use
% the sprintf c conventions: s=sprintf(fmtstr,number);
% e.g.  Pftoa('%7.2f',value) is identical to sprintf('%7.2f',value).

% Optional format modifiers are allowed between the % sign and the field width.
% An optional modifier is one of the characters "+-jJkL". The + and - modifiers
% control padding the output with blanks and the other four modifiers allow
% the conversion of complex numbers. These modifies are described fully in
% the prin.pdf help file.

% BACKGROUND: ---------------------------------------------------------------
% Pftoa() expands on the number conversion capabilities of sprintf's d,f,e,g
% conversion codes (which are identical to the c language conventions). These
% formats are ideal for many situations, however the fact that these formats
% will sometimes output more characters than the specified field width make
% them inappropriate when used to generate a number that is displayed in a
% fixed sized GUI element (such as an edit box) or in a table of numbers
% arranged in fixed width columns. This motivated the invention of Pftoa's new
% V and W formats. With the e & g formats one is often forced to specify a very
% small number of signifcant digits since otherwise on the possibly rare
% occations when the numbers are very big or very small an unintelligable
% display is produced in the GUI, or the generated table becomes hopelessly
% misallined. For example, suppose a column of numbers of width 8 characters
% normally contains numbers that look something like 1.234567 but could
% occationally contain a number such as 7.654321E-100. The best you could
      % do with a g format would be '%8.2g' which would produce the strings
% 1.2 and 7.6E-100 which means the numbers we see most often are truncated
% far more than necessary. Essentially with the e and g formats, you specify
% the precision you want and you accept whatever size string is produced.
% With the V and W formats, this is turned around. You specify the length
% of the string you want, and Pftoa supplies as much precision as possible
% with this constraint.
%
% For displaying columns of numbers (using a fixed spaced font) the V format
% is best since it always outputs a character string of the specified length.
% For example, the format string '%7V' will output seven characters. Never
% less and never more.
%
% For displaying a number in an edit box the W format is best. For example
% '%7W' will output at most seven characters, although it will output fewer
% than 7 characters if this does not reduce the precision of the output.
% For example, the integer "34" will be displayed with just two characters
% (instead of padding with blanks like the V format does). Since the text
% in an edit box is most often center aligned, this produces a more pleasing
% result. Using a lower case w (i.e. the '%7w' format) behaves similarly
% except that periods are not counted in the character count. This means
% that if a decimal point is needed to represent the number, 8 characters
% will be output and if a decimal point is not included in the representation
% then only 7 characters are output. This is most useful when using
% proportional width fonts. The period is not counted because the character
% width of the period is small compared with the '0-9' digits. Actually
% since most GUIs are rendered using proportially spaced fonts, the w format
% is used more often than the W format.
%
if nargin~=2 disp('Calling sequence: resultString = Pftoa(formatString,val)'); return; end;

fmtstr = deblank(fmtstr);                 % make sure format code is the last character
fcode = fmtstr(end);  fc = upper(fcode);  % extract format code. Convert to upper case
fw = fmtstr(2:end-1);                     % extract field width
pad = 0;                                  % no final padding (+1/-1 = pad on left/right)
if isempty(fw) fw = '7'; end;             % if field width is omited, use the default (7)
mf = fw(1);  sp = ' '; lmf = lower(mf);   % get possible modifier
if lmf=='k' mf=mf-1; lmf=lmf-1; sp=''; end; % k/K modifiers are more "Kompact" than j/J (no spaces)
if lmf=='j'                               % is there a complex modifier?
  fw(1) = '%';  fw = [fw fcode];          % yes, create the format string without the modifier
  ival = imag(val);  rval = real(val);
  if mf=='J' | (rval*ival)~=0
              if ival<0 pp = '-'; ival = -ival; else pp = '+'; end;
              if rval==0 rval=abs(rval); end;
              s = [Pftoa(fw,rval) sp pp sp Pftoa(fw,ival) 'i'];  % both real/imag parts
  elseif ival s = [Pftoa(fw,ival) 'i']; else s = Pftoa(fw,rval); % only need one part
  end;
  return;
end;
if fc~='W' & fc~='V'  s = sprintf(fmtstr,val); return; end; % use sprintf if format isn't v,V,w,W
uc = fc==fcode;                           % upper case code (i.e. V or W)
val = real(val);                          % ignore imaginary part
if mf=='+' pad=1; elseif mf=='-' pad=-1; end;
if pad fw = fw(2:end); if isempty(fw) fw = '7'; end; end;
w = sscanf(fw,'%d'); v = w;               % get field width
if ~w s = ''; return; end;                % zero field width returns an empty string
if fc=='V' s = [blanks(v-1) '*']; else s = '*'; end;  ss = s;  neg = [];
if     val==0     s = strrep(s,'*','0');

elseif isnan(val) s = [blanks(length(s)-3) 'NaN']; if v<3 s=s(1:v); end;
elseif isinf(val) if val>0 s = [blanks(length(s)-3) 'Inf'];  if v<3 s=s(1:v); end;
                  else     s = [blanks(length(s)-4) '-Inf']; if v<4 s=s(1:v); end;
                  end;
else neg = val<0;
end;
if isempty(neg)  % special cases (0,Inf,Nan) come here
  if pad
    if fc=='W'    p = v-length(s);  if p<1 return; end;
                  if pad>0 s = [blanks(p) s]; else s = [s blanks(p)]; end;
    elseif pad<0  while s(1)==' ' s = [s(2:end) ' ']; end;
    elseif val==0 s = ['0.' repmat('0',1,v-2)]; 
    end;
  end;
  return;
end;
q = [6 0 1 1; 5 1 1 2; 4 0 3 3; 0 0 0 0; 3 0 4 4; 4 1 2 3; 5 2 0 2;   % v,w formats
     7 1 0 0; 6 2 0 1; 5 1 2 2; 0 0 0 0; 4 1 3 3; 5 2 1 2; 6 3 -1 1]; % V,W formats
q = q(7*uc+min(find(abs(val) < [1e-99 1e-9 .01 10^(v-neg) 1e10 1e100 inf])),:);
fp = v - q(1) - neg;                         % compute fp, the format precision
if fp==-1 & uc fp=0; v=v+1; end;
if fp<0 return; end;                         % not enough digits available
if q(1) fmt = 'e'; else fmt = 'g'; end;      % select the e or g format
if ~fp  q = q + [0,1,-1,-1]; end;            % e format sometimes removes the "."
s = sprintf(sprintf('%%1.%d%c',fp,fmt),val); % convert to decimal string
n = length(s);                               % length of result
if n>3 & s(n-3)=='e'                         % is it a 2 digit exponent (for MAC)
  s = [s(1:n-2) '0' s(n-1:n)];               % change it to 3 digits
  n = n + 1;
end;
if q(1) q = [1:v-q(2) v+q(3):v+q(4)];        % here for e format
else                                         % here for g format
  fdot = findstr('.',s);
  if length(fdot)
    i = uc;  lz = 0;
    if fdot==2 & s(1)=='0' | length(findstr('-0.',s))
       i = i + 1;
       lz = length(findstr('0.0',s));
    end;
    if i s = sprintf(sprintf('%%1.%dg',fp-i),val); % use one or two fewer digits
         n = length(s);
    end;
    if lz s = strrep(s,'0.0','.0'); n=n-1; end;
  end;
  q = 1:min(~uc+v,n);
end; % end if q(1)
if max(q)>n s = ss; return; end;             % don't go over array bounds
s = s(q);  n = length(s);
if length([findstr(s,'0') findstr(s,'.') findstr(s,'-')]) == n % is there at least one nonzero digit?
  s = ss; return;                                              % return if not
end;
if fc=='V'
  p = w-length(s);                           % number of padding characters required
  isp = length(findstr('.',s));              % true if there is a period
  if ~uc p=p+isp; end;
  if p<=0 return; end;                       % no padding required
  if fmt=='e' s = [' ' s]; return; end;      % pad with blanks on left (p will be 1)
  if ~isp s=[s '.']; if uc p=p-1; end; end;  % if there is no period, add one before padding
  s = [s repmat('0',1,p)];                   % pad with zeros on the right
elseif pad & fc=='W'
  p = v-length(s);  if p<1 return; end;
  if pad>0 s = [blanks(p) s]; else s = [s blanks(p)]; end;
end;
end
% end function Pftoa

function [Res] = divisible(n,x)
if(rem(n,x)==0)
Res = 1;
else
Res = 0;
end
end

