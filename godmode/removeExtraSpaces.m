function [strOut]=removeExtraSpaces(strIn);
temp2=strIn;
temp1='';
%Replaces double spaces with single spaces until the string doesn't change for an iteration.
while ~strcmp(temp1,temp2)
  temp1=temp2;
  temp2=regexprep(temp1,'  ',' ');
end
strOut=temp2;
