% for st=2:1:678
%     if abs(velocity_matlab7(st))>abs(velocity_matlab8(st))
%         temp=velocity_matlab7(st);
%        velocity_matlab7(st)=velocity_matlab8(st);
%         velocity_matlab8(st)=temp;
%     end
% end

% for st=518:1:678
%     velocity_matlab5(st)= velocity_matlab5(558);
%     velocity_matlab6(st)= velocity_matlab6(558);
%     velocity_matlab7(st)= velocity_matlab7(558);
%     velocity_matlab8(st)= velocity_matlab8(558);
% end

temp=Ti(30:677,2);
Ti(30:677,2)=Ti(30:677,1);
Ti(30:677,1)=temp;

temp=Ti(30:677,4);
Ti(30:677,4)=Ti(30:677,3);
Ti(30:677,3)=temp;