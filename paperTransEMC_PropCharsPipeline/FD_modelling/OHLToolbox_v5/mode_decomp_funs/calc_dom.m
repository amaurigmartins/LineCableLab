function [f,vel_dom_OC,vel_dom_SC]=calc_dom(freq,length)

f=0:1:1000000;

vel_dom_OC_dis1=4*freq*length/1;
vel_dom_OC_dis2=4*freq*length/3;
vel_dom_OC_dis3=4*freq*length/5;

vel_dom_SC_dis1=4*freq*length/2;
vel_dom_SC_dis2=4*freq*length/4;
vel_dom_SC_dis3=4*freq*length/6;

vel_dom_OC_temp1=spline(freq,vel_dom_OC_dis1);
vel_dom_OC_temp2=spline(freq,vel_dom_OC_dis2);
vel_dom_OC_temp3=spline(freq,vel_dom_OC_dis3);

vel_dom_SC_temp1=spline(freq,vel_dom_SC_dis1);
vel_dom_SC_temp2=spline(freq,vel_dom_SC_dis2);
vel_dom_SC_temp3=spline(freq,vel_dom_SC_dis3);

vel_dom_OC(:,1)=ppval(vel_dom_OC_temp1,f);
vel_dom_OC(:,2)=ppval(vel_dom_OC_temp2,f);
vel_dom_OC(:,3)=ppval(vel_dom_OC_temp3,f);

vel_dom_SC(:,1)=ppval(vel_dom_SC_temp1,f);
vel_dom_SC(:,2)=ppval(vel_dom_SC_temp2,f);
vel_dom_SC(:,3)=ppval(vel_dom_SC_temp3,f);