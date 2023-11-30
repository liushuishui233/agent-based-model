save("C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal","sum_all","sum_elder","sum_job","inf_spacetype","intup_all","intlow_all","intup_job","intlow_job","intup_elder","intlow_elder","intup_spacetype","intlow_spacetype","Rt_mc","cumulative_job_vac","infectionspace_job","inf_pair")
%space_type
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',inf_spacetype,"space_type","B2:CM9");
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',intup_spacetype,"space_type","B13:CM20");
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',intlow_spacetype,"space_type","B24:CM31");

%sum_all
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',sum_all,"sum_all","B2:CM9");
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',intup_all,"sum_all","B13:CM20");
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',intlow_all,"sum_all","B24:CM31");

%sum_elder
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',sum_elder,"sum_elder","C2:CN25");
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',intup_elder,"sum_elder","C29:CN52");
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',intlow_elder,"sum_elder","C56:CN79");

%sum_job
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',sum_job,"sum_job","C2:CN33");
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',intup_job,"sum_job","C37:CN68");
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',intlow_job,"sum_job","C72:CN103");

%Rt_mc
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',Rt_mc,"rt");

%cumulative_job_vac
xlswrite('C:/Users/lxy/OneDrive/benyuan/sci writing/anew simulation/520/520_normal.xlsx',cumulative_job_vac,"cumulative job  vac","C2:E17");










% 响一声
 
%sound(sin(2*pi*25*(1:4000)/100));
 
% 鸟声
 
load chirp
sound(y,Fs);
 
% 锣声
 
%load gong
%sound(y,Fs);
 
% 哈里路亚
 
load handel
sound(y,Fs);
 
% 笑声
 
load laughter
sound(y,Fs);
 
% 啪哒声
 
load splat
sound(y,Fs);
 
% 火车
 
load train
sound(y,Fs);