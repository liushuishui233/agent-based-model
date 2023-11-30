function [col_up,col_low] = confint(col,MC_run,conf_int)
low_index=floor((1-conf_int)/2*MC_run+0.01+1);
up_index=MC_run-floor((1-conf_int)/2*MC_run+0.01+1)+1;
row=size(col,1);
column=size(col,2);
row_unit=floor(row/MC_run+0.01);
col_up=zeros(row_unit,column);
col_low=zeros(row_unit,column);
for j=1:column
    for k=1:row_unit
       A=zeros(MC_run,1);
       for i=1:MC_run
           A(i,1)=col(k+row_unit*(i-1),j);
       end
       B=sort(A);
       col_low(k,j)=B(low_index,1);
       col_up(k,j)=B(up_index,1);
    end
end

end

