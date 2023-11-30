function [sample] = alias_matrix(Wgt,N)
L=length(Wgt);
Swgt=sum(Wgt);
% initialize the two array
accept_prob=zeros(L,1); % save the probability
alias_index=zeros(L,1); % save the sequence number
small_queue=[]; % save index with area less than 1.0
large_queue=[]; % save index with area larger than 1.0
for i=1:L
    accept_prob(i,1)=L*Wgt(i)/Swgt;
    if accept_prob(i,1)<1
        small_queue=[small_queue;i];
    else
        large_queue=[large_queue;i];
    end
end
while not(isempty(small_queue)) && not(isempty(large_queue))
    small_index=small_queue(1,1);
    large_index=large_queue(1,1);
    small_queue(1,:)=[];
    large_queue(1,:)=[];
    alias_index(small_index,1)=large_index;
    accept_prob(large_index,1)=accept_prob(large_index,1)+accept_prob(small_index,1)-1;
    if accept_prob(large_index,1)<1
        small_queue=[small_queue;large_index];
    else
        large_queue=[large_queue;large_index];
    end
end
sample=zeros(N,1);
for i=1:N
random_num1=floor(rand()*L+1);
random_num2=rand();
if random_num2<accept_prob(random_num1,1)
    sample(i,1)=random_num1;
else
    sample(i,1)=alias_index(random_num1,1);
end
end
end

