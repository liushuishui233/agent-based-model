clear
clc
% Intervention (1): increasing of vaccine rate
% As we focus on the protection of elder people, we first set a background
% vaccine rate of 0,1,2,3 vaccine for employee 1, employee 2, elder and
% student, respectively
vac_ratio=[0.014,0.027,0.281,0.678;0.014,0.027,0.281,0.678;0.089,0.113,0.159,0.639;0.238,0.122,0.481,0.158];
hospital_ratio =[0.0057,0.0057,0.0057,0.0057;0.0057,0.0057,0.0057,0.0057;0.0377,0.0377,0.0377,0.0377;0.0101,0.0101,0.0101,0.0101]; % hospitalization rate for residents with various vaccine
critical_ratio=[0.1217,0.1217,0.1217, 0.1217;0.1217,0.1217,0.1217,0.1217;0.5276,0.5276,0.5276,0.5276;0.0288,0.0288,0.0288,0.0288]; % Ratio for residents transferring from hospitalization to critical states with various vaccine
death_ratio=[0.2,0.2,0.2,0.2;0.2,0.2,0.2,0.2;0.4,0.4,0.4,0.4;0.1,0.1,0.1,0.1]; % Ratio for residents transferring from critical states to recovery with various vaccine

% Intervention (2): reduction of going out and playing chess and cards for elder
% original frequency of activities
freq_act=[2/7,2/7,1/7;2/7,2/7,1/7;1/28,1/7,1/28;2/7,1/14,1/28];
% column no. for retail, restaurant, leisure in sequence. row no. for
% employee 1, employee2, elder, student
freq_elder_out=1; % frequency of going to retail buildings of the elder people as well as the playing chess

% Intervention (3): reverse protection for elder people with periodic
% antigen test of companion
antigen_test=0; % logical test for antigen test, 1 for execution, 0 for none

% Intervention (4) for validation: school closure
school_closure=0; % proportion of school closure, 0 for go to school, 1 for school closure

% Intervention (5) for validation: work from home
wfh=0; % proportion of working from home, 0 for go to office, 1 for work from home

% simulation duration of the pandemic
day_dua=90; % day, epidemic duration
MC_run=100; % runs of Monte Carlo simulation
conf_int=0.9; % -, confidence interval


load 287.mat
% Rd represents residence including center coordinate (km) and population
% Rr represents restaurant including center coordinat  e (km) and area (m2)
% Rt represents retail including center coordinate (km) and area (m2)
% Ls represents leisure including center coordinate (km) and area (m2)
% Fs represents family structure, 1st column family No., 2nd column job
% no., 3rd column family number, 4th column family type
% job no. 1 employee 1, 2 employee 2, 3 elder, 4 child
% Construct the motion trail for residents
n_emp=sum(Rd(:,3)); % -, total number of residents
Rv=rand(n_emp,1); % Immunity characteristics of every resident also the initial sample No.
hospitalization_index=rand(n_emp,1); % hospitalization characteristics of every resident
critical_index=rand(n_emp,1); % critical characteristics of every resident
death_index=rand(n_emp,1); % death characteristics of every resident
n_Rd=size(Rd,1); % number of residences
n_Rr=size(Rr,1); % number of restaurants
n_Rt=size(Rt,1); % number of retail buildings
n_Ls=size(Ls,1); % number of leisure places
n_Fs=size(Fs,1); % population number of unit family structure 
n_family=max(Fs(:,1)); % family number in family structure unit
num_family=floor(n_emp/n_Fs*n_family+0.1); % family number in the whole community

% construct characters of every resident with no., family no., family type,
% job, and compulsory activity location (office for employee, classroom for
% students, and chess room for elder), vaccine status
height=3;  % m, height of space
house_area=30; % m2/per, house area per person
ach_house=1; % h-1, air change rate for houses
office_scale=10; % people, average officer number in a single office reflecting the contact during office hour
office_area=30; % m2, office area
ven_office=30*4; % m3/h.per, required fresh air volume for each occupant in offices
class_scale=50; % people, average school classroom scale
class_area=70; % m2, classroom area
ven_class=30*4; % m3/h.per, required fresh air volume for each occupant in classrooms
ach_class=4.5; % h-1, air change rate of classrooms
chess_scale=100; % people, member of the chess room not the exact occupant number at real time
chess_area=60; % m2, chess room area
chess_proportion=0.2; % -, attendence rate for chess and card activity for elder people
ach_chess=0.5; % -, m3/h.per, required fresh air volume for each occupant in chess room
job_pro=[18;16;10;16]./n_Fs; % occupant number with different jobs in unit family structure
job_dis=job_pro.*n_emp; % job distribution among residents
num_occloc=[office_scale;office_scale;chess_scale;class_scale]; % average occupant number in each locations
num_comloc=ceil(job_dis./num_occloc); % location number for compulsory activities
num_chessloc=ceil(Rd(:,3).*job_pro(3,1)./chess_scale); % number of chess location in this residence
num_comloc(3,1)=sum(num_chessloc);
area_rtcri=500; % m2, critical area of retail building
area_rrcri=300; % m2, critical area of restaurants
ach_rt=1; % h-1, air change rate of retail buildings, 1 for small shops and 0.7 for large supermarket
ach_rr=1; % h-1, air change rate of restaurants, 1 for small shops and 0.7 for large restaurants
ach_ls=1; % h-1, air change rate of leisure places
% construct the area and air change rate for compulsory activity locations
info_space=zeros(num_family+sum(num_comloc)+n_Rt+n_Rr+n_Ls+1,2); % no. of spaces with volume (m3) and fresh air volume (m3/h)
info_space(num_family+1:num_family+sum(num_comloc),:)=[height*office_area*ones(num_comloc(1)+num_comloc(2),1),office_scale*ven_office*ones(num_comloc(1)+num_comloc(2),1);height*chess_area*ones(num_comloc(3),1),chess_area*height*ach_chess*ones(num_comloc(3),1);height*class_area*ones(num_comloc(4),1),class_scale*ven_class*ones(num_comloc(4),1)];
info_space(num_family+sum(num_comloc)+1:num_family+sum(num_comloc)+n_Rt+n_Rr+n_Ls,:)=[height*Rt(:,3),(0.7+(ach_rt-0.7)*(Rt(:,3)<area_rtcri))*height.*Rt(:,3);height*Rr(:,3),(0.7+(ach_rr-0.7)*(Rr(:,3)<area_rrcri))*height.*Rr(:,3);height*Ls(:,3),ach_ls*height*Ls(:,3)];
info_space(num_family+sum(num_comloc)+n_Rt+n_Rr+n_Ls+1,:)=[100000,100000]; % transportation
trans_no=num_family+sum(num_comloc)+n_Rt+n_Rr+n_Ls+1; % transportation number

% construct the weighted random sample for each residence
Wgt_Rr=zeros(n_Rr,n_Rd); % -, probability from i_th residence to j_th Restaurant
Wgt_Rt=zeros(n_Rt,n_Rd); % -, probability from i_th residence to j_th retail
Wgt_Ls=zeros(n_Ls,n_Rd); % -, probability from i_th residence to j_th leisure
for i=1:n_Rd
    for j=1:n_Rr
        dis_Rr=6371.004*acos(1-(power((sin((90-Rr(j,1))*pi/180)*cos(Rr(j,2)*pi/180)-sin((90-Rd(i,2))*pi/180)*cos(Rd(i,1)*pi/180)),2)+power((sin((90-Rr(j,2))*pi/180)*sin(Rr(j,1)*pi/180)-sin((90-Rd(i,2))*pi/180)*sin(Rd(i,1)*pi/180)),2)+power((cos((90-Rr(j,2))*pi/180)-cos((90-Rd(i,2))*pi/180)),2))/2);
        Wgt_Rr(j,i)=Rr(j,3)/(dis_Rr^2);
    end
    for j=1:n_Rt
        dis_Rt=6371.004*acos(1-(power((sin((90-Rt(j,1))*pi/180)*cos(Rt(j,2)*pi/180)-sin((90-Rd(i,2))*pi/180)*cos(Rd(i,1)*pi/180)),2)+power((sin((90-Rt(j,2))*pi/180)*sin(Rt(j,1)*pi/180)-sin((90-Rd(i,2))*pi/180)*sin(Rd(i,1)*pi/180)),2)+power((cos((90-Rt(j,2))*pi/180)-cos((90-Rd(i,2))*pi/180)),2))/2);
        Wgt_Rt(j,i)=Rt(j,3)/(dis_Rt^2);
    end
    for j=1:n_Ls
        dis_Ls=6371.004*acos(1-(power((sin((90-Ls(j,1))*pi/180)*cos(Ls(j,2)*pi/180)-sin((90-Rd(i,2))*pi/180)*cos(Rd(i,1)*pi/180)),2)+power((sin((90-Ls(j,2))*pi/180)*sin(Ls(j,1)*pi/180)-sin((90-Rd(i,2))*pi/180)*sin(Rd(i,1)*pi/180)),2)+power((cos((90-Ls(j,2))*pi/180)-cos((90-Rd(i,2))*pi/180)),2))/2);
        Wgt_Ls(j,i)=Ls(j,3)/(dis_Ls^2);
    end
end

% construct the random matrix for retail, restaurants and leisure places
N_Rs=10000; % length of the random array
RS_Rr=zeros(N_Rs,n_Rd); % random sample matrix for restaurants
RS_Rt=zeros(N_Rs,n_Rd); % random sample matrix for retail buildings
RS_Ls=zeros(N_Rs,n_Rd); % random sample matrix for leisure places
for i=1:n_Rd
    [RS_Rr(:,i)] = alias_matrix(Wgt_Rr(:,i),N_Rs);
    [RS_Rt(:,i)] = alias_matrix(Wgt_Rt(:,i),N_Rs);
    [RS_Ls(:,i)] = alias_matrix(Wgt_Ls(:,i),N_Rs);
end

% urban-community mobility
dt_int=0.1; % h, time interval for the time schedule
dt_trans=0.5; % h, transporation duration from one place to another
dt_ls=0.5; % h, dwell time in leisure places
% elder purchasing time 8am-10am or 4pm-5pm, and the elder purchasing every day
t_purch=[80:100,160:170]'; % time window for purchasement of the elder
% Infection transmission in communities
dt_ex_inf=1; % days, days from exposed to infected
dt_inf_qua=5; % days, days from infected to quarantined or hospitalized
dt_qua_cov=14; % days, days from quarantined to recovery or death
dt_qua_hos=3; % days, days from quarantined to hospitalization
dt_hos_cri=10; % days, days from hospitalization to recovery or critical state
dt_cri_cov=7; % days, days from critical state to recovery or death
% quanta/h could be variable in infection period
quan=[60;100;110;125;120]; % quanta/h, quanta emission rate from 1st to 5th day
k_vac=[1;1.179;1.333;7.8]; % anti-infection effect with vaccine, 0,1,2,3 vaccine
p_inhale=0.48; % m3/h, inhalation rate
antigen_period=3; % antigen test per three days for residents living with elder
antigen_sen=[0.884;0.95;0.96;0.90;0.75;0.63]; % antigen test sensitivity from 0th to 5th day

% For simplicity, the initial infector number is 1 for several communities,
% we could change the job of the infector

% Data collection
sum_all=zeros(8,day_dua); % temporal variation for all people with various states
sum_job=zeros(32,day_dua); % temporal variation for all people of different jobs with various states
sum_elder=zeros(24,day_dua); % temporal variation for elder people in differnt family types with various states
inf_space=zeros(trans_no,day_dua); % temporal variation of infector number in each space
inf_spacetype=zeros(8,day_dua); % temporal variation of infector number in different space types
% confidence interval
col_all=zeros(8*MC_run,day_dua); %temporal variation for all people with various status of all runs
col_job=zeros(32*MC_run, day_dua);
col_elder=zeros(24*MC_run,day_dua);
col_spacetype=zeros(8*MC_run,day_dua);
intup_all=zeros(8,day_dua); 
intlow_all=zeros(8,day_dua); 
intup_job=zeros(32,day_dua);
intlow_job=zeros(32,day_dua);
intup_elder=zeros(24,day_dua); 
intlow_elder=zeros(24,day_dua);
intup_spacetype=zeros(8,day_dua);
intlow_spacetype=zeros(8,day_dua);

inf_pair=zeros(16,16); % infector-infectee pairs divided by jobs and vaccines
infectionspace_job=zeros(32,day_dua); % temporal infection space divided by jobs
cumulative_job_vac=zeros(16,4); % row: 4 jobs x 4 vaccine; column: infection, hospitalization, ICU, death

% map the family type of elder people
map_family=[1;1;2;3;3;3;3];
% space type mapping: 1 family; 2 office; 3 chess; 4 classroom; 5 retail; 6
% restaurant; 7 Leisure; 8 transportation
map_space=[1*ones(num_family,1);2*ones(num_comloc(1,1)+num_comloc(2,1),1);3*ones(num_comloc(3,1),1);4*ones(num_comloc(4,1),1);5*ones(n_Rt,1);6*ones(n_Rr,1);7*ones(n_Ls,1);8];
Rt_mc=zeros(day_dua+1,1); % Rt in Monte Carlo simulation, including the initial date
Rt_deno=zeros(day_dua+1,1); % Monte Carlo runs Rt nonzero denominator
% disp_mc=zeros(1,4); % dispersion parameter of superspreading, k_value, p, p0, t10
% disp_num=0; % effective dispersion paramerter in Monte Carlo
% Monte Carlo run
for run=1:MC_run
    run
% elder people would play chess and card activity with neighbors
id_tot=zeros(n_emp,7); % characters of every residents
n_rest=0; % no. for the residents
Rt_cal=zeros(day_dua+1,3); % including the original date; infectee number at date t, newly infected by 
for i=1:n_Rd
    for j=1:Rd(i,3)
        n_rest=n_rest+1;
        n_familystruc=floor((n_rest-1)/n_Fs); % n_th family structure unit
        seq=mod(j-1,n_Fs)+1; % sequence of residents
        id_tot(n_rest,1)=Fs(seq,1)+n_family*n_familystruc; % family number
        id_tot(n_rest,2)=i; % Residence no.
        id_tot(n_rest,3)=Fs(seq,2); % job no.
        if Fs(seq,2)<2.5 
            % for employees
            if wfh<0.5
                % go to office
                id_tot(n_rest,4)=floor(rand*num_comloc(Fs(seq,2),1))+1+sum(num_comloc(1:Fs(seq,2),1))-num_comloc(Fs(seq,2),1)+num_family; % compulsory activity location no.
            else
                % work from home
                id_tot(n_rest,4)=Fs(seq,1)+n_family*n_familystruc; % compulsory activity location no.
            end
        elseif Fs(seq,2)<3.5 
            % for elder
            if wfh<0.5
                % they could play chess
                if rand<chess_proportion
                    id_tot(n_rest,4)=floor(rand*num_chessloc(i,1))+1+sum(num_chessloc(1:i,1))-num_chessloc(i,1)+sum(num_comloc(1:2,1))+num_family; % compulsory activity location no.
                else
                    id_tot(n_rest,4)=Fs(seq,1)+n_family*n_familystruc; % compulsory activity location no.
                end
            else
                % they are forbidden playing chess
                id_tot(n_rest,4)=Fs(seq,1)+n_family*n_familystruc; % compulsory activity location no.
            end
        else
            % for students
            if school_closure<0.5
                % go to school
                id_tot(n_rest,4)=floor(rand*num_comloc(Fs(seq,2),1))+1+sum(num_comloc(1:Fs(seq,2),1))-num_comloc(Fs(seq,2),1)+num_family; % compulsory activity location no.
            else
                % school closure
                id_tot(n_rest,4)=Fs(seq,1)+n_family*n_familystruc; % compulsory activity location no.
            end            
        end
        
        id_tot(n_rest,5)=Fs(seq,3); % family member number
        id_tot(n_rest,6)=Fs(seq,4); % family type
        info_space(Fs(seq,1)+n_family*n_familystruc,1)=Fs(seq,3)*house_area*height; % m3, space volume
        info_space(Fs(seq,1)+n_family*n_familystruc,2)=Fs(seq,3)*house_area*height*ach_house; % m3/h, space fresh air volume
        % vaccine index
        vaccine_index=rand;
        if vaccine_index<vac_ratio(Fs(seq,2),1)
            id_tot(n_rest,7)=0; % number of vaccine
        elseif vaccine_index<vac_ratio(Fs(seq,2),1)+vac_ratio(Fs(seq,2),2)
            id_tot(n_rest,7)=1; % number of vaccine
        elseif vaccine_index<vac_ratio(Fs(seq,2),1)+vac_ratio(Fs(seq,2),2)+vac_ratio(Fs(seq,2),3)
            id_tot(n_rest,7)=2; % number of vaccine
        else
            id_tot(n_rest,7)=3; % number of vaccine
        end
    end
end

% status of the individuals
status_occ=zeros(n_emp,2); % -, first column, 0: susceptible; 1: exposed; 2: 
% infected; 3: quarantine; 4: recovery; 5: hospitalization; 6: critical; 7: death; second column, days
ori_infector_no=floor(rand*n_emp)+1;
status_occ(ori_infector_no,1)=2;
Rout_0=id_tot(:,1); % original location no.
dose_0=zeros(n_emp,1); % dose exposure at yesterday 24:00
Infector_0=zeros(n_emp,1); % yesterday 24:00 infector no. 
Inf_info=zeros(n_emp,5); % -, status (0 for susceptible, 1 for infected), infector no., infected date of infectee, infected date of infector, newly infected number by the resident
Inf_info(ori_infector_no,1)=1;
for time=1:day_dua
    time
    highrisk_loc=zeros(num_family+sum(num_comloc)+n_Rt+n_Rr+n_Ls+1,1); % mark the compulsory activity locations with high risk
    Rout=zeros(floor(24/dt_int+0.1),n_emp); % time schedule for every person in a day
    infector_info=0;
% construct the routine of every individual
occ_no=0; % resident number
crowd_no=zeros(floor(24/dt_int+0.1),n_Rt+n_Rr+n_Ls); % temporal occupant number in locations
Num_infector=0; % present infector number
for i=1:n_emp
    if status_occ(i,1)>1.5 && status_occ(i,1)<2.5
        Num_infector=Num_infector+1;
    end
end
if Num_infector>0.5
for i=1:n_Rd
    for j=1:Rd(i,3)
        occ_no=occ_no+1;
        if status_occ(occ_no,1)<2.5
            house_no=id_tot(occ_no,1); % house number
            job=id_tot(occ_no,3);
            Rank_rd=floor(rand*N_Rs)+1; % select in the alias array
            nightlife_start=[175:195]'; % the starting moment of nightlife
            rout_no=nightlife_start(floor(rand*size(nightlife_start,1))+1,1);
            act_rand=rand(3,1); % random value for three activities
            act_array=[];
            if act_rand(1,1)<freq_act(job,1)
                dt_rt=0.3+0.1*(Rt(RS_Rt(Rank_rd,i),3)>area_rtcri); % h,dwell time in retail buildings, 0.1 h for small shops (<500 m2) and 0.5 h for large supermarkets (>500 m2)
                act_array=[act_array;[RS_Rt(Rank_rd,i)+num_family+sum(num_comloc),floor(dt_rt/dt_int+0.1)]];
            end
            if act_rand(2,1)<freq_act(job,2)
                dt_rr=0.4+0.6*(Rr(RS_Rr(Rank_rd,i),3)>area_rrcri); % h, dwell time in restaurants, 0.5 h for small restaurants (<300 m2) and 1 h for large restaurants (>300 m2)
                act_array=[act_array;[RS_Rr(Rank_rd,i)+num_family+sum(num_comloc)+n_Rt,floor(dt_rr/dt_int+0.1)]];
            end
            if act_rand(3,1)<freq_act(job,3)
                act_array=[act_array;[RS_Ls(Rank_rd,i)+num_family+sum(num_comloc)+n_Rt+n_Rr,floor(dt_ls/dt_int+0.1)]];
            end
            [act_arr] = KS(act_array);
            Rout(:,occ_no)=house_no*ones(floor(24/dt_int+0.1),1);
            if job<3.5 && job>2.5
                % for elder
                if rand<freq_elder_out
                Rout(125:150,occ_no)=id_tot(occ_no,4)*ones(26,1); % playing chess and card
                time_purch=floor(rand*size(t_purch,1))+1;
                rt_rank=floor(rand*N_Rs)+1; % select a retail building
                dt_rt=0.3+0.1*(Rt(RS_Rt(rt_rank,i),3)>area_rtcri); 
                Rout(t_purch(time_purch,1)-floor(dt_trans/dt_int+0.1):t_purch(time_purch,1)-1,occ_no)=trans_no*ones(floor(dt_trans/dt_int+0.1),1); % leave the house
                Rout(t_purch(time_purch,1):t_purch(time_purch,1)+floor(dt_rt/dt_int+0.1)-1,occ_no)=(RS_Rt(rt_rank,i)+num_family+sum(num_comloc))*ones(floor(dt_rt/dt_int+0.1),1); % arrive at a retail building
                crowd_no(t_purch(time_purch,1):t_purch(time_purch,1)+floor(dt_rt/dt_int+0.1)-1,RS_Rt(rt_rank,i))=crowd_no(t_purch(time_purch,1):t_purch(time_purch,1)+floor(dt_rt/dt_int+0.1)-1,RS_Rt(rt_rank,i))+ones(floor(dt_rt/dt_int+0.1),1);
                Rout(t_purch(time_purch,1)+floor(dt_rt/dt_int+0.1):t_purch(time_purch,1)+floor(dt_rt/dt_int+0.1)+floor(dt_trans/dt_int+0.1)-1,occ_no)=trans_no*ones(floor(dt_trans/dt_int+0.1),1); % arrive at the house
                end           
                if not(isempty(act_arr))
                    Rout(175:rout_no,occ_no)=trans_no*ones(rout_no-175+1,1); % leave the house
                end
            else
                % for employees and students
                Rout(95:170,occ_no)=id_tot(occ_no,4)*ones(170-95+1,1);
                Rout((95-floor(dt_trans/dt_int+0.1)):(95-1),occ_no)=trans_no*ones(floor(dt_trans/dt_int+0.1),1); % arrive at office or classroom from house
                Rout(171:rout_no,occ_no)=trans_no*ones(rout_no-171+1,1); % leave the house          
            end
            if not(isempty(act_arr))
                for k=1:size(act_arr,1)
                    Rout(rout_no+1:rout_no+act_arr(k,2),occ_no)=act_arr(k,1)*ones(act_arr(k,2),1);
                    crowd_no(rout_no+1:rout_no+act_arr(k,2),act_arr(k,1)-num_family-sum(num_comloc))=crowd_no(rout_no+1:rout_no+act_arr(k,2),act_arr(k,1)-num_family-sum(num_comloc))+ones(act_arr(k,2),1);
                    Rout(rout_no+act_arr(k,2)+1:rout_no+act_arr(k,2)+floor(dt_trans/dt_int+0.1),occ_no)=trans_no*ones(floor(dt_trans/dt_int+0.1),1);
                    rout_no=rout_no+act_arr(k,2)+1;
                end
            end
        end
     end    
end
end

% calculate the infection probability
status_occ(:,2)=status_occ(:,2)+ones(n_emp,1);
con=zeros(floor(24/dt_int+0.1),num_family+sum(num_comloc)+n_Rt+n_Rr+n_Ls+1); % quanta/m3, concentration in locations
infector_loc=zeros(floor(24/dt_int+0.1),num_family+sum(num_comloc)+n_Rt+n_Rr+n_Ls+1); % -, dominant infector no. 
quan_loc=zeros(floor(24/dt_int+0.1),num_family+sum(num_comloc)+n_Rt+n_Rr+n_Ls+1);  % quanta, quanta emission by infector in locations
if Num_infector>0.5
for i=1:n_emp
    if status_occ(i,1)>1.5 && status_occ(i,1)<2.5
        % for infectors
        for j=1:floor(24/dt_int+0.1)
            location=Rout(j,i);
            con(j,location)=con(j,location)+quan(status_occ(i,2),1)/info_space(location,2);
            if quan(status_occ(i,2),1)>quan_loc(j,location)
                quan_loc(j,location)=quan(status_occ(i,2),1);
                infector_loc(j,location)=i; % ith resident is the dominant infector in location at j moment
            else
            end
            con(j,trans_no)=0;
        end
    end
end
end

for i=1:n_emp
if Num_infector>0.5
    if status_occ(i,1)<0.5
        dose=0;
        for j=1:floor(24/dt_int+0.1)
            location=Rout(j,i);
            if j<1.5
                dose=dose_0(i,1);
                if abs(Rout(j,i)-Rout_0(i,1))>0.5
                    prob_inf=1-exp(-dose/k_vac(id_tot(i,7)+1,1));
                    dose=0;
                    if Rv(i,1)<prob_inf
                        status_occ(i,1)=1; % transfer from susceptible to exposed
                        status_occ(i,2)=0; % the first day of exposed
                        Inf_info(i,1)=1; % infector status index
                        Inf_info(i,2)=Infector_0(i); % infector no.
                        Inf_info(i,3)=time; % infected date of infectee
                        Inf_info(i,4)=Inf_info(Infector_0(i),3); % infected date of infector
                        inf_spacetype(map_space(Rout_0(i,1),1),time)=inf_spacetype(map_space(Rout_0(i,1),1),time)+1/MC_run;
                        infectionspace_job(map_space(Rout_0(i,1),1)+8*(id_tot(i,3)-1),time)=infectionspace_job(map_space(Rout_0(i,1),1)+8*(id_tot(i,3)-1),time)+1/MC_run;
                        col_spacetype(map_space(Rout_0(i,1),1)+8*(run-1),time)=col_spacetype(map_space(Rout_0(i,1),1)+8*(run-1),time)+1;
                        if Inf_info(i,2)<0.5
                            Infector_00=i;
                        else
                            Infector_00=Inf_info(i,2);
                        end
                        inf_pair(id_tot(i,7)+1+4*(id_tot(i,3)-1),id_tot(Infector_00,7)+1+4*(id_tot(Infector_00,3)-1))=inf_pair(id_tot(i,7)+1+4*(id_tot(i,3)-1),id_tot(Infector_00,7)+1+4*(id_tot(Infector_00,3)-1))+1/MC_run; % infector-infectee pair
                        
                    end
                end
            else
                if abs(Rout(j,i)-Rout(j-1,i))>0.5
                % resident transfers from one place to another
                    prob_inf=1-exp(-dose/k_vac(id_tot(i,7)+1,1));
                    dose=0;
                    if Rv(i,1)<prob_inf
                        status_occ(i,1)=1; % transfer from susceptible to exposed
                        status_occ(i,2)=0; % the first day of exposed
                        Inf_info(i,1)=1; % infector status index
                        Inf_info(i,3)=time; % infected date of infectee
                        if infector_info>0.5
                            Inf_info(i,2)=infector_info; % infector no.
                            Inf_info(i,4)=Inf_info(infector_info,3); % infected date of infector
                        else
                            Inf_info(i,2)=Infector_0(i); % infector no.
                            Inf_info(i,4)=Inf_info(Infector_0(i),3); % infected date of infector                            
                        end
                        inf_spacetype(map_space(Rout(j-1,i),1),time)=inf_spacetype(map_space(Rout(j-1,i),1),time)+1/MC_run;
                        infectionspace_job(map_space(Rout(j-1,i),1)+8*(id_tot(i,3)-1),time)=infectionspace_job(map_space(Rout(j-1,i),1)+8*(id_tot(i,3)-1),time)+1/MC_run;
                        col_spacetype(map_space(Rout(j-1,i),1)+8*(run-1),time)=col_spacetype(map_space(Rout(j-1,i),1)+8*(run-1),time)+1;
                        if Inf_info(i,2)<0.5
                            Infector_00=i;
                        else
                            Infector_00=Inf_info(i,2);
                        end
                        inf_pair(id_tot(i,7)+1+4*(id_tot(i,3)-1),id_tot(Infector_00,7)+1+4*(id_tot(Infector_00,3)-1))=inf_pair(id_tot(i,7)+1+4*(id_tot(i,3)-1),id_tot(Infector_00,7)+1+4*(id_tot(Infector_00,3)-1))+1/MC_run; % infector-infectee pair
                    end
                end
            end
            dose=dose+con(j,location)*p_inhale*dt_int;
            if con(j,location)>0
                infector_info=infector_loc(j,location);
            end
        end
        dose_0(i,1)=dose;
        Rout_0(i,1)=Rout(floor(24/dt_int+0.1),i);
        if infector_info>0.5
        Infector_0(i)=infector_info;
        end
        infector_info=0;
    end
end
    if (status_occ(i,1)-0.5)*(status_occ(i,1)-1.5)<0 && status_occ(i,2)>dt_ex_inf-0.1
            status_occ(i,1)=2; % transfer from exposed to infected
            status_occ(i,2)=0; % the first day of infector
            cumulative_job_vac(id_tot(i,7)+1+4*(id_tot(i,3)-1),1)=cumulative_job_vac(id_tot(i,7)+1+4*(id_tot(i,3)-1),1)+1/MC_run; % infection divided by jobs and vaccines
    elseif (status_occ(i,1)-1.5)*(status_occ(i,1)-2.5)<0 && status_occ(i,2)>dt_inf_qua-0.1
            status_occ(i,1)=3; % transfer from infected to quarantined
            status_occ(i,2)=0; % the first day of quarantined
            % intervention strategy of check the companion in compulsory
            % activities, i.e., house, office, chess, classroom
            highrisk_loc(id_tot(i,1),1)=1; % mark the house no. with high risk
            highrisk_loc(id_tot(i,4),1)=1; % mark the compulsory activity location no. with high risk
    else
        if (status_occ(i,1)-2.5)*(status_occ(i,1)-3.5)<0
            if (status_occ(i,2)-dt_qua_hos-0.1)*(status_occ(i,2)-dt_qua_hos+0.1)<0
                if hospitalization_index(i,1)<hospital_ratio(id_tot(i,3),id_tot(i,7)+1)
                status_occ(i,1)=5; % transfer from quarantined to hospitalization
                status_occ(i,2)=0; % the first day of hospitalization
                cumulative_job_vac(id_tot(i,7)+1+4*(id_tot(i,3)-1),2)=cumulative_job_vac(id_tot(i,7)+1+4*(id_tot(i,3)-1),2)+1/MC_run; % hospitalization divided by jobs and vaccines
                end
            elseif status_occ(i,2)>dt_qua_cov-0.1
                status_occ(i,1)=4; % transfer from quarantined to recovery
                status_occ(i,2)=0; % the first day of recovery
            end
        elseif (status_occ(i,1)-4.5)*(status_occ(i,1)-5.5)<0
            if status_occ(i,2)>dt_hos_cri-0.1
                if critical_index(i,1)<critical_ratio(id_tot(i,3),id_tot(i,7)+1)
                    status_occ(i,1)=6; % transfer from hospitalization to critical state
                    status_occ(i,2)=0; % the first day of critical state
                    cumulative_job_vac(id_tot(i,7)+1+4*(id_tot(i,3)-1),3)=cumulative_job_vac(id_tot(i,7)+1+4*(id_tot(i,3)-1),3)+1/MC_run; % ICU divided by jobs and vaccines

                else
                    status_occ(i,1)=4; % transfer from hospitalization to recovered
                    status_occ(i,2)=0; % the first day of recovery
                end
            end
        elseif (status_occ(i,1)-5.5)*(status_occ(i,1)-6.5)<0
            if status_occ(i,2)>dt_cri_cov-0.1
                if death_index(i,1)<death_ratio(id_tot(i,3),id_tot(i,7)+1)
                    status_occ(i,1)=7; % transfer from critical state to death
                    status_occ(i,2)=0; % the first day of death
                    cumulative_job_vac(id_tot(i,7)+1+4*(id_tot(i,3)-1),4)=cumulative_job_vac(id_tot(i,7)+1+4*(id_tot(i,3)-1),4)+1/MC_run; % Deaths divided by jobs and vaccines

                else
                    status_occ(i,1)=4; % transfer from critical state to recovered
                    status_occ(i,2)=0; % the first day of recovery
                end
            end            
        end
     end
end

for i=1:n_emp
% Intervention: PCR test for the susceptible close contacting with newly quarantined
% (confirmed) in compulsory activity locations
    if highrisk_loc(id_tot(i,1),1)>0.5 || highrisk_loc(id_tot(i,4),1)>0.5
        if (status_occ(i,1)- 0.5)*(status_occ(i,1)-2.5)<0 
            % i is an infector or exposed then he/she will be quarantined
            status_occ(i,1)=3; % transfer from infected to quarantined
            status_occ(i,2)=0; % the first day of quarantined            
        end
    end
    if antigen_test>0.5
    % Intervention: Periodic antigen test for companions of the elder
        if id_tot(i,6)>6.5 && mod(time,antigen_period)<0.5
            if status_occ(i,1)<2.5 && status_occ(i,1)>1.5
                if rand<antigen_sen(status_occ(i,2)+1,1)
                    status_occ(i,1)=3; % transfer from infected to quarantined
                    status_occ(i,2)=0; % the first day of quarantined     
                end
            end
        end
    end
end
for i=1:n_emp
    sum_all(status_occ(i,1)+1,time)=sum_all(status_occ(i,1)+1,time)+1/MC_run;
    sum_job(status_occ(i,1)+1+8*(id_tot(i,3)-1),time)=sum_job(status_occ(i,1)+1+8*(id_tot(i,3)-1),time)+1/MC_run;
    col_all(status_occ(i,1)+1+8*(run-1),time)=col_all(status_occ(i,1)+1+8*(run-1),time)+1;
    col_job(status_occ(i,1)+1+8*(id_tot(i,3)-1)+32*(run-1),time)=col_job(status_occ(i,1)+1+8*(id_tot(i,3)-1)+32*(run-1),time)+1;
    if (id_tot(i,3)-2.5)*(id_tot(i,3)-3.5)<0
        % elder people
        sum_elder(status_occ(i,1)+1+8*(map_family(id_tot(i,6),1)-1),time)=sum_elder(status_occ(i,1)+1+8*(map_family(id_tot(i,6),1)-1),time)+1/MC_run;
        col_elder(status_occ(i,1)+1+8*(map_family(id_tot(i,6),1)-1)+24*(run-1),time)=col_elder(status_occ(i,1)+1+8*(map_family(id_tot(i,6),1)-1)+24*(run-1),time)+1;
    end
end
end
for i=1:n_emp
    if Inf_info(i,1)>0.5
        Rt_cal(Inf_info(i,3)+1,1)=Rt_cal(Inf_info(i,3)+1,1)+1; % denominator of Rt 
        Rt_cal(Inf_info(i,4)+1,2)=Rt_cal(Inf_info(i,4)+1,2)+1; % numerator of Rt
        if Inf_info(i,2)>0.5
            Inf_info(Inf_info(i,2),5)=Inf_info(Inf_info(i,2),5)+1;
        end
    end
end
Rt_cal(1,2)=Rt_cal(1,2)-1;
for time=0:day_dua
    if Rt_cal(time+1,1)>0.5
        Rt_deno(time+1,1)=Rt_deno(time+1,1)+1;
        Rt_cal(time+1,3)=Rt_cal(time+1,2)/Rt_cal(time+1,1);
    end
end
Rt_mc=Rt_mc+Rt_cal(:,3);
%{
infector_tot=sum(Inf_info(:,1)); % total number of infectors during day_dua
Inf_newinf=zeros(infector_tot,1); % number of newly infected
inf_seq=0;
for i=1:n_emp
    if Inf_info(i,1)>0.5
        inf_seq=inf_seq+1;
        Inf_newinf(inf_seq,1)=Inf_info(i,5);
    end
end
disp=kvalue(Inf_newinf);
if not(isnan(disp(1,1))) && not(isnan(disp(1,2)))
    disp_num=disp_num+1;
    disp_mc=disp_mc+disp;
end
%}
end
% disp_mc=disp_mc/disp_num;

% Confidence interval calculation
% 90% confidence interval
[intup_all,intlow_all]=confint(col_all,MC_run,conf_int);
[intup_job,intlow_job]=confint(col_job,MC_run,conf_int);
[intup_elder,intlow_elder]=confint(col_elder,MC_run,conf_int);
[intup_spacetype,intlow_spacetype]=confint(col_spacetype,MC_run,conf_int);

for time=0:day_dua
    if Rt_deno(time+1,1)>0.5
        Rt_mc(time+1,1)=Rt_mc(time+1,1)/Rt_deno(time+1,1);
    else
        Rt_mc(time+1,1)=-1;
    end
end

%鸟
%load chirp
%sound(y,Fs)
