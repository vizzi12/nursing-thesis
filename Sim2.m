%%Simulation code for the bachelor thesis: Improvement of Care Management
%%at Nursing homes.

clear all
close all
clc
rng(123)


%Loading the poisson rates
lambda=zeros(24,3);
lambda(:,1)=table2array(readtable("PostServiceDelay\afd1\delay_rates_afd1"));
lambda(:,2)=table2array(readtable("PostServiceDelay\afd2\delay_rates_afd2"));
lambda(:,3)=table2array(readtable("PostServiceDelay\afd3\delay_rates_afd3"));

%Extracting the maximum value of lambda.
maxlambda=max(lambda);


%Loading the phase type distributions.
ServiceTime=zeros(4,5,6,3);
for i=1:3
    for j=1:6
        ServiceTime(:,:,j,i)=table2array(readtable(strcat("ServiceTime\afd",num2str(i),"\phases_",num2str((j-1)*4),"_",num2str(j*4-1),"_afd",num2str(i))));
    end
end


%Intializing the base case
residents=[14,15,16];
no_sections=length(residents);

all_caregivers=[repmat([11,12,13,21,22,23,31,32,33],8,1);repmat([11,12,0,21,22,0,31,32,0],8,1)];
section_caregivers=zeros(16,3,3);

section_caregivers(:,:,1)=[repmat([11,12,13],8,1);repmat([11,12,0,],8,1)];
section_caregivers(:,:,2)=[repmat([21,22,23],8,1);repmat([21,22,0,],8,1)];
section_caregivers(:,:,3)=[repmat([31,32,33],8,1);repmat([31,32,0,],8,1)];

hours=16; %7 to 23 
Maxgroup1=2;
primary=cell(3,1);
for i=1:no_sections
    primary{i}=zeros(hours,Maxgroup1,residents(i));
    % Adjust primary caregivers here
    n=round(residents(i)/3);
    primary{i}(1:8,1:2,1:n)=repmat(10*i+(1:2),8,1,n);
    primary{i}(1:8,1:2,n+1:2*n)=repmat(10*i+(2:3),8,1,n);
    primary{i}(1:8,1:2,2*n+1:end)=repmat(10*i+[1,3],8,1,residents(i)-2*n);
    
    primary{i}(9:16,1:2,:)=repmat(10*i+(1:2),8,1,residents(i));
end



secondary=cell(3,1);
Maxgroup2=3;
for i=1:no_sections
    secondary{i}=zeros(hours,Maxgroup2,residents(i));
    % Adjust secondary caregivers here
    n=round(residents(i)/3);
    secondary{i}(1:8,1,1:n)=(i*10+3)*ones(8,1,n);
    secondary{i}(1:8,1,n+1:2*n)=(i*10+1)*ones(8,1,n);
    secondary{i}(1:8,1,2*n+1:end)=(i*10+2)*ones(8,1,residents(i)-2*n);
    
    index=1:3;
    index(i)=[];
    secondary{i}(1:8,2:3,1:n)=(index*10+1).*ones(8,1,n);
    secondary{i}(1:8,2:3,n+1:2*n)=(index*10+2).*ones(8,1,n);
    secondary{i}(1:8,2:3,2*n+1:end)=(index*10+3).*ones(8,1,residents(i)-2*n); 
    
    n=round(residents(i)/2);
    secondary{i}(9:16,1:2,1:n)=(index*10+1).*ones(8,1,n);
    secondary{i}(9:16,1:2,n+1:end)=(index*10+2).*ones(8,1,residents(i)-n);
end


tertiary=cell(3,1);
Maxgroup3=2;

for i=1:no_sections
    tertiary{i}=zeros(hours,Maxgroup3,residents(i));
    % Adjust tertiary caregivers here
    n=round(residents(i)/3);
    
    index=1:3;
    index(i)=[];
    
    tertiary{i}(1:8,1:2,1:n)=(index*10+3).*ones(8,1,n);
    tertiary{i}(1:8,1:2,n+1:2*n)=(index*10+2).*ones(8,1,n);
    tertiary{i}(1:8,1:2,2*n+1:end)=(index*10+1).*ones(8,1,residents(i)-2*n);

    n=round(residents(i)/2);
    tertiary{i}(9:16,1,1:n)=((mod(i,3)+1)*10+1)*ones(8,1,n);
    tertiary{i}(9:16,1,n+1:end)=((mod(i,3)+1)*10+2)*ones(8,1,residents(i)-n);
end


%Randomizing the order of the caregivers to prevent section 1 caregivers
%always accepting care first.
for i=1:no_sections
    for j=1:residents(i)
        for k=1:16
            primary{i}(k,:,j)=primary{i}(k,randperm(2),j);
            secondary{i}(k,:,j)=secondary{i}(k,randperm(3),j);
            tertiary{i}(k,:,j)=tertiary{i}(k,randperm(2),j);
        end

    end
    
    for k=1:16
        section_caregivers(k,:,i)=section_caregivers(k,randperm(3),i);
    end
end


for i=1:16
   all_caregivers(i,:)=all_caregivers(i,randperm(9));
end



Priority_Probabilites=[1/10,8/10,1/10];
a=[0,cumsum(Priority_Probabilites)];
tau=180/3600;



%The scheduled care
Scheduled=[linspace(7+5*eps,8.75,residents(1))', 20/60*ones(residents(1),1), ones(residents(1),1), 2*ones(residents(1),1);
           linspace(7+5*eps,8.75,residents(2))', 20/60*ones(residents(2),1), 2*ones(residents(2),1), 2*ones(residents(2),1);
           linspace(7+5*eps,8.75,residents(3))', 20/60*ones(residents(3),1), 3*ones(residents(3),1), 2*ones(residents(3),1);
           
           linspace(9,11,3)', 60/60*ones(3,1), ones(3,1), 2*ones(3,1);
           linspace(9,11.5,6)', 30/60*ones(6,1), ones(6,1), 3*ones(6,1);
           linspace(9,11,3)', 60/60*ones(3,1), 2*ones(3,1), 2*ones(3,1);
           linspace(9,11.5,6)', 30/60*ones(6,1), 2*ones(6,1), 3*ones(6,1);
           linspace(9,11,3)', 60/60*ones(3,1), 3*ones(3,1), 2*ones(3,1);
           linspace(9,11.5,6)', 30/60*ones(6,1), 3*ones(6,1), 3*ones(6,1);
           
           12,1,1,1;
           12,1,2,1;
           12,1,3,1;
           linspace(12,12.75,residents(1))', 2/60*ones(residents(1),1), ones(residents(1),1), 2*ones(residents(1),1);
           linspace(12,12.75,residents(2))', 2/60*ones(residents(2),1), 2*ones(residents(2),1), 2*ones(residents(2),1);
           linspace(12,12.75,residents(3))', 2/60*ones(residents(3),1), 3*ones(residents(3),1), 2*ones(residents(3),1);
           
           linspace(13,14,2)', 60/60*ones(2,1), ones(2,1), 2*ones(2,1);
           linspace(13,14.5,3)', 30/60*ones(3,1), ones(3,1), 3*ones(3,1);
           linspace(13,14,2)', 60/60*ones(2,1), 2*ones(2,1), 2*ones(2,1);
           linspace(13,14.5,3)', 30/60*ones(3,1), 2*ones(3,1), 3*ones(3,1);
           linspace(13,14,2)', 60/60*ones(2,1), 3*ones(2,1), 2*ones(2,1);
           linspace(13,14.5,3)', 30/60*ones(3,1), 3*ones(3,1), 3*ones(3,1);
           
           linspace(15,16,2)', 60/60*ones(2,1), ones(2,1), 2*ones(2,1);
           linspace(15,16,2)', 60/60*ones(2,1), 2*ones(2,1), 2*ones(2,1);
           linspace(15,16,2)', 60/60*ones(2,1), 3*ones(2,1), 2*ones(2,1);

           17,1,1,1;
           17,1,2,1;
           17,1,3,1;
           linspace(17,17.75,residents(1))', 2/60*ones(residents(1),1), ones(residents(1),1), 2*ones(residents(1),1);
           linspace(17,17.75,residents(2))', 2/60*ones(residents(2),1), 2*ones(residents(2),1), 2*ones(residents(2),1);
           linspace(17,17.75,residents(3))', 2/60*ones(residents(3),1), 3*ones(residents(3),1), 2*ones(residents(3),1);
           
           linspace(18,20,3)', 60/60*ones(3,1), ones(3,1), 2*ones(3,1);
           linspace(18,20,3)', 60/60*ones(3,1), 2*ones(3,1), 2*ones(3,1);
           linspace(18,20,3)', 60/60*ones(3,1), 3*ones(3,1), 2*ones(3,1);
        
           linspace(21,23-10/60,residents(1))', 10/60*ones(residents(1),1), ones(residents(1),1), 2*ones(residents(1),1);
           linspace(21,23-10/60,residents(2))', 10/60*ones(residents(2),1), 2*ones(residents(2),1), 2*ones(residents(2),1);
           linspace(21,23-10/60,residents(3))', 10/60*ones(residents(3),1), 3*ones(residents(3),1), 2*ones(residents(3),1);
           
           ];
       
       
Scheduled=Scheduled(randperm(size(Scheduled,1)),:);
Scheduled(:,2)=Scheduled(:,2)*0.99999;

%%

days=5000;


%Arrays to store information
max_jobs=35;
finished_jobs=zeros(10,days*max_jobs);

no_scheduled=size(Scheduled,1);
finished_scheduled=zeros(10,days*no_scheduled);
    c=0; % counter
scheduled_c=0;


for d=1:days %Looping over all days
 
    t=7;
    call_times=cell(length(residents),1);
    for i=1:no_sections
        call_times{i}=zeros(residents(i),1);
        for j=1:residents(i)
            call_times{i}(j)=generate(lambda(:,i),maxlambda(i),t);
        end
    end

    queue=[];
    jobs=[0;0;100;100;0;0];
    jobs_server_id=0;
    jobs_priority=0;

    times=[call_times; Scheduled(:,1); 100;100;100];

    while true % Looping over one day
        
        %Extracting next event
        [x,idx]=cellfun(@min,times,'UniformOutput',0);
        if any(cellfun(@isempty,x))
            idx{cellfun(@isempty,x)}=1;
            x{cellfun(@isempty,x)}=100;
        end

        x=cell2mat(x);
        idx=cell2mat(idx);
        
        [t,idx2]=min(x);
        

        if t > 23
            break
        end
        
        if idx2 <= no_sections %Call
            c=c+1;
            resident_id=idx(idx2);
            p = find(rand(1) > a ,1,'last');
            T=simphasetype(ServiceTime(:,:,ceil(t/4),idx2))/3600;

            times{idx2}(resident_id)=100;
            
            finished_jobs(1,c)=idx2;
            finished_jobs(2,c)=t;
            finished_jobs(4,c)=T;
            finished_jobs(7,c)=p;
            finished_jobs(8,c)=d;
            
            if ismember(primary{idx2}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority <= p )) %Add to queue if primary caregivers are busy doing more important or equally important
                queue=[queue,[idx2;resident_id;t;T;1;c;p]]; % 1 for call (escalation level)
                times{no_sections+3}=queue(3,:)+queue(5,:)*tau;
                times{no_sections+3}(queue(5,:) == 4)= 100;


            elseif ismember(primary{idx2}(ceil(t)-7,:,resident_id),jobs_server_id) %If a primary caregiver needs to reprioritize
                for i=3:-1:p+1
                   index=find(ismember(primary{idx2}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority == i)),1,'first');
                   if ~isempty(index)  
                       job_idx=find(primary{idx2}(ceil(t)-7,index,resident_id)== jobs_server_id,1);
                       break
                   end
                end

                queue=[queue, [jobs(:,job_idx);i]];
                times{no_sections+4}=t-1;

                server_id=jobs_server_id(job_idx);
                jobs_server_id(job_idx)=[];
                jobs_priority(job_idx)=[];
                jobs(:,job_idx)=[];

                jobs_server_id=[jobs_server_id,server_id];
                jobs_priority=[jobs_priority,p];
                jobs=[jobs,[idx2;resident_id;t;T;1;c]];
                times{no_sections+2}=jobs(3,:)+jobs(4,:);
                
                finished_jobs(3,c)=t;


            else %If the call can be attended right away
                idx3=find(~ismember(primary{idx2}(ceil(t)-7,:,resident_id),jobs_server_id),1,'first');

                jobs_server_id=[jobs_server_id,primary{idx2}(ceil(t)-7,idx3,resident_id)];
                jobs_priority=[jobs_priority,p];
                jobs=[jobs,[idx2;resident_id;t;T;1;c]];
                times{no_sections+2}=jobs(3,:)+jobs(4,:);
                
                finished_jobs(3,c)=t;
            end


        elseif idx2 == no_sections+1 %Scheduled
            scheduled_c=scheduled_c+1;
            
            T=Scheduled(idx(idx2),2);
            section_id=Scheduled(idx(idx2),3);
            resident_id=0;
            p = Scheduled(idx(idx2),4);

            times{idx2}(idx(idx2))=100;
            
            finished_scheduled(1,scheduled_c)=section_id;
            finished_scheduled(2,scheduled_c)=t;
 
            finished_scheduled(4,scheduled_c)=T;
            finished_scheduled(7,scheduled_c)=p;
            finished_scheduled(8,scheduled_c)=d;

            if ismember(section_caregivers(ceil(t)-7,:,section_id),jobs_server_id(jobs_priority <= p | jobs(2,:) > 0)) %Add to queue if caregivers are busy
                queue=[queue,[section_id;resident_id;t;T;10^3;scheduled_c;p]]; % 10^3 for scheduled care
                times{no_sections+3}=queue(3,:)+queue(5,:)*tau;
                times{no_sections+3}(queue(5,:) == 4)= 100;
                
            elseif ismember(section_caregivers(ceil(t)-7,:,section_id),jobs_server_id) %reprio scheduled care
                for i=3:-1:p+1
                   index=find(ismember(section_caregivers(ceil(t)-7,:,section_id),jobs_server_id(jobs(2,:) == 0 & jobs_priority == i)),1,'first');
                   if ~isempty(index)  
                       job_idx=find(section_caregivers(ceil(t)-7,index,section_id)== jobs_server_id,1);
                       break
                   end
                end

                queue=[queue, [jobs(:,job_idx);i]];
                times{no_sections+4}=t-1;

                server_id=jobs_server_id(job_idx);
                jobs_server_id(job_idx)=[];
                jobs_priority(job_idx)=[];
                jobs(:,job_idx)=[];

                jobs_server_id=[jobs_server_id,server_id];
                jobs_priority=[jobs_priority,p];
                jobs=[jobs,[section_id;resident_id;t;T;10^3;scheduled_c]];
                times{no_sections+2}=jobs(3,:)+jobs(4,:);
                
                finished_scheduled(3,scheduled_c)=t; 
            
            

            else %If scheduled can be performed right away
                idx3=find(~ismember(section_caregivers(ceil(t)-7,:,section_id),jobs_server_id),1,'first');

                jobs_server_id=[jobs_server_id,section_caregivers(ceil(t)-7,idx3,section_id)];
                jobs_priority=[jobs_priority,p];
                jobs=[jobs,[section_id;resident_id;t;T;10^3;scheduled_c]];
                times{no_sections+2}=jobs(3,:)+jobs(4,:);
                
                finished_scheduled(3,scheduled_c)=t;
            end


        elseif idx2 == no_sections+2 %Job finish
            job_idx=idx(idx2);

            if jobs(2,job_idx) % If the job is a call generate a new call
               times{jobs(1,job_idx)}(jobs(2,job_idx))=generate(lambda(:,jobs(1,job_idx)),maxlambda(jobs(1,job_idx)),t); 
               
               finished_jobs(5,jobs(6,job_idx))=t;
               
               if ismember(jobs_server_id(job_idx),primary{jobs(1,job_idx)}(ceil(jobs(3,job_idx))-7,:,jobs(2,job_idx)))
                finished_jobs(6,jobs(6,job_idx))=1;
               elseif ismember(jobs_server_id(job_idx),secondary{jobs(1,job_idx)}(ceil(jobs(3,job_idx))-7,:,jobs(2,job_idx)))
                   finished_jobs(6,jobs(6,job_idx))=2;
               elseif ismember(jobs_server_id(job_idx),tertiary{jobs(1,job_idx)}(ceil(jobs(3,job_idx))-7,:,jobs(2,job_idx)))
                   finished_jobs(6,jobs(6,job_idx))=3;
               else
                   finished_jobs(6,jobs(6,job_idx))=4;
               end
               
            else
               finished_scheduled(5,jobs(6,job_idx))=t;
            end

            free_server=jobs_server_id(job_idx);
            jobs_server_id(job_idx)=[];
            jobs_priority(job_idx)=[];
            jobs(:,job_idx)=[];

            times{no_sections+2}(job_idx)=[];

            section_id=round(free_server-6,-1)/10;
            flag=0;
  
            if ~isempty(queue)      
                for k=1:3
                    idx_array=queue(7,:)==k;
                    priorityq=queue(:,idx_array);

                    for i=1:size(priorityq,2)
                        if priorityq(2,i) && ismember(free_server,primary{priorityq(1,i)}(ceil(t)-7,:,priorityq(2,i)))
                            
                           jobs=[jobs,priorityq(1:6,i)];
                           jobs(3,end)=t;
                           jobs_server_id=[jobs_server_id,free_server];
                           jobs_priority=[jobs_priority,priorityq(7,i)];
                          
                           if finished_jobs(3,jobs(6,end))== 0
                                finished_jobs(3,jobs(6,end))=t;
                           end

                           ii=find(idx_array,i);
                           queue(:,ii(end))=[];

                           times{no_sections+3}(ii(end))=[];
                           times{no_sections+2}=jobs(3,:)+jobs(4,:);

                           flag=1;
                           break
                        end
                    end
                    if flag==1
                        break
                    end
                    for i=1:size(priorityq,2)
                            if priorityq(5,i) > 1  && priorityq(2,i) && ismember(free_server,secondary{priorityq(1,i)}(ceil(t)-7,:,priorityq(2,i)))
                                jobs=[jobs,priorityq(1:6,i)];
                                jobs(3,end)=t;
                                jobs_server_id=[jobs_server_id,free_server];
                                jobs_priority=[jobs_priority,priorityq(7,i)];
                                
                                if finished_jobs(3,jobs(6,end))== 0
                                    finished_jobs(3,jobs(6,end))=t;
                                end
                                
                                ii=find(idx_array,i);
                                queue(:,ii(end))=[];

                                times{no_sections+3}(ii(end))=[];
                                times{no_sections+2}=jobs(3,:)+jobs(4,:);

                                flag=1;
                                break
                            end
                    end
                    if flag==1
                        break
                    end
                    for i=1:size(priorityq,2)
                            if priorityq(5,i) > 2 && priorityq(2,i) && ismember(free_server,tertiary{priorityq(1,i)}(ceil(t)-7,:,priorityq(2,i)))
                                jobs=[jobs,priorityq(1:6,i)];
                                jobs(3,end)=t;
                                jobs_server_id=[jobs_server_id,free_server];
                                jobs_priority=[jobs_priority,priorityq(7,i)];
                                
                                if finished_jobs(3,jobs(6,end))== 0
                                    finished_jobs(3,jobs(6,end))=t;
                                end

                                ii=find(idx_array,i);
                                queue(:,ii(end))=[];

                                times{no_sections+3}(ii(end))=[];
                                times{no_sections+2}=jobs(3,:)+jobs(4,:);                    
                                flag=1;
                                break
                            end
                    end
                    if flag==1
                        break
                    end
                    for i=1:size(priorityq,2)
                            if priorityq(2,i) && priorityq(5,i) ==4
                                jobs=[jobs,priorityq(1:6,i)];
                                jobs(3,end)=t;
                                jobs_server_id=[jobs_server_id,free_server];
                                jobs_priority=[jobs_priority,priorityq(7,i)];
                                
                                if finished_jobs(3,jobs(6,end))== 0
                                    finished_jobs(3,jobs(6,end))=t;
                                end

                                ii=find(idx_array,i);
                                queue(:,ii(end))=[];

                                times{no_sections+3}(ii(end))=[];
                                times{no_sections+2}=jobs(3,:)+jobs(4,:);   

                                flag=1;
                                break
                            end
                    end
                    if flag==1
                        break
                    end
                    for i=1:size(priorityq,2)
                        if ~priorityq(2,i) && priorityq(1,i)== section_id
                            jobs=[jobs,priorityq(1:6,i)];
                            jobs(3,end)=t;
                            jobs_server_id=[jobs_server_id,free_server];
                            jobs_priority=[jobs_priority,priorityq(7,i)];
                            
                            if finished_scheduled(3,jobs(6,end))== 0
                                    finished_scheduled(3,jobs(6,end))=t;
                            end
                            
                            ii=find(idx_array,i);
                            
                            times{no_sections+3}(ii(end))=[];
                            times{no_sections+2}=jobs(3,:)+jobs(4,:);
                            ii=find(idx_array,i);
                            queue(:,ii(end))=[];
                            flag=1;
                            break
                        end
                    end
                    if flag==1
                        break
                    end
                end
            end

        elseif idx2 == no_sections+3 %Queue escalation
            queue_idx=idx(idx2);

            queue(5,queue_idx)=queue(5,queue_idx)+1;
            section_id=queue(1,queue_idx);
            resident_id=queue(2,queue_idx);
            p=queue(7,queue_idx);

            switch queue(5,queue_idx)
                case 2 % Check secondary caregivers

                    if any(~ismember(secondary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id)) %secondary available
                         idx3=find(~ismember(secondary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id),1,'first');

                         jobs_server_id=[jobs_server_id,secondary{section_id}(ceil(t)-7,idx3,resident_id)];
                         jobs_priority=[jobs_priority,p];

                         jobs=[jobs,queue(1:6,queue_idx)];
                         jobs(3,end)=t;
                         
                          if finished_jobs(3,jobs(6,end))== 0
                                    finished_jobs(3,jobs(6,end))=t;
                          end

                         times{no_sections+2}=jobs(3,:)+jobs(4,:);

                         queue(:,queue_idx)=[];
                         times{no_sections+3}(queue_idx)=[];
                    elseif any(~ismember(secondary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority <= p))) %secondary needs to prioritize
                        for i=3:-1:p+1
                            index=find(ismember(secondary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority == i)),1,'first');
                            if ~isempty(index)  
                                job_idx=find(secondary{section_id}(ceil(t)-7,index,resident_id)== jobs_server_id,1);
                                break    
                            end
                        end


                        queue=[queue, [jobs(:,job_idx);i]];
                        times{no_sections+4}=t-1;

                        server_id=jobs_server_id(job_idx);
                        jobs_server_id(job_idx)=[];
                        jobs_priority(job_idx)=[];
                        jobs(:,job_idx)=[];
                        

                        jobs_server_id=[jobs_server_id,server_id];
                        jobs_priority=[jobs_priority,p];
                        jobs=[jobs,queue(1:6,queue_idx)];
                        jobs(3,end)=t;
                        times{no_sections+2}=jobs(3,:)+jobs(4,:); 
                        
                        if finished_jobs(3,jobs(6,end))== 0
                                finished_jobs(3,jobs(6,end))=t;
                        end

                        queue(:,queue_idx)=[];
                        times{no_sections+3}(queue_idx)=[];

                    else
                        times{no_sections+3}=queue(3,:)+queue(5,:)*tau;
                        times{no_sections+3}(queue(5,:) == 4)= 100;
                    end

                case 3 % Check tertiary caregivers
                    if any(~ismember(tertiary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id)) %tertiary available
                         idx3=find(~ismember(tertiary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id),1,'first');

                         jobs_server_id=[jobs_server_id,tertiary{section_id}(ceil(t)-7,idx3,resident_id)];
                         jobs_priority=[jobs_priority,p];
                         jobs=[jobs,queue(1:6,queue_idx)];
                         jobs(3,end)=t;
                         
                         if finished_jobs(3,jobs(6,end))== 0
                                finished_jobs(3,jobs(6,end))=t;
                         end
                         times{no_sections+2}=jobs(3,:)+jobs(4,:);
                         queue(:,queue_idx)=[];
                         times{no_sections+3}(queue_idx)=[];
                    elseif any(~ismember(tertiary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority <= p))) % tertiary needs to prioritize
                        for i=3:-1:p+1
                            index=find(ismember(tertiary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority == i)),1,'first');
                            if ~isempty(index)  
                                job_idx=find(tertiary{section_id}(ceil(t)-7,index,resident_id)== jobs_server_id,1);
                                break %%%%%%%%%%%%%%%%%
                            end
                        end

                         times{no_sections+4}=t-1;
                         queue=[queue, [jobs(:,job_idx);i]];


                        server_id=jobs_server_id(job_idx);
                        jobs_server_id(job_idx)=[];
                        jobs_priority(job_idx)=[];
                        jobs(:,job_idx)=[];

                        jobs_server_id=[jobs_server_id,server_id];
                        jobs_priority=[jobs_priority,p];
                        jobs=[jobs,queue(1:6,queue_idx)];
                        jobs(3,end)=t; 
                        
                        if finished_jobs(3,jobs(6,end))== 0
                                finished_jobs(3,jobs(6,end))=t;
                        end

                        queue(:,queue_idx)=[];
                        times{no_sections+3}(queue_idx)=[];


                    else
                        times{no_sections+3}=queue(3,:)+queue(5,:)*tau;
                        times{no_sections+3}(queue(5,:) == 4)= 100;
                    end

                case 4 % Check all caregivers
                    if any(~ismember(all_caregivers(ceil(t)-7,:),jobs_server_id))
                         idx3=find(~ismember(all_caregivers(ceil(t)-7,:),jobs_server_id),1,'first');

                         jobs_server_id=[jobs_server_id,all_caregivers(ceil(t)-7,idx3)];
                         jobs_priority=[jobs_priority,p];
                         jobs=[jobs,queue(1:6,queue_idx)];
                         jobs(3,end)=t;
                         times{no_sections+2}=jobs(3,:)+jobs(4,:);
                         
                         if finished_jobs(3,jobs(6,end))== 0
                                finished_jobs(3,jobs(6,end))=t;
                         end

                         queue(:,queue_idx)=[];
                         times{no_sections+3}(queue_idx)=[];

                    elseif any(~ismember(all_caregivers(ceil(t)-7,:),jobs_server_id(jobs_priority <= p))) %some needs to reprio
                        for i=3:-1:p+1
                            index=find(ismember(all_caregivers(ceil(t)-7,:),jobs_server_id(jobs_priority == i)),1,'first');
                            if ~isempty(index)  
                                job_idx=find(all_caregivers(ceil(t)-7,index)== jobs_server_id,1);
                                break %%%%%%%%%%%%%%%%%
                            end
                        end


                        times{no_sections+4}=t-1;
                        queue=[queue, [jobs(:,job_idx);i]];


                        server_id=jobs_server_id(job_idx);
                        jobs_server_id(job_idx)=[];
                        jobs_priority(job_idx)=[];
                        jobs(:,job_idx)=[];

                        jobs_server_id=[jobs_server_id,server_id];
                        jobs_priority=[jobs_priority,p];
                        jobs=[jobs,queue(1:6,queue_idx)];
                        jobs(3,end)=t; 
                        times{no_sections+2}=jobs(3,:)+jobs(4,:);
                        
                        if finished_jobs(3,jobs(6,end))== 0
                                finished_jobs(3,jobs(6,end))=t;
                        end

                        queue(:,queue_idx)=[];
                        times{no_sections+3}(queue_idx)=[];

                    else
                        times{no_sections+3}(queue_idx)=100;
                    end
            end
        else % reprio job
            t=t+1;
            queue_idx=size(queue,2);
            section_id=queue(1,queue_idx);
            resident_id=queue(2,queue_idx);
            p=queue(7,queue_idx);
            
            cc=queue(6,queue_idx);

            T=queue(4,queue_idx)-(t-queue(3,queue_idx));

           times{no_sections+4}=100;

            if resident_id %if it is a call
                
                finished_jobs(9,cc)=finished_jobs(9,cc)+1;
                
                queue(5,queue_idx)=min(4,ceil((t-queue(3,queue_idx))/tau));
                queue(3:4,end)=[t;T];

                    if any(~ismember(primary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id)) %primary avalable available
                         idx3=find(~ismember(primary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id),1,'first');

                         jobs_server_id=[jobs_server_id,primary{section_id}(ceil(t)-7,idx3,resident_id)];
                         jobs_priority=[jobs_priority,p];

                         jobs=[jobs,queue(1:6,queue_idx)];
                         jobs(3,end)=t;

                         times{no_sections+2}=jobs(3,:)+jobs(4,:);

                         queue(:,queue_idx)=[];

                    elseif queue(5,queue_idx) > 1 &&  any(~ismember(secondary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id))
                         idx3=find(~ismember(secondary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id),1,'first');

                         jobs_server_id=[jobs_server_id,secondary{section_id}(ceil(t)-7,idx3,resident_id)];
                         jobs_priority=[jobs_priority,p];

                         jobs=[jobs,queue(1:6,queue_idx)];
                         jobs(3,end)=t;

                         times{no_sections+2}=jobs(3,:)+jobs(4,:);

                         queue(:,queue_idx)=[];

                    elseif queue(5,queue_idx) > 2 &&  any(~ismember(tertiary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id)) 
                         idx3=find(~ismember(tertiary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id),1,'first');

                         jobs_server_id=[jobs_server_id,tertiary{section_id}(ceil(t)-7,idx3,resident_id)];
                         jobs_priority=[jobs_priority,p];

                         jobs=[jobs,queue(1:6,queue_idx)];
                         jobs(3,end)=t;

                         times{no_sections+2}=jobs(3,:)+jobs(4,:);

                         queue(:,queue_idx)=[];

                    elseif queue(5,queue_idx) > 3 &&  any(~ismember(all_caregivers(ceil(t)-7,:),jobs_server_id)) 
                         idx3=find(~ismember(all_caregivers(ceil(t)-7,:),jobs_server_id),1,'first');

                         jobs_server_id=[jobs_server_id,all_caregivers(ceil(t)-7,idx3)];
                         jobs_priority=[jobs_priority,p];

                         jobs=[jobs,queue(1:6,queue_idx)];
                         jobs(3,end)=t;

                         times{no_sections+2}=jobs(3,:)+jobs(4,:);

                         queue(:,queue_idx)=[];

                    elseif  any(~ismember(primary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority <= p)))%primary reprio
                        for i=3:-1:p+1
                            index=find(ismember(primary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority == i)),1,'first');
                            if ~isempty(index)  
                                job_idx=find(primary{section_id}(ceil(t)-7,index,resident_id)== jobs_server_id,1);
                                break %%%%%%%%%%%%%
                            end
                        end

                        queue=[queue, [jobs(:,job_idx);i]];
                        times{no_sections+4}=t-1;

                        server_id=jobs_server_id(job_idx);
                        jobs_server_id(job_idx)=[];
                        jobs_priority(job_idx)=[];
                        jobs(:,job_idx)=[];

                        jobs_server_id=[jobs_server_id,server_id];
                        jobs_priority=[jobs_priority,p];
                        jobs=[jobs,queue(1:6,queue_idx)];
                        jobs(3,end)=t;
                        times{no_sections+2}=jobs(3,:)+jobs(4,:); 

                        queue(:,queue_idx)=[];

                    elseif queue(5,queue_idx) > 1 && any(~ismember(secondary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority <= p))) %sec reprio
                        for i=3:-1:p+1
                            index=find(ismember(secondary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority == i)),1,'first');
                            if ~isempty(index)  
                                job_idx=find(secondary{section_id}(ceil(t)-7,index,resident_id)== jobs_server_id,1);
                                break %%%%%%%%%%%%%%%%%%%%%%%%%%
                            end
                        end


                        queue=[queue, [jobs(:,job_idx);i]];
                        times{no_sections+4}=t-1;

                        server_id=jobs_server_id(job_idx);
                        jobs_server_id(job_idx)=[];
                        jobs_priority(job_idx)=[];
                        jobs(:,job_idx)=[];

                        jobs_server_id=[jobs_server_id,server_id];
                        jobs_priority=[jobs_priority,p];
                        jobs=[jobs,queue(1:6,queue_idx)];
                        jobs(3,end)=t;
                        times{no_sections+2}=jobs(3,:)+jobs(4,:); 

                        queue(:,queue_idx)=[];

                    elseif queue(5,queue_idx) > 2 && any(~ismember(tertiary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority <= p))) %tert reprio
                        for i=3:-1:p+1
                            index=find(ismember(tertiary{section_id}(ceil(t)-7,:,resident_id),jobs_server_id(jobs_priority == i)),1,'first');
                            if ~isempty(index)  
                                job_idx=find(tertiary{section_id}(ceil(t)-7,index,resident_id)== jobs_server_id,1);
                                break %%%%%%%%%%%%%%%%%%%%%%%
                            end
                        end     

                        queue=[queue, [jobs(:,job_idx);i]];
                        times{no_sections+4}=t-1;

                        server_id=jobs_server_id(job_idx);
                        jobs_server_id(job_idx)=[];
                        jobs_priority(job_idx)=[];
                        jobs(:,job_idx)=[];

                        jobs_server_id=[jobs_server_id,server_id];
                        jobs_priority=[jobs_priority,p];
                        jobs=[jobs,queue(1:6,queue_idx)];
                        jobs(3,end)=t;
                        times{no_sections+2}=jobs(3,:)+jobs(4,:); 

                        queue(:,queue_idx)=[];

                    elseif queue(5,queue_idx) > 3 && any(~ismember(all_caregivers(ceil(t)-7,:),jobs_server_id(jobs_priority <= p))) %some needs to reprio
                        for i=3:-1:p+1
                            index=find(ismember(all_caregivers(ceil(t)-7,:),jobs_server_id(jobs_priority == i)),1,'first');
                            if ~isempty(index)  
                                job_idx=find(all_caregivers(ceil(t)-7,index)== jobs_server_id,1);
                                break %%%%%%%%%%%%%%%
                            end
                        end
                        times{no_sections+4}=t-1;
                        queue=[queue, [jobs(:,job_idx);i]];


                        server_id=jobs_server_id(job_idx);
                        jobs_server_id(job_idx)=[];
                        jobs_priority(job_idx)=[];
                        jobs(:,job_idx)=[];

                        jobs_server_id=[jobs_server_id,server_id];
                        jobs_priority=[jobs_priority,p];
                        jobs=[jobs,queue(1:6,queue_idx)];
                        jobs(3,end)=t; 

                        queue(:,queue_idx)=[];

                    else
                        finished_jobs(10,cc)=finished_jobs(10,cc)+1;
                        if queue(5,end) <= 4
                            times{no_sections+3}=queue(3,:)+queue(5,:)*tau;
                            times{no_sections+3}(queue(5,:) == 4)= 100;
                        elseif isempty(times{no_sections+3})
                            times{no_sections+3}=100;
                        end
                    end

            else % scheduled care
                finished_scheduled(9,cc)=finished_scheduled(9,cc)+1;
                if ismember(section_caregivers(ceil(t)-7,:,section_id),jobs_server_id) %remain in queue if caregivers are busy
                    queue(3:4,end)=[t;T];
                    times{no_sections+3}=queue(3,:)+queue(5,:)*tau;
                    times{no_sections+3}(queue(5,:) == 4)= 100;
                    finished_scheduled(10,cc)=finished_scheduled(10,cc)+1;
                    
                    
                else %If scheduled can be performed right away
                    idx3=find(~ismember(section_caregivers(ceil(t)-7,:,section_id),jobs_server_id),1,'first');

                    jobs_server_id=[jobs_server_id,section_caregivers(ceil(t)-7,idx3,section_id)];
                    jobs_priority=[jobs_priority,p];
                    jobs=[jobs,[section_id;resident_id;t;T;10^3; queue(6,end)]];
                    times{no_sections+2}=jobs(3,:)+jobs(4,:);

                    queue(:,queue_idx)=[];
                end
            end 
                     
        end        
    end


end


%% Plotting and exporting data as the plots in the thesis are done in Tikz

finished_jobs(:,finished_jobs(5,:)==0)=[];
finished_scheduled(:,finished_scheduled(5,:)==0)=[];


%waiting times for each section

wait=finished_jobs(3,:)-finished_jobs(2,:);

wait_conf=zeros(16,2);

wait_q=zeros(16,5);

k=1;
for t=7:22
    wait_conf(k,1)=mean(wait(finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600;
    wait_conf(k,2)=var(wait(finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1)*3600);
   
    
    wait_q(k,:)=quantile(wait(finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),[0.25 0.5 0.75 0.95, 0.99])*3600;
    
    
     k=k+1;
end

figure
plot(7:22,wait_conf(:,1))
hold on
plot(7:22,wait_conf(:,1)+norminv(0.975)*sqrt(wait_conf(:,2))/(sqrt(days)))
plot(7:22,wait_conf(:,1)-norminv(0.975)*sqrt(wait_conf(:,2))/(sqrt(days)))

fprintf(fopen('base.dat','w'), '%f %f %f %f\n', [7:22;wait_conf(:,1)';(wait_conf(:,1)-norminv(0.975)*sqrt(wait_conf(:,2))/(sqrt(days)))';(wait_conf(:,1)-norminv(0.975)*sqrt(wait_conf(:,2))/(sqrt(days)))']);


fprintf(fopen('basemin.dat','w'), '%f %f %f %f\n', [7:22;(wait_conf(:,1)/60)';((wait_conf(:,1)-norminv(0.975)*sqrt(wait_conf(:,2))/(sqrt(days)))/60)';((wait_conf(:,1)-norminv(0.975)*sqrt(wait_conf(:,2))/(sqrt(days)))/60)']);


fprintf(fopen('baseqmin.dat','w'), '%f %f %f %f %f %f\n', [7:22;(wait_q/60)']);

fprintf(fopen('baseq.dat','w'), '%f %f %f %f %f %f\n', [7:22;wait_q']);


wait1_q=zeros(16,5);
wait2_q=zeros(16,5);
wait3_q=zeros(16,5);


wait_times=zeros(16,3);
wait_sd=zeros(16,3);
k=1;
for t=7:22
    wait_times(k,1)=mean(wait(finished_jobs(1,:) == 1 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600;
    wait_times(k,2)=mean(wait(finished_jobs(1,:) == 2 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600;
    wait_times(k,3)=mean(wait(finished_jobs(1,:) == 3 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600;
    
    
    wait1_q(k,:)=quantile(wait(finished_jobs(1,:) == 1 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),[0.25 0.5 0.75 0.95, 0.99])*3600;
    wait2_q(k,:)=quantile(wait(finished_jobs(1,:) == 2 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),[0.25 0.5 0.75 0.95, 0.99])*3600;
    wait3_q(k,:)=quantile(wait(finished_jobs(1,:) == 3 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),[0.25 0.5 0.75 0.95, 0.99])*3600;

    N=size(wait(finished_jobs(1,:) == 1 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),2);
    wait_sd(k,1)=sqrt(var((wait(finished_jobs(1,:) == 1 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600)/N);
    
    N=size(wait(finished_jobs(1,:) == 2 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),2);
    wait_sd(k,2)=sqrt(var((wait(finished_jobs(1,:) == 2 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600)/N);
    
    N=size(wait(finished_jobs(1,:) == 3 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),2);
    wait_sd(k,3)=sqrt(var((wait(finished_jobs(1,:) == 3 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600)/N);
    
    k=k+1;
end

fprintf(fopen('wait1.dat','w'), '%f %f %f %f %f %f\n', [7:22; wait1_q']);
fprintf(fopen('wait2.dat','w'), '%f %f %f %f %f %f\n', [7:22; wait2_q']);
fprintf(fopen('wait3.dat','w'), '%f %f %f %f %f %f\n', [7:22; wait3_q']);


fprintf(fopen('waitsec1.dat','w'), '%f %f %f %f\n', [7:22;wait_times(:,1)';wait_times(:,1)'-norminv(0.995)*wait_sd(:,1)';wait_times(:,1)'+norminv(0.975)*wait_sd(:,1)']);

fprintf(fopen('waitsec2.dat','w'), '%f %f %f %f\n', [7:22;wait_times(:,2)';wait_times(:,2)'-norminv(0.995)*wait_sd(:,2)';wait_times(:,2)'+norminv(0.975)*wait_sd(:,2)']);

fprintf(fopen('waitsec3.dat','w'), '%f %f %f %f\n', [7:22;wait_times(:,3)';wait_times(:,3)'-norminv(0.995)*wait_sd(:,3)';wait_times(:,3)'+norminv(0.975)*wait_sd(:,3)']);






figure
plot(7:22,wait_times(:,1))
hold on
plot(7:22,wait_times(:,2))
plot(7:22,wait_times(:,3))

legend('Section 1','Section 2','Section 3','location', 'best')
xlabel('$t$ (hour of day)','Interpreter', 'Latex','FontSize', 10)
ylabel('Waiting time (seconds)', 'FontSize', 10)
title('Simulated average waiting time over the day for the three sections')


% Waiting times for each priority
wait1_p=zeros(16,5);
wait2_p=zeros(16,5);
wait3_p=zeros(16,5);


wait_prio=zeros(16,3);
prio_sd=zeros(16,3);
k=1;
for t=7:22
    wait_prio(k,1)=mean(wait(finished_jobs(7,:) == 1 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600;
    wait_prio(k,2)=mean(wait(finished_jobs(7,:) == 2 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600;
    wait_prio(k,3)=mean(wait(finished_jobs(7,:) == 3 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600;
    
    wait1_p(k,:)=quantile(wait(finished_jobs(7,:) == 1 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),[0.25 0.5 0.75 0.95, 0.99])*3600;
    wait2_p(k,:)=quantile(wait(finished_jobs(7,:) == 2 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),[0.25 0.5 0.75 0.95, 0.99])*3600;
    wait3_p(k,:)=quantile(wait(finished_jobs(7,:) == 3 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),[0.25 0.5 0.75 0.95, 0.99])*3600;
    
    
    N=size(wait(finished_jobs(7,:) == 1 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),2);
    prio_sd(k,1)=sqrt(var((wait(finished_jobs(7,:) == 1 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600)/N);
    
    N=size(wait(finished_jobs(7,:) == 2 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),2);
    prio_sd(k,2)=sqrt(var((wait(finished_jobs(7,:) == 2 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600)/N);
    
    N=size(wait(finished_jobs(7,:) == 3 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1),2);
    prio_sd(k,3)=sqrt(var((wait(finished_jobs(7,:) == 3 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1))*3600)/N);
    
    k=k+1;
end

p1=[7:22;wait_prio(:,1)';wait_prio(:,1)'-norminv(0.995)*prio_sd(:,1)';wait_prio(:,1)'+norminv(0.995)*prio_sd(:,1)'];

p1( p1 < 0) =0;

fprintf(fopen('prioq1.dat','w'), '%f %f %f %f %f %f\n', [7:22; wait1_p']);
fprintf(fopen('prioq2.dat','w'), '%f %f %f %f %f %f\n', [7:22; wait2_p']);
fprintf(fopen('prioq3.dat','w'), '%f %f %f %f %f %f\n', [7:22; wait3_p']);

fprintf(fopen('waitprio1.dat','w'), '%f %f %f %f\n', p1);

fprintf(fopen('waitprio2.dat','w'), '%f %f %f %f\n', [7:22;wait_prio(:,2)';wait_prio(:,2)'-norminv(0.995)*prio_sd(:,2)';wait_prio(:,2)'+norminv(0.995)*prio_sd(:,2)']);

fprintf(fopen('waitprio3.dat','w'), '%f %f %f %f\n', [7:22;wait_prio(:,3)';wait_prio(:,3)'-norminv(0.995)*prio_sd(:,3)';wait_prio(:,3)'+norminv(0.995)*prio_sd(:,3)']);


figure
plot(7:22,wait_prio(:,1))
hold on
plot(7:22,wait_prio(:,2))
plot(7:22,wait_prio(:,3))

legend('Priority 1','Priority 2','Priority 3','location', 'best')
xlabel('$t$ (hour of day)','Interpreter', 'Latex','FontSize', 10)
ylabel('Waiting time (seconds)', 'FontSize', 10)
title('Simulated average waiting time over the day for the three priorities','FontSize', 10)

wait1=wait(finished_jobs(7,:) == 1);
wait2=wait(finished_jobs(7,:) == 2);
wait3=wait(finished_jobs(7,:) == 3);

figure
[f,x]=ecdf(wait1*3600,'Function','survivor');
[f2,x2]=ecdf(wait2*3600,'Function','survivor');
[f3,x3]=ecdf(wait3*3600,'Function','survivor');


figure
semilogy(x,f)
hold on
semilogy(x2,f2)
semilogy(x3,f3)

f2=f2([1:60:end-400,end-399:5:end]);
x2=x2([1:60:end-400,end-399:5:end]);

f3=f3([1:10:end-200,end-199:2:end]);
x3=x3([1:10:end-200,end-199:2:end]);

fprintf(fopen('survivor1.dat','w'), '%f %f \n', [x';f']);
fprintf(fopen('survivor2.dat','w'), '%f %f \n', [x2';f2']);
fprintf(fopen('survivor3.dat','w'), '%f %f \n', [x3';f3']);



figure
histogram(wait1*3600,26,'Normalization','pdf')
figure
histogram(wait2*3600,26,'Normalization','pdf')
figure
histogram(wait3*3600,26,'Normalization','pdf')

[N2,edges2]=histcounts(wait2*3600,26,'Normalization','pdf');
[N3,edges3]=histcounts(wait3*3600,edges2,'Normalization','pdf');

[N1,edges1]=histcounts(wait1*3600,edges2,'Normalization','pdf');


fprintf(fopen('distprio.dat','w'), '%f %f %f %f\n', [edges3(1:end-1);N1;N2;N3]);



figure
[N,edges]=histcounts(wait*3600,26,'Normalization','pdf');
fprintf(fopen('dist.dat','w'), '%f %f\n', [edges(1:end-1);N]);
histogram(wait*3600,26,'Normalization','pdf')


job1=finished_jobs(:,finished_jobs(7,:) == 1);
job2=finished_jobs(:,finished_jobs(7,:) == 2);
job3=finished_jobs(:,finished_jobs(7,:) == 3);

%number of primary
distri=zeros(4,4);
distri(1,1)=length(job1(6,job1(6,:)==1))/size(job1,2);
distri(2,1)=length(job1(6,job1(6,:)==2))/size(job1,2);
distri(3,1)=length(job1(6,job1(6,:)==3))/size(job1,2);
distri(4,1)=length(job1(6,job1(6,:)==4))/size(job1,2);

distri(1,2)=length(job2(6,job2(6,:)==1))/size(job2,2);
distri(2,2)=length(job2(6,job2(6,:)==2))/size(job2,2);
distri(3,2)=length(job2(6,job2(6,:)==3))/size(job2,2);
distri(4,2)=length(job2(6,job2(6,:)==4))/size(job2,2);

distri(1,3)=length(job3(6,job3(6,:)==1))/size(job3,2);
distri(2,3)=length(job3(6,job3(6,:)==2))/size(job3,2);
distri(3,3)=length(job3(6,job3(6,:)==3))/size(job3,2);
distri(4,3)=length(job3(6,job3(6,:)==4))/size(job3,2);

distri(1,4)=length(finished_jobs(6,finished_jobs(6,:)==1))/size(finished_jobs,2);
distri(2,4)=length(finished_jobs(6,finished_jobs(6,:)==2))/size(finished_jobs,2);
distri(3,4)=length(finished_jobs(6,finished_jobs(6,:)==3))/size(finished_jobs,2);
distri(4,4)=length(finished_jobs(6,finished_jobs(6,:)==4))/size(finished_jobs,2);



NN=size(job2,2);
NN2=size(job3,2);
%number of interruped
interrupt=zeros(6,4);
interrupt(1,1)=nnz(job2(9,:)==1)/NN;
interrupt(2,1)=nnz(job2(9,:)==2)/NN;
interrupt(3,1)=nnz(job2(9,:)==3)/NN;
interrupt(4,1)=nnz(job2(9,:)==4)/NN;
interrupt(5,1)=nnz(job2(9,:)==5)/NN;
interrupt(6,1)=nnz(job2(9,:)==6)/NN;

interrupt(1,2)=nnz(job2(10,:)==1)/NN;
interrupt(2,2)=nnz(job2(10,:)==2)/NN;
interrupt(3,2)=nnz(job2(10,:)==3)/NN;
interrupt(4,2)=nnz(job2(10,:)==4)/NN;
interrupt(5,2)=nnz(job2(10,:)==5)/NN;
interrupt(6,2)=nnz(job2(10,:)==6)/NN;

interrupt(1,3)=nnz(job3(9,:)==1)/NN2;
interrupt(2,3)=nnz(job3(9,:)==2)/NN2;
interrupt(3,3)=nnz(job3(9,:)==3)/NN2;
interrupt(4,3)=nnz(job3(9,:)==4)/NN2;
interrupt(5,3)=nnz(job3(9,:)==5)/NN2;
interrupt(6,3)=nnz(job3(9,:)==6)/NN2;

interrupt(1,4)=nnz(job3(10,:)==1)/NN2;
interrupt(2,4)=nnz(job3(10,:)==2)/NN2;
interrupt(3,4)=nnz(job3(10,:)==3)/NN2;
interrupt(4,4)=nnz(job3(10,:)==4)/NN2;
interrupt(5,4)=nnz(job3(10,:)==5)/NN2;
interrupt(6,4)=nnz(job3(10,:)==6)/NN2;

interrupt=interrupt*100;



sum(wait > 0)/length(wait)


inter_wait=finished_jobs(5,:)-(finished_jobs(3,:)+finished_jobs(4,:));

inter_wait3=inter_wait(finished_jobs(7,:) == 3);

hej=mean(inter_wait3(job3(10,:)>=1))*3600

nz_inter=inter_wait(:,inter_wait > 10e-12);

mean(nz_inter)*3600

%time spent
job_time=sum(finished_jobs(4,:));

scheduled_time=sum(finished_scheduled(4,:));

total=3*3*8*days+3*2*8*days;

job_time/scheduled_time
job_time/total
scheduled_time/total
1-job_time/total-scheduled_time/total


%Scheduled care
waitq_sched=zeros(16,5);


wait_scheduled=finished_scheduled(3,:)-finished_scheduled(2,:);

wait_sched=zeros(16,4);
sched_sd=zeros(16,4);
k=1;
for t=7:22
    wait_sched(k,1)=mean(wait_scheduled(finished_scheduled(7,:) == 1 & finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1))*3600;
    wait_sched(k,2)=mean(wait_scheduled(finished_scheduled(7,:) == 2 & finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1))*3600;
    wait_sched(k,3)=mean(wait_scheduled(finished_scheduled(7,:) == 3 & finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1))*3600;
    
    wait_sched(k,4)=mean(wait_scheduled(finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1))*3600;
    
    
    waitq_sched(k,:)=quantile(wait_scheduled(finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1),[0.25 0.5 0.75 0.95, 0.99])*3600;    
    
    N=size(wait_scheduled(finished_scheduled(7,:) == 1 & finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1),2);
    sched_sd(k,1)=sqrt(var((wait_scheduled(finished_scheduled(7,:) == 1 & finished_scheduled(2,:) > t  & finished_scheduled(2,:) <= t+1))*3600)/N);
    
    N=size(wait_scheduled(finished_scheduled(7,:) == 2 & finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1),2);
    sched_sd(k,2)=sqrt(var((wait_scheduled(finished_scheduled(7,:) == 2 & finished_scheduled(2,:) > t  & finished_scheduled(2,:) <= t+1))*3600)/N);
    
    N=size(wait_scheduled(finished_scheduled(7,:) == 3 & finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1),2);
    sched_sd(k,3)=sqrt(var((wait_scheduled(finished_scheduled(7,:) == 3 & finished_scheduled(2,:) > t  & finished_scheduled(2,:) <= t+1))*3600)/N);
    
    N=size(wait_scheduled(finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1),2);
    sched_sd(k,4)=sqrt(var((wait_scheduled(finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1))*3600)/N);
    
    k=k+1;
end

p2=[7:22;wait_sched(:,1)';wait_sched(:,1)'-norminv(0.995)*sched_sd(:,1)';wait_sched(:,1)'+norminv(0.995)*sched_sd(:,1)'];

p2( p2 < 0) =0;

fprintf(fopen('schedq1.dat','w'), '%f %f %f %f %f %f\n', [7:22; (waitq_sched/60)']);
fprintf(fopen('waitsched1.dat','w'), '%f %f %f %f\n', p2);

fprintf(fopen('waitsched2.dat','w'), '%f %f %f %f\n', [7:22;wait_sched(:,2)';wait_sched(:,2)'-norminv(0.995)*sched_sd(:,2)';wait_sched(:,2)'+norminv(0.995)*sched_sd(:,2)']);

fprintf(fopen('waitsched3.dat','w'), '%f %f %f %f\n', [7:22;wait_sched(:,3)';wait_sched(:,3)'-norminv(0.995)*sched_sd(:,3)';wait_sched(:,3)'+norminv(0.995)*sched_sd(:,3)']);

waitsched=[7:22;(wait_sched(:,4)');(wait_sched(:,4)'-norminv(0.995)*sched_sd(:,4)');(wait_sched(:,4)'+norminv(0.995)*sched_sd(:,4)')];

waitschedmin=[7:22;(wait_sched(:,4)')/60;(wait_sched(:,4)'-norminv(0.995)*sched_sd(:,4)')/60;(wait_sched(:,4)'+norminv(0.995)*sched_sd(:,4)')/60];

% waitsched(waitsched < 0) = 0;

fprintf(fopen('waitsched.dat','w'), '%f %f %f %f\n', waitsched);


fprintf(fopen('waitschedmin.dat','w'), '%f %f %f %f\n', waitschedmin);


sched1=wait_scheduled(finished_scheduled(7,:) == 2 & finished_scheduled(2,:) >= 8  & finished_scheduled(2,:) < 10)*3600;
sched2=wait_scheduled(finished_scheduled(7,:) == 3 & finished_scheduled(2,:) >= 8  & finished_scheduled(2,:) < 10)*3600;


[fs,xs]=ecdf(sched1/60,'Function','survivor');
[fs2,xs2]=ecdf(sched2/60,'Function','survivor');


xs=xs([1:20:end-400,end-399:4:end]);
fs=fs([1:20:end-400,end-399:4:end]);

fs2=fs2([1:10:end-200,end-199:2:end]);
xs2=xs2([1:10:end-200,end-199:2:end]);

fprintf(fopen('survsched1.dat','w'), '%f %f \n', [xs';fs']);
fprintf(fopen('survsched2.dat','w'), '%f %f \n', [xs2';fs2']);



[fs3,xs3]=ecdf(wait_scheduled*60,'Function','survivor');
xs3=xs3([1:80:end-400,end-399:4:end]);
fs3=fs3([1:80:end-400,end-399:4:end]);


fprintf(fopen('basesurvsched.dat','w'), '%f %f \n', [xs3';fs3']);




NN2=size(finished_scheduled,2);
%number of interruped
interrupt2=zeros(6,2);
interrupt2(1,1)=nnz(finished_scheduled(9,:)>=1)/NN2;
interrupt2(2,1)=nnz(finished_scheduled(9,:)>=2)/NN2;
interrupt2(3,1)=nnz(finished_scheduled(9,:)>=3)/NN2;
interrupt2(4,1)=nnz(finished_scheduled(9,:)>=4)/NN2;
interrupt2(5,1)=nnz(finished_scheduled(9,:)>=5)/NN2;
interrupt2(6,1)=nnz(finished_scheduled(9,:)>=6)/NN2;

interrupt2(1,2)=nnz(finished_scheduled(10,:)>=1)/NN2;
interrupt2(2,2)=nnz(finished_scheduled(10,:)>=2)/NN2;
interrupt2(3,2)=nnz(finished_scheduled(10,:)>=3)/NN2;
interrupt2(4,2)=nnz(finished_scheduled(10,:)>=4)/NN2;
interrupt2(5,1)=nnz(finished_scheduled(10,:)>=6)/NN2;
interrupt2(6,1)=nnz(finished_scheduled(10,:)>=6)/NN2;



busy=zeros(16,2);

k=1;
for t=7:22
    busy(k,1)=sum(finished_jobs(4,finished_jobs(2,:) >= t  & finished_jobs(2,:) < t+1));
    busy(k,2)=sum(finished_scheduled(4,finished_scheduled(2,:) >= t  & finished_scheduled(2,:) < t+1));
    k=k+1;
end

busy(1:8,:)=busy(1:8,:)/(3*3*days);
busy(9:16,:)=busy(9:16,:)/(3*2*days);


fprintf(fopen('busy.dat','w'), '%f %f %f\n', [7:22;busy']);


inter_wait2=finished_scheduled(5,:)-(finished_scheduled(3,:)+finished_scheduled(4,:));

nz_inter2=inter_wait2(:,inter_wait2 > 10e-12);

mean(nz_inter2)

test=finished_scheduled(:,inter_wait2 > 10e-12);

w1=(length(wait(wait <= tau + 10e-12 & wait >= tau - 10e-12))/length(wait))*100
w2=(length(wait(wait <= 10e-12))/length(wait))*100



finished_scheduled2=finished_scheduled(:,finished_scheduled(2,:) > 7  & finished_scheduled(2,:) <= 10);

NN2=size(finished_scheduled2,2);
%number of interruped
interrupt2=zeros(6,2);
interrupt2(1,1)=nnz(finished_scheduled2(9,:)>=1)/NN2;
interrupt2(2,1)=nnz(finished_scheduled2(9,:)>=2)/NN2;
interrupt2(3,1)=nnz(finished_scheduled2(9,:)>=3)/NN2;
interrupt2(4,1)=nnz(finished_scheduled2(9,:)>=4)/NN2;
interrupt2(5,1)=nnz(finished_scheduled2(9,:)>=5)/NN2;
interrupt2(6,1)=nnz(finished_scheduled2(9,:)>=6)/NN2;

interrupt2(1,2)=nnz(finished_scheduled2(10,:)>=1)/NN2;
interrupt2(2,2)=nnz(finished_scheduled2(10,:)>=2)/NN2;
interrupt2(3,2)=nnz(finished_scheduled2(10,:)>=3)/NN2;
interrupt2(4,2)=nnz(finished_scheduled2(10,:)>=4)/NN2;
interrupt2(5,1)=nnz(finished_scheduled2(10,:)>=6)/NN2;
interrupt2(6,1)=nnz(finished_scheduled2(10,:)>=6)/NN2;


%%

nocall=zeros(16,3);



k=1;
for t=7:22
    nocall(k,1)=sum(finished_jobs(finished_jobs(1,:)==1 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1));
    nocall(k,2)=sum(finished_jobs(finished_jobs(1,:)==2 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1));
    nocall(k,3)=sum(finished_jobs(finished_jobs(1,:)==3 & finished_jobs(2,:) > t  & finished_jobs(2,:) <= t+1));
    k=k+1;
end

figure
plot(7:22,nocall(:,1))
hold on
plot(7:22,nocall(:,2))
plot(7:22,nocall(:,3))


%% functions
function t = simphasetype(P)
%Function to simulate phase type distribution
    n=size(P,1);
    t=0;
    a=[0;cumsum(P(:,1))];
    S=find(rand(1) > a ,1,'last');
    
    Q=[-P(:,2:end)*ones(n,1),P(:,2:end)];
 
    P=bsxfun(@rdivide,Q,-diag(Q,1)); 
    P(n+1:n+1:end)=0;
 
     while S ~= 0 
        a=[0,cumsum(P(S,:))];  
        t=t+exprnd(-1/Q((n+1)*S));
        
        S=find(rand(1) > a,1,'last')-1;   
     end
end


function t = generate(lambda,maxLambda,t)
%Function to simulate the time-dependent poisson process.
    p=0;
    while rand(1) > p
         t=t+exprnd(1/maxLambda);
         p=lambda(ceil(mod(t,24)))/maxLambda;
    end
end
