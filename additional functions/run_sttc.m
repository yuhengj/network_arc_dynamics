% STCC derivation, based on Cutts and Eglen (2014)

function index = run_sttc(N1, N2, dt, Time, spike_times_1, spike_times_2)

if N1==0 || N2==0
    index=NaN;
    
else
    T=Time(2)-Time(1);
    TA=run_T(N1,dt,Time(1),Time(2), spike_times_1);
    TA=TA/T;
    TB=run_T(N2,dt,Time(1),Time(2), spike_times_2);
    TB=TB/T;
    PA=run_P(N1,N2,dt, spike_times_1, spike_times_2);
    PA=PA/N1;
    PB=run_P(N2,N1,dt, spike_times_2, spike_times_1);
    PB=PB/N2;
    index=0.5*(PA-TB)/(1-TB*PA)+0.5*(PB-TA)/(1-TA*PB);
    
end



    function Nab = run_P(N1,N2,dt,spike_times_1,spike_times_2)
        Nab=0;
        j = 1;
        for i=1:N1
            while j<=N2
                %check every spike in train 1 to see if there's a spike in train 2 within dt  (don't count spike pairs)
                %don't need to search all j each iteration
                
                if abs(spike_times_1(i)-spike_times_2(j))<= dt
                    Nab=Nab+1;
                    break;
                    
                elseif spike_times_2(j)>spike_times_1(i)
                    break;
                else
                    j=j+1;
                    
                end
            end
        end
        
    end

    function time_A = run_T(N1,dt, startv, endv, spike_times_1)
        i=1;
        time_A=2*N1*dt;
        
        if N1==1
            if(spike_times_1(1)-startv)<dt
                time_A=time_A-startv+spike_times_1(1)-dt;
                
            elseif(spike_times_1(1)+dt)>endv
                time_A=time_A-spike_times_1(1)-dt+endv;
            end
            
        else
            while i<(N1)
                
                diff=spike_times_1(i+1)-spike_times_1(i);
                
                if(diff<2*dt)
                    %//subtract overlap
                    time_A = time_A-2*dt+diff;
                end
                i = i +1;
            end
            
            %check if spikes are within dt of the startv and/or end, if so just need to subract
            %overlap of first and/or last spike as all within-train overlaps have been accounted for
            
            if (spike_times_1(1)-startv)<dt
                
                time_A=time_A-startv+spike_times_1(1)-dt;
                
            end
            
            if(endv-spike_times_1(N1))<dt
                
                time_A = time_A-spike_times_1(N1)-dt+endv;
            end
            
            
        end
        
    end
end