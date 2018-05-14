function [Statematrix,t] = MarkovModel0(N,edges,Hurricane,neighborWeight,maxtimesteps,startingstate)
%MARKOVMODEL by Molly Rose Kelly-Gorham given a network and stresses on
%each point Hurricane occurs for 3 timesteps and then is over the program
%goes through a while loop until the maxtimesteps is reached or all the
%nodes have recovered.

statematrix = zeros(maxtimesteps + 1,N);
statematrix(1,:) = startingstate;
stress = 1-exp(-Hurricane);
t = 1;
while t < 4 || ((sum(statematrix(t,:)) ~=0)  &&  t<=maxtimesteps)
    state = statematrix(t,:);
    for n = 1:N
        neighborsdown = sum(state(edges{n}))/length(edges{n});
        % markov because probability depends on current state
        if ~state(n)
            probfail = stress(n)*(1-neighborWeight) + neighborWeight*(neighborsdown);
            statematrix(t+1,n) = (rand(1) <= probfail);
        else
            if sum(stress) == 0
                probrecover = 1-((stress(n))*(1-neighborWeight) + neighborWeight*(neighborsdown));
            else 
                probrecover = 0;
            end
            statematrix(t+1,n) = (rand >= probrecover); % false means return to operational state
        end
    end
    if t == 3
        stress = zeros(length(Hurricane),1);
    end
    t = t+1;
end
Statematrix = statematrix(1:t,:);
end

