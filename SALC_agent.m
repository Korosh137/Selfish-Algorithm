% The Selfish Algorithm (SA)
% Korosh Mahmoodi 111618
% Cite: Mahmoodi, Korosh, Bruce J. West, and Cleotilde Gonzalez. "Selfish algorithm and emergence of collective intelligence." Journal of Complex Networks 8.3 (2020): cnaa019.
tic
clc ;
clear all ;
close all ;

% In SALC, there are two decisions for each agent to make:
% not-Connect or Connect (1 or 2)
% Defect or Cooperate (1 or 2)
% For each choice there is a corresponding adaptive propensity/probability

Trial = 1e6 ; % Simulaton length
Size = 20 ; % Number of the agents

% Constants relating payoff (as feedback) to the change of the propensity/probability of the corresponding decision made
Chi_Connection = 0.1 ; % SAC mechanism (connect decision)
Chi_RL = 0.1 ; %  SAL mechanism (Cooperation or Defection decision)

% Parameters of the Prisoner's Dilema game's
s = 0 ; % s is the payoff of the cooperator agent if another agent defected
p = 0 ; % p is the payoff of the agents if both defected
tc = 0.5 ; % temptation to cheat. (1+ tc) is the payoff of the defector agent if the other agent cooperated (tc < s + 1 and  tc < 1)

Ratio_CC = zeros(Trial, 1) ; CC = 0 ; % Ratio of Cooperation-Cooperation decisions

Out1 = zeros(Size, 1) ; % Previous payoff of agents
Out2 = zeros(Size, 1) ; % Current payoff of agents

D_C_decision = zeros(Size, 1) ; % D (defection) as 1 and C (cooperation) as 2

D_C_decisionColor = zeros(Trial, Size) ; % Assigns color for the nodes at each trial

P_Connection = zeros(Size, Size, 2) ; % Propensity of the pari of agents to connect/play
P_RL = zeros(Size, Size, 2) ; % For agent m paired with agent n the P_RL(m, n, 1) and P_RL(m, n, 2) are the propensites of decison D and C, respectively

% Initial conditions
for jjj = 1 : Size
    r = rand ;
    if r < 0.5
        D_C_decision(jjj) = 1 ;
    else
        D_C_decision(jjj) = 2 ;
    end

end
for i = 1 : Size
    for j = 1 : Size
        for k = 1 : 2 % number of choices/propensities
            if i ~= j
                P_Connection(i, j, k) = 0.5 ;
                P_RL(i, j, k) = 0.5 ;
            end
        end
    end
end

Connection_decision = ones(Size, Size) ; % Not-connct as 1 and connect as 2

gh =  1 ;

for ti = 2 : Trial

    P_Connection_Used = ones(Size, 1) ; % Tracks if the Connection decision mechanism has used (2) in the final decision making or not (1)
    P_RL_Used = ones(Size, 1) ; % Tracks if the RL decision mechanism has used (2) in the final decision making or not (1)

    % Selecting two agents
    m = floor(1 + Size * rand) ;
    n = floor(1 + Size * rand) ;

    while m == n
        n = floor(1 + Size * rand) ;
    end

    % Connection decisions of the two agents
    Connection_decision(m, n) = Decision(P_Connection(m, n, :)) ;
    Connection_decision(n, m) = Decision(P_Connection(n, m, :)) ;

    while Connection_decision(m,n) == 1 || Connection_decision(n, m) == 1

        % Selecting two agents
        m = floor(1 + Size * rand) ;
        n = floor(1 + Size * rand) ;

        while m == n
            n = floor(1 + Size * rand) ;
        end

        % Connection decisions
        Connection_decision(m, n) = Decision(P_Connection(m, n, :)) ;
        Connection_decision(n, m) = Decision(P_Connection(n, m, :)) ;
    end

    P_Connection_Used(m) = 2 ; % i.e. connection decision mechanism has used
    P_Connection_Used(n) = 2 ;

    % Decision of C or D
    D_C_decision(m) = Decision(P_RL(m, n, :)) ;
    D_C_decision(n) = Decision(P_RL(n, m, :)) ;

    P_RL_Used(m) = 2 ; % i.e. RL decision mechanism has used
    P_RL_Used(n) = 2 ;

    % update the color of the nodes in the demo
    D_C_decisionColor(ti, :) =  D_C_decisionColor(ti-1, :) ;
    D_C_decisionColor(ti, m) =  D_C_decision(m) ;
    D_C_decisionColor(ti, n) =  D_C_decision(n) ;

    % Payoffs from the Prisoner's Dilemma game
    Out2(m) = Payoff_PD(D_C_decision(m), D_C_decision(n), s, p, tc) ;
    Out2(n) = Payoff_PD(D_C_decision(n), D_C_decision(m), s, p, tc) ;

    % Updating the propensities
    PConnectionm = zeros(1, 2) ;
    PConnectionn = zeros(1, 2) ;
    PRLm = zeros(1, 2) ;
    PRLn = zeros(1, 2) ;

    for k = 1 : 2
        PConnectionm(k) = P_Connection(m, n, k) ;
        PConnectionn(k) = P_Connection(n, m, k) ;

        PRLm(k) = P_RL(m, n, k) ;
        PRLn(k) = P_RL(n, m, k) ;
    end

    P_Connection(m, n, :) = SA_Update(P_Connection_Used(m), PConnectionm, Connection_decision(m,n), Out2(m) , Out1(m) , Chi_Connection) ;
    P_Connection(n, m, :) = SA_Update(P_Connection_Used(n), PConnectionn, Connection_decision(n,m), Out2(n) , Out1(n) , Chi_Connection) ;

    P_RL(m, n, :) = SA_Update(P_RL_Used(m), PRLm , D_C_decision(m), Out2(m), Out1(m), Chi_RL) ;
    P_RL(n, m, :) = SA_Update(P_RL_Used(n), PRLn , D_C_decision(n), Out2(n), Out1(n), Chi_RL) ;

    % Updating the previous payoffs
    Out1(m) = Out2(m) ;
    Out1(n) = Out2(n) ;

    if  D_C_decision(m) == 2  &&  D_C_decision(n) == 2
        CC = CC + 1 ;
    end
    Ratio_CC(ti) = CC/ ti ;

end   % end trial


plot(Ratio_CC)

hold off

toc

% Decision making
function D = Decision(Prop)

g = length(Prop) ;
CS =   cumsum(Prop) ;

r = rand ;
for i = 1 : g
    if r < CS(i)
        D = i ;
        break
    end
end
end

% Payoff of the Prisoner's Dilemma game
function Payoff = Payoff_PD(DC1, DC2, s, p, tc)

if DC2 == 2
    NC = 1 ;
else
    NC = 0 ;
end

if DC1 == 2
    Payoff = NC * (1) + (1 - NC) * (-s) ;
else
    Payoff = NC * (1 + tc) + (1 - NC) * (p) ;
end

end

% Updating the propensities
function P = SA_Update(MechanismUsed, pp, Dec, Pay, Paybefore, Chi)

% This function updates the propensities of the decision mechanisms which contributed to the final decision of the agent
% MechanismUsed: 2 if the decision mechanism is used and 1 otherwise
% pp: The input propensities
% Dec: Indicates the index of the propensity that resulted in the decision
% Pay: Payoff of the agent at the current trial
% Paybefore: Payoff of the agent at the previous trial
% Chi: Constants relating payoff to the change of the propensity/probability

if MechanismUsed == 2 % If the decision mechanism has used (2) in the final decision making or not (1)
    g = length(pp) ;
    P = zeros(1 , g) ;

    Delta =  Chi * (Pay - Paybefore) / (abs(Pay) + abs(Paybefore)) ;

    if   Pay == Paybefore
        P = pp ;
    else

        P(Dec) = pp(Dec) + Delta ;

        for i = 1 : g
            if i ~= Dec
                P(i) = pp(i) - Delta/(g-1) ;
            end
        end

        % % % Boundary condition
        for j = 1 : g
            if  P(j) <= 0
                P(j) =  1e-10 ;
            end
            if  P(j) >= 1
                P(j) =  1-1e-10 ;
            end
        end

        % % % Normalization
        R = cumsum(P);
        P = P ./ R(g) ;
    end

else
    P = pp ;
end

end





