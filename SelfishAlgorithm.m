% The Selfish Algorithm (SA)
% Korosh Mahmoodi 111618
% Cite: Mahmoodi, Korosh, Bruce J. West, and Cleotilde Gonzalez. "Selfish algorithm and emergence of collective intelligence." Journal of Complex Networks 8.3 (2020): cnaa019.
tic
clc ;
clear all ;
close all ;
% In SA, there are three decisions for each agent to make:
% Connect or not-Connect
% Defect or Cooperate
% not-Trust or Trsut
% For each decision, there is a moving threshold that divides the interval [0 1] into two parts, P and 1-P, where P is the propensity of the corresponding decision 
                                        Trial = 1e5 ; % Simulaton length
                                        Size = 10 ; % Number of the agents

% Constant relating payoff to the change of the propensity/probability of
Chi_Connection = 0.1 ; % SAC mechanism (connect/play decision)
Chi_RL = 0.1 ; %  SAL mechanism (Cooperation or Defection decision)
Chi_Trust = 0.1 ; % SAT mechanism (trust decision)
% You can deactivate SAT, SAC, or both mechanisms by setting the corresponding Chi constant(s) to zero (for deactivating SAT must set its initial condition to zero)

DemoStart = Trial - 1 ; % Time for the start of the demo
Demoend = Trial ;

% Prisoner's Dilema game's payoffs
s = 0 ; % s is the payoff of the cooperator agent if another agent defected 
p = 0 ; % p is the payoff of the agents if both defected
tt = 0.9 ; % (1+ tt) is the payoff of the defector agent if the other agent cooperated (tt < s + 1)
Cost = 0 ; % cost of cooperation

mov = VideoWriter('SA.avi') ;
open(mov)
Ratio_CC = zeros(Trial, 1) ; CC = 0 ; % Ration of mutual cooperation

Out1 = zeros(Size, 1) ; % Previous payoff
Out2 = zeros(Size, 1) ; % Current payoff

D_C_decision = zeros(Size, 1) ; % D (defection) as 1 and C (cooperation) as 2
D_C_decisionColor = zeros(Trial, Size) ; % Assigns color for the nodes at each trial

Connection =  zeros(Size, Size) ; % Tendensy of the agents to connect/play
P_RL = zeros(Size, Size) ; % P_RL and (1-P_RL) are the propensity of decison D and C, respectively
P_Trust = zeros(Size, Size) ; % P_Trust and (1-P_Trust) are the propensity of decison "not to trust" and "trust" the decison of another agent, respectively

Conect2ty = zeros(Trial, Size, Size) ; % Collects the Connection matrix of each trial

% Initial conditions
for jjj = 1 : Size
    r = rand ;
    if r < 0.5
        D_C_decision(jjj) = 1 ; % 1 as defect (D) decision
    else
        D_C_decision(jjj) = 2 ; % 2 as cooperation (C) decision
    end
end
for jjj = 1 : Size
    for iii = 1 : Size
        if iii ~= jjj
            P_RL(jjj , iii) = 0.5 ;
            P_Trust(jjj, iii) =  0.5 ;
            Connection(jjj, iii) = 1 ;
        end
    end
end

gh =  1 ;
for ti = 2 : Trial

    % Randomly selecting two agents to connect/play
    m = floor(1 + Size * rand) ;
    n = floor(1 + Size * rand) ;

    while m == n
        n = floor(1 + Size * rand) ;
    end

    Suum1 = cumsum(Connection(m , :)) ;
    P1_Connection = Connection(m , n) / Suum1(Size) ; % Propensity of agent m to connect/play with agent n

    Suum2 = cumsum(Connection(n , :));
    P2_Connection = Connection(n , m) / Suum2(Size) ;

    r1 = rand ;
    r2 = rand ;

    while r1  >  P1_Connection  ||  r2  >  P2_Connection % If the two agents decided not to connect/play, two other agents get selectd

        m = floor(1 + Size * rand) ;
        n = floor(1 + Size * rand) ;

        while m == n
            n = floor(1 + Size * rand) ;
        end

        Suum1 = cumsum(Connection(m , :)) ;
        P1_Connection = Connection(m , n) / Suum1(Size) ;

        Suum2 = cumsum(Connection(n , :)) ;
        P2_Connection = Connection(n , m) / Suum2(Size) ;

        r1 = rand ;
        r2 = rand ;
    end

    r = rand ; 
    if r      >     P_RL(m , n)  % Decision of C or D
        D_C_decision(m) = 2 ;
    else
        D_C_decision(m) = 1 ;
    end

    r = rand;
    if r       >       P_RL(n, m)
        D_C_decision(n) = 2 ;
    else
        D_C_decision(n) = 1 ;
    end

    %  Decision to trust or not the decision of the other agent 
    D_C_m = D_C_decision(m) ;
    D_C_n = D_C_decision(n) ;

    Trustm = 1 ;
    r = rand ;
    if r > P_Trust(m , n)
        Trustm = 2 ;
        D_C_decision(m) = D_C_n ;
    end

    Trustn = 1 ;
    r = rand ;
    if r  > P_Trust(n , m)
        Trustn = 2 ;
        D_C_decision(n) = D_C_m ;
    end
    % update the color of the nodes
    D_C_decisionColor(ti, :) =  D_C_decisionColor(ti-1, :) ;
    D_C_decisionColor(ti, m) =  D_C_decision(m) ;
    D_C_decisionColor(ti, n) =  D_C_decision(n) ;

    % Payoffs from the Prisoner's Dilemma game
    if D_C_decision(m) == 2 
        ggn = 1 ;
    else
        ggn = 0 ;
    end

    if D_C_decision(n) == 2
        ggm = 1 ;
    else
        ggm = 0 ;
    end

    if D_C_decision(m) == 2
        Out2(m) = ggm * (1 - Cost) + (1 - ggm) * (-s) ;
    else
        Out2(m) = ggm * (1 + tt) + (1 - ggm) * (p) ;
    end

    if D_C_decision(n) == 2
        Out2(n) = ggn * (1 - Cost) + (1 - ggn) * (-s) ;
    else
        Out2(n) = ggn * (1 + tt) + (1 - ggn) * (p) ;
    end

    %   Updating the thresholds
    Connection(m, n) = UpdateSA(Connection(m, n) , 1 , Out2(m) , Out1(m) , Chi_Connection) ;
    Connection(n , m) = UpdateSA(Connection(n , m) , 1 , Out2(n) , Out1(n) , Chi_Connection) ;
   
    if Trustm == 1  % Only update P_RL if the agent didn't Trust
        P_RL(m, n) = UpdateSA(P_RL(m, n) , D_C_decision(m), Out2(m), Out1(m), Chi_RL) ;
    end
    if Trustn == 1
        P_RL(n, m) = UpdateSA(P_RL(n, m) , D_C_decision(n)  , Out2(n) , Out1(n), Chi_RL) ;
    end
 
    P_Trust(m, n) = UpdateSA(P_Trust(m, n), Trustm , Out2(m), Out1(m), Chi_Trust) ;
    P_Trust(n, m) = UpdateSA(P_Trust(n, m), Trustn , Out2(n), Out1(n), Chi_Trust) ;

    Out1(m) = Out2(m) ;
    Out1(n) = Out2(n) ;

    if  D_C_decision(m) == 2  &&  D_C_decision(n) == 2
        CC = CC + 1 ;
    end
    Ratio_CC(ti) = CC/ ti ;

    gh = gh + 1 ;
    Conect2ty(gh, :, :) = Conect2ty(gh-1, :, :) ;
    Conect2ty(gh, m, n) = Connection(m , n) ;
    Conect2ty(gh, n, m) = Connection(n , m) ;

end   % end trial

%  Demo
for ty = DemoStart : Demoend

            weights(:, :) = Conect2ty(ty, :, :) ; % The intensity of the lines represents the propensity of the agents to connect/play with one another

    for tn  = 1 : Size
           Summ = cumsum(weights(tn, :)) ;
           weights(tn, :) = weights(tn, :) / Summ(Size) ;
    end

    G = digraph(weights) ;
    LWidths = 1*G.Edges.Weight ;

    plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths)
    intensityValue = LWidths  ;
    OOOO =   intensityValue  ;
    hhh = plot(G,'LineWidth',LWidths) ;
    set(gcf,'color','w') ;
    cccc = 0 ;
    for uuu = 1 : Size
        axis off
        for vvv = 1 : Size
            if uuu ~= vvv
                cccc =  cccc + 1 ;
                if   D_C_decisionColor(ty ,vvv) == 2
                    Col = 'g' ; % Cooperator
                else
                    Col = 'r' ; % Defector
                end
                highlight(hhh, uuu, vvv, 'EdgeColor', [ (1 -OOOO(cccc, 1))  (1-OOOO(cccc, 1))  (1-OOOO(cccc, 1)) ], 'LineWidth', 2, 'MarkerSize', 8)  ;  
                highlight(hhh,  vvv, 'NodeColor', Col)  ;
            end
        end
    end
    frame = getframe(gcf) ;
    writeVideo(mov,frame)
end
hold off

close(mov)
toc

function aa = UpdateSA(pi0, Decision, Pay, Paybefore, Chi)

if Decision == 2   % Decision (Cooperated, Trusted, or not-Connected/Played )
    dp = -1 * Chi ;
else
    dp = Chi ;
end

if   Pay == Paybefore
    aa = pi0 ;
else
    aa = pi0 + dp * (Pay - Paybefore) / (abs(Pay) + abs(Paybefore)) ;
end

% Boundary condition
if  aa < 0
    aa = 1e-2 ;
end
if  aa > 1
    aa = 1-1e-2 ;
end
end

