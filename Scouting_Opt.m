clear, clc, close all
load("VineData.mat")
load("EnvironmentalForcing.mat")
findSwitch = 0;
cost = 0;
loops = 20;
minCostAll = 1000000*ones(1,loops);
optAmtAll = zeros(1,loops);
optSpdAll = zeros(1,loops);
optTFoundAll = zeros(1,loops);
for i = 1:loops
    disp(i)
    for amt = 12:-1:1
        for speed =.005:.005:.5
            for t = 1:length(tspan) 
                if mod(t,24) == 0 && findSwitch ==0
                    [infects,infectsFound] = Scouting(speed,amt,vine,t,2);
                    cost = cost + amt*100;
                    if (t-1)/24 > 10
                        cost = cost + 1000;
                    end
                    if infectsFound == 1 && findSwitch == 0
                        tFound = t;
                        findSwitch = 1;
    %                     disp('Infection Found')
                    end
                end
            end
            if cost<minCostAll(i)
                minCostAll(i) = cost;
                optAmtAll(i) = amt;
                optSpdAll(i) = speed;
                optTFoundAll(i) = tFound;
            end
            cost = 0;
            findSwitch = 0;
        end
    end
end
% disp(optTFound/24)
meanMinCost = mean(minCostAll)
meanOptAmt = mean(optAmtAll)
meanOptSpd = mean(optSpdAll)
meanTFound = mean(optTFoundAll)
% figure
% scatter(optAmtAll,optSpdAll)


function [infects,infectsFound] = Scouting(speed,amt,vine,t,opts)
    NpX = 50;
    NpY = 50;
    infects = zeros(NpX,NpY);
    infectsFound = 0;
    DetectSize = (20*speed/10)^2/4*pi/5000;
    distMax = speed*3600;
    distUsed = 0;
%     if opts == 1
%         gridSize = floor(sqrt(Npnts));
%         for i = 1:gridSize:NpX
%             for j = 1:gridSize:NpY
    if opts == 2
        for a = 1:amt
            currLoc = [0,0];
            while distUsed < distMax && infectsFound ~= 1
                RandSearch = randi(NpX*NpY);
                %fprintf('day:%i (%i,%i)\n',round(t/24),vine(RandSearch).X+.5,vine(RandSearch).Y+.5)
                distUsed = distUsed + sqrt((vine(RandSearch).X - currLoc(1))^2 + (vine(RandSearch).Y - currLoc(2))^2);
                if distUsed > distMax
                    break
                end
                if vine(RandSearch).I(t) >= DetectSize
                    infects(vine(RandSearch).X+0.5,vine(RandSearch).Y+0.5) = 1;
                    infectsFound = 1;
                    return
                end
                currLoc = [vine(RandSearch).X,vine(RandSearch).Y];

            end
        end
    end
end
