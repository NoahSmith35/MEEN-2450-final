clear, clc, close all
load("VineData.mat")
load("EnvironmentalForcing.mat")
findSwitch = 0;
cost = 0;
minCost = 1000000;
for amt = 1:12
    for speed =.0004:.0001:.5
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
                    disp('Infection Found')
                end
            end
        end
        if cost<minCost
            minCost = cost
            optAmt = amt
            optSpeed = speed
        end
        cost = 0;
    end
end
disp(tFound/24)


function [infects,infectsFound] = Scouting(speed,amt,vine,t,opts)
    NpX =10;
    NpY = 10;
    infects = zeros(NpX,NpY);
    infectsFound = 0;
    DetectSize = (20*speed/10)^2/4*pi/5000;
    Npnts = floor(speed*3600);
    distMax = speed*3600;
    distUsed = 0;
    currLoc = [0,0];
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
