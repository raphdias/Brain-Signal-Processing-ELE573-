
distanceSum = 0; 

for sumIndx=1:numElectrodesWithinRad
    distanceSum = distanceSum + 1/e_within_rad(sumIndx,2); 
end

gij = []; 

for i=1:numElectrodesWithinRad
    gijSingle = (1/e_within_rad(i,2))/(distanceSum);
    gij = [gij gijSingle]; 
end