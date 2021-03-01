function [oversupply, pshort] = calculate_oversupply_and_pshort(actual, forecast)
% oversupply: MWh of total FRP oversupply
% pshort: Probability of FRP shortage

if size(actual, 1) ~= size(forecast, 1)
    error('Height does not match!');
    oversupply = nan;
    pshort = nan;
else
    forecast(forecast<0) = 0;
    difference = forecast - actual;
    nshort_15min = difference<0;
    nshort_hourly = sum(nshort_15min, 2);
    nshort = sum(nshort_hourly); 
    pshort = nshort/numel(actual);
    
    need_hourly = max(actual, [], 2);
    need_hourly(need_hourly<0) = 0;
    oversupply_hourly = forecast - need_hourly;
    oversupply_hourly(oversupply_hourly<0) = 0;
    oversupply = sum(oversupply_hourly);
end


end