function array_history= return_history_days(thisday, ndays)
% Return ndays weekends/weekdays before YYYY/MM/DD
if nargin==1
    ndays = 30;
end
% thisday = datetime(YYYY, MM, DD);
if mod(weekday(thisday), 6) == 1
    todayisweekday = false;
else
    todayisweekday = true;
end
% array_history = nan(ndays, 1);
% 
% validnumber = 0;
% historyday  = thisday;
% while validnumber < ndays
%     historyday = historyday - days(1);
%     if mod(weekday(historyday), 6) == 1
%         historyisweekday = false;
%     else
%         historyisweekday = true;
%     end
%     if historyisweekday == todayisweekday
%         validnumber = validnumber + 1;
%         array_history(validnumber) = historyday;
%     end
% end
nweeks = ceil(ndays/2);

thisday_utc = datetime(thisday.Year, thisday.Month, thisday.Day, 'TimeZone', 'UTC'); % Damn daylight saving time! Convert everything to UTC

candidate_days = thisday_utc-days(nweeks*7):days(1):thisday_utc-days(1);
candidate_days = flipud(candidate_days');
candidate_isweekday = ~(mod(weekday(candidate_days(:)), 6) == 1);
array_history = candidate_days(candidate_isweekday==todayisweekday);
array_history = array_history(1:ndays);

% Now convert back to the original timezone
array_history = datetime(array_history.Year, array_history.Month, array_history.Day, 'TimeZone', thisday.TimeZone);
end
