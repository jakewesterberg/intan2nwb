function I = photodiodeTrigger(photodiode_events, coded_events, expected_window_size)

if nargin < 3
    expected_window_size = [-1/60*3 1/60*3]; % assuming refresh rate 60Hz
end

if size(coded_events,2) > 1
    coded_events = coded_events.';
end

if size(photodiode_events,2) == 1
    photodiode_events = photodiode_events.';
end

d = nan(numel(coded_events), 1);
idx = nan(numel(coded_events), 1);

edges = [-Inf, mean([photodiode_events(2:end); photodiode_events(1:end-1)]), +Inf];
I = discretize(coded_events, edges);

mistriggers = sum((photodiode_events(I) - coded_events) < expected_window_size(1)) | ...
      (photodiode_events(I) - coded_events) > expected_window_size(2);

if mistriggers > 0
    warning([num2str(mistriggers) ' photodiode alignments exceed the concern window.'])
end

end