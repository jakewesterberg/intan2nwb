%%
% this script shows how to use the microsaccade detection method
% published in Otero-Millan et al. Journal of Vision 2014

% Set up variables --------------------------------------------------------
folder = fileparts(mfilename('fullpath'));
if ( isempty( folder) )
    folder = pwd;
end
folder = [folder '\data'];


% =========================================================================
% EDIT THIS ===============================================================
% Make sure to fill the variables samples with the eye movement data in the
% proper colums. In blinks mark all the samples that are not good data. In
% most VOG systems it is best to remove 100 ms before and after the blink
% to remove all artifacts.
% =========================================================================
    session = 'test';
    samplerate = 500;

%     samples = [];
%     %  samples(:,1)     timestamps of the recording in miliseconds
%     %  samples(:,2)     horizontal position of the left eye in degrees
%     %  samples(:,3)     vertical position of the left eye in degrees  
%     %  samples(:,4)     horizontal position of the right eye in degrees     
%     %  samples(:,5)     vertical position of the right eye in degrees  
% 
%     blinks = [];
%     %  blinks           binary vector indicating for each sample if it
    %                   belons to a blink (1) or not (0)
% =========================================================================
% END EDIT THIS ===========================================================
% =========================================================================

% Loads the recording and prepares it por processing
recording = ClusterDetection.EyeMovRecording.Create(folder, session, samples, blinks, samplerate);

% Runs the saccade detection
[saccades stats] = recording.FindSaccades();

% Plots a main sequence
enum = ClusterDetection.SaccadeDetector.GetEnum;
figure
subplot(2,2,1)
plot(saccades(:,enum.amplitude),saccades(:,enum.peakVelocity),'o')
set(gca,'xlim',[0 1],'ylim',[0 100]);
xlabel('Saccade amplitude (deg)');
ylabel('Saccade peak velocity (deg/s)');


% Plots the traces with the labeled microsaccades
subplot(2,2,[3:4])
plot(samples(:,1), samples(:,2:end));
hold
yl = get(gca,'ylim');
u1= zeros(size(samples(:,1)))+yl(1);
u2= zeros(size(samples(:,1)))+yl(1);
u1((saccades(:,enum.startIndex))) = yl(2);
u2(saccades(:,enum.endIndex)) = yl(2);
u = cumsum(u1)-cumsum(u2);
plot(samples(:,1), u,'k')

xlabel('Time (ms)');
ylabel('Eye Position (deg)');

legend({'Left Horiz', 'Left Vert', 'Right Horiz' , 'Right Vert', 'Microsaccades'})
