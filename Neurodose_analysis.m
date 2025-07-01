% Analysis start for Neurodose SocSciSurvey data: SDT Type I and II,
% empiricist index on intelligence test data!

% add the analysis script folder to MATLAB path, so it knows where to
% find it. Change this to your folder!
clear; restoredefaultpath
scriptsfolder = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/Teaching/24-25SS/Fopra/1DecisionStyle/analysis';
addpath(scriptsfolder)

% go to the folder where the data is at. Change this to your folder!
datafolder = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/Teaching/24-25SS/Fopra/1DecisionStyle/data/Neurodose/';
cd(datafolder) 
mkdir('preproc') % make a folder for the preprocessed data
cd raw

t = readtable('NEU Neurodose ScoSci Survey (101 VP)'); % read raw data into table TODO get correct answers!
missedtrials = string(t.key_resp1_keys) == ""; % empty string for no response
t = t(not(missedtrials), :); % remove missed trials from raw data

signalpresence = string(t.condition_type) == "with_circles"; % 0 for target absent, 1 for target present
response = string(t.key_resp1_keys) == "d"; % d was pressed for reporting target presence

out_table = table(); % make a new table with the output variables of interest
out_table.VPcode = unique(string(t.VP_Code_bitteEINGEBEN___1_ErsteBuchstabeDesVornamensDerMutter_z));
out_table.nmissedresponses = sum(missedtrials); % how often did the subject not press? Keep track
if out_table.nmissedresponses > 30
  warning('More than 30 missed trials! Reomve this subject?')
end
out_table.accuracy = mean(signalpresence == response); % proportion correct

% compute dprime and criterion Type I
hitrate = sum(signalpresence == true & response == true) / sum(signalpresence == true); % N hits / N signal presence trials
farate = sum(signalpresence == false & response == true) / sum(signalpresence == false);
hitrate = max(min(hitrate, 1 - eps), eps); % Avoid divide by zero or infinity due to norminv computation  
farate = max(min(farate, 1 - eps), eps);

out_table.TypeI_dprime = norminv(hitrate) - norminv(farate); % Zscored H – FA rates
out_table.TypeI_criterion = -0.5 * (norminv(hitrate) + norminv(farate)); % Zscored H + FA times -0.5

% compute dprime Type II and confidence
% Step 1: Create variables for accuracy and confidence
t.accuracy = signalpresence == response; % response is in line with stimulus shown
t.confidence(string(t.key_resp_5_keys) == "d") = "High"; % Assuming button d was high confidence
t.confidence(string(t.key_resp_5_keys) == "k") = "Low";
% Step 2: Create contingency table for high confidence
high_conf_correct = sum(t.confidence == "High" & t.accuracy); % High confidence, correct responses
high_conf_incorrect = sum(t.confidence == "High" & not(t.accuracy)); % High confidence, incorrect responses
% Step 3: Create contingency table for low confidence
low_conf_correct = sum(t.confidence == "Low" & t.accuracy); % Low confidence, correct responses
low_conf_incorrect = sum(t.confidence == "Low" & not(t.accuracy)); % Low confidence, incorrect responses
% Step 4: Calculate hit and false alarm rates for Type II
type_II_hit_rate = high_conf_correct / (high_conf_correct + low_conf_correct); % Type II hit rate
type_II_false_alarm_rate = high_conf_incorrect / (high_conf_incorrect + low_conf_incorrect); % Type II false alarm rate
% Step 5: Avoid divide by zero or infinity due to norminv computation  
type_II_hit_rate = max(min(type_II_hit_rate, 1 - eps), eps);
type_II_false_alarm_rate = max(min(type_II_false_alarm_rate, 1 - eps), eps);

out_table.TypeII_dprime = norminv(type_II_hit_rate) - norminv(type_II_false_alarm_rate);
out_table.confidence = sum(t.confidence == "High") / height(t); % confidence 0 to 1, 1 means only high confidence responses
out_table.empiricist_index = (1 - out_table.confidence) * out_table.TypeII_dprime; 
% Explanation: Empiricist will have high TypeII dprime and also low confidence. The 1 - confidence will boost the score
% for a subject who does just that. Therefore higher EI score means more
% empiricist
disp(out_table)
% save table as csv to the preproc folder, use VPcode as the filename
writetable(out_table, fullfile(datafolder, 'preproc', append(out_table.VPcode, ".csv")))