%% Vorwort
% Es wurden mehr Daten ausgelesen als wir am Ende verwendet haben, das hat
% aber beim Verständnis der Ausreißer geholfen.

%% Analyse
% Analysis start for Neurodose: SDT Type I and II, empiricist index
% Detection task, intelligence test
% add the analysis script folder to MATLAB path, so it knows where to
% find it. Change this to your folder!
clear; restoredefaultpath
scriptsfolder = '';
addpath(scriptsfolder)

% go to the folder where the data is at. Change this to your folder!
datafolder = '';
cd(datafolder)
mkdir('preproc') % make a folder for the preprocessed data
files = dir(fullfile(datafolder, '*.csv'));
big_table = table(); % make a new table with the output variables of interest

% loop over data files
for i = 1:length(files)
    disp(files(i).name);

    t = readtable(fullfile(datafolder, files(i).name)); % read raw data into table
    missedtrials = string(t.key_resp1_keys) == ""; % empty string for no response
    missedtrials2 = string(t.gdms) == "";

    t1 = t(not(missedtrials), :); % remove missed trials from raw data
    t1 = t1(11:90, :);          % data table of the 80 sdt trials without the training trials
    t2 = t(not(missedtrials2), :);      % data table of the gdms questionnaire

    signalpresence = string(t1.condition_type) == "with_circles"; % 0 for target absent, 1 for target present
    response = string(t1.key_resp1_keys) == "d"; % d was pressed for reporting target presence

    
    out_table = table(); % make a new table with the output variables of interest

    % read out the socioeconomic information
    out_table.Geschlecht = t.key_resp_12_keys(3);
    out_table.Bildung = t.key_resp_13_keys(4);
    out_table.Psych = t.key_resp_14_keys(5);

    % read out the mean gdms subscale scores 
    out_table.Rational = mean(t2.slider_response([1,6,11,16,21]));
    out_table.Intuitive = mean(t2.slider_response([2,7,12,17,22]));
    out_table.Dependent = mean(t2.slider_response([3,8,13,18,23]));
    out_table.Spontaneous = mean(t2.slider_response([5,10,15, 20, 25]));
    out_table.Avoidant = mean(t2.slider_response([4,9,14,19,24]));

    % read out the VPCodes being the VPCode or a random number if missing
    vp_col_idx = find(contains(t1.Properties.VariableNames, 'VP_Code'));

    if ismissing(t1.VP_Code_eingeben_FallsGew_nscht_1_DerErsteBuchstabeDesVornamens(1))
        out_table.VPCode = randi(999999);
    else
        out_table.VPCode = string(t1.VP_Code_eingeben_FallsGew_nscht_1_DerErsteBuchstabeDesVornamens(1));
    end

    % read out the age being the age or 0 if missing
    if ismissing(t.Alter_bitteEingeben__(1))
        out_table.Alter = 0;
    else
        out_table.Alter = t.Alter_bitteEingeben__(1);
    end
   

    out_table.filename = string(files(i).name);
    out_table.nmissedresponses = sum(missedtrials); % how often did the subject not press? Keep track
    if out_table.nmissedresponses > 54
        warning('More than 10 missed trials! Reomve this subject?')
    end
    out_table.accuracy = mean(signalpresence == response); % proportion correct
    
    

    
    % compute dprime and criterion Type I using log-linear correction for
    % extreme proportions see Hautus (1995)
    hitrate = (sum(signalpresence == true & response == true) + 0.5) / (sum(signalpresence == true) + 1); % N hits / N signal presence trials
    farate = (sum(signalpresence == false & response == true) + 0.5) / (sum(signalpresence == false) + 1);
    
    
    hitrate = max(min(hitrate, 1 - eps), eps); % Avoid divide by zero or infinity due to norminv computation
    farate = max(min(farate, 1 - eps), eps);
    



    % read out the TypeI sdt measures
    out_table.TypeI_dprime = norminv(hitrate) - norminv(farate); % Zscored H – FA rates
    out_table.TypeI_criterion = -0.5 * (norminv(hitrate) + norminv(farate)); % Zscored H + FA times -0.5

    % read out Reaction time
    out_table.StimulusRT = mean(t1.key_resp1_rt);
    out_table.ConfidenceRT = mean(t1.key_resp_5_rt);
    out_table.SliderRT = mean(t2.slider_rt);



    % compute dprime Type II and confidence
    % Step 1: Create variables for accuracy and confidence
    t1.accuracy = signalpresence == response; % response is in line with stimulus shown
    t1.confidence(string(t1.key_resp_5_keys) == "d") = "High"; % Assuming button d was high confidence
    t1.confidence(string(t1.key_resp_5_keys) == "k") = "Low";
    % Step 2: Create contingency table for high confidence
    high_conf_correct = sum(t1.confidence == "High" & t1.accuracy); % High confidence, correct responses
    high_conf_incorrect = sum(t1.confidence == "High" & not(t1.accuracy)); % High confidence, incorrect responses
    % Step 3: Create contingency table for low confidence
    low_conf_correct = sum(t1.confidence == "Low" & t1.accuracy); % Low confidence, correct responses
    low_conf_incorrect = sum(t1.confidence == "Low" & not(t1.accuracy)); % Low confidence, incorrect responses
    % Step 4: Calculate hit and false alarm rates for Type II
    type_II_hit_rate = (high_conf_correct + 0.5)/ (high_conf_correct + low_conf_correct + 1); % Type II hit rate
    type_II_false_alarm_rate = (high_conf_incorrect + 0.5) / (high_conf_incorrect + low_conf_incorrect + 1); % Type II false alarm rate
    
    
    % Step 5: Avoid divide by zero or infinity due to norminv computation
    type_II_hit_rate = max(min(type_II_hit_rate, 1 - eps), eps);
    type_II_false_alarm_rate = max(min(type_II_false_alarm_rate, 1 - eps), eps);



    out_table.TypeII_dprime = norminv(type_II_hit_rate) - norminv(type_II_false_alarm_rate); % AKA typeII-dprime, see Fleming and Lau (2014)More actions
    out_table.TypeII_criterion = -0.5 * (norminv(type_II_hit_rate) + norminv(type_II_false_alarm_rate));
    out_table.metacog_efficiency = out_table.TypeII_dprime / out_table.TypeI_dprime; % Meta-dprime normalized by objective performance, see Fleming and Lau (2014)
    out_table.confidence = sum(t1.confidence == "High") / height(t1); % confidence 0 to 1, 1 means only high confidence responses
    out_table.empiricist_index = (1 - out_table.confidence) * out_table.TypeII_dprime;
    % Explanation: Empiricist will have high TypeII dprime and also low confidence. The 1 - confidence will boost the score
    % for a subject who does just that. Therefore higher EI score means more
    % empiricist
    % for a subject who does just that. Therefore higher EI score means more empiricistMore actions
    out_table.empiricist_index_efficiency = (1 - out_table.confidence) * out_table.metacog_efficiency;
    % Explanation: metacog_efficiency instead of TypeII_dprime, to control for
    % objective sensitivity differences between subjects.
    out_table.empiricist_index_efficiency_log = (1 - out_table.confidence) * -log(out_table.metacog_efficiency);
    % Explanation:  when the denominator (d′) is small, meta-d′/d′ can give rather extreme values which may undermine
    % power in a groupstatistical analysis. However, this problem can also be addressedby taking log of meta- d′/d′,
    % as is often done to correct for the non-normality of ratio measures
    % (Howell, 2009). we take -log to make the values positive.
    


    % Read out the gdms items
    out_table.Rational1 = mean(t2.slider_response(1));
    out_table.Rational2 = mean(t2.slider_response(6));
    out_table.Rational3 = mean(t2.slider_response(11));
    out_table.Rational4 = mean(t2.slider_response(16));
    out_table.Rational5 = mean(t2.slider_response(21));

    out_table.Spontaneous1 = mean(t2.slider_response(5));
    out_table.Spontaneous2 = mean(t2.slider_response(10));
    out_table.Spontaneous3 = mean(t2.slider_response(15));
    out_table.Spontaneous4 = mean(t2.slider_response(20));
    out_table.Spontaneous5 = mean(t2.slider_response(25));

    out_table.Intuitive1 = mean(t2.slider_response(2));
    out_table.Intuitive2 = mean(t2.slider_response(7));
    out_table.Intuitive3 = mean(t2.slider_response(12));
    out_table.Intuitive4 = mean(t2.slider_response(17));
    out_table.Intuitive5 = mean(t2.slider_response(22));

    out_table.Dependent1 = mean(t2.slider_response(3));
    out_table.Dependent2 = mean(t2.slider_response(8));
    out_table.Dependent3 = mean(t2.slider_response(13));
    out_table.Dependent4 = mean(t2.slider_response(18));
    out_table.Dependent5 = mean(t2.slider_response(23));

    out_table.Avoidant1 = mean(t2.slider_response(4));
    out_table.Avoidant2 = mean(t2.slider_response(9));
    out_table.Avoidant3 = mean(t2.slider_response(14));
    out_table.Avoidant4 = mean(t2.slider_response(19));
    out_table.Avoidant5 = mean(t2.slider_response(24));

    big_table = [big_table; out_table];
end

disp(big_table);
% save table as csv to the preproc folder, use date and time as filename
tablename = string(datetime('now', 'Format', 'd-MMM-y-HH-mm-ss')) + '.csv';
writetable(big_table, fullfile(datafolder, 'preproc', tablename));

