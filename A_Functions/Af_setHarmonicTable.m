function HarmonicTable = Af_setHarmonicTable(nB)
% map FAST channel to my harmonic naming scheme
% nB is number of blades

HarmonicTable = ...
    {'RootMyc1',    'Flap', 1;
     'RootMxb1',    'Edge', 1;
     'LSSTipMya',   'Msy',  1;
     'LSSTipMza',   'Msz',  1;
     'LSSTipMys',   'Mmy',  nB;
     'LSSTipMzs',   'Mmz',  nB;
     'YawBrMyp',    'Myy',  nB;
     'YawBrMzp',    'Myz',  nB;
     'TwrBsMxt',    'Mtx',  nB;
     'TwrBsMyt',    'Mty',  nB};
 
 