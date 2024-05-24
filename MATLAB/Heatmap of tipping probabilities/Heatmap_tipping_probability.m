%% Heatmap_Tipping
%RRC Sept 2023

%input data as matrix, in rows from top column
probability = readmatrix('C:\Users\rc686\OneDrive - University of Exeter\RuthPeteRichard GRL paper\Scripts\RRC\MC simulations\tipping_probabilities_decadal_all.txt');

%reshape the matrix of probabilities
probability = reshape(probability,9,8);
probability = flipud(probability);

%write out the array of probabilities
xvalues = {'Zero','u01 100yrs','u01 1000yrs','u03 20yrs','u03 50yrs','u03 70yrs','u03 100yrs','u03 1000yrs'};
yvalues = {'20 \times HadGEM3MM','15 \times HadGEM3MM','10 \times HadGEM3MM','5 \times HadGEM3MM','MPI','CanESM5','HadGEM3MM','HadGEM3LL','0'};

figure(1);
h = heatmap(xvalues, yvalues, probability);
