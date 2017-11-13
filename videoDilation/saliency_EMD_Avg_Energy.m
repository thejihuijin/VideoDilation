function [ energy_graph ] = saliency_EMD_Avg_Energy( saliencyFrames )
%SALIENCYENERGY 
% Takes in a series of saliency maps and computes the
% energy as the earth mover's distance from the average saliency map. It 
% then returns the normalized energy graph from 0 to 1

% Read sizes
[~, ~, nFrames] = size(saliencyFrames);

% Generate CDF of Saliency Frames
edges = .025:.05:1;

avgsalmap = mean(saliencyFrames,3);
avg_hist = hist(avgsalmap(:),edges)/numel(avgsalmap);
avg_cdf = cumsum(avg_hist);
energy_graph = zeros([1 nFrames]);
for ii=1:nFrames
    salIm = saliencyFrames(:,:,ii);
    N = hist(salIm(:),edges)/numel(salIm);
    
    X = cumsum(N);
    energy_graph(ii) = sum(abs(X-avg_cdf));
end
mov_avg_window = 15;
energy_graph = smooth(energy_graph,mov_avg_window);
energy_graph = energy_graph - min(energy_graph);
energy_graph = energy_graph/max(energy_graph);
end

