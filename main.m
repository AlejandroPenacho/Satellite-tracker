
addpath(genpath("."));

gs_data = [0; 0];
particles = [10000; 1000; 0];
time = 0;

data = get_observations(gs_data, particles, time)