# LFP_CSD
A few MATLAB scripts for analysing multi-channel LFP data recorded using TDT systems. 

In order for looping to work, data folders should be organised as: **AnimalTank/PositionNumber/RecordingBlocks**.

The analysis pipeline involves the following steps:

## 1) LFP_CSD.m
Organises and saves the data as LFP.mat

## 2) CSD_Plotting.m 
Plots the CSD and allows manual picking of granular layer limits.

## 3) LFP_Plotting_LayerAv.m 
Uses the granular layer limits and LFP data to plot LFPs for each layer.

## 4) Oscillatory_LFP.m 
Uses the LFP data to look at oscillatory activity across frequency bands. 

To-do list:
- [ ] Check all required functions are also in the repository
