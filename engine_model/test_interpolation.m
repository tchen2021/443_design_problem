clc;clear;close all
%torque_percent = interpolateTorque(31000, 0)
%[eta, D] = interpolate_eta(280, 2000, 25000)
[percent_power_avi, percent_power_total] = engine_model(25000,280,0,0.0311,224)
