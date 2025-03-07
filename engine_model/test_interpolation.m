torque_percent = interpolateTorque(31000, 0)
[eta, D] = interpolate_eta(280, 2000, 25000)
[percent_power_req] = engine_model(31000,280,0,0.035,224)
