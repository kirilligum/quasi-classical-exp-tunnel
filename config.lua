
----------------------------------------
-- total system parameters
m=1.0;
nn_coup = 5.0;
n1_shift=0.5;
n2_shift=0.5;
energy_threshold_factor_relative=2.0;
energy_threshold_factor_absolute=1e-1;

----------------------------------------
-- initial conditions
pp=3.0;
qp=-4.0;
n1=1.0;
n2=1.0-n1;
n1_bin_start = 0.05;
n2_bin_start = 0.05;
n1_bin_end = 0.5;
n2_bin_end = 0.5;

----------------------------------------
-- propagation parameters
end_time = 1.5e1;
--end_time = 2e1*1e-15/atomic_time; -- sec to atomic_time
n_times = 1000; 
time_step = end_time/n_times;

----------------------------------------
-- method parameters
seed = 1;
n_trajs = 1;
--n_trajs = 500;
--n_trajs = 2000;

-- dofile "config.lua"
-- print(sys_k2,"\n",sys_ki,"\n",sys_k4,"\n",sys_m,"\n",sys_hb,"\n",s0)
-- sed -i "s@./data/ki_........e...@$(find ./data/ki* -type d)@g" *.plt
-- rm -rf data; time ./cmqhd nano.lua; sed -i "s@./data/ki_........e...@$(find ./data/ki* -type d)@g" *.plt;killall gnuplot; gnuplot -p -e 'load "pqpss.plt"'
--lua -e 'dofile "config.lua"; print( "\n bohr_radius = ",bohr_radius, "\n hartree_energy = ",hartree_energy, "\n ev = ", ev, "\n atomic_velocity = ", atomic_velocity, "\n speed_of_light = ",speed_of_light, "\n w0 = ", w0, "\n end_time = ", end_time,"\n temperature = ",temperature )'
