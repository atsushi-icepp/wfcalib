Simulation codes for SiPM gain calibration study.
"phFunc.h" and "phrec.cpp" are code for generating generalized Poisson distribution.

"Waveform.h" is a header file for waveform simulation. Waveform generation method is defined and implemented in this file.

"StatGainCalib.cpp"(caluculation) and "StatConfig.xml"(configration) are the code for simulating "statistical method" by generating waveform. After execution, "fout.root" file is generated and the simulation result is written in the file. "stat_draw.cpp" can be used to analyze the result in detail.

"wf_simulation.cpp" and "SimConfig.xml" are the codes for simulating waveform method. These codes generate waveform and simulated results are written into "fout.root". "draw_wfsim.cpp" can be used to analyze the result in detail.
