@echo OFF
set CONDAPATH=C:\Users\westerja\Anaconda3
set ENVNAME=ecephys
if %ENVNAME%==base (set ENVPATH=%CONDAPATH%) else (set ENVPATH=%CONDAPATH%\envs\%ENVNAME%)
call %CONDAPATH%\Scripts\activate.bat %ENVPATH%

set GIT_PYTHON_REFRESH=quiet
set PYTHONIOENCODING=utf-8

cd C:\Users\westerja\Documents\GitHub\intan2nwb\forked_toolboxes\ecephys_spike_sorting

python -m ecephys_spike_sorting.modules.kilosort_postprocessing --input_json C:\Users\westerja\Desktop\ecephys_test\ecephys_spike_sorting_input.json --output_json C:\Users\westerja\Desktop\ecephys_test\ecephys_spike_sorting_output.json

python -m ecephys_spike_sorting.modules.mean_waveforms --input_json C:\Users\westerja\Desktop\ecephys_test\ecephys_spike_sorting_input.json --output_json C:\Users\westerja\Desktop\ecephys_test\ecephys_spike_sorting_output.json

python -m ecephys_spike_sorting.modules.noise_templates --input_json C:\Users\westerja\Desktop\ecephys_test\ecephys_spike_sorting_input.json --output_json C:\Users\westerja\Desktop\ecephys_test\ecephys_spike_sorting_noise_output.json

python -m ecephys_spike_sorting.modules.quality_metrics --input_json C:\Users\westerja\Desktop\ecephys_test\ecephys_spike_sorting_input.json --output_json C:\Users\westerja\Desktop\ecephys_test\ecephys_spike_sorting_quality_output.json

call conda deactivate