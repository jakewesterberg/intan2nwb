# intan2nwb
 Transforms intan recording session to distributable NWB files.
 1. Extracts recording device(s)
 2. Extracts Analog I/O data
 3. Extracts Digital I/O data
 4. Binarizes raw neural data
 5. Filters raw data to LFP (1-250 Hz) and MUA (|500-5000 Hz|)
 6. Kilosorts binary data

# requirements
 1. matnwb (forked version from jake)
 2. npy-matlab
 3. Kilosort (2.0, preferred atm, jake's forked version has config files)
  a. Visual Studio C++ compiler (Community 2019 confirmed working)
 4. ecephys_spike_sorting (from Allen Institute) 
  a. requires py env with v3.7
  b. argschema==1.17.5
  c. marshmallow==2.19.2
  d. Anaconda 3.0 works well for managing this...

