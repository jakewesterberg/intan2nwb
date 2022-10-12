Noise Templates
==============
Identifies "noise" units based on template shape

[Kilosort2](https://github.com/MouseLand/Kilosort2) generates templates of a fixed length (2 ms) that matches the time coures of an extracellularly detected spike waveform. However, there are no constraints on template shape, which means that the algorithm often fits templates to voltage fluctuations that could not physically result from the current flow associated with an action potential. The units associated with these templates are considered "noise," and must be filtered out prior to analysis. This is true for other spike sorters as well, but the characteristics of the noise waveforms may be highly algorithm-dependent.

This module contains code for two different approaches to noise template identification:

(1) `id_noise_templates()` uses a variety of heuristics to find units with abnormal spatial spread (single channel or >25 channels), abnormal shape (no peak and trough), or multiple spatial peaks. These are based on many observations of typical noise template shapes from Neuropixels recordings in cortex, hippocampus, thalamus, and midbrain. The appropriate heuristics will likely need to be updated for different types of electrodes or different brain regions.

(2)  `id_noise_templates_rf()` uses a [random forest classifier](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html) trained on manually annotated templates. A pickle file containing the classifier object is included in this repository. A PyQt-based app (`template_classifier_app.py`) is available if you'd like train your own classifier.

Because there's so much variation in the shape of noise templates, we've found it hard to get the false negative rate down to zero with either approach (i.e., there are always some obvious noise units that pass through). Therefore, we still need a manual curation step to remove the remaining noise units. Any suggestions for how to improve the classifier's performance are welcome.

Running
-------
```
python -m ecephys_spike_sorting.modules.noise_templates --input_json <path to input json> --output_json <path to output json>
```
Two arguments must be included:
1. The location of an existing file in JSON format containing a list of paths and parameters.
2. The location to write a file in JSON format containing information generated by the module while it was run.

See the `_schemas.py` file for detailed information about the contents of the input JSON.


Input data
----------
- **Kilosort outputs** : includes spike times, spike clusters, templates, etc.


Output data
-----------
- **cluster_group.tsv** : labels for each cluster in spike_clusters.npy