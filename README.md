# bromecast
Data processing and visualization for landscape genomics project of cheatgrass (*Bromus tectorum*). 

# Directory Stucture

```txt
finish updating this structure!!
bromecast/
│   
├── 01_data/
│   │   ├── data
                             # Test set 
│   └──  processing_scripts/
│       ├── raw_train_data               
│      ├── clean                        # Train set (trimmed 1-second recordings)
│      ├── spects                       # Mel-spectrograms of cleaned samples                          
│      ├── crop_wav_files.ipynb         # Crop categorized .wav files into 1-second clips
│      └── get_high_rated_recordings.R  # Extract recordings with a 3-5 star rating
│   
│   └── covariates/                      
│       ├── grid_code.R                  # Create a regular grid across region of interest with extracted covariate information for predictions of bioacoustic compositions
│       ├── make_datasets.R              # Make datasets of observed calls and associated covariate information
│       └── make_seasonal_data.R         # Make datasets corresponding to the spring of each year analyzed (2020-2023)
│
├── 02_train_and_predict/              
│   ├── train.ipynb                      # Train the model on classified samples
│   └── predict.ipynb                    # Run predictions using trained models on out-of-sample recordings
│
└── 03_summarize/                      
    ├── eta_plot.jpeg                            # Plot of eta coefficients
    ├── plot_etas.R                              # R script to create eta_plot.jpeg
    ├── spect_samples.jpg                        # Example spectrograms used in modeling (Fig. 2)
    ├── spect_fig_and_sample_nums.ipynb          # Plot random sample of spectrograms for each call type (Fig. 2) 
    ├── spatial_prediction_maps.jpeg             # Visual summary of spatial predictions (Fig. 5)
    ├── prediction_map.R                         # R script to generate prediction maps (Fig. 5)
    ├── covs_across_years.jpeg                   # Precipitation and temperature trends from 2020-2023 (Fig. 6) 
    └── plot_avg_temp_and_precip_across_years.R  # R script for temp/precip trends (Fig. 6)

```
