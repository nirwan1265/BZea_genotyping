cd /Users/nirwantandukar/Documents/Github/BZea_genotyping/scripts/HMM_introgression                                                                                                                      
                                                                                                                                                                                                           
  # Copy scripts to server, then run:                                                                                                                                                                      
  source /opt/anaconda3/etc/profile.d/conda.sh                                                                                                                                                             
  conda activate hmm                                                                                                                                                                                       
                                                                                                                                                                                                           
  export HMM_PY="$(pwd)/hmm_viterbi_bc2s3.py"                                                                                                                                                              
  export MAPDIR="$(pwd)/genetic_maps"                                                                                                                                                                      
  export POST_R="$(pwd)/batch_introgression_analysis.R"                                                                                                                                                    
  export OUTROOT="$(pwd)/optuna_tuning_runs"                                                                                                                                                               
  export MAX_GAP_KB=10000                                                                                                                                                                                  
  export MIN_BLOCK_KB=5000                                                                                                                                                                                 
  export DROP_FRAC=0.10                                                                                                                                                                                    
  export K_REPS=3                                                                                                                                                                                          
                                                                                                                                                                                                           
  # Run with your calibration samples                                                                                                                                                                      
  python tune_hmm_optuna.py calib_samples.txt 100 


#   Input Files:                                                                                                                                                                                             
                                                                                                                                                                                                           
#   1. GL files (genotype likelihoods) - e.g., PN9_SID849.chr1-10.GL.filtered.tsv.gz                                                                                                                         
#   2. Genetic maps - chr1.map through chr10.map (in genetic_maps/ folder)                                                                                                                                   
#   3. Calibration sample list - calib_samples.txt (paths to GL files)                                                                                                                                       
                                                                                                                                                                                                           
#   Scripts Used:                                                                                                                                                                                            
                                                                                                                                                                                                           
#   1. tune_hmm_optuna.py - Main optimization script                                                                                                                                                         
#   2. hmm_viterbi_bc2s3.py - HMM that produces .statepath.tsv.gz files                                                                                                                                      
#   3. batch_introgression_analysis.R - Post-processing that produces .complete_blocks.bed files                                                                                                             
                                                                                                                                                                                                           
#   Yes, the introgression blocks are used!                                                                                                                                                                  
                                                                                                                                                                                                           
#   The optimization pipeline:                                                                                                                                                                               
#   GL file → HMM → statepath.tsv.gz → R script → complete_blocks.bed                                                                                                                                        
#                                                         ↓                                                                                                                                                  
#                                                 Used to compute:                                                                                                                                           
#                                                 - Jaccard similarity                                                                                                                                       
#                                                 - Breakpoints per Mb                                                                                                                                       
#                                                 - Het% guardrail                                                                                                                                           
                                                                                                                                                                                                           
#   The .complete_blocks.bed files are compared between baseline and thinned replicates to measure stability (Jaccard similarity of donor intervals).                                                        
                                                                                                                                                                                                           
#   Summary of Required Files:                                                                                                                                                                               
                                                                                                                                                                                                           
#   scripts/HMM_introgression/                                                                                                                                                                               
#   ├── tune_hmm_optuna.py          # Optimizer                                                                                                                                                              
#   ├── hmm_viterbi_bc2s3.py        # HMM script                                                                                                                                                             
#   ├── batch_introgression_analysis.R  # Post-analysis                                                                                                                                                      
#   ├── genetic_maps/                                                                                                                                                                                        
#   │   ├── chr1.map ... chr10.map  # Genetic maps                                                                                                                                                           
#   ├── calib_samples.txt           # List of GL file paths                                                                                                                                                  
#   └── (your GL files or paths to them)               