Code for "Global translation speed affects native folding of multi-domain proteins" Tippmann et al 2025, journal

List of scripts:

0_assignment.py (footprint assignment)<br />
1_CI_calc.py (confidence interval calculation for SeRP data sets)<br />
2_onset.py (chaperone engagement onset determination)<br />
3_metagene.py (SeRP metagene calculation)<br />
4_get_AF_data.py (read residue contact information from mmcif files)<br />
5_structure_assisted_domain_annotation.py (combines Af structure derived contact information with CATH server domain search to define protein domains on AF structures)<br />
6_sasa_calc.py (calculate solvent accessible surface areas of AF structures)<br />
7_unsatisfied_residue_calculation.py (determine unsatisfied residues based on AF structures)<br />
8_unsatisfied_inter_intra_domain_residue_calculation.py (determine unsatisfied residues within and between domains based on AF structures)<br />

Files used:

E. coli MG1655 CDS reference (NC_000913.3.gff3)<br />
Alphafold files, https://alphafold.ebi.ac.uk/download (UP000000625_83333_ECOLI_v4)<br />
E. coli MG1655 amino acid sequences<br />
CATH files are obtained using the FunFHMMer web server<br />
PAE files are obtained via web scraping with the python “urllib.request” module<br />
System Requirements / Hardware Requirements: Normal desktop computer, at least 8 GB of RAM. For faster cluster analysis, a multi-core processor is helpful.<br />
Software Requirements: Python 3.9 with the following additional modules: numpy, statsmodels, scipy, matplotlib, freesasa<br />

Installation guide via conda: Install your conda environment following: https://conda.io/projects/conda/en/latest/user-guide/getting-started.html(https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)<br />
Install python and all modules, for example scipy following: https://anaconda.org/anaconda/scipy(https://anaconda.org/anaconda/scipy) Installation usually take 2-3 hours.

Instructions for use How to run the software on your data Activate the conda environment where python and all modules are installed. To test all pythonscripts, the “path” or “main_path” variable should be a string with the path to the demo folder.<br />
Reproductive instructions CATH and Alphafold library queries might change over time due to newly deposited data sets.
