This software is developed to calculate the activity history of transposable elements based on the UPGMA algorithm.
This can be used for all kinds of repetitive sequences.
It is recommended to use repetitive sequences with <1kb.

There will be four results files;
1. A blast result file (xml format)
2. A parsed blast result (csv format)
3. A UPGMA phylogenetic tree (nwk format, it can be opened with MEGA software)
4. Raw data for activity history (csv format, tagged with "_UPGMA_distance")
   => The raw data for activity history can be used to draw a histogram to visualize the activity history.


Prerequisites;
1. Any version of local blast search software.
2. Python 3.x
All the prerequisite software has to be in PATH.

Usage;
To see options;
python TE_activity_calculator.py.txt -help

Example;
python TE_activity_calculator.py.txt -fasta test_sequence.txt -score 60


This software is developed by Minkyu Park, University of Georgia Athens.
If there is any question, please send an email to minkju3003@gmail.com
1/13/2021