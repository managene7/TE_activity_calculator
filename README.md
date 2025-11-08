<h4>This software is developed to calculate the activity history of transposable elements based on the UPGMA algorithm.</h4>
<h4>This can be used for all kinds of repetitive sequences.</h4>
<h4>It is recommended to use repetitive sequences with <1kb.</h4>

<h4>There will be four results files;</h4>
1. A blast result file (xml format)<br>
2. A parsed blast result (csv format)<br>
3. A UPGMA phylogenetic tree (nwk format, it can be opened with MEGA software)<br>
4. Raw data for activity history (csv format, tagged with "_UPGMA_distance")<br>
   => The raw data for activity history can be used to draw a histogram to visualize the activity history.<br>


<h4>Prerequisites;</h4>
1. Any version of local blast search software.<br>
2. Python 3.x<br>
All the prerequisite software has to be in PATH.<br>

<h4>Usage;</h4>
To see options;

```
python TE_activity_calculator.py.txt -help
```

<h4>Example;</h4>

```
python TE_activity_calculator.py.txt -fasta test_sequence.txt -score 60
```

This software is developed by Minkyu Park, University of Georgia Athens.<br>
If there is any question, please send an email to minkju3003@gmail.com<br>
1/13/2021
