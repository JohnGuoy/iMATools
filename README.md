
## iMATools
&emsp;&emsp;We developed an integrated DNA methylation pattern region identification and annotation platform iMATools (integrated Methylation Analysis Tools) based on long-read sequencing data, which is used for encoding processing, format conversion, methylation in ultra-large-scale long-read methylation information It aims to provide methylation researchers with professional methylation analysis and visualization tools for long-read sequencing, and accurately reveal DNA methylation patterns at the cellular and read levels.

&emsp;&emsp;iMATools has three modules:
* **towig** - methylation level file to Wiggle format
* **pattern** - identification of methylation patterns  of genomic reigons
* **mrv** - visualization of CpG sites in specific regions of the genome

### Workflow

### Install
You first need to install Python v3.8+ and Perl v5.16+, then iMATools can be used directly after decompressing. 
```
unzip iMATools-main.zip
```

For the mrv tool software, you need to install some dependency packages first:
```
pip install portion
pip install tqdm
pip install matplotlib
```

### Manual

* These are simple examples.

__Usage:__ Convert "H1_bismark.cov" into "wig" format. Methy counts is in col 5,unmethy counts is in col 6.
```shell
towig -t bismark -n H1 -r hg19.fa -i H1_bt2.bismark.cov -o H1_wig
```
__Usage:__ Identify methylation patterns from "H1_wig/".
```shell
pattern -i H1_wig/ -o H1_pattern/ -n H1
```
__Usage:__ Visualizing methylation of CpG sites in specific regions of Y chromosome third-generation sequencing data
```shell
python mrv.py --data-file ./Y10895.txt --chromosome Y --cpg-range [10084283,10090100]
```
Output `Y_10084283_10087018_visualization.svg` file, you can open this file with Chrome browser, where black dots represent methylated CpG sites and white dots represent unmethylated CpG sites:
![ CpG ranges [10084283,10087018] of Y chromosome ](https://github.com/JohnGuoy/iMATools/blob/main/test_data/Y_10084283_10087018_visualization.jpg)

### Using Tips

1 If you use PBS(Portable Batch System) in your cluster server, **avoid to appoint relative path** for `-o,--outdir` and other parameters which need to assign path because workspace will be changed when pbs file is submitted. 

2 wiggle format
More detail information in [UCSC Genome Browser: Wiggle Track Format (WIG)](http://genome.ucsc.edu/goldenPath/help/wiggle.html).

3 construction information
iMATools is contructed in `Python v3.8` and `Perl v5.16.3`. 
and iMATools is tested in Python v3.8 and Perl v5.16.3. 

4 dependence relationship
towig is independent. Input could come from `BSMAP`,`Bismark` or ENCODE, Roadmap, TCGA.
pattern is independent. 
mrv requires libraries of shutil, portion, matplotlib and tqdm. 

### A work we did with iMATools platform
We used the towig and pattern tool software in the iMATools platform to mine the intermediate methylation pattern regions in the ONT sequencing data of normal breast cells. The average methylation level of this methylation pattern region is about 0.5, which can be found on the IGV Genome Browser. It is intuitive to see that the methylation level of this region is at an intermediate level.

We also used the mrv tool software in the iMATools platform to visualize and analyze the intermediate methylation pattern regions of imprinted genes and non-imprinted genes at the read level, and found the CpG sites of the reads in the pattern region. Whether or not they are methylated varies greatly in overall form, and we speculate that imprinted and non-imprinted genes differ in the formation of intermediate methylation pattern regions.

Visualization of the intermediate methylation pattern region (IMR) and its reads adjacent to the imprinted gene **FAM50B**:
![FAM50B](https://github.com/JohnGuoy/iMATools/blob/main/test_data/FAM50B.png)

It can be seen intuitively that the multiple reads in this region can be clearly divided into two categories. One category of reads has almost all methylated CpG sites (the black dots in the figure represent methylated reads). CpG sites, white dots represent unmethylated CpG sites), while another class of reads has almost no CpG sites unmethylated. This is the obvious feature of allele-specific methylation regions, and we speculate that the intermediate methylation pattern regions of imprinted genes are caused by their allelic differential methylation. As for the reason why a CpG site on a read in the visualization is not fully methylated or not methylated at all, it may be due to errors in third-generation sequencing.

Visualization of intermediate methylation pattern regions (IMRs) and their reads adjacent to the non-imprinted gene **FRMD6**:
![FRMD6](https://github.com/JohnGuoy/iMATools/blob/main/test_data/FRMD6.png)

The intermediate methylation pattern regions of non-imprinted genes cannot be clearly divided into two categories. We therefore speculate that the intermediate methylation pattern regions of common genes are phenotypes of intermediate states in their demethylation processes, rather than caused by allelic differential methylation.

Studies have shown that FAM50B (family with sequence similarity 50 member B) is a protein-coding gene. Diseases associated with the FAM50B gene include Temple syndrome and gestational trophoblastic tumor. The function of the FAM50B gene may be regulated by allele-specific methylation, and the method used in this section provides an alternative tool for the analysis and validation of this hypothesis.

The ONT sequencing data file of normal breast cells used in this work, R143-0107.filt.fq.gz, can be obtained from the person in charge of the Institute of Biomedical Big Data, Eye Hospital Affiliated to Wenzhou Medical University.
