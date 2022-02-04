
## iMATools
&emsp;&emsp;We developed an integrated DNA methylation pattern region identification and annotation platform iMATools (integrated Methylation Analysis Tools) based on long-read sequencing data, which is used for encoding processing, format conversion, methylation in ultra-large-scale long-read methylation information It aims to provide methylation researchers with professional methylation analysis and visualization tools for long-read sequencing, and accurately reveal DNA methylation patterns at the cellular and read levels.

&emsp;&emsp;iMATools has four modules:
* **towig** - methylation level file to Wiggle format
* **pattern** - identification of methylation patterns  of genomic reigons
* **mrv** - visualization of CpG sites in specific regions of the genome

--
### Workflow

--
### Install
You first need to install Ptyon v3.8+ and Perl v5.16+, then iMATools can be used directly after decompressing. 
```
unzip iMATools-master.zip
```

--
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

--
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
