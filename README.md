
## iBSTools
&emsp;&emsp;iBStools(integrated Bisulfite Sequencing Tools) is an integrated tools for comprehensive analysis of bisulfite sequencing reads including whole genome bisulfite sequencing(WGBS) and reduced representation bisulfite sequencing (RRBS). 
>iBSTools is unpublished.

&emsp;&emsp;iBSTools has four modules:
* **towig** - methylation level file to Wiggle format
* **pattern** - identification of methylation patterns  of genomic reigons
* **refumr** - identification of reference methyalted regions among mutiple methylomes.
* **dmr** - identification of differentially methylated regions among two groups

![workflow](https://github.com/methylation/iBSTools/blob/master/imgs/sketch_map.png "foo")

--
### Install
iBSTools can be used directly after decompressing. 
```
unzip iBSTools-master.zip
```

--
### Manual

* These are simple examples, more details please read the [iBSTools wiki](https://github.com/methylation/iBSTools/wiki)

__Usage:__ Convert "H1_bismark.cov" into "wig" format. Methy counts is in col 5,unmethy counts is in col 6.
```shell
towig -t bismark -n H1 -r hg19.fa -i H1_bt2.bismark.cov -o H1_wig
```
__Usage:__ Identify methylation patterns from "H1_wig/".
```shell
pattern -i H1_wig/ -o H1_pattern/ -n H1
```
__Usage:__ Identify reference methylation patterns regions from mutiple methylomes in "wig_list.txt".
```shell
refumr -p UM -path ./software/iBSTools_v1.1.0/ -w wig_list.txt -o ref_UM
```
__Usage:__ Identify differentially methylated regions for a specific genomic regions between two groups of methylomes.
```shell
dmr -r ref_UM/ref_UM.bed -rh 1 -w1 file_list_1.txt  -w2 file_list_2.txt -o diff
```

--
### Using Tips

1. If you use PBS(Portable Batch System) in your cluster server, **avoid to appoint relative path** for `-o,--outdir` and other parameters which need to assign path because workspace will be changed when pbs file is submitted. 

2. wiggle format
More detail information in [UCSC Genome Browser: Wiggle Track Format (WIG)](http://genome.ucsc.edu/goldenPath/help/wiggle.html).

3. construction information
iBSTools is contructed in `R3.1.3` and `perl v5.16.3`. 
and iBSTools is tested in R2.x and perl v5.10.x... 
iBSTools just employs basic funtions in `R` and `perl`. So almost all of versions of R and perl is available.

4. dependence relationship
towig is independent. Input could come from `BSMAP`,`Bismark` or ENCODE, Roadmap, TCGA.
pattern is independent. 
refumr requires pattern. 

--
#### Realease history
More details please read the [Realease History](https://github.com/methylation/iBSTools/blob/master/REALEASE_HISTORY.md)

