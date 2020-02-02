# autodada2
Script suite to analyse amplicon sequencing data using the dada2 pipeline. Parameters for certain steps of the pipeline can be edited for specific genes. The default parameters are currently tailored towards ITS2 data.

Authored by Alex Veglia and Lauren Howe-Kerr

This is the master script to perform dada2 anlaysis with amplicon sequencing data

### Setting everything up

Programs needed to run script ->

	1. AdapterRemoval2
	
	2. blast suite(specifically blastn and makeblastdb)
	
	3. reformat.sh from the bbtools kit
	
	4. R (packages will be installed when script is run
	
	5. FastQC

**NOTE:** PATHs to these programs must be set in the autodada2_v1.1.sh script 

```bash 
nano autodada2_v1.1.sh 
```

### General execution

Run the autodada2v1.1.sh executable in the same directory as the read libraries and make sure the supp_scripts file is in the same directory as script, if not, you have to set the path to the directory in the script file 

``` 
autodada2_v1.1.sh [-hrbplaguf;] 2>&1 | tee outputfile.out
```

Options ->

   **[ -h ]** help
   
   **[ -r ]**        		
   Rename read library files. Expected to have *_S* and *_001* -> must be true to have the rename work ok, can be    customized in script though
   
   **[ -a ]**  Remove specified adapters with AdapterRemoval2
   
   **[ -b ]**           		
   Remove primers with BBDuk defaults
   
   **[ -p <f,r|n>]**    		
   Remove primers by choppin off specified amount of basepairs: f (forward) and r (reverse), or n (both)
   
   **[ -l ]**
   Run LULU on asvs, will not run if specified mainly for ITS2
   
   **[ -u ]**
   Set this option to only run LULU, but make sure to have files in LULUdir
   
   **[ -g ]**
   Run the auto phyloseq script to produce phyloseq object and relative abundance plots
   
   **[ -f ]**
   Run ONLY the phyloseq script to produce phyloseq object and realtive abundance plots, need files within $"workdir"/phylo

**Quick note:** make sure to always run from designated working directory; the scripts anchors itself to "'$workdir'"
  

**Examples:**

```autodada2_v1.1.sh 2>&1 | tee outputfile.out```   -> This will just run the dada2 analysis on the fastq files present in workdir; will produce labeled ASV_fasta file, counts table and files needed for LULU and phyloseq.

```autodada2_v1.1.sh -lab 2>&1 | tee outputfile.out```  -> This will remove adapters, remove primers with bbduk, and run LULU analysis on asvs genrated by dada2

```autodada2_v1.1.sh -labg 2>&1 | tee outputfile.out```  -> This will remove adapters, remove primers with bbduk, run LULU AND run phyloseq to generate relative abundance plots

