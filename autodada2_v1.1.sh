#!/bin/bash

#title: "Autodada master script"
#author: "Alex Veglia and Lauren Howe-Kerr"
#date: "January 2020"
##
workdir="`pwd`"
SCRIPTDIR=""$workdir"/supp_scripts"
sampleset=""				#basename for files created (renamed at end)
db=""		#Path to directory holding blastn database
dbname=""					#Name of blastn database
AdapterRemoval=""					#Patho to AdapterRemoval exicutable
reformat=""				#Path to reformat.sh exicutable from bbtools
blastn=""						#Path to blastn script
makeblastdb=""					#Path to makeblastdb
python="python##"					#version of python to call 
bbduk=""
renamerr=""$SCRIPTDIR"/rename_rr.sh"			#this should not have to be changed
meta=""							##path to "metadata" file with /PATH/TO/FILE/FILENAME
fastqc=""
usage () {
echo "

This is the master script to perform dada2 anlaysis with amplicon sequencing data

General exicution:

autodada2_v1.1.sh [-hrbplagu;] 2>&1 | tee outputfile.out

||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   [ -h ]           		help
   [ -r ]           		Rename read library files. Expected to have *_S* and *_001* -> must be true to have the rename work ok, can be customized in script though
   [ -a ]                       Remove specified adapters with AdapterRemoval2
   [ -b ]           		Remove primers with BBDuk defaults
   [ -p <f,r|n>]    		Remove primers by choppin off specified amount of basepairs: f (forward) and r (reverse), or n (both)
   [ -l ]			Run LULU on asvs, will not run if specified mainly for ITS2
   [ -u ]                       Set this option to only run LULU, but make sure to have files in LULUdir
   [ -g ]			Run the auto phylo script to produce phyloseq object and relative abundance plots
   [ -f ]			Run ONLY the pyhlo script to produce phyloseq object and realtive abundance plots, need files within $"workdir"/phylo

Quick note: make sure to always run from designated working directory; the scripts anchors itself to "'$workdir'"
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  
  
Example scripts:

autodada2_v1.1.sh 2>&1 | tee outputfile.out   -> This will just run the dada2 analysis on the fastq files present in workdir; will produce labeled ASV_fasta file, counts table and files needed for LULU and phyloseq.

autodada2_v1.1.sh -lab 2>&1 | tee outputfile.out  -> This will remove adapters, remove primers with bbduk, and run LULU analysis on asvs genrated by dada2

autodada2_v1.1.sh -labg 2>&1 | tee outputfile.out  -> This will remove adapters, remove primers with bbduk, run LULU AND run phyloseq to generate relative abundance plots


"

}

while getopts "hrbp:laguf" OPTION; do
     case $OPTION in
         h) usage; exit;;
         r) RENAME="yes";;
         b) BBDUKPRIMERTRIM="yes";;
         p) MANPRIMERTRIM=${OPTARG:-"20,22"};;  ##still need to add bbduk if/then 1_28
	 l) LULU="yes";;
	 a) ADAPTERREM="yes";;
	 g) PHYLO="yes";;
	 u) LULUONLY="yes";;
	 f) PHYONLY="yes";;
     esac
done
shift $((OPTIND-1))      # required, to "eat" the options that have been processed

echo "Working directory = $workdir"
echo "Supplemental scripts held in = $SCRIPTDIR"
echo "Defined sample set ID = $sampleset"
echo "Database directory = $db"
echo "Database to be used for blastn $workdir"
echo "AdapterRemoval script held in = $AdapterRemoval"
echo "Reformat script held in = $reformat"
echo "Location of blastn script $blastn"
echo "Location of makeblastdb script $makeblastdb"

##set adapter and primer sequences below
echo "setting primer and adapter sequences"
adapt1="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
adapt2="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
fwd="GAATTGCAGAACTCCGTGAACC"
rev="CGGGTTCWCTTGTYTGACTTCATGC"
echo "Forward adapter = $adapt1"
echo "Reverse adapter = $adapt2"
echo "Forward primer = $fwd"
echo "Reverse primer = $rev"

echo "Options set currently (blank = not set):"

echo "Renaming: $RENAME"
echo "Running AdapterRemoval: $ADAPTERREM"
echo "Primer trim w/ specified sequences and bbduk: $BBDUKPRIMERTRIM"
echo "Manual primer trim with bbduk: $MANPRIMERTRIM "
echo "Running LULU: $LULU"
echo "Running strictly LULU: $LULUONLY"
echo "Running Phyloseq script: $PHYLO"
echo "Running strictly LULU: $PHYONLY"
echo -n "Review information above, if good, type y and hit [Enter] "
read res
if [ "$res" != "y" ]; then exit 1; fi

if [[ "$PHYONLY" == "yes" ]]
then	echo "Only Phylo set"
	echo "There should be a /phylo directory with necessary files"
	echo "Moving there now"
	cd ./phylo/
        cp $meta ./variable_phyloseq.csvs
	echo "Phyloseq .txt files should be present and meta file too"
        echo "Files currently in phylo:"
        ls
        echo "OK, now I'll check for correct input for LULU..."
        if [[ $(ls | grep -c "variable_phyloseq.csvs") -ge 1 ]]
        then
                echo "variable_phyloseq.csvs present..."
        else
                echo "No variable_phyloseq.csvs found, program stopped."
                exit 1
        fi
	if [[ $(ls | grep -c "phyloseqfile.txt") -ge 1 ]]
        then
                echo "phyloseqfile.txt present..."
        else
                echo "No phyloseqfile.txt found, program stopped."
                exit 1
        fi
	echo "Everything here, lets ball."
	cp $SCRIPTDIR/autophy.r .
	Rscript autophy.r
	echo "Phyloseq complete, check the phylo dir for plots"
fi


if [[ "$LULUONLY" == "yes" ]]
then	
	echo "Only LULU option set"
	echo "Running LULU with files in LULUdir"
	cd ./LULUdir/
	echo " But first, YOU check if you labeled the fasta file and csv correctly"
	echo "ASV fasta should be named ASV_lulu_labs.fasta"
	echo "csv file should be named seqtab.nochim_forLULU.csv"
	echo "Files currently in LULUdir:"
	ls
	echo -n "Review information above, if good, type y and hit [Enter] "
	read res
	if [ "$res" != "y" ]; then exit 1; fi
	echo "OK, now I'll check for correct input for LULU..."
	if [[ $(ls | grep -c "ASV_lulu_labs.fasta") -ge 1 ]]
	then
		echo "ASV_lulu_labs.fasta present..."
	else 
		echo "No fasta file found, program stopped. Check to see everything ran well above"
		exit 1
	fi
	if [[ $(ls | grep -c "seqtab.nochim_forLULU.csv") -ge 1 ]]
	then
		echo "seqtab.nochim_forLULU.csv present, moving forward..."
	else 
		echo "No csv file, program stopped. Check to see everything ran well above"
		exit 1
	fi
	echo "Ok, now performing ASVxASV blastn"
	$makeblastdb -in ASV_lulu_labs.fasta -parse_seqids -dbtype nucl
	$blastn -db ASV_lulu_labs.fasta -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 90 -perc_identity 84 -query ASV_lulu_labs.fasta
	cp $SCRIPTDIR/autoluluv1.r .
	Rscript autoluluv1.r
	echo "Ok LULU completed, should have created a new fasta with confirmed ASVs"
	echo "Checking now..."
	if [[ $(ls | grep -c "LULUasv.fasta") -ge 1 ]]
	then
		echo "LULUasv.fasta present, moving forward..."
	else 
		echo "No LULUasv.fasta file, program stopped. Check to see everything ran well above"
		exit 1
	fi
	echo "Asigning taxonomy to asvs now"
	mkdir ../blast
	mv LULUasv.fasta ../blast
	cd ./blast
	cp $db/"$dbname"* .
	echo "Assesing state of blastdb"
	if [ $(ls | grep -c $dbname) -ge 1 ]
	then 
		echo "Blastdb present. Now checking if built..."
		if [ $(ls | grep -c $dbname) -ge 4 ]
		then
			echo "Uniprot Built moving to run blastn"
		else 
			$makeblastdb -in $dbname -parse_seqids -dbtype nucl &&
			echo "Blast db built." ||
			echo "error in building uniprot db"
		fi
	else
		echo "Error: Blastn db not found in db directly"
		exit 1
	fi
	echo "Running blastn"
	$blastn -db $dbname -query LULUasv.fasta -evalue 1E-4 -out asvlulublast.out -outfmt 6 -num_alignments 1 
	num_asvs=$( grep -c ">" LULUasv.fasta)
        num_hits=$( wc -l asvlulublast.out | awk '{print $1}' )
        if [[ "$num_hits" == "$num_asvs" ]]
	then	
		echo "we good, but you should still check the blast out real quick ;)"
                echo "Number of ASVs = "$num_asvs""
                echo "Number of blast hits = "$num_hits""
                echo -n "check names of libraries and make sure R scripts are available. If all good above, type y and hit [Enter] "
                read res
                if [ "$res" != "y" ]; then exit 1; fi

	else 
		echo "Something wrong with blast output, edit manually and rerun blast script in suppscripts"
		echo "Number of ASVs = "$num_asvs""
		echo "Number of blast hits = "$num_hits""
		echo "check or edit blast output and rerun command with -blast option"
		exit 1
	fi
	awk '{print $2}' asvlulublast.out | sed 's/;/_/g' > newasv_names.txt
	j=1
	for i in $(cat newasv_names.txt)
	do	echo ">seq${j}_$i" >> newnew_asvnames.txt
		j=$(($x+1))
	done
	rm newasv_names.txt
	$python $renamerr LULUasv.fasta newnew_asvnames.txt The_labeled_asvs.fasta
	$reformat in=The_labeled_asvs.fasta out=tmpssasv.fasta fastawrap=500
	mv The_labeled_asvs.fasta ../Rdoe
	touch sequencelist.txt
	echo "      " >> sequencelist.txt
	grep -E '^[ACGT]+$' tmpssasv.fasta >> sequencelist.txt
	touch phylumlist.txt
	echo "Phylum" >> phylumlist.txt
	touch orderlist.txt
	echo "Order" >> orderlist.txt
	touch classlist.txt
	echo "Class" >> classlist.txt
	grep ">" tmpssasv.fasta | awk -F "_" '{print $1}' | awk -F ">" '{print $2}' >> phylumlist.txt
	grep ">" tmpssasv.fasta | awk -F "_" '{print $2}' >> orderlist.txt
	grep ">" tmpssasv.fasta | awk -F "_" '{print $3}' | awk -F "." '{print $1}' >> classlist.txt	
	paste -d' ' sequencelist.txt phylumlist.txt orderlist.txt classlist.txt >> phyloseqfile.txt
	rm *list* 
	rm tmpssasv.fasta
	echo "ASVs renamed and tentative phyloseq input file produced"
	cd $workdir
	if [[ "PHYLO" == "yes" ]]
        then
                mkdir phylo
                mv phyloseqfile.txt ./phylo
		cd ./phylo/
                cp $meta ./variables_phyloseq.csv
                mv ../LULUdir/lulu.RData .
                mv $SCRIPTDIR/autophy_lulu.v1.r .
                Rscript autophy_lulu.v1.r
                echo "done with phylo"
                exit 1
        fi
	echo "Statement ran"
	exit 1
fi


if [[ "$RENAME" == "yes" ]]; then
	echo "renaming files"
	for x in *fastq.gz;
	do	newname=$(ls $x | awk -F "_S" '{print $1}')
        red=$(ls $x | awk -F "_L001_" '{print $2}' | awk -F "_001" '{print $1}')
		mv $x ./"$newname"_"$red".fastq.gz
	done
	echo "files renamed"
else
	echo "Not renaming, moving to FastQC analysis"
fi

if [ $(ls | grep -c "fastqc1_out") == 1 ]
then 	
	echo "FastQC ran"
else
	mkdir fastqc1_out
	$fastqc *fastq.gz -t 15 -o ./fastqc1_out
fi

ls 

echo -n "Check names of libraries for _R*.fastq.gz and all if good above, type y and hit [Enter] "
read res
if [ "$res" != "y" ]; then exit 1; fi


echo "Starting adapter removal with AdapterRemoval2"
echo "Beginning loop"
echo "Removing adapters:"
if [[ "$ADAPTERREM" == "yes" ]]; then
	mkdir originals
	for R1 in  *R1.fastq.gz ;
	do
		echo "$R1"
		directoryname=$(ls $R1 | awk -F "_R" '{print $1}')
		mkdir $directoryname
		filename="${R1%%_R1.fastq.gz}"
		mv $R1 ./$directoryname
		mv ${filename}_R2.fastq.gz ./$directoryname
		cd ./$directoryname
		$AdapterRemoval --file1 $R1 --file2 ${filename}_R2.fastq.gz --threads 15 --adapter1 $adapt1 --adapter2 $adapt2 --minlength 35 --gzip --basename $directoryname
		mv *pair*.truncated.gz ..
		mv *fastq.gz ../originals
		cd ..
		rm -r $directoryname/
		echo "$directoryname is done"
	done
else
	echo "Adapters not being removed"
fi

if [[ "$BBDUKPRIMERTRIM" == "yes" ]]
then
	echo "Primer removal using primer sequences indicated:"
	echo "Forward = $fwd"
	echo "Reverse = $rev"
	echo "making new directory"
	mkdir prime_rem
	mv *truncated.gz ./prime_rem
	cd ./prime_rem
	echo "beginning loop doe"
	for R1 in  *.pair1.truncated.gz ;
	do
			directoryname=$(ls $R1 | awk -F "." '{print $1}')
			mkdir $directoryname
			filename="${R1%%.pair1.truncated.gz}"
			mv $R1 ./$directoryname
			mv ${filename}.pair2.truncated.gz ./$directoryname
			cd ./$directoryname
		    $bbduk in=$R1 in2=${filename}.pair2.truncated.gz out="$directoryname"_R1.fastq.gz out2="$directoryname"_R2.fastq.gz literal=$fwd,$rev copyundefined=t restrictleft=25 k=10 ordered=t mink=2 ktrim=l ecco=f rcomp=t minlength=70 tbo tpe
			mv *.fastq.gz ../../
			echo "$directoryname is done"
			cd ..
	done
	cd ..
	rm -r ./prime_rem
	echo "primer removal completed"
fi

if [ "$MANPRIMERTRIM" != "" ]
then
	echo "Manual trimming by cutting sepcified number of bases from reads, parameters: $MANPRIMERTRIM"
	mkdir prime_rem
	mv *truncated.gz ./prime_rem
	cd ./prime_rem
	FTRIM=`echo "$MANPRIMERTRIM" | cut -d ',' -f 1`
	RTRIM=`echo "$MANPRIMERTRIM" | cut -d ',' -f 2`
	echo "Forward primer length = $FTRIM ; Reverse primer length = $RTRIM"
	for R1 in  *.pair1.truncated.gz ;
	do
			directoryname=$(ls $R1 | awk -F "." '{print $1}')
			mkdir $directoryname
			filename="${R1%%.pair1.truncated.gz}"
			mv $R1 ./$directoryname
			mv ${filename}.pair2.truncated.gz ./$directoryname
			cd ./$directoryname
	    		$bbduk in=$R1 in2=${filename}.pair2.truncated.gz out="$directoryname"_R1.fastq.gz ftl="$FTRIM"
			$bbduk in=${filename}.pair2.truncated.gz out2="$directoryname"_R2.fastq.gz ftl="$RTRIM"
			echo "$directoryname is done"
			mv *fastq.gz ../../
			cd ..
	done
	cd ..
	rm -r ./prime_rem
	echo "primer removal completed"
fi

if [ $(ls | grep -c "truncated.gz") -gt 0 ]
then
	echo "Files still in adapter removal output form, converting now"
	for R1 in  *.pair1.truncated.gz ;
	do		basename=$(ls $R1 | awk -F "." '{print $1}')
			filename="${R1%%.pair1.truncated.gz}"
			mv $R1 ./"$basename"_R1.fastq.gz
			mv ${filename}.pair2.truncated.gz ./"$basename"_R2.fastq.gz
			echo "$basename renamed"
	done
fi
	
if [ $(ls | grep -c "fastq.gz") -gt 0 ] 
then
	echo "Now gunzipping doe"
	gunzip *.fastq.gz
	echo "Gunzipped; moving to R"
else
	echo "Files decompressed already, moving to R"
fi

mkdir Rdir
cp $SCRIPTDIR/autoDADA2Rscript1* ./Rdir
mv *fastq ./Rdir
cd ./Rdir
echo "Before moving to R check location and file names, should be in Rdir SAMPLENAME_RX.fastq"
echo "..."
pwd
ls
echo -n "check names of libraries and make sure R scripts are available. If all good above, type y and hit [Enter] "
read res
if [ "$res" != "y" ]; then exit 1; fi

if [ $(ls | grep -c "autoDADA2Rscript1v1.r" ) -ge 1 ]
        then
            echo "Rscript1 available, will run now"
			Rscript autoDADA2Rscript1v1.r
			echo "Rscript1 ran"
        else
			echo "ERROR: Rscript1 not present in current working directory, please make sure Rscripts present and run again"
			exit 1
fi

echo "Moving outputs around"
if [ $(ls | grep -c "ASV_og_seqs.fasta" ) != 1 ]
then 	echo "Fasta with asvs not present"
		exit 1 
fi

mkdir ../LULUdir
mkdir ../plots
mkdir ../csvs_r

mv *.pdf ../plots
mv *.plot ../plots

mv seqtab.nochim_forLULU.csv ../LULUdir/ 
mv *.csv ../csvs_r/
echo "Editing fasta headers for LULU"
cut -d ';' -f 1 ASV_og_seqs.fasta > ASV_lulu_labs.fasta
mv ASV_lulu_labs.fasta ../LULUdir/


if [[ "$LULU" == "yes" ]]; then
	echo "Running LULU with files in LULUdir"
	cd ../LULUdir/
	echo "Checking for correct input for LULU"
	if [[ $(ls | grep -c "ASV_lulu_labs.fasta") -ge 1 ]]
	then
		echo "ASV_lulu_labs.fasta present..."
	else 
		echo "No fasta file found, program stopped. Check to see everything ran well above"
		exit 1
	fi
	if [[ $(ls | grep -c "seqtab.nochim_forLULU.csv") -ge 1 ]]
	then
		echo "seqtab.nochim_forLULU.csv present, moving forward..."
	else 
		echo "No csv file, program stopped. Check to see everything ran well above"
		exit 1
	fi
	echo "Ok, now performing ASVxASV blastn"
	$makeblastdb -in ASV_lulu_labs.fasta -parse_seqids -dbtype nucl
	$blastn -db ASV_lulu_labs.fasta -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 90 -perc_identity 84 -query ASV_lulu_labs.fasta
	cp $SCRIPTDIR/autoluluv1.r .
	if [ $(ls | grep -c autoluluv1.r) == 1 ]
	then
		Rscript autoluluv1.r
	else
		echo "Auto LULU script not present, stopping"
		exit 1
	fi
	echo "Ok LULU completed, should have created a new fasta with confirmed ASVs"
	echo "Checking now..."
	if [[ $(ls | grep -c "LULUasv.fasta") -ge 1 ]]
	then
		echo "LULUasv.fasta present, moving forward..."
	else 
		echo "No LULUasv.fasta file, program stopped. Check to see everything ran well above"
		exit 1
	fi
	echo "Asigning taxonomy to asvs now"
	mkdir ./blast
	mv LULUasv.fasta ./blast
	cd ./blast
	cp $db/"$dbname"* .
	echo "Assesing state of blastdb"
	if [ $(ls | grep -c $dbname) -ge 1 ]
		then 
			echo "Blastdb present. Now checking if built..."
			if [ $(ls | grep -c $dbname) -ge 4 ]
				then
					echo "Uniprot Built moving to run blastn"
				else 
					$makeblastdb -in $dbname -parse_seqids -dbtype nucl &&
					echo "Blast db built." ||
					echo "error in building uniprot db"
			fi
		else
			echo "Error: Blastn db not found in db directly"
			exit 1
	fi
	echo "Running blastn"
	$blastn -db $dbname -query LULUasv.fasta -evalue 1E-4 -out asvlulublast.out -outfmt 6 -num_alignments 1 
	num_asvs="$( grep -c ">" LULUasv.fasta)"
	num_hits="$( wc -l asvlulublast.out | awk '{print $1}')" 
	if [[ "$num_asv" == "$num_hits" ]]    
	then    
		echo "we good, but you should still check the blast out real quick ;)"
                echo "Number of ASVs = "$num_asvs""
                echo "Number of blast hits = "$num_hits""
                echo -n "check names of libraries and make sure R scripts are available. If all good above, type y and hit [Enter] "
                read res
                if [ "$res" != "y" ]; then exit 1; fi
	else
		echo "Something wrong with blast output, edit manually and rerun blast script in suppscripts"
		echo "Number of ASVs = "$num_asvs""
                echo "Number of blast hits = "$num_hits""
		echo "Unfotunealtely, you must edit the bast output manyallu or use a different db and try again"		exit 1
	fi
	awk '{print $2}' asvlulublast.out | sed 's/;/_/g' > newasv_names.txt
	j=1
	for i in $(cat newasv_names.txt)
	do	echo ">seq${j}_$i" >> newnew_asvnames.txt
        	j=$(($x+1))
	done
	rm newasv_names.txt
	$python $renamerr LULUasv.fasta newnew_asvnames.txt The_labeled_asvs.fasta
	$reformat in=The_labeled_asvs.fasta out=tmpssasv.fasta fastawrap=500
	mv The_labeled_asvs.fasta $workdir/"$sampleset"_luasvs_tax.fasta
	mv LULUasv.fasta ..
	touch sequencelist.txt
	echo "      " >> sequencelist.txt
	grep -E '^[ACGT]+$' tmpssasv.fasta >> sequencelist.txt
	touch phylumlist.txt
	echo "Phylum" >> phylumlist.txt
	touch orderlist.txt
	echo "Order" >> orderlist.txt
	touch classlist.txt
	echo "Class" >> classlist.txt
	grep ">" tmpssasv.fasta | awk -F "_" '{print $1}' | awk -F ">" '{print $2}' >> phylumlist.txt
	grep ">" tmpssasv.fasta | awk -F "_" '{print $2}' >> orderlist.txt
    	grep ">" tmpssasv.fasta | awk -F "_" '{print $3}' | awk -F "." '{print $1}' >> classlist.txt	
	paste -d' ' sequencelist.txt phylumlist.txt orderlist.txt classlist.txt >> phyloseqfile.txt
	rm *list* 
	rm tmpssasv.fasta
	mkdir $workdir/phylo
	mv phyloseqfile.txt $workdir/phylo/
	cd ..
	echo "ASVs renamed and tentative phyloseq input file produced"
else
	echo "LULU not selected to run, files still sitting in the LULUdir though"
	echo "We in bash now for taxonomy of asv fasta file"
	cd ..
	mkdir blast
	mv ./Rdoe/ASV_og_its2.fasta ./blast
	cd ./blast
	cp $db/"$dbname"* .
	echo "Assesing state of blastdb"
	if [ $(ls | grep -c $dbname) -ge 1 ]
		then 
			echo "Blastdb present. Now checking if built..."
			if [ $(ls | grep -c $dbname) -ge 4 ]
			then
				echo "Blastdb Built moving to run blastn"
			else 
				$makeblastdb -in $dbname -dbtype nucl &&
				echo "Blast db built." ||
				echo "error in building uniprot db"
			fi
		else
			echo "Error: Blastn db not found in db directly"
			exit 1
	fi
	echo "Running blastn now"
	$blastn -db $dbname -query ASV_og_seqs.fasta -evalue 1E-4 -out asvblast.out -outfmt 6 -num_alignments 1 
        echo "Number of ASVs = "$num_asvs""
        echo "Number of blast hits = "$num_hits""
  	if [[ "$num_asv" == "$num_hits" ]]
        then
                echo "we good, but you should still check the blast out real quick ;)"
                echo "Number of ASVs = "$num_asvs""
                echo "Number of blast hits = "$num_hits""
                echo -n "check names of libraries and make sure R scripts are available. If all good above, type y and hit [Enter] "
                read res
                if [ "$res" != "y" ]; then exit 1; fi
        else
                echo "Something wrong with blast output, edit manually and rerun blast script in suppscripts"
                echo "Number of ASVs = "$num_asvs""
                echo "Number of blast hits = "$num_hits""
                echo "Unfotunealtely, you must edit the bast output manyallu or use a different db and try again"               exit 1
        fi
	awk '{print $2}' asvblast.out | sed 's/;/_/g' > newasv_names.txt
	j=1
	for i in $(cat newasv_names.txt)
	do	echo ">seq${j}_$i" >> newnew_asvnames.txt
        j=$(($x+1))
	done
	rm newasv_names.txt
	$python $renamerr ASV_og_seqs.fasta newnew_asvnames.txt The_labeled_asvs.fasta
	$reformat in=The_labeled_asvs.fasta out=tmpssasv.fasta fastawrap=500
	mv The_labeled_asvs.fasta $workdir/"$sampleset"_dadaasvs_tax.fasta
	touch sequencelist.txt
	echo "      " >> sequencelist.txt
	grep -E '^[ACGT]+$' tmpssasv.fasta >> sequencelist.txt
	touch phylumlist.txt
	echo "Phylum" >> phylumlist.txt
	touch orderlist.txt
	echo "Order" >> orderlist.txt
	touch classlist.txt
	echo "Class" >> classlist.txt
	grep ">" tmpssasv.fasta | awk -F "_" '{print $1}' | awk -F ">" '{print $2}' >> phylumlist.txt
	grep ">" tmpssasv.fasta | awk -F "_" '{print $2}' >> orderlist.txt
    	grep ">" tmpssasv.fasta | awk -F "_" '{print $3}' | awk -F "." '{print $1}' >> classlist.txt	
	paste -d' ' sequencelist.txt phylumlist.txt orderlist.txt classlist.txt >> phyloseqfile.txt
	rm *list* 
	rm tmpssasv.fasta
	mkdir $workdir/phylo
	mv phyloseqfile.txt $workdir/phylo/
	echo "ASVs renamed and tentative phyloseq input file produced"
fi

pwd
cd $workdir
pwd
echo "Everything done."

if [[ "$PHYLO" == "yes" ]]
then 
	cd ./phylo/
	cp $meta ./variables_phyloseq.csv
	if [[  "$LULU" == "yes" ]]
        then
                mv ../LULUdir/lulu.RData .
                cp $SCRIPTDIR/autophy_lulu.v1.r .
                Rscript autophy_lulu.v1.r
                echo "done with phylo"
        else
                mv ../Rdir/amptmp.RData .
                cp $SCRIPTDIR/autophy_nolulu.v1.r .
                Rscript autophy_nolulu.v1.r
                echo "done with phylo"
        fi
else
	echo "Phyloseq option not set"

fi




echo "Completed analysis"
