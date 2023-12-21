# Drug resistance research

#update to process December 13.
```sh
db=/public/home/2022122/xugang/project/maize_genome/maize_db2
db=/public/home/2022122/xugang/project/alfalfa_genome/alfalfa.zm1
db=/public/home/2022122/xugang/project/antidrug/reference/Ecoli.v1
db2=/public/home/2022122/xugang/project/antidrug/reference/abm18396.puc57
ref=/public/home/2022122/xugang/project/antidrug/reference/Ecoli.v1.fa
ref2=/public/home/2022122/xugang/project/antidrug/reference/abm18396.puc57.fa
datapath2=`pwd`/rawdata/231219/
output=`pwd`/output2
adapt1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
adapt2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
thread=56
#partnum=28
node=Fnode2
#node=Fnode1
#node=Cnode
#node=Gnode
node=Cnode2
#################


readapter(){
name=$1

file='a0-clean'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log



echo -e "#!/bin/bash
#SBATCH -o ${output}/${file}/log/$name.cutadaptor.%j.out
#SBATCH -e ${output}/${file}/log/$name.cutadaptor.%j.error
#SBATCH --partition=${node}
#SBATCH -J 0${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc
echo 'start remove adapter'
#there are some errors.
cutadapt -a ${adapt1} -A ${adapt2} -j ${thread} -m 70 -o ${output}/$file/${name}.1.fq -p ${output}/$file/${name}.2.fq rawdata/${name}_R1.fq.gz rawdata/${name}_R2.fq.gz
echo 'start compressed file'
gzip --force ${output}/$file/$name.1.fq
gzip --force ${output}/${file}/$name.2.fq
">a0.remove.adap.$name.sh

}

readapter2(){
name=$1
name2=$2
name3=$3
file='a0-clean'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log



echo -e "#!/bin/bash
#SBATCH -o ${output}/${file}/log/$name.cutadaptor.%j.out
#SBATCH -e ${output}/${file}/log/$name.cutadaptor.%j.error
#SBATCH --partition=${node}
#SBATCH -J 0${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc
echo 'start remove adapter'
cutadapt -a ${adapt1} -A ${adapt2} -j ${thread} -m 70 -o ${output}/$file/${name}.1.fq -p ${output}/$file/${name}.2.fq ${datapath2}/${name2} ${datapath2}/${name3}
echo 'start compressed file'
gzip --force ${output}/$file/$name.1.fq
gzip --force ${output}/${file}/$name.2.fq
">a0.remove.adap.$name.sh

}
readapter2 EC_15 EC_15.cleaned_1.fastq.gz EC_15.cleaned_2.fastq.gz
readapter2 EC_20 EC_20.cleaned_1.fastq.gz EC_20.cleaned_2.fastq.gz
readapter2 EC_5 EC_5.cleaned_1.fastq.gz EC_5.cleaned_2.fastq.gz
mapf(){
name=$1
file='a1-map'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log



echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 1${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc
bowtie2 -p  ${thread} --end-to-end --sensitive-local -x ${db} -1 $output/a0-clean/$name.1.fq.gz -2 $output/a0-clean/$name.2.fq.gz  --un-conc-gz  $output/$file/${name}.decontaminate.fq.gz -S $output/$file/${name}.sam >  $output/$file/${name}.txt 2>  $output/$file/${name}.txt

" > a1.mapf.${name}.sh
}
mapplasmid(){
name=$1
file='a1-map'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log



echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 1${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc
bowtie2 -p  ${thread} --end-to-end --sensitive-local -x ${db2} -1 $output/$file/${name}.decontaminate.fq.1.gz -2  $output/$file/${name}.decontaminate.fq.2.gz  --un-conc-gz  $output/$file/${name}.other.fq.gz -S $output/$file/${name}.plasmid.sam >  $output/$file/${name}.plasmid.txt 2>  $output/$file/${name}.plasmid.txt

" > a1.mapf.plasmid.${name}.sh


}
samtobam(){
name=$1
file='a2-bam'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 2${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

samtools view -F 4 --threads ${thread} -bS $output/a1-map/${name}.sam > $output/$file/${name}.bam
samtools sort --threads ${thread} $output/$file/${name}.bam > $output/$file/${name}.sort.bam 
samtools index $output/$file/${name}.sort.bam
rm  $output/$file/${name}.bam
" > a2.bam.${name}.sh

}

bamtobg(){
name=$1
file='a3-bigwig'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 3${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc


conda run -n deeptool bamCoverage -b ${output}/a2-bam/${name}.sort.bam -of bigwig -o ${output}/${file}/${name}.bw -p $((thread)) --ignoreDuplicates --binSize 1000 --normalizeUsing RPKM
"> a3.bam2gw.${name}.sh
}

bamtofq(){
name=$1
file='a4-fastq'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 4${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

samtools sort -n -o ${output}/${file}/${name}.bam $output/a2-bam/${name}.sort.bam
bedtools bamtofastq -i ${output}/${file}/${name}.bam -fq ${output}/${file}/${name}.1.fq -fq2 ${output}/${file}/${name}.2.fq
"> a4.bam2fq.${name}.sh
}

fqtofa(){
name=$1
file='a5-fasta'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

cat ${output}/a4-fastq/${name}.1.fq | fastq_to_fasta -n | fastx_trimmer -f 1 -l 50 | fastx_collapser >  ${output}/${file}/${name}.1.50.fa
cat ${output}/a4-fastq/${name}.2.fq | fastq_to_fasta -n | fastx_trimmer -f 1 -l 50 | fastx_collapser >  ${output}/${file}/${name}.2.50.fa

cat ${output}/a4-fastq/${name}.1.fq | fastq_to_fasta -n | fastx_trimmer -f 1 -l 100 | fastx_collapser >  ${output}/${file}/${name}.1.100.fa
cat ${output}/a4-fastq/${name}.2.fq | fastq_to_fasta -n | fastx_trimmer -f 1 -l 100 | fastx_collapser >  ${output}/${file}/${name}.2.100.fa

"> a5.fq2fa.${name}.sh

}
findmotif(){
name=$1
len=$2
file='a6-motif'

[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
[[ -d $output/${file}/${name}.1.$len ]] || mkdir -p $output/${file}/${name}.1.$len
[[ -d $output/${file}/${name}.2.$len ]] || mkdir -p $output/${file}/${name}.2.$len
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 6${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

findMotifs.pl ${output}/a5-fasta/${name}.1.${len}.fa fasta $output/${file}/${name}.1.$len -fasta ${ref} -p ${thread}
findMotifs.pl ${output}/a5-fasta/${name}.2.${len}.fa fasta $output/${file}/${name}.2.$len -fasta ${ref} -p ${thread}
"> a6.motif.${name}.${len}.sh
}
findmotif2(){
name=$1
len=$2
file='a6-motif'

[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
[[ -d $output/${file}/${name}.1.$len ]] || mkdir -p $output/${file}/${name}.1.$len
[[ -d $output/${file}/${name}.2.$len ]] || mkdir -p $output/${file}/${name}.2.$len
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 6${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

findMotifs.pl ${output}/a5-fasta/${name}.1.${len}.fa fasta $output/${file}/${name}.1.$len -fasta ${ref2} -p ${thread}
findMotifs.pl ${output}/a5-fasta/${name}.2.${len}.fa fasta $output/${file}/${name}.2.$len -fasta ${ref}2 -p ${thread}
"> a6.motif.${name}.${len}.sh
}
memef(){
name=$1
len=$2
file='a7-meme'

[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
[[ -d $output/${file}/${name}.1.$len ]] || mkdir -p $output/${file}/${name}.1.$len
[[ -d $output/${file}/${name}.2.$len ]] || mkdir -p $output/${file}/${name}.2.$len
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 7${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc


conda run -n meme meme ${output}/a5-fasta/${name}.1.${len}.fa -oc $output/${file}/${name}.1.$len -p 8 -dna -mod zoops -nmotifs 3 -revcomp
conda run -n meme meme ${output}/a5-fasta/${name}.2.${len}.fa -oc $output/${file}/${name}.2.$len -p 8 -dna -mod zoops -nmotifs 3 -revcomp
"> a7.meme.${name}.${len}.sh


}

#readapter EC10_1
#readapter EC10_2

#mapf EC10_1
#mapf EC10_2
#samtobam EC10_1
#samtobam EC10_2
#mapplasmid EC10_1
#mapplasmid EC10_2
#samtobam EC10_1.plasmid
#samtobam EC10_2.plasmid

#bamtobg EC10_1 
#bamtobg EC10_2
#bamtobg EC10_1.plasmid
#bamtobg EC10_2.plasmid

#bamtofq  EC10_1 
#bamtofq  EC10_2
#bamtofq EC10_1.plasmid
#bamtofq EC10_2.plasmid

#fqtofa  EC10_1 
#fqtofa  EC10_2
#fqtofa EC10_1.plasmid
#fqtofa EC10_2.plasmid

#findmotif EC10_1 50
#findmotif EC10_1 100
#findmotif EC10_2 50
#findmotif EC10_2 100
#findmotif2 EC10_1.plasmid 50 
#findmotif2 EC10_1.plasmid 100
#findmotif2 EC10_2.plasmid 50
#findmotif2 EC10_2.plasmid 100


#memef EC10_1 50 
#memef EC10_1 100
#memef EC10_2 50
#memef EC10_2 100
#memef EC10_1.plasmid 50
#memef EC10_1.plasmid 100
#memef EC10_2.plasmid 50
#memef EC10_2.plasmid 100




```
#update to process September.
```sh
db=/public/home/2022122/xugang/project/maize_genome/maize_db2
db=/public/home/2022122/xugang/project/alfalfa_genome/alfalfa.zm1
db=/public/home/2022122/xugang/project/antidrug/reference/Ecoli.v1
db2=/public/home/2022122/xugang/project/antidrug/reference/abm18396.puc57
ref=/public/home/2022122/xugang/project/antidrug/reference/Ecoli.v1.fa
ref2=/public/home/2022122/xugang/project/antidrug/reference/abm18396.puc57.fa
datapath2=`pwd`/rawdata/
output=`pwd`/output
adapt1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
adapt2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
thread=56
#partnum=28
node=Fnode2
#node=Fnode1
#node=Cnode
#node=Gnode
node=Cnode2
#################


readapter(){
name=$1

file='a0-clean'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log



echo -e "#!/bin/bash
#SBATCH -o ${output}/${file}/log/$name.cutadaptor.%j.out
#SBATCH -e ${output}/${file}/log/$name.cutadaptor.%j.error
#SBATCH --partition=${node}
#SBATCH -J 0${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc
echo 'start remove adapter'
#there are some errors.
cutadapt -a ${adapt1} -A ${adapt2} -j ${thread} -m 70 -o ${output}/$file/${name}.1.fq -p ${output}/$file/${name}.2.fq rawdata/${name}_R1.fq.gz rawdata/${name}_R2.fq.gz
echo 'start compressed file'
gzip --force ${output}/$file/$name.1.fq
gzip --force ${output}/${file}/$name.2.fq
">a0.remove.adap.$name.sh

}

readapter2(){
name=$1
name2=$2
name3=$3
file='a0-clean'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log



echo -e "#!/bin/bash
#SBATCH -o ${output}/${file}/log/$name.cutadaptor.%j.out
#SBATCH -e ${output}/${file}/log/$name.cutadaptor.%j.error
#SBATCH --partition=${node}
#SBATCH -J 0${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc
echo 'start remove adapter'
cutadapt -a ${adapt1} -A ${adapt2} -j ${thread} -m 70 -o ${output}/$file/${name}.1.fq -p ${output}/$file/${name}.2.fq ${datapath2}/${name2} ${datapath2}/${name3}
echo 'start compressed file'
gzip --force ${output}/$file/$name.1.fq
gzip --force ${output}/${file}/$name.2.fq
">a0.remove.adap.$name.sh

}



mapf(){
name=$1
file='a1-map'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log



echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 1${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc
bowtie2 -p  ${thread} --end-to-end --sensitive-local -x ${db} -1 $output/a0-clean/$name.1.fq.gz -2 $output/a0-clean/$name.2.fq.gz  --un-conc-gz  $output/$file/${name}.decontaminate.fq.gz -S $output/$file/${name}.sam >  $output/$file/${name}.txt 2>  $output/$file/${name}.txt

" > a1.mapf.${name}.sh
}
mapplasmid(){
name=$1
file='a1-map'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log



echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 1${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc
bowtie2 -p  ${thread} --end-to-end --sensitive-local -x ${db2} -1 $output/$file/${name}.decontaminate.fq.1.gz -2  $output/$file/${name}.decontaminate.fq.2.gz  --un-conc-gz  $output/$file/${name}.other.fq.gz -S $output/$file/${name}.plasmid.sam >  $output/$file/${name}.plasmid.txt 2>  $output/$file/${name}.plasmid.txt

" > a1.mapf.plasmid.${name}.sh


}
samtobam(){
name=$1
file='a2-bam'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 2${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

samtools view -F 4 --threads ${thread} -bS $output/a1-map/${name}.sam > $output/$file/${name}.bam
samtools sort --threads ${thread} $output/$file/${name}.bam > $output/$file/${name}.sort.bam 
samtools index $output/$file/${name}.sort.bam
rm  $output/$file/${name}.bam
" > a2.bam.${name}.sh

}

bamtobg(){
name=$1
file='a3-bigwig'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 3${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc


conda run -n deeptool bamCoverage -b ${output}/a2-bam/${name}.sort.bam -of bigwig -o ${output}/${file}/${name}.bw -p $((thread)) --ignoreDuplicates --binSize 1000 --normalizeUsing RPKM
"> a3.bam2gw.${name}.sh
}

bamtofq(){
name=$1
file='a4-fastq'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 4${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

samtools sort -n -o ${output}/${file}/${name}.bam $output/a2-bam/${name}.sort.bam
bedtools bamtofastq -i ${output}/${file}/${name}.bam -fq ${output}/${file}/${name}.1.fq -fq2 ${output}/${file}/${name}.2.fq
"> a4.bam2fq.${name}.sh
}

fqtofa(){
name=$1
file='a5-fasta'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

cat ${output}/a4-fastq/${name}.1.fq | fastq_to_fasta -n | fastx_trimmer -f 1 -l 50 | fastx_collapser >  ${output}/${file}/${name}.1.50.fa
cat ${output}/a4-fastq/${name}.2.fq | fastq_to_fasta -n | fastx_trimmer -f 1 -l 50 | fastx_collapser >  ${output}/${file}/${name}.2.50.fa

cat ${output}/a4-fastq/${name}.1.fq | fastq_to_fasta -n | fastx_trimmer -f 1 -l 100 | fastx_collapser >  ${output}/${file}/${name}.1.100.fa
cat ${output}/a4-fastq/${name}.2.fq | fastq_to_fasta -n | fastx_trimmer -f 1 -l 100 | fastx_collapser >  ${output}/${file}/${name}.2.100.fa

"> a5.fq2fa.${name}.sh

}
findmotif(){
name=$1
len=$2
file='a6-motif'

[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
[[ -d $output/${file}/${name}.1.$len ]] || mkdir -p $output/${file}/${name}.1.$len
[[ -d $output/${file}/${name}.2.$len ]] || mkdir -p $output/${file}/${name}.2.$len
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 6${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

findMotifs.pl ${output}/a5-fasta/${name}.1.${len}.fa fasta $output/${file}/${name}.1.$len -fasta ${ref} -p ${thread}
findMotifs.pl ${output}/a5-fasta/${name}.2.${len}.fa fasta $output/${file}/${name}.2.$len -fasta ${ref} -p ${thread}
"> a6.motif.${name}.${len}.sh
}
findmotif2(){
name=$1
len=$2
file='a6-motif'

[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
[[ -d $output/${file}/${name}.1.$len ]] || mkdir -p $output/${file}/${name}.1.$len
[[ -d $output/${file}/${name}.2.$len ]] || mkdir -p $output/${file}/${name}.2.$len
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 6${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc

findMotifs.pl ${output}/a5-fasta/${name}.1.${len}.fa fasta $output/${file}/${name}.1.$len -fasta ${ref2} -p ${thread}
findMotifs.pl ${output}/a5-fasta/${name}.2.${len}.fa fasta $output/${file}/${name}.2.$len -fasta ${ref}2 -p ${thread}
"> a6.motif.${name}.${len}.sh
}
memef(){
name=$1
len=$2
file='a7-meme'

[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
[[ -d $output/${file}/${name}.1.$len ]] || mkdir -p $output/${file}/${name}.1.$len
[[ -d $output/${file}/${name}.2.$len ]] || mkdir -p $output/${file}/${name}.2.$len
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 7${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date

source /public/home/2022122/xugang/bashrc


conda run -n meme meme ${output}/a5-fasta/${name}.1.${len}.fa -oc $output/${file}/${name}.1.$len -p 8 -dna -mod zoops -nmotifs 3 -revcomp
conda run -n meme meme ${output}/a5-fasta/${name}.2.${len}.fa -oc $output/${file}/${name}.2.$len -p 8 -dna -mod zoops -nmotifs 3 -revcomp
"> a7.meme.${name}.${len}.sh


}

#readapter EC10_1
#readapter EC10_2

#mapf EC10_1
#mapf EC10_2
#samtobam EC10_1
#samtobam EC10_2
#mapplasmid EC10_1
#mapplasmid EC10_2
#samtobam EC10_1.plasmid
#samtobam EC10_2.plasmid

#bamtobg EC10_1 
#bamtobg EC10_2
#bamtobg EC10_1.plasmid
#bamtobg EC10_2.plasmid

#bamtofq  EC10_1 
#bamtofq  EC10_2
#bamtofq EC10_1.plasmid
#bamtofq EC10_2.plasmid

#fqtofa  EC10_1 
#fqtofa  EC10_2
#fqtofa EC10_1.plasmid
#fqtofa EC10_2.plasmid

#findmotif EC10_1 50
#findmotif EC10_1 100
#findmotif EC10_2 50
#findmotif EC10_2 100
#findmotif2 EC10_1.plasmid 50 
#findmotif2 EC10_1.plasmid 100
#findmotif2 EC10_2.plasmid 50
#findmotif2 EC10_2.plasmid 100


#memef EC10_1 50 
#memef EC10_1 100
#memef EC10_2 50
#memef EC10_2 100
#memef EC10_1.plasmid 50
#memef EC10_1.plasmid 100
#memef EC10_2.plasmid 50
#memef EC10_2.plasmid 100




```