# 软件说明
capTEs测序技术自动化分析软件


# 安装
tar zxvf capTEtools.tar.gz
cd capTEtools
tar xvf citeTools.tar
conda env create -f capTE_env.yml
conda activate capTE


# Examples
python capTEtools.py -id Test -i example/test.fastq -o test -r GRCh38.p13.genome.fa -g gencode.v37.chr_patch_hapl_scaff.annotation.gtf -b database/GRCh38.p13.genome.fa.TE.bed 2> test.log


