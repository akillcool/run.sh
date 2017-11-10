
# fastq_quality_filter -h
# 	usage: fastq_quality_filter [-h] [-v] [-q N] [-p N] [-z] [-i INFILE] [-o OUTFILE]

# 	version 0.0.6
# 	   [-h]         = This helpful help screen.
# 	   [-q N]       = Minimum quality score to keep.
# 	   [-p N]       = Minimum percent of bases that must have [-q] quality.
# 	   [-z]         = Compress output with GZIP.
# 	   [-i INFILE]  = FASTA/Q input file. default is STDIN.
# 	   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.
# 	   [-v]         = Verbose - report number of sequences.
# 			  If [-o] is specified,  report will be printed to STDOUT.
# 			  If [-o] is not specified (and output goes to STDOUT),
# 			  report will be printed to STDERR.

#   fastq_quality_filter作用是去除测序数据中的低质序列
#   [-q N]       = 最小的需要留下的质量值
#   [-p N]       = 每个reads中最少有百分之多少的碱基需要有-q的质量值
#   [-v]         = 详细-报告序列编号，如果使用了-o则报告会直接在STDOUT，如果没有则输入到STDERR
#   [-i INFILE]  = FASTA/Q输入文件，默认为STDIN。
# 	[-o OUTFILE] = FASTA/Q输出文件，默认为STDOUT。
#   >   将运行日志覆盖写入log文件
fastq_quality_filter -q 20 -p 50 -v -i data/E6.step3.R1.fq -o output/1_filter/E6-1.fq > output/1_filter/E6-1.qc.log

#   BWA 是一种能够将差异度较小的序列比对到一个较大的参考基因组上的软件包。它有三个不同的算法：
#   BWA-backtrack: 是用来比对 Illumina 的序列的，reads 长度最长能到 100bp。
#   BWA-SW: 用于比对 long-read ，支持的长度为 70bp-1Mbp；同时支持剪接性比对。
#   BWA-MEM: 都支持较长的read长度，同时都支持剪接性比对（split alignments)，
#   但是BWA-MEM是更新的算法，也更 快，更准确，且 BWA-MEM 对于 70bp-100bp 的 
#   Illumina 数据来说，效果也更好些。

#   在运用这三种算法之前，需要利用 BWA 的 index 命令，构建出参考基因组的 FM-index，
#   而对与上述的三种不同的算法而言，又有不同的命令：
#   aln/samse/sampe ----> BWA-backtrack   
#   (samse 中的 se 是 single-end 的简写，而 sampe 中的 pe 是 paired-end 的简写）。
#   bwasw ----> BWA-SW
#   mem ----> BWA-MEM

#   index   Usage：bwa index [ –p prefix ] [ –a algoType ] <in.db.fasta>
#           Index database sequence in the FASTA format.
#           OPTIONS:
#           -P STR  输出数据库的前缀；
#               【默认和输入的文件名一致，输出的数据库在其输入文件所在的文件夹，并以该文件名为前缀。】
#           -a [is|bwtsw]   构建index的算法，有两个算法：
#                           is  是默认的算法，虽然相对较快，但是需要较大的内存，
#                               当构建的数据库大于2GB的时候就不能正常工作了。
#                           bwtsw   对于短的参考序列式不工作的，必须要大于等于10MB, 
#                                   但能用于较大的基因组数据，比如人的全基因组。


#   bwa index对fasta文件构建FM-index索引
#   >   将运行日志覆盖写入log文件
bwa index reference/hsa_hpv18.fa > reference/hsa_hpv18.fa.index.log

#   parallel
#       -k  由于多个任务是并行执行的， -k 参数可以让输出的顺序与输入的顺序一致。
#       --pipe  管道前面的文件不是作为参数，而是标准输入传给后面的命令
#       -L N   GNU parallel把输入看作是record，默认一条record只有一行，通过-L N参数可以把record设定为N行
#       参考：https://coyee.com/article/compare/10721-running-bash-commands-in-parallel
#            http://www.gnu.org/software/parallel/parallel_tutorial.html

#   读取上一步去除低质序列后的序列
#   通过python命令对读取出的序列进行并行处理并将序列保存，同时通过 2> 生成错误信息日志文件
cat output/1_filter/E6-1.fq | parallel -k --pipe -L 4 "python3 /Users/agosto/SciWorkflow/bio/ngs/rules/tag/fastx-rm_tag.py -i - -t TTAATTGAGTTGTCATATGTTAATAACGGT -f fastq -a 5,-4,-15,-10 -c 100" > output/1_filter/E6-1.filtered.fq 2> output/1_filter/E6-1.filtered.log

#   mem比对   bwa mem [options] ref.fa reads.fq [mates.fq]
#   -t INT  线程数，默认是1。
#   -M      将 shorter split hits 标记为次优，以兼容 Picard’s markDuplicates 软件。
#   -p      若无此参数：输入文件只有1个，则进行单端比对；若输入文件有2个，则作为paired reads进行比对。
#           若加入此参数：则仅以第1个文件作为输入(输入的文件若有2个，则忽略之)，该文件必须是read1.fq
#           和read2.fa进行reads交叉的数据。
#   -R STR  完整的read group的头部，可以用 '\t' 作为分隔符， 在输出的SAM文件中被解释为制表符TAB. 
#           read group 的ID，会被添加到输出文件的每一个read的头部。
#   -T INT  当比对的分值比 INT 小时，不输出该比对结果，这个参数只影响输出的结果，不影响比对的过程。
#   -a      将所有的比对结果都输出，包括 single-end 和 unpaired paired-end的 reads，
#           但是这些比对的结果会被标记为次优。
#   参考：http://starsyi.github.io/2016/05/24/BWA-%E5%91%BD%E4%BB%A4%E8%AF%A6%E8%A7%A3/#3-2-mem比对

#   进行mem比对，将比对信息存储到sam文件中，并将错误信息输出到日志文件
bwa mem -t 4 -M reference/hsa_hpv18.fa output/1_filter/E6-1.filtered.fq > output/2_align/E6-1.sam.0 2> reference/hsa_hpv18.fa.index.log

#   参考：https://www.ibm.com/developerworks/cn/education/aix/au-gawk/
#   http://starsyi.github.io/2016/05/24/SAM%E6%96%87%E4%BB%B6%E6%A0%BC%E5%BC%8F%E8%AF%B4%E6%98%8E/

#   如果记录E6-1.sam.0第一个字段的不包含@（即sam文件中的非注释行）或者2^8和第二个字段按位与（FLAG位标识）为假，输出记录到E6-1.sam
gawk '{if ($1 ~ /^@/ || !and(256, $2)) print $0}' output/2_align/E6-1.sam.0 > output/2_align/E6-1.sam


#   Usage: samtools view [options] <in.bam>|<in.sam> [region1 [...]]
#   默认情况下不加 region，则是输出所有的 region.
# Options: -b       output BAM
#                   默认下输出是 SAM 格式文件，该参数设置输出 BAM 格式
#          -h       print header for the SAM output
#                   默认下输出的 sam 格式文件不带 header，该参数设定输出sam文件时带 header 信息
#          -H       print header only (no alignments)
#          -S       input is SAM
#                   默认下输入是 BAM 文件，若是输入是 SAM 文件，则最好加该参数，否则有时候会报错。
#          -u       uncompressed BAM output (force -b)
#                   该参数的使用需要有-b参数，能节约时间，但是需要更多磁盘空间。
#          -c       Instead of printing the alignments, only count them and print the 
#                   total number. All filter options, such as ‘-f’, ‘-F’ and ‘-q’ , 
#                   are taken into account.
#          -1       fast compression (force -b)
#          -x       output FLAG in HEX (samtools-C specific)
#          -X       output FLAG in string (samtools-C specific)
#          -c       print only the count of matching records
#          -L FILE  output alignments overlapping the input BED FILE [null]
#          -t FILE  list of reference names and lengths (force -S) [null]
#                   使用一个list文件来作为header的输入
#          -T FILE  reference sequence file (force -S) [null]
#                   使用序列fasta文件作为header的输入
#          -o FILE  output file name [stdout]
#          -R FILE  list of read groups to be outputted [null]
#          -f INT   required flag, 0 for unset [0]
#          -F INT   filtering flag, 0 for unset [0] 
#                   Skip alignments with bits present in INT [0]
#                   数字4代表该序列没有比对到参考序列上
#                   数字8代表该序列的mate序列没有比对到参考序列上
#          -q INT   minimum mapping quality [0]
#          -l STR   only output reads in library STR [null]
#          -r STR   only output reads in read group STR [null]
#          -s FLOAT fraction of templates to subsample; integer part as seed [-1]
#          -?       longer help

#   Usage: samtools sort [-n] [-m <maxMem>] <in.bam> <out.prefix>  
#       -m 参数默认下是 500,000,000 即500M（不支持K，M，G等缩写）。对于处理大数据时，如果内存够用，则设置大点的值，以节约时间。
#       -n 设定排序方式按short reads的ID排序。默认下是按序列在fasta文件中的顺序（即header）和序列从左往右的位点排序。

#   输出BAM格式文件并排序
samtools view -b output/2_align/E6-1.sam | samtools sort > output/2_align/E6-1.bam


#   参考：http://starsyi.github.io/2016/05/25/%E5%8F%98%E5%BC%82%E6%A3%80%E6%B5%8B%EF%BC%88BWA-SAMtools-picard-GATK%EF%BC%89/#7-Duplicates-Marking

#   去除由PCR扩增所形成的duplicates，若失败，程序退出
picard MarkDuplicates REMOVE_DUPLICATES=true I=output/2_align/E6-1.bam O=output/2_align/E6-1.md.bam M=output/2_align/E6-1.md.txt 2> output/2_align/E6-1.index.log || exit 0


#   参考：http://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html

#   将BAM文件转换为BED文件。
bedtools bamtobed -i output/2_align/E6-1.md.bam > output/3_parse/E6-1.bed

#   参考：http://blog.sina.com.cn/s/blog_46d703950101fqi5.html

#   把碱基对数量多于50的记录筛选出来，保存到E6-1.bed.0
awk '$3 - $2 >= 50' output/3_parse/E6-1.bed > output/3_parse/E6-1.bed.0


#   参考：http://bedtools.readthedocs.io/en/latest/content/tools/merge.html

#   合并位于同一个bed文件中的重叠区域
bedtools merge -i output/3_parse/E6-1.bed.0 > output/3_parse/E6-1.bed.1

#   -c	Specify columns from the input file to operate upon (see -o option, below). 
#       Multiple columns can be specified in a comma-delimited list.
# -o	
# Specify the operation that should be applied to -c.
# Valid operations:
# sum, min, max, absmin, absmax,
# mean, median,
# collapse (i.e., print a delimited list (duplicates allowed)),
# distinct (i.e., print a delimited list (NO duplicates allowed)),
# count
# count_distinct (i.e., a count of the unique values in the column),
# Default: sum
# Multiple operations can be specified in a comma-delimited list.
# If there is only column, but multiple operations, all operations will be
# applied on that column. Likewise, if there is only one operation, but
# multiple columns, that operation will be applied to all columns.
# Otherwise, the number of columns must match the the number of operations,
# and will be applied in respective order.

# E.g., -c 5,4,6 -o sum,mean,count will give the sum of column 5,
# the mean of column 4, and the count of column 6.
# The order of output columns will match the ordering given in the command.

#   指定输入文件中的4，6，1列，对4，6列执行分隔打印，对1列（计算重叠区间的个数）计数
bedtools merge -i output/3_parse/E6-1.bed.0 -c 4,6,1 -o collapse,collapse,count > output/3_parse/E6-1.bed.2


#   参考：http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
#   Usage:
#       bedtools coverage [OPTIONS] -a <FILE> \
#                              -b <FILE1, FILE2, ..., FILEN>
#   -d	Report the depth at each position in each A feature. 
#       Positions reported are one based. Each position and depth follow the complete A feature.

#   计算覆盖度，并将结果中的第一列、第五列合并输出到E6-1.stat
bedtools coverage -d -a output/3_parse/E6-1.bed.1 -b output/3_parse/E6-1.bed.0 | bedtools merge -c 1,5 -o count,mean > output/3_parse/E6-1.stat

#   若上一步记录的第5个字段和第4个字段都大于等于1，则输出到E6-1.bed.3
awk '$5 >= 1 && $4 >= 1' output/3_parse/E6-1.stat > output/3_parse/E6-1.bed.3

#   参考http://blog.csdn.net/sunnyyoona/article/details/52981239

#   将文件E6-1.bed.3的记录的第1，2，3个字段以“_”和缩进隔开，并再输出一遍整条记录，并排序
#   对文件E6-1.bed.2的记录进行类似的操作
#   将两个文件按照第一列进行匹配，然后把第一列删除，结果保存到E6-1.stat.filtered
#   注意输入重定向<前面有空格
join -t $'\t' -j 1 <(gawk '{print $1 "_" $2 "_" $3 "\t" $0}' output/3_parse/E6-1.bed.3 | sort) <(gawk '{print $1 "_" $2 "_" $3 "\t" $4 "\t" $5 "\t" $6}' output/3_parse/E6-1.bed.2 | sort) | cut -f 2- > output/3_parse/E6-1.stat.filtered

#   删除E6-1.bed.?文件，‘？’匹配所有字符
rm output/3_parse/E6-1.bed.?

#   去除换行，将各列以制表符分隔，$f字段以“，”分隔=>array1；$g字段也以“，”分隔=>array2
#   将array1[$_]、$array2[$_]、$a、$b、$c、$d、$e以制表符分隔输出到E6-1.stat.filtered
#   排序后的文件为E6-1.lst
perl -ne 'chomp; ($a, $b, $c, $d, $e, $f, $g, $h) = split /\t/; @array1 = split(/,/, $f); @array2 = split(/,/, $g); print "$array1[$_]\t$array2[$_]\t$a\t$b\t$c\t$d\t$e\n" for (0..$#array1);' output/3_parse/E6-1.stat.filtered | sort > output/4_validate/E6-1.lst

#   将E6-1.stat.filtered的第六列的每个记录用“，”分割，将分割后的结果随机排序
#   选择随机后的前十个记录排序后输出到E6-1.site_10.lst
cut -f 6 output/3_parse/E6-1.stat.filtered | perl -MList::Util=shuffle -ne 'chomp; @array = split /,/; @shuffled_indexes = shuffle(0..$#array); @pick_indexes = @shuffled_indexes[0..(10-1)]; @picks = @array[ @pick_indexes ]; print "$_\n" for @picks;' | sort > output/4_validate/E6-1.site_10.lst


#   参考：https://github.com/lh3/seqtk

#   extract subsequences from FASTA/Q
seqtk subseq data/E6.step3.R1.fq output/4_validate/E6-1.site_10.lst > output/4_validate/E6-1.site_10.fq

python3 /Users/agosto/SciWorkflow/bio/ngs/rules/tag/fastx-view_tag.py -i output/4_validate/E6-1.site_10.fq -t TTAATTGAGTTGTCATATGTTAATAACGGT -c 100 > output/4_validate/E6-1.site_10.align.txt


#   参考shell中join的用法，作用就是将两个表格文件整合到一起
#   通过perl正则表达式处理后再排序输出
#   sed命令，参考http://linux.51yip.com/search/sed
#   最后输出为表格
join -t $'\t' -j 1 -o 1.3,1.4,1.5,1.6,1.7,1.2,1.1,2.2 <(join -t $'\t' -j 1 output/4_validate/E6-1.site_10.lst output/4_validate/E6-1.lst | sort) <(seqtk seq -a output/4_validate/E6-1.site_10.fq | perl -pe 's/^>(\S+).*\n/$1\t/' | sort) | sort -k 1,1g -k 2,2n | sed '1s/^/Chrom\tStart\tEnd\tLength\tDepth\tStrand\tReadName\tReadSeq\n/' > output/4_validate/E6-1.site_10.fq.xls

if [ -s output/3_parse/E6-1.stat.filtered ]; then
    cut -f 1-3 output/3_parse/E6-1.stat.filtered | awk '{print $1 "\t" $2-500 "\t" $3+1000}' | perl -pe 's/\t-\d+\t/\t0\t/' > `dirname output/5_visualization/E6-1_offtargets.svg`/E6-1.bed

    seqtk subseq reference/hsa_hpv18.fa `dirname output/5_visualization/E6-1_offtargets.svg`/E6-1.bed > `dirname output/5_visualization/E6-1_offtargets.svg`/E6-1.fa

    python3 /Users/agosto/SciWorkflow/bio/ngs/rules/tag/fastx-find_tag.py -i `dirname output/5_visualization/E6-1_offtargets.svg`/E6-1.fa -t GAAGCTACCTGATCTGTGCANGG -f fasta -a 5,1,-150,-100 -c 52 | sort > `dirname output/5_visualization/E6-1_offtargets.svg`/E6-1.align.txt

    join -t $'\t' -j 1 -o 1.1,1.2,1.3,2.2 `dirname output/5_visualization/E6-1_offtargets.svg`/E6-1.align.txt <(awk '{print $1 ":" $2-500+1 "-" $3+1000 "\t" $8}' output/3_parse/E6-1.stat.filtered | perl -pe 's/:-\d+-/:1-/' | sort) > `dirname output/5_visualization/E6-1_offtargets.svg`/E6-1_identifiedOfftargets.xls

    sed -i '1s/^/#BED Chromosome\tOff-Target Sequence\tTarget Sequence\tbi.sum.mi\n/' `dirname output/5_visualization/E6-1_offtargets.svg`/E6-1_identifiedOfftargets.xls

    python /Users/agosto/SciWorkflow/bio/ngs/rules/visualization/guide-seq.py `dirname output/5_visualization/E6-1_offtargets.svg`/E6-1_identifiedOfftargets.xls output/5_visualization/E6-1_offtargets.svg E6-1
else
    touch output/5_visualization/E6-1_offtargets.svg
fi

