# Scissors

Пайплайн для обработки экзомных данных.

## Протокол

### Конвертация (опционально)

1. *SOLiD* в *Illumina*:

```bash
cutadapt -j $core_number -c --format=sra-fastq --bwa --action=none -o $output_fastq $input_fastq
```

### Оценка качества входных данных 

```bash
fastqc -o $output_dir -t $core_number $input_fastq_s
```
* Размер последовательностей: >=8 букв
* Сверхрепрезентированные последовательности: отсутствуют

### Обрезка адаптеров

```bash
cutadapt -j $core_number -m 8 -a $adapter_R1 -A $adapter_R2 -o $output_fastq_R1 -p $output_fastq_R2 $input_fastq_R1 $input_fastq_R2
```

### Выравнивание, сортировка и слияние

* Убрать неканонические хромосомы из bed-файла
* Присвоить `@RG`

Выравнивание и сортировка:

```bash
bwa mem -R $RGHeader -t $core_number -v 1 $ref $input_fastq_R1 [$input_fastq_R2] | gatk SortSam --VERBOSITY ERROR -SO queryname -I /dev/stdin -O $output_bam
```

Присваивание `@RG` вручную (если используются сторонние *bam*):

```bash
PicardCommandLine AddOrReplaceReadGroups I=$strandless_bam O=$readgroups_bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
```

Слияние:

```bash
gatk MergeSamFiles --USE_THREADING true -SO queryname -I $input_bam_1 -I $input_bam_2 -O $output_bam
```

Анализ:

```bash
samtools flagstat $bam > $stats_txt
```

### Удаление дубликатов

```bash
gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=2" MarkDuplicates --REMOVE_DUPLICATES true --VERBOSITY ERROR --ASSUME_SORT_ORDER queryname -M $metrics_txt -I $input_bam -O /dev/stdout | gatk SortSam --VERBOSITY ERROR -SO coordinate -I /dev/stdin -O $output_bam
```

### Рекалибровка qual'ов

Для обучения модели требуются вариации в VCF формате (для человеческого генома - [dbSNP >132](https://ftp.ncbi.nih.gov/snp/organisms/)).
Нативная база данных с NCBI требует перепарсинг - другие контиги, а также удаление точек в Ref/Alt.
Обжать базу нужно *bgzip*.

Далее выполняется индексирование (и одновременно проверка на пригодность):

```bash
gatk IndexFeatureFile -I $bgzipped_vcf_gz
```

Рекалибровка:

```bash
gatk --java-options "-Xmx3G -XX:+UseParallelGC -XX:ParallelGCThreads=2" BaseRecalibrator --gatk-config-file $GATK_config -I $input_bam --known-sites $dbSNP -O $BQSRTable -R $ref;
gatk --java-options "-Xmx3G" ApplyBQSR --gatk-config-file $GATK_config -RF MappedReadFilter -RF MappingQualityAvailableReadFilter -RF MappingQualityNotZeroReadFilter -RF NotDuplicateReadFilter -RF NotSecondaryAlignmentReadFilter -RF PassesVendorQualityCheckReadFilter -OBI false -bqsr $BQSRTable -I $input_bam -O $output_bam
```

Отчёт (требует R-пакет *gsalib*):

```bash
gatk AnalyzeCovariates -bqsr $BQSRTable -plots $recalibration_plots_pdf
```

### Просмотр bam-файлов глазом

Используется браузер *IGV*.

### Анализ bam-файлов

```bash
bedtools coverage -hist -sorted -g $ref_index -a $capture_bed -b $bam | grep -P "^all.*$" > $coverage_capture;
bedtools coverage -hist -sorted -g $ref_index -a $NOT_capture_bed -b $bam | grep -P "^all.*$" > $coverage_NOT_capture;
```

* Доля выравненных последовательностей: >80%
* Доля дубликатов в экзонах: <30%
* Обогащение (среднее покрытие capture / среднее покрытие вне capture): >10
* Покрытие 95%: >=10
* Среднее покрытие (арифметическое): >= 70

Это сильно завышенные показатели, но к ним надо стремиться.

### Идентификация вариантов

```bash
gatk --java-options "-Xmx3G" HaplotypeCaller --native-pair-hmm-threads 2 -OVI false --dont-use-soft-clipped-bases true -I $input_bam -O $output_vcf -R $ref
```

* Skip-limit (skip regions of extremely high coverage): стартовое значение 200 (уточнить по ходу работы).
* DP для идеальных данных (см. контроль качества выше): >12-15.
* DP при среднем покрытим <70: >8-10.
* DP при анализе конкретного гена: не учитывать.

### Аннотация вариантов (*Annovar*, *Ensembl VEP*)

Базы данных *ANNOVAR*:

| Protocol | Operation | Comment | Database ID |
|:-----|:---:|:-----|:-----|
| refGene | g | FASTA sequences for all annotated transcripts in RefSeq Gene | refGene |
| knownGene | g | FASTA sequences for all annotated transcripts in UCSC Known Gene | knownGene |
| ensGene | g | FASTA sequences for all annotated transcripts in Gencode v31 Basic collection | ensGene |
| dbnsfp35c | f | whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, PROVEAN, MetaSVM, MetaLR, VEST, M-CAP, CADD, GERP++, DANN, fathmm-MKL, Eigen, GenoCanyon, fitCons, PhyloP and SiPhy scores from dbNSFP version 3.5c, suitable for commercial use | dbnsfp35c |
| dbscsnv11 | f | dbscSNV version 1.1 for splice site prediction by AdaBoost and Random Forest | dbscsnv11 |
| intervar_20180118 | f | InterVar: clinical interpretation of missense variants (indels not supported) | intervar_20180118 |
| gnomad211_genome | f | gnomAD exome collection (v2.1.1), with \"AF AF_popmax AF_male AF_female AF_raw AF_afr AF_sas AF_amr AF_eas AF_nfe AF_fin AF_asj AF_oth non_topmed_AF_popmax non_neuro_AF_popmax non_cancer_AF_popmax controls_AF_popmax\" header | gnomad211_genome |
| gene4denovo201907 | f | gene4denovo database | gene4denovo201907 |
| kaviar_20150923 | f | 170 million Known VARiants from 13K genomes and 64K exomes in 34 projects | kaviar_20150923 |
| hrcr1 | f | 40 million variants from 32K samples in haplotype reference consortium | hrcr1 |
| abraom | f | 2.3 million Brazilian genomic variants abraom |
| ALL.sites.2015_08 | f | The 1000G team fixed a bug in chrX frequency calculation. Based on 201508 collection v5b (based on 201305 alignment) | 1000g2015aug |
| AFR.sites.2015_08 | f | The 1000G team fixed a bug in chrX frequency calculation. Based on 201508 collection v5b (based on 201305 alignment) | 1000g2015aug |
| AMR.sites.2015_08 | f | The 1000G team fixed a bug in chrX frequency calculation. Based on 201508 collection v5b (based on 201305 alignment) | 1000g2015aug |
| EUR.sites.2015_08 | f | The 1000G team fixed a bug in chrX frequency calculation. Based on 201508 collection v5b (based on 201305 alignment) | 1000g2015aug |
| EAS.sites.2015_08 | f | The 1000G team fixed a bug in chrX frequency calculation. Based on 201508 collection v5b (based on 201305 alignment) | 1000g2015aug |
| SAS.sites.2015_08 | f | The 1000G team fixed a bug in chrX frequency calculation. Based on 201508 collection v5b (based on 201305 alignment) | 1000g2015aug |
| gme | f | Great Middle East allele frequency including NWA (northwest Africa), NEA (northeast Africa), AP (Arabian peninsula), Israel, SD (Syrian desert), TP (Turkish peninsula) and CA (Central Asia) | gme |
| esp6500siv2_all | f | alternative allele frequency in All subjects in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls. This is lifted over from hg19 by myself. | esp6500siv2_all |
| avsnp150 | f | dbSNP150 with allelic splitting and left-normalization | avsnp150 |
| clinvar_20200316 | f | Clinvar version 20200316 with separate columns (CLNALLELEID CLNDN CLNDISDB CLNREVSTAT CLNSIG) | clinvar_20200316 |
| regsnpintron | f | prioritize the disease-causing probability of intronic SNVs | regsnpintron |
| revel | f | REVEL scores for non-synonymous variants | revel |
| gwava | f | whole genome GWAVA_region_score and GWAVA_tss_score (GWAVA_unmatched_score has bug in file) | gwava |

```bash
perl table_annovar.pl $output_vcf ./humandb -buildver $genome_assembly -protocol $comma_separated_protocols --operation $comma_separated_operations --remove --vcfinput --thread $threads;
```
Все описанные БД надо доставлять отдельно:

```bash
perl annotate_variation.pl -downdb -buildver $genome_assembly -webfrom annovar $database_ID humandb/
```

### Поиск по таблице

* Отбор по частоте в популяционных данных (берём максимальную из всех БД, но в процессе анализа смотрим данные по всем).
	* Если частое аутосомно-рецессивное заболевание: частота по БД <5%.
	* Если редкое: <0.05-0.5%. 
	* Нижняя граница частоты: 0.01%. При анализе всегда берётся чуть выше.

* Поиск вариантов
	* Поиск loss-of-function variants: frameshift, stop codon, splice-site.
	* Если ничего не нашлось, ищем missence.

### Проверка найденного

* Проверка сырых данных в bam-файле (не последний ли/первый/альтернативный экзон).
* Перепроверка частот в популяции по разным БД.
* GnomAD/ExAC loss-of-function intolerant (pLI) score для гена. Сверять с мышиными фенотипами, литературой и проч.
* HGMD. Описана ли в этой позиции другая замена и какой её эффект.
* ClinVar. Ей доверия меньше, поэтому после остальных.
* OMIM. Сравнение известного фенотипа с фенотипом по БД.
* SIFT, PolyPhen и пр.: предсказание эффекта.

### Прочие особенности анализа

* Если есть подозрение на ограниченный список генов, сначала анализируются варианты в этих генах.
* Если нет: полный анализ.

## Требуемая информация по мутации

* Референсный аллель
* Замена
* Плоидность
* QUAL (*freebayes*)
* Имя гена
* Импакт
* Максимальная частота встречаемости из всех БД, частоты по разным БД

**Примечание.** Отдельным пунктом надо отмечать присутствие мутации в разных БД и разброс по частотам.

* pLI-score
* OMIM
* HGMD
* ClinVar
* Номер экзона, общее число экзонов
* DP: суммарное + для каждого из вариантов, если гетерозигота
* Координата (сортировка)
* Все остальные данные, в т.ч. данные программ, предсказывающих эффекты сплайсинга, замен АК, etc.

## Ссылки

1. [Калькулятор патогенности вариантов нуклеотидной последовательности](http://calc.generesearch.ru/)
