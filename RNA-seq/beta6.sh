#!/bin/sh

#SBATCH -J RNAseq
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000

source ~/.bash_profile

usage(){
    echo "R N A - S E Q   W O R K F L O W - @bixBeta"
    echo ""
    echo ""
    echo "Usage: bash" $0 "[-h arg] [-p arg] [-d args] [-t arg] [-g arg] [-r arg] [-s arg] [-c arg]"
    echo
    echo "---------------------------------------------------------------------------------------------------------------------------"
    echo "[-h] --> Display Help"
    echo "[-p] --> Project Identifier Number"
    echo "[-d] --> Comma Spearated Values for Delimiter and Field <delim,field or default> default: _,5 "
    echo "[-t] --> Fastq Trimming <nextSE, nextPE, 4colorSE, miSeqPE, novaPE >"
    echo "[-g] --> Reference Genome <hg38, GRCh38, mm10, GRCm38, etc.>"
    echo "[-r] --> <SE> <SES> <PE> <PES> <SEBS> <PEBS> <UNMS> or <UNMP> "
    echo "[-s] --> Library Strandedness <0, 1, 2> where 1 = first strand, 2 = reverse strand, 0 for unstranded counts "
    echo "[-c] --> GeneBody Coverage <yes, no> "
    echo "---------------------------------------------------------------------------------------------------------------------------"
    echo ""
    echo "******************************************** "
    echo "**********/ EXTENDED GENOME LIST /********** "
    echo "[hg38]=/workdir/genomes/Homo_sapiens/hg38/UCSC/hg38.star "
    echo "[mm10]=/workdir/genomes/Mus_musculus/mm10/UCSC/mm10.star "
    echo "[GRCh38]=/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/GRCh38.star "
    echo "[GRCm38]=/workdir/genomes/Mus_musculus/mm10/ENSEMBL/GRCm38.star "
    echo "[cat]=/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/genomeDir "
    echo "[chicken]=/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/galgal5.star "
    echo "[chicken6]=/workdir/genomes/Gallus_gallus/Galgal6/ENSEMBL/gtf.102/star.index.102 "
    echo "[horse]=/workdir/genomes/Equus_caballus/EquCab3/ENSEMBL/gtf.102/star.index.102 "
    echo "[horse2]=/workdir/genomes/Equus_caballus/EquCab2/ENSEMBL/EquCab2.star.index "
    echo "[rat]=/workdir/genomes/Rattus_norvegicus/rn6/ENSEMBL/rat.star "
    echo "[ercc]=/workdir/genomes/contaminants/ERCC_spikeIns/ercc.star "
    echo "[lonchura]=/workdir/genomes/Lonchura_striata/LonStrDom1/ENSEMBL/lonchura.star "
    echo "[goose]=/workdir/genomes/Anser_brachyrhynchus/ASM259213v1/ENSEMBL/goose.star "
    echo "[ehv8]=/workdir/genomes/FastQ_Screen_Genomes/EHV8/ehv8.star "
    echo "[erdman]=/workdir/genomes/Mycobacterium_tuberculosis/Ensembl_GCA_000668235/GCA_000668235.star "
    echo "[TB]=/workdir/genomes/Mycobacterium_tuberculosis/CDC1551_Ensembl/cdc1551.star "
    echo "[maize]=/workdir/genomes/Zea_mays/B73_RefGen_v4/ENSEMBL/star.maize "
    echo "[finch]=/workdir/genomes/Taeniopygia_guttata/taeGut3.2.4/ENSEMBL/UPDATED.ANNOTS/star.index.updated "
    echo "[finch2]=/workdir/genomes/Geospiza_fortis_ground_finch/GeoFor_1.0/NCBI/star.index "
    echo "[dog]=/workdir/genomes/Canis_familiaris/canFam3/ENSEMBL/CanFam3.1.101/genomeDir "
    echo "[green]=/workdir/genomes/Chlorocebus_sabaeus/CHlSab1/ENSEMBL/genomeDir "
    echo "[BG8]=/workdir/genomes/Methylomicrobium_album_BG8/ASM21427v3/NCBI/custom3-CDS2exon/genomeDir_fixSpaces "
    echo "[faba]=/workdir/genomes/Vicia_faba/VfEP_Reference-Unigene/NCBI/genomeDir-wo-gff"
    echo "[aphid]=/workdir/genomes/Acyrthosiphon_pisum/pea_aphid_22Mar2018_4r6ur/NCBI/genomeDir"
    echo "[cholera]=/workdir/genomes/Vibrio_cholerae/N16961/NCBI/genomeDir"
    echo "[salmonella]=/workdir/genomes/Salmonella_enterica/ASM21085v2/NCBI/genomeDir"
    echo "[aedes]=/workdir/genomes/Aedes_aegypti/AaegL5.0/NCBI/AaegL5.0.star"
    echo "[DC3000]=/workdir/genomes/Pseudomonas_syringae/Tomato_DC3000/NCBI/star.index"
    echo "[Theileria]=/workdir/genomes/Theileria/Tannulata_ASM322v1/NCBI/starIndex"
    echo "[cow]=/workdir/genomes/Bos_taurus/ARS-UCD1.2_GCA_002263795.2/ENSEMBL/star.index"
    echo "[EA273]=/workdir/genomes/Erwinia_amylovora/GCF_000091565.1/ncbi/genomeDir"
    echo "[ddSmed]=/workdir/genomes/Schmidtea_mediterranea/dd_Smed_v6/NCBI/genomeDir"

}



trimSE(){

        mkdir TrimQC_stats fastQC trimmed_fastqs
        for i in fastqs/*.gz
        do
            trim_galore -j 8  --nextseq 20 --gzip --length 50  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
        done
        mv *_trimming_report.txt TrimQC_stats
        mv *trimmed.fq.gz trimmed_fastqs

}

trimHiSeq(){

        mkdir TrimQC_stats fastQC trimmed_fastqs
        for i in fastqs/*.gz
        do
            trim_galore -j 8  --quality 20 --gzip --length 50  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
        done
        mv *_trimming_report.txt TrimQC_stats
        mv *trimmed.fq.gz trimmed_fastqs

}

trimPE(){

                cd fastqs
                ls -1 *_R1.fastq* > .R1
                ls -1 *_R2.fastq* > .R2
                paste -d " " .R1 .R2 > Reads.list

                readarray fastqs < Reads.list
                mkdir fastQC

                for i in "${fastqs[@]}"
                do
                        trim_galore -j 8  --nextseq 20 --gzip --length 50  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
                done

                mkdir TrimQC_stats trimmed_fastqs
                mv *_trimming_report.txt TrimQC_stats
                mv *_val* trimmed_fastqs
                mv TrimQC_stats fastQC trimmed_fastqs ..

                cd ..
}

trimMiSeqPE(){

                cd fastqs
                ls -1 *_R1_* > .R1
                ls -1 *_R2_* > .R2
                paste -d " " .R1 .R2 > Reads.list

                readarray fastqs < Reads.list
                mkdir fastQC

                for i in "${fastqs[@]}"
                do
                        trim_galore -j 8  --quality 20 --gzip --length 50  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
                done

                mkdir TrimQC_stats trimmed_fastqs
                mv *_trimming_report.txt TrimQC_stats
                mv *_val* trimmed_fastqs
                mv TrimQC_stats fastQC trimmed_fastqs ..

                cd ..
}

trimHiSeqPE(){

                cd fastqs
                ls -1 *_1.fq* > .R1
                ls -1 *_2.fq* > .R2
                paste -d " " .R1 .R2 > Reads.list

                readarray fastqs < Reads.list
                mkdir fastQC

                for i in "${fastqs[@]}"
                do
                        trim_galore -j 8  --quality 20 --gzip --length 50  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
                done

                mkdir TrimQC_stats trimmed_fastqs
                mv *_trimming_report.txt TrimQC_stats
                mv *_val* trimmed_fastqs
                mv TrimQC_stats fastQC trimmed_fastqs ..

                cd ..
}

declare -A genomeDir

genomeDir=( ["hg38"]="/workdir/genomes/Homo_sapiens/hg38/UCSC/hg38.star" \
["mm10"]="/workdir/genomes/Mus_musculus/mm10/UCSC/mm10.star" \
["GRCh38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/GRCh38.star" \
["GRCm38"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/GRCm38.star" \
["cat"]="/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/genomeDir" \
["chicken"]="/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/galgal5.star" \
["chicken6"]="/workdir/genomes/Gallus_gallus/Galgal6/ENSEMBL/gtf.102/star.index.102" \
["horse"]="/workdir/genomes/Equus_caballus/EquCab3/ENSEMBL/gtf.102/star.index.102" \
["horse2"]="/workdir/genomes/Equus_caballus/EquCab2/ENSEMBL/EquCab2.star.index" \
["ATCC_13047"]="/workdir/genomes/Enterobacter_cloacae/ATCC_13047/custom/ATCC_13047.GTF" \
["grape"]="/workdir/genomes/Vitis_vinifera/GCA_000003745.2/ENSEMBL/Vitis_vinifera.12X.43.bed12" \
["rat"]="/workdir/genomes/Rattus_norvegicus/rn6/ENSEMBL/rat.star" \
["ercc"]="/workdir/genomes/contaminants/ERCC_spikeIns/ercc.star" \
["lonchura"]="/workdir/genomes/Lonchura_striata/LonStrDom1/ENSEMBL/lonchura.star" \
["goose"]="/workdir/genomes/Anser_brachyrhynchus/ASM259213v1/ENSEMBL/goose.star" \
["ehv8"]="/workdir/genomes/FastQ_Screen_Genomes/EHV8/ehv8.star" \
["erdman"]="/workdir/genomes/Mycobacterium_tuberculosis/Ensembl_GCA_000668235/GCA_000668235.star" \
["TB"]="/workdir/genomes/Mycobacterium_tuberculosis/CDC1551_Ensembl/cdc1551.star" \
["maize"]="/workdir/genomes/Zea_mays/B73_RefGen_v4/ENSEMBL/star.maize" \
["finch"]="/workdir/genomes/Taeniopygia_guttata/taeGut3.2.4/ENSEMBL/UPDATED.ANNOTS/star.index.updated" \
["finch2"]="/workdir/genomes/Geospiza_fortis_ground_finch/GeoFor_1.0/NCBI/star.index" \
["green"]="/workdir/genomes/Chlorocebus_sabaeus/CHlSab1/ENSEMBL/genomeDir" \
["dog"]="/workdir/genomes/Canis_familiaris/canFam3/ENSEMBL/CanFam3.1.101/genomeDir" \
["faba"]="/workdir/genomes/Vicia_faba/VfEP_Reference-Unigene/NCBI/genomeDir-wo-gff" \
["aphid"]="/workdir/genomes/Acyrthosiphon_pisum/pea_aphid_22Mar2018_4r6ur/NCBI/genomeDir" \
["cholera"]="/workdir/genomes/Vibrio_cholerae/N16961/NCBI/genomeDir" \
["BG8"]="/workdir/genomes/Methylomicrobium_album_BG8/ASM21427v3/NCBI/custom3-CDS2exon/genomeDir_fixSpaces" \
["aedes"]="/workdir/genomes/Aedes_aegypti/AaegL5.0/NCBI/AaegL5.0.star" \
["aedesVB"]="/workdir/genomes/Aedes_aegypti/AaegL5.0/Vectorbase/AaegyptiLVP_AGWG_release49/genomeDir.vectorBase" \
["DC3000"]="/workdir/genomes/Pseudomonas_syringae/Tomato_DC3000/NCBI/star.index" \
["Theileria"]="/workdir/genomes/Theileria/Tannulata_ASM322v1/NCBI/starIndex" \
["cow"]="/workdir/genomes/Bos_taurus/ARS-UCD1.2_GCA_002263795.2/ENSEMBL/star.index" \
["salmonella"]="/workdir/genomes/Salmonella_enterica/ASM21085v2/NCBI/genomeDir" \
["EA273"]="/workdir/genomes/Erwinia_amylovora/GCF_000091565.1/ncbi/genomeDir" \
["ddSmed"]="/workdir/genomes/Schmidtea_mediterranea/dd_Smed_v6/NCBI/genomeDir" \
["yeast"]="/workdir/genomes/Saccharomyces_cerevisiae/R64-1-1_GCA_000146045.2/ENSEMBL/star.index" )

declare -A bed12

bed12=(	["hg38"]="/workdir/genomes/Homo_sapiens/hg38/UCSC/genes.bed12" \
["mm10"]="/workdir/genomes/Mus_musculus/mm10/UCSC/BED12/mm10.ucsc.bed12" \
["GRCh38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.bed12" \
["GRCm38"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.bed12" \
["cat"]="/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/Felis_catus.Felis_catus_9.0.95.bed12" \
["chicken"]="/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/Gallus_gallus.Gallus_gallus-5.0.bed12" \
["chicken6"]="/workdir/genomes/Gallus_gallus/Galgal6/ENSEMBL/gtf.102/Gallus_gallus.GRCg6a.102.bed12" \
["horse"]="/workdir/genomes/Equus_caballus/EquCab3/ENSEMBL/gtf.102/Equus_caballus.EquCab3.0.102.bed12" \
["horse2"]="/workdir/genomes/Equus_caballus/EquCab2/ENSEMBL/Equus_caballus.EquCab2.94.BED12" \
["ATCC_13047"]="/workdir/genomes/Enterobacter_cloacae/ATCC_13047/GCF_000025565.1_ASM2556v1_genomic.bed12" \
["grape"]="/workdir/genomes/Vitis_vinifera/GCA_000003745.2/ENSEMBL/Vitis_vinifera.star" \
["rat"]="/workdir/genomes/Rattus_norvegicus/rn6/ENSEMBL/Rattus_norvegicus.Rnor_6.0.bed12" \
["lonchura"]="/workdir/genomes/Lonchura_striata/LonStrDom1/ENSEMBL/Lonchura_striata_domestica.LonStrDom1.bed12" \
["goose"]="/workdir/genomes/Anser_brachyrhynchus/ASM259213v1/ENSEMBL/Anser_brachyrhynchus.ASM259213v1.bed12" \
["maize"]="/workdir/genomes/Zea_mays/B73_RefGen_v4/ENSEMBL/Zea_mays.B73_RefGen_v4.bed12" \
["finch"]="/workdir/genomes/Taeniopygia_guttata/taeGut3.2.4/ENSEMBL/UPDATED.ANNOTS/Taeniopygia_guttata.bTaeGut1_v1.p.bed12" \
["finch2"]="/workdir/genomes/Geospiza_fortis_ground_finch/GeoFor_1.0/NCBI/GCF_000277835.1_GeoFor_1.0_genomic.bed12" \
["dog"]="/workdir/genomes/Canis_familiaris/canFam3/ENSEMBL/CanFam3.1.101/Canis_lupus_familiaris.CanFam3.1.101.bed12" \
["aedes"]="/workdir/genomes/Aedes_aegypti/AaegL5.0/NCBI/GCF_002204515.2_AaegL5.0_genomic.bed12" \
["aedesVB"]="/workdir/genomes/Aedes_aegypti/AaegL5.0/Vectorbase/AaegyptiLVP_AGWG_release49/VectorBase-49_AaegyptiLVP_AGWG.BED12" \
["cow"]="/workdir/genomes/Bos_taurus/ARS-UCD1.2_GCA_002263795.2/ENSEMBL/Bos_taurus.ARS-UCD1.2.BED12" \
["BG8"]="/workdir/genomes/Methylomicrobium_album_BG8/ASM21427v3/NCBI/test/ncbi-genomes-2020-08-27/GCF_000214275.2_ASM21427v3_genomic.bed12" \
["ddSmed"]="/workdir/genomes/Schmidtea_mediterranea/dd_Smed_v6/NCBI/dd_Smed_v6.pcf.bed12" \
["yeast"]="/workdir/genomes/Saccharomyces_cerevisiae/R64-1-1_GCA_000146045.2/ENSEMBL/Saccharomyces_cerevisiae.R64-1-1.bed12" )




align(){

    cd trimmed_fastqs

    for i in *_trimmed.fq.gz

        do

        iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

        STAR \
        --runThreadN 12 \
        --genomeDir ${genomeDir[${DIR}]} \
        --readFilesIn $i \
        --readFilesCommand gunzip -c \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $iSUB. \
        --limitBAMsortRAM 61675612266 \
        --quantMode GeneCounts


        done


        # source activate
        multiqc -f -n ${PIN}.star.multiqc.report .
        mkdir STAR.COUNTS STAR.BAMS STAR.LOGS
        mv *.ReadsPerGene.out.tab STAR.COUNTS
        mv *.bam STAR.BAMS
        mv *.out *.tab *_STARtmp *.list *star.multiqc.report_data STAR.LOGS
        mkdir STAR
        mv STAR.* *.html STAR

        mv STAR ..
    cd ..

}

se_split(){

        cd trimmed_fastqs

        for i in *_trimmed.fq.gz

        do

          iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

          STAR \
          --runThreadN 12 \
          --genomeDir ${genomeDir[${DIR}]} \
          --readFilesIn $i \
          --readFilesCommand gunzip -c \
          --outSAMstrandField intronMotif \
          --outFilterIntronMotifs RemoveNoncanonical \
          --outReadsUnmapped Fastx \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix $iSUB. \
          --limitBAMsortRAM 61675612266 \
          --quantMode GeneCounts

        done

    #			source activate RSC
                multiqc -f -n ${PIN}.starSPLIT.multiqc.report .
                mkdir STAR.SPLIT.COUNTS STAR.SPLIT.BAMS STAR.SPLIT.LOGS STAR.SPLIT.Unmapped
                mv *Unmapped.out.mate* STAR.SPLIT.Unmapped
                cd STAR.SPLIT.Unmapped
                    for i in *mate*
                        do
                            mv $i `echo $i | sed "s/Unmapped/not.$DIR/g"`
                        done
                cd ..
                mv *.ReadsPerGene.out.tab STAR.SPLIT.COUNTS
                mv *.bam STAR.SPLIT.BAMS
                mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SPLIT.LOGS
                mkdir STAR.SPLIT
                mv STAR.* *.html STAR.SPLIT
                mv STAR.SPLIT ..
                cd ..

}



alignPE(){

    cd trimmed_fastqs
    ls -1 *_1.fq.gz > .trR1
    ls -1 *_2.fq.gz > .trR2
    paste -d " " .trR1 .trR2 > Trimmed.list

    readarray trimmedFastqs < Trimmed.list

    for i in "${trimmedFastqs[@]}"

        do

            iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

            STAR \
            --runThreadN 12 \
            --genomeDir ${genomeDir[${DIR}]} \
            --readFilesIn $i \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${iSUB}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

        done

    #    source activate RSC
        multiqc -f -n ${PIN}.star.multiqc.report .
        mkdir STAR.COUNTS STAR.BAMS STAR.LOGS
        mv *.ReadsPerGene.out.tab STAR.COUNTS
        mv *.bam STAR.BAMS
        mv *.out *.tab *_STARtmp *.list *star.multiqc.report_data STAR.LOGS
        mkdir STAR
        mv STAR.* *.html STAR

        mv STAR ..
    cd ..


    # cd STAR/STAR.BAMS
    # 	for i in *.bam
    # 	do
    # 		/programs/bin/samtools/samtools index -b $i
    # 	done
    # cd ..
    # echo
    # echo
    # pwd
    # echo
    # echo
    # source activate RSeQC
    # geneBody_coverage.py -r ${bed12[${DIR}]} -i STAR.BAMS/ -o ${PIN}
    # mkdir geneBodyCov
    # mv *geneBodyCoverage.* log.txt geneBodyCov
    # cd ..


}

pe_split(){

          cd trimmed_fastqs
          ls -1 *_1.fq.gz > .trR1
          ls -1 *_2.fq.gz > .trR2
          paste -d " " .trR1 .trR2 > Trimmed.list

          readarray trimmedFastqs < Trimmed.list

          for i in "${trimmedFastqs[@]}"

          do

            iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

            STAR \
            --runThreadN 12 \
            --genomeDir ${genomeDir[${DIR}]} \
            --readFilesIn $i \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${iSUB}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

          done

                    # source activate RSC
                    multiqc -f -n ${PIN}.starSPLIT.multiqc.report .
                    mkdir STAR.SPLIT.COUNTS STAR.SPLIT.BAMS STAR.SPLIT.LOGS STAR.SPLIT.Unmapped
                    mv *Unmapped.out.mate* STAR.SPLIT.Unmapped
                    cd STAR.SPLIT.Unmapped
                        for i in *mate*
                            do
                                mv $i `echo $i | sed "s/Unmapped/not.$DIR/g"`
                            done
                    cd ..
                    mv *.ReadsPerGene.out.tab STAR.SPLIT.COUNTS
                    mv *.bam STAR.SPLIT.BAMS
                    mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SPLIT.LOGS
                    mkdir STAR.SPLIT
                    mv STAR.* *.html STAR.SPLIT
                    mv STAR.SPLIT ..
                    cd ..
}


pe_bacteria_split(){

          cd trimmed_fastqs
          ls -1 *_1.fq.gz > .trR1
          ls -1 *_2.fq.gz > .trR2
          paste -d " " .trR1 .trR2 > Trimmed.list

          readarray trimmedFastqs < Trimmed.list

          for i in "${trimmedFastqs[@]}"

          do

            iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

            STAR \
            --runThreadN 12 \
            --genomeDir ${genomeDir[${DIR}]} \
            --readFilesIn $i \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${iSUB}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

          done

                    # source activate RSC
                    multiqc -f -n ${PIN}.star.multiqc.report .
                    mkdir STAR.SPLIT.COUNTS STAR.SPLIT.BAMS STAR.SPLIT.LOGS STAR.SPLIT.Unmapped
                    mv *Unmapped.out.mate* STAR.SPLIT.Unmapped
                    cd STAR.SPLIT.Unmapped
                        for i in *mate*
                            do
                                mv $i `echo $i | sed "s/Unmapped/not.$DIR/g"`
                            done
                    cd ..
                    mv *.ReadsPerGene.out.tab STAR.SPLIT.COUNTS
                    mv *.bam STAR.SPLIT.BAMS
                    mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SPLIT.LOGS
                    mkdir STAR.SPLIT
                    mv STAR.* *.html STAR.SPLIT
                    mv STAR.SPLIT ..
                    cd ..
}


se_bacteria_split(){
          cd trimmed_fastqs

          for i in *_trimmed.fq.gz

          do

            iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

            STAR \
            --runThreadN 12 \
            --genomeDir ${genomeDir[${DIR}]} \
            --readFilesIn $i \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${iSUB}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

          done

  #			source activate RSC
          multiqc -f -n ${PIN}.starSPLIT.multiqc.report .
          mkdir STAR.SPLIT.COUNTS STAR.SPLIT.BAMS STAR.SPLIT.LOGS STAR.SPLIT.Unmapped
          mv *Unmapped.out.mate* STAR.SPLIT.Unmapped
          cd STAR.SPLIT.Unmapped
              for i in *mate*
                  do
                      mv $i `echo $i | sed "s/Unmapped/not.$DIR/g"`
                  done
          cd ..
          mv *.ReadsPerGene.out.tab STAR.SPLIT.COUNTS
          mv *.bam STAR.SPLIT.BAMS
          mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SPLIT.LOGS
          mkdir STAR.SPLIT
          mv STAR.* *.html STAR.SPLIT
          mv STAR.SPLIT ..
          cd ..

}



UNMSE() {
          cd STAR.SPLIT/STAR.SPLIT.Unmapped/

          ls -1 *.mate* > unmapped.list

          grepOldRef=`head -1 unmapped.list | grep -o -P '(?<=.not.).*(?=.out.mate*)'`

          readarray unmappedSE < unmapped.list

          for i in "${unmappedSE[@]}"

          do

            iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

            STAR \
            --runThreadN 12 \
            --genomeDir ${genomeDir[${DIR}]} \
            --readFilesIn $i \
            --outSAMstrandField intronMotif \
            --outReadsUnmapped Fastx \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${iSUB}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

          done

          multiqc -f -n ${PIN}.star.unmapped.${grepOldRef}.mapped.to.${DIR}.multiqc.report .

          mkdir STAR.SPLIT.COUNTS STAR.SPLIT.BAMS STAR.SPLIT.LOGS STAR.SPLIT.Unmapped
          mv *Unmapped.out.mate* STAR.SPLIT.Unmapped
          cd STAR.SPLIT.Unmapped
              for i in *Unmapped.out.mate*
                  do
                      mv $i `echo $i | sed "s/Unmapped/not.$DIR/g"`
                  done
          cd ..
          mv *.ReadsPerGene.out.tab STAR.SPLIT.COUNTS
          mv *.bam STAR.SPLIT.BAMS
          mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SPLIT.LOGS
          mkdir STAR.SPLIT.unmapped.${grepOldRef}.mapped.to.${DIR}
          mv STAR.* *.html STAR.SPLIT.unmapped.${grepOldRef}.mapped.to.${DIR}
          # mv STAR.SPLIT ..
          cd ../../


}

UNMPE() {
          cd STAR.SPLIT/STAR.SPLIT.Unmapped/

          ls -1 *mate1* > .unR1
          ls -1 *mate2 > .unR2
          paste -d " " .unR1 .unR2 > unmapped.list

          grepOldRef=`head -1 .unR1 | grep -o -P '(?<=.not.).*(?=.out.mate*)'`

          readarray unmappedPE < unmapped.list

          for i in "${unmappedPE[@]}"

          do

            iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

            STAR \
            --runThreadN 12 \
            --genomeDir ${genomeDir[${DIR}]} \
            --readFilesIn $i \
            --outSAMstrandField intronMotif \
            --outReadsUnmapped Fastx \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${iSUB}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

          done

          multiqc -f -n ${PIN}.star.unmapped.${grepOldRef}.mapped.to.${DIR}.multiqc.report .

          mkdir STAR.SPLIT.COUNTS STAR.SPLIT.BAMS STAR.SPLIT.LOGS STAR.SPLIT.Unmapped
          mv *Unmapped.out.mate* STAR.SPLIT.Unmapped
          cd STAR.SPLIT.Unmapped
              for i in *Unmapped.out.mate*
                  do
                      mv $i `echo $i | sed "s/Unmapped/not.$DIR/g"`
                  done
          cd ..
          mv *.ReadsPerGene.out.tab STAR.SPLIT.COUNTS
          mv *.bam STAR.SPLIT.BAMS
          mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SPLIT.LOGS
          mkdir STAR.SPLIT.unmapped.${grepOldRef}.mapped.to.${DIR}
          mv STAR.* *.html STAR.SPLIT.unmapped.${grepOldRef}.mapped.to.${DIR}
          # mv STAR.SPLIT ..
          cd ../../


}


geneBodyCov(){
        cd STAR*/*.BAMS

        for i in *.bam
        do
          BASE=`basename $(echo $i) .Aligned.sortedByCoord.out.bam `
          mv $i ${BASE}.bam
        done

        for i in *.bam
        do
            /programs/bin/samtools/samtools index -b $i
        done
        cd ..
        echo
        echo
        pwd
        echo
        echo
        source activate RSeQC
        geneBody_coverage.py -r ${bed12[${DIR}]} -i *.BAMS/ -o ${PIN}
        mkdir geneBodyCov
        mv *geneBodyCoverage.* log.txt geneBodyCov
        cd ..
}

while getopts "hp:t:g:s:r:c:d:" opt; do
    case ${opt} in

    h)
        echo
        echo
        echo
        usage
        echo
        echo
        exit 1

    ;;

    p )

        PIN=$OPTARG
        echo "Project Identifier = " $PIN
    ;;

    t )

        T=$OPTARG

    ;;

    g)

        DIR=$OPTARG

    ;;

    s)

        STRAND=$OPTARG

    ;;

    r)

        RUN=$OPTARG

    ;;

    c)

    GBCOV=$OPTARG

    ;;

    d)
    DELIM=$OPTARG

    ;;

    \?)
        echo
        echo
        echo
        usage

    ;;

    esac

done
# shift $((OPTIND -1))
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if PIN is provided

if [[ -z "${PIN+x}" ]]; then

    PIN="PIN_Null"
fi


#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if delimiter parameter exists
if [[ ! -z "${DELIM+x}" ]]; then
    #statements
    if [[ $DELIM == default ]]; then

    DELIMITER="_"
    FIELD="5"
    echo "file naming will be done using the default delimiter settings"
  else

    DELIMITER=`echo $DELIM | cut -d , -f1`
    FIELD=`echo $DELIM | cut -d , -f2-`
    echo "file naming will be done using the delim = $DELIMITER and field = $FIELD settings"

  fi

fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if trimming parameter exists

if [[ ! -z "${T+x}" ]]; then
    #statements

    if   [[ $T == nextSE ]]; then
        trimSE
    elif [[ $T == nextPE ]]; then
        trimPE
    elif [[ $T == 4colorSE ]]; then
        trimHiSeq
    elif [[ $T == miSeqPE ]]; then
        trimMiSeqPE
  elif [[ $T == novaPE ]]; then
    trimHiSeqPE
    else
        echo "-t only accepts nextSE, nextPE, 4colorSE, miSeqPE, novaPE as arguments"
        exit 1
    fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if genomeDir provided

if [[ ! -z "${DIR+x}" ]]; then
    if [ ${genomeDir[${DIR}]+_} ]; then
        echo Reference genome selected = $DIR
        echo

        if [[ ! -z "${RUN+x}" ]] && [[ $RUN == "PE" ]]; then
                alignPE

            elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "PES" ]]; then
                pe_split

            elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "SE" ]]; then
                align

            elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "SES" ]]; then
                se_split

            elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "UNMS" ]]; then
                UNMSE

            elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "UNMP" ]]; then
                UNMPE

            elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "PEBS" ]]; then
              pe_bacteria_split

            elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "SEBS" ]]; then
              se_bacteria_split
            else
                echo "missing -r option "
                usage
                exit 1
        fi


    else
        echo "The reference genome provided '"$DIR"' is not available"
        echo " OR missing -r "
        exit 1

    fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if strandedness info provided

if [[ ! -z "${STRAND+x}" ]]; then
    if [[ $STRAND = "1" ]]; then
        echo first strand selected

        for i in STAR/STAR.COUNTS/*.ReadsPerGene.out.tab
            do
            awk 'NR > 4 {print $1 "\t" $3}' $i > $i.rawCounts
            done
            mkdir STAR/STAR.COUNTS/rawCounts
            mv STAR/STAR.COUNTS/*.rawCounts STAR/STAR.COUNTS/rawCounts

    elif [[ $STRAND = "2" ]]; then
        echo reverse strand selected

        for i in STAR/STAR.COUNTS/*.ReadsPerGene.out.tab
            do
            awk 'NR > 4 {print $1 "\t" $4}' $i > $i.rawCounts
            done
            mkdir STAR/STAR.COUNTS/rawCounts
            mv STAR/STAR.COUNTS/*.rawCounts STAR/STAR.COUNTS/rawCounts

    else
        echo unstranded selected

        for i in STAR/STAR.COUNTS/*.ReadsPerGene.out.tab
            do
            awk 'NR > 4 {print $1 "\t" $2}' $i > $i.rawCounts
            done
            mkdir STAR/STAR.COUNTS/rawCounts
            mv STAR/STAR.COUNTS/*.rawCounts STAR/STAR.COUNTS/rawCounts
    fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if run geneBodyCov True
if [[ ! -z "${GBCOV+x}" ]]; then
    if [[ $GBCOV = "yes" ]]; then
    geneBodyCov
  else
    echo "geneBodyCov not True"
  fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------

if [[ -z $1 ]] || [[  $1 = "help"  ]] ; then
    #statements
    echo
    echo
    usage
    echo
    echo
    exit 1

else
    echo >> beta6.run.log
    echo `date` >> beta6.run.log
    echo "Project Identifier Specified = " $PIN >> beta6.run.log
    echo "Reference Genome Specified   = " $DIR >> beta6.run.log
    echo "Trimming for smRNA seq       = " $T >> beta6.run.log
    echo "SE or PE                     = " $RUN >> beta6.run.log
    echo "Strandedness specified       = " $STRAND >> beta6.run.log
    echo "GeneBody Coverage            = " $GBCOV >> beta6.run.log
    echo >> beta6.run.log

    echo "ENV INFO: " >> beta6.run.log
    echo >> beta6.run.log
    echo "trim_galore -j 8  version:"`trim_galore -j 8  --version | grep 'version' | cut -d "n" -f2` >> beta6.run.log
    echo "STAR version:" `~/bin/STAR-2.7.0e/bin/Linux_x86_64/STAR --version` >> beta6.run.log
    echo "multiqc version:" `~/miniconda2/bin/multiqc --version` >> beta6.run.log
    echo "samtools version:" `/programs/bin/samtools/samtools --version` >> beta6.run.log
    echo "rseqc version: rseqc=2.6.4 " >> beta6.run.log
    echo -------------------------------------------------------------------------------------------------- >> beta6.run.log

fi
