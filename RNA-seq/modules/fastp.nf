runmode = params.mode

process FASTPSEM {
    
    tag "$id"
    label 'process_medium'

    input:
        tuple val(id), path(reads)
    
    output:
        tuple val(id), path("*gz")         
        
    script:

    if ( runmode == "SE" || runmode == "SES" || runmode == "SEBS" ){
        
        """
        fastp \
        -z 4 -w 16 \
        --length_required 50 --qualified_quality_phred 20 \
        --trim_poly_g \
        -i ${reads} \
        -o ${id}_val_1.fq.gz \
        -h ${id}.fastp.html \
        -j ${id}.fastp.json
    
        """

    }

    else if ( runmode == "PE" || runmode == "PES" || runmode == "PEBS" ){

        """
            fastp \
            -z 4 -w 16 \
            --length_required 50 --qualified_quality_phred 20 \
            --trim_poly_g \
            -i ${reads[0]} \
            -I ${reads[1]} \
            -o ${id}_val_1.fq.gz \
            -O ${id}_val_2.fq.gz \
            -h ${id}.fastp.html \
            -j ${id}.fastp.json
        
        """


    } else {

        error "Runmode ${runmode} is not supported"
        exit 0
    }



}

