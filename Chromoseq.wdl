import "ChromoseqAnalysis.wdl" as subWF

workflow ChromoSeq {

  # This is the path of the illumina run directory on staging drive of dragen node after rsync from storage0
  String? RunDir

  # Illumina-compatible fastq_list.csv, take fastq list as input instead of RunDir
  String? FastqList
   
  # Illumina-compatible samplesheet for demuxing, prepared from launcher script using excel spreadsheet as input
  String? SampleSheet

  # lane, if one is specified
  String? Lane

  String? DragenEnv

  # array of genders, DOBs, sample_names, exceptions. Will have to be prepared from excel spreadsheet
  Array[String] Genders
  Array[String] DOBs
  Array[String] Samples
  Array[String] Exceptions

  # name of batch
  String Batch

  # The root directory where all the batches are stored
  String SeqDir
  
  # The batch output directory
  String BatchDir = SeqDir + '/' + Batch
  
  String DBSNP       = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_align_inputs/hg38/dbsnp.vcf.gz"
  String DOCM        = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_align_inputs/hg38/docm.vcf.gz"
  String POPAllele   = "/staging/runs/Chromoseq/dragen_align_inputs/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  String NoiseFile   = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_align_inputs/hg38/dragen_v1.0_systematic_noise.nextera_wgs.120920.bed.gz"
  String SvNoiseFile = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_align_inputs/hg38/WGS_v1.0.0_hg38_sv_systematic_noise.bedpe.gz"
  String NirvanaDB   = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_align_inputs/hg38/nirvana_annotation_data"

  String Translocations
  String GenesBed
  
  String Cytobands
  String SVDB

  String CustomAnnotationVcf 
  String CustomAnnotationIndex
  String CustomAnnotationParameters
  String? GeneFilterString
  
  String HotspotVCF
  String MantaConfig
  String MantaRegionConfig
  
  String HaplotectBed
  
  String Reference
  String ReferenceDict
  String ReferenceIndex
  String ReferenceBED

  String DragenReference    = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen393_hg38"
  String VEP

  String adapter1
  String adapter2

  String refWig
  String gcWig
  String mapWig
  String ponRds
  String centromeres
  String genomeStyle
  String genome

  String RefRangeJSON
  String RunInfoString
  String tmp
  
  Float minVarFreq
  Int MinReads
  Float varscanPvalindel
  Float varscanPvalsnv

  Int CNAbinsize = 500000
  Int MinCNASize = 5000000
  Float MinCNAabund = 10.0
  Int CNVmergedist
  Int CNVfilterlength

  Int MinValidatedReads
  Float MinValidatedVAF

  Int MinCovFraction
  Int MinGeneCov
  Int MinRegionCov
  
  String JobGroup
  String Queue
  String DragenQueue

  String chromoseq_docker
  String DragenDockerImage
  String DragenMEM
  Int DragenCPU

  String DemuxFastqDir = "/scratch1/fs1/gtac-mgi/CLE/chromoseq/demux_fastq"

  call dragen_demux {
    input: rundir=RunDir,
    FastqList=FastqList,
    BatchDir=BatchDir,
    Batch=Batch,
    sheet=SampleSheet,
    lane=Lane,
    jobGroup=JobGroup,
    queue=DragenQueue,
    DragenCPU=DragenCPU,
    DragenMEM=DragenMEM,
    DragenEnv=DragenEnv,
    DragenDockerImage=DragenDockerImage
  }

  scatter(i in range(length(Samples))){
    call dragen_align {
      input: BatchDir=BatchDir,
      Batch=Batch,
      fastqfile=dragen_demux.fastqfile,
      sample=Samples[i],
      gender=Genders[i],
      Reference=DragenReference,
      adapter1=adapter1,
      adapter2=adapter2,
      refWig=refWig,
      GeneBed=GenesBed,
      SVBed=Translocations,
      CNVmergedist=CNVmergedist,
      CNVfilterlength=CNVfilterlength,
      POPAllele=POPAllele,
      NirvanaDB=NirvanaDB,
      DBSNP=DBSNP,
      DOCM=DOCM,
      NoiseFile=NoiseFile,
      SvNoiseFile=SvNoiseFile,
      jobGroup=JobGroup,
      queue=DragenQueue,
      DragenCPU=DragenCPU,
      DragenMEM=DragenMEM,
      DragenEnv=DragenEnv,
      DragenDockerImage=DragenDockerImage
    }

    call subWF.ChromoseqAnalysis {
      input: Cram=dragen_align.cram,
      CramIndex=dragen_align.index,
      Name=Samples[i],
      Gender=Genders[i],
      DOB=DOBs[i],
      Exception=Exceptions[i],
      MappingSummary=dragen_align.mapping_summary,
      CoverageSummary=dragen_align.coverage_summary,
      TumorCounts=dragen_align.counts,
      Translocations=Translocations,
      GenesBed=GenesBed,
      Cytobands=Cytobands,
      SVDB=SVDB,
      CustomAnnotationVcf=CustomAnnotationVcf,
      CustomAnnotationIndex=CustomAnnotationIndex,
      CustomAnnotationParameters=CustomAnnotationParameters,
      GeneFilterString=GeneFilterString,
      HotspotVCF=HotspotVCF,
      MantaConfig=MantaConfig,
      MantaRegionConfig=MantaRegionConfig,
      HaplotectBed=HaplotectBed,
      Reference=Reference,
      ReferenceDict=ReferenceDict,
      ReferenceIndex=ReferenceIndex,
      ReferenceBED=ReferenceBED,
      VEP=VEP,
      gcWig=gcWig,
      mapWig=mapWig,
      ponRds=ponRds,
      centromeres=centromeres,
      genomeStyle=genomeStyle,
      genome=genome,
      minVarFreq=minVarFreq,
      MinReads=MinReads,
      varscanPvalindel=varscanPvalindel,
      varscanPvalsnv=varscanPvalsnv,
      CNAbinsize=CNAbinsize,
      MinCNASize=MinCNASize,
      MinCNAabund=MinCNAabund,
      MinValidatedReads=MinValidatedReads,
      MinValidatedVAF=MinValidatedVAF,
      MinCovFraction=MinCovFraction,
      MinGeneCov=MinGeneCov,
      MinRegionCov=MinRegionCov,
      RefRangeJSON=RefRangeJSON,
      RunInfoString=RunInfoString,
      tmp=tmp,
      JobGroup=JobGroup,
      Queue=Queue,
      chromoseq_docker=chromoseq_docker,
      OutputDir=BatchDir + "/" + Samples[i]
    }
  }
  
  call batch_qc {
    input: order_by=ChromoseqAnalysis.all_done,
    BatchDir=BatchDir,
    queue=Queue,
    jobGroup=JobGroup,
    docker=chromoseq_docker
  }

  call move_demux_fastq {
    input: order_by=ChromoseqAnalysis.all_done,
    Batch=Batch,
    DemuxFastqDir=DemuxFastqDir,
    queue=DragenQueue,
    jobGroup=JobGroup
  }

  call remove_rundir {
    input: order_by=ChromoseqAnalysis.all_done,
    rundir=RunDir,
    queue=DragenQueue,
    jobGroup=JobGroup
  }
}

task dragen_demux {
  String Batch
  String BatchDir
  String DemuxReportDir = BatchDir + "/dragen_demux_reports"

  String rootdir = "/staging/runs/Chromoseq/"
  String LocalFastqDir = rootdir + "demux_fastq/" + Batch
  String LocalDemuxReportDir = LocalFastqDir + "/Reports"
  String LocalFastqList = rootdir + "sample_sheet/" + Batch + '_fastq_list.csv'
  String LocalSampleSheet = rootdir + "sample_sheet/" + Batch + '.csv'
  String log = rootdir + "log/" + Batch + "_demux.log"

  String? rundir
  String? FastqList
  String? sheet
  String? lane

  String? DragenEnv

  String DragenDockerImage
  String DragenMEM
  Int DragenCPU

  String queue
  String jobGroup
  
  command {
    if [ -n "${FastqList}" ]; then
      /bin/cp ${FastqList} ${LocalFastqList}
    else
      /bin/cp ${sheet} ${LocalSampleSheet}

      if [ -n "${lane}" ]; then
        /opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --strict-mode true --sample-sheet ${LocalSampleSheet} --bcl-input-directory ${rundir} --output-directory ${LocalFastqDir} --bcl-only-lane ${lane} &> ${log}
      else
        /opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --strict-mode true --sample-sheet ${LocalSampleSheet} --bcl-input-directory ${rundir} --output-directory ${LocalFastqDir} &> ${log}
      fi

      /bin/mv ${log} ./ && \
      /bin/rm -f ${LocalSampleSheet} && \
      /bin/cp -r ${LocalDemuxReportDir} ${DemuxReportDir} && \
      /bin/cp "${LocalFastqDir}/Reports/fastq_list.csv" ${LocalFastqList}
    fi
  }

  runtime {
    docker_image: DragenDockerImage
    dragen_env: DragenEnv
    cpu: DragenCPU
    memory: DragenMEM
    queue: queue
    job_group: jobGroup
  }

  output {
    String fastqfile = "${LocalFastqList}"
  }
}

task dragen_align {
  String Batch
  String BatchDir

  String rootdir = "/staging/runs/Chromoseq/"
  String LocalAlignDir = rootdir + "align/" + Batch

  String fastqfile
  String sample
  String gender

  String adapter1
  String adapter2

  Int CNVmergedist
  Int CNVfilterlength

  String Reference
  
  String refWig
  String GeneBed
  String SVBed
  String DBSNP
  String DOCM
  String NirvanaDB
  String NoiseFile
  String POPAllele
  String SvNoiseFile

  String CovLevels = '10,20,40,60,80'

  String? DragenEnv
  String DragenDockerImage
  String queue
  String jobGroup
  String DragenMEM
  Int DragenCPU

  String outdir = BatchDir + "/" + sample
  String dragen_outdir = outdir + "/dragen"
  String LocalSampleDir = LocalAlignDir + "/" + sample
  String log = rootdir + "log/" + sample + "_align.log"
  
  command <<<
    if [ ! -d "${LocalAlignDir}" ]; then
      /bin/mkdir ${LocalAlignDir}
    fi

    if [ ! -d "${BatchDir}" ]; then
      /bin/mkdir ${BatchDir}
    fi

    awk -v OFS="\t" '{ split($7,a,"_"); print $1,$2,$3,a[1],".",$9; print $4,$5,$6,a[2],".",$10; }' ${SVBed} | sort -u -k 1,1V -k 2,2n > sv.bed

    /bin/mkdir ${LocalSampleDir} && \
    /bin/mkdir ${outdir} && \
    /opt/edico/bin/dragen -r ${Reference} --sample-sex ${gender} \
    --read-trimmers adapter --trim-adapter-read1 ${adapter1} --trim-adapter-read2 ${adapter2} \
    --tumor-fastq-list ${fastqfile} --tumor-fastq-list-sample-id ${sample} \
    --qc-coverage-region-1 ${refWig} --qc-coverage-filters-1 'mapq<0,bq<0' \
    --qc-coverage-region-2 ${GeneBed} --qc-coverage-reports-2 cov_report --qc-coverage-region-2-thresholds ${CovLevels} \
    --qc-coverage-region-3 sv.bed --qc-coverage-reports-3 cov_report --qc-coverage-region-3-thresholds ${CovLevels} \
    --qc-coverage-ignore-overlaps=true \
    --gc-metrics-enable true --enable-map-align-output true --enable-bam-indexing true --enable-duplicate-marking true \
    --enable-variant-caller true --dbsnp ${DBSNP} --vc-somatic-hotspots ${DOCM} --vc-systematic-noise ${NoiseFile} --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 \
    --enable-variant-annotation true --variant-annotation-assembly GRCh38 --variant-annotation-data ${NirvanaDB} \
    --enable-cnv true --cnv-somatic-enable-het-calling true --cnv-enable-ref-calls false --cnv-merge-distance ${CNVmergedist} --cnv-filter-length ${CNVfilterlength} --cnv-population-b-allele-vcf ${POPAllele} \
    --enable-sv true --sv-output-contigs true --sv-hyper-sensitivity true --sv-min-edge-observations 2 --sv-min-candidate-spanning-count 1 \
    --sv-use-overlap-pair-evidence true --sv-enable-somatic-ins-tandup-hotspot-regions true --sv-systematic-noise ${SvNoiseFile} \
    --output-format CRAM --output-directory ${LocalSampleDir} --output-file-prefix ${sample} &> ${log} && \
    /bin/mv ${log} ./ && \
    /bin/mv ${LocalSampleDir} ${dragen_outdir}
    
  >>>

  runtime {
    docker_image: DragenDockerImage
    dragen_env: DragenEnv
    cpu: DragenCPU
    memory: DragenMEM
    queue: queue
    job_group: jobGroup
  }

  output {
    File cram = "${dragen_outdir}/${sample}_tumor.cram"
    File index = "${dragen_outdir}/${sample}_tumor.cram.crai"
    File counts = "${dragen_outdir}/${sample}.qc-coverage-region-1_read_cov_report.bed"
    File mapping_summary = "${dragen_outdir}/${sample}.mapping_metrics.csv"
    File coverage_summary = "${dragen_outdir}/${sample}.wgs_coverage_metrics.csv"
  }
}

task batch_qc {
  Array[String] order_by
  String BatchDir
  String queue
  String jobGroup
  String docker

  String qcOut = BatchDir + "/QC_info.txt"

  command {
    if [ -n "$(/bin/ls -d ${BatchDir}/TWDY-*)" ]; then
        /bin/chmod -R 777 ${BatchDir}/TWDY-*
    fi

    /usr/bin/perl /usr/local/bin/QC_info.pl "${BatchDir}/*/*.chromoseq.txt" > ${qcOut}
  }
  runtime {
    docker_image: docker
    memory: "4 G"
    queue: queue
    job_group: jobGroup
  }
  output {
    String done = stdout()
  }
}

task move_demux_fastq {
  Array[String] order_by
  String Batch
  String DemuxFastqDir
  String queue
  String jobGroup

  String LocalDemuxFastqDir = "/staging/runs/Chromoseq/demux_fastq/" + Batch

  command {
    if [ -n "${LocalDemuxFastqDir}" ]; then
      /bin/mv ${LocalDemuxFastqDir} ${DemuxFastqDir}
    fi
  }
  runtime {
    docker_image: "docker1(ubuntu:xenial)"
    queue: queue
    job_group: jobGroup
  }
  output {
    String done = stdout()
  }
}

task remove_rundir {
  Array[String] order_by
  String? rundir
  String queue
  String jobGroup
  
  command {
    if [ -n "${rundir}" ]; then 
      /bin/rm -Rf ${rundir}
    fi
  }
  runtime {
    docker_image: "docker1(ubuntu:xenial)"
    queue: queue
    job_group: jobGroup
  }
  output {
    String done = stdout()
  }
}
