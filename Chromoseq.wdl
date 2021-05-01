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

  # array of genders. Will have to be prepared from excel spreadsheet
  Array[String] Genders

  # sample names
  Array[String] Samples

  # sample exceptions
  Array[String] Exceptions

  # name of batch
  String Batch

  # The root directory where all the batches are stored
  String SeqDir
  
  # The batch output directory
  String BatchDir = SeqDir + '/' + Batch
  
  String DBSNP     = "/staging/runs/Chromoseq/dragen_align_inputs/hg38/dbsnp.vcf.gz"
  String DOCM      = "/staging/runs/Chromoseq/dragen_align_inputs/hg38/docm.vcf.gz"
  String NoiseFile = "/staging/runs/Chromoseq/dragen_align_inputs/hg38/dragen_v1.0_systematic_noise.nextera_wgs.120920.bed.gz"

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

  String DragenReference    = "/staging/runs/Chromoseq/refdata/dragen_hg38"
  String DragenReferenceBED = "/staging/runs/Chromoseq/refdata/dragen_hg38/all_sequences.fa.bed.gz"
  String VEP

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

  Int MinValidatedReads
  Float MinValidatedVAF

  Int MinCovFraction
  Int MinGeneCov
  Int MinRegionCov
  
  String JobGroup
  String Queue
  String DragenQueue

  String chromoseq_docker
  String DragenDocker

  call dragen_demux {
    input: rundir=RunDir,
    FastqList=FastqList,
    BatchDir=BatchDir,
    Batch=Batch,
    sheet=SampleSheet,
    lane=Lane,
    jobGroup=JobGroup,
    queue=DragenQueue,
    docker=DragenDocker
  }

  scatter(i in range(length(Samples))){
    call dragen_align {
      input: BatchDir=BatchDir,
      Batch=Batch,
      fastqfile=dragen_demux.fastqfile,
      sample=Samples[i],
      gender=Genders[i],
      Reference=DragenReference,
      ReferenceBed=DragenReferenceBED,
      CNAbinsize=CNAbinsize,
      DBSNP=DBSNP,
      DOCM=DOCM,
      NoiseFile=NoiseFile,
      jobGroup=JobGroup,
      queue=DragenQueue,
      docker=DragenDocker
    }

    call subWF.ChromoseqAnalysis {
      input: Cram=dragen_align.cram,
      CramIndex=dragen_align.index,
      Name=Samples[i],
      Gender=Genders[i],
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

  String queue
  String docker
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
    docker_image: docker
    cpu: "20"
    memory: "200 G"
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

  String Reference
  String ReferenceBed
  Int CNAbinsize
  String DBSNP
  String DOCM
  String NoiseFile
  String queue
  String docker
  String jobGroup

  String outdir = BatchDir + "/" + sample
  String dragen_outdir = outdir + "/dragen"
  String LocalSampleDir = LocalAlignDir + "/" + sample
  String log = rootdir + "log/" + sample + "_align.log"
  
  command {
    if [ ! -d "${LocalAlignDir}" ]; then
      /bin/mkdir ${LocalAlignDir}
    fi

    if [ ! -d "${BatchDir}" ]; then
      /bin/mkdir ${BatchDir}
    fi

    /bin/mkdir ${LocalSampleDir} && \
    /bin/mkdir ${outdir} && \
    /opt/edico/bin/dragen -r ${Reference} --sample-sex ${gender} \
    --tumor-fastq-list ${fastqfile} --tumor-fastq-list-sample-id ${sample} \
    --enable-map-align-output true --enable-bam-indexing true --enable-duplicate-marking true \
    --enable-variant-caller true --dbsnp ${DBSNP}  --vc-somatic-hotspots ${DOCM} --vc-systematic-noise ${NoiseFile} \
    --enable-cnv true --cnv-target-bed ${ReferenceBed} --cnv-interval-width ${CNAbinsize} \
    --enable-sv true --sv-exome true --sv-output-contigs true --sv-hyper-sensitivity true \
    --output-format CRAM --output-directory ${LocalSampleDir} --output-file-prefix ${sample} &> ${log} && \
    /bin/mv ${log} ./ && \
    /bin/mv ${LocalSampleDir} ${dragen_outdir}
  }

  runtime {
    docker_image: docker
    cpu: "20"
    memory: "200 G"
    queue: queue
    job_group: jobGroup
  }

  output {
    File cram = "${dragen_outdir}/${sample}_tumor.cram"
    File index = "${dragen_outdir}/${sample}_tumor.cram.crai"
    File counts = "${dragen_outdir}/${sample}.target.counts.gz"
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
    docker_image: "ubuntu:xenial"
    queue: queue
    job_group: jobGroup
  }
  output {
    String done = stdout()
  }
}
