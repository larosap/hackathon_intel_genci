#!/usr/bin/env nextflow

/*
kate: syntax groovy; space-indent on; indent-width 2;
================================================================================
=                                 S  A  R  E  K                                =
================================================================================
 New Germline (+ Somatic) Analysis Workflow. Started March 2016.
--------------------------------------------------------------------------------
 @Authors
 Sebastian DiLorenzo <sebastian.dilorenzo@bils.se> [@Sebastian-D]
 Jesper Eisfeldt <jesper.eisfeldt@scilifelab.se> [@J35P312]
 Phil Ewels <phil.ewels@scilifelab.se> [@ewels]
 Maxime Garcia <maxime.garcia@scilifelab.se> [@MaxUlysse]
 Szilveszter Juhos <szilveszter.juhos@scilifelab.se> [@szilvajuhos]
 Max Käller <max.kaller@scilifelab.se> [@gulfshores]
 Malin Larsson <malin.larsson@scilifelab.se> [@malinlarsson]
 Marcel Martin <marcel.martin@scilifelab.se> [@marcelm]
 Björn Nystedt <bjorn.nystedt@scilifelab.se> [@bjornnystedt]
 Pall Olason <pall.olason@scilifelab.se> [@pallolason]
--------------------------------------------------------------------------------
 @Homepage
 http://opensource.scilifelab.se/projects/sarek/
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/SciLifeLab/Sarek/README.md
--------------------------------------------------------------------------------
 Processes overview
 - GetVersionAll - Get version of tools
 - RunMultiQC - Run MultiQC on reports
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

if (params.help) exit 0, helpMessage()
if (!SarekUtils.isAllowedParams(params)) exit 1, "params unknown, see --help for more information"
if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"
step = params.step.toLowerCase()
tools = params.tools

// Check for awsbatch profile configuration
// make sure queue is defined
if (workflow.profile == 'awsbatch') {
    if (!params.awsqueue) exit 1, "Provide the job queue for aws batch!"
}

params.sampleDir = false
params.name = false
params.monochrome_logs = false
//params.tools = false
params.skipQC = false
//params.outdir = '$baseDir/Reports'

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


/*
================================================================================
                                PRINTING SUMMARY
================================================================================
*/

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision)          summary['Pipeline Release']    = workflow.revision
summary['Run Name']          = custom_runName ?: workflow.runName
if (params.sample)           summary['Sample File']        = params.sample
if (params.sampleDir)           summary['Sample Dir']        = params.sampleDir
summary['Genome']       = params.genome
summary['Genome Base']  = params.genome_base
summary['Sequencing Center']  = params.sequencing_center
if (params.targetBED)           summary['Target BED']        = params.targetBED
if (params.step)                summary['Step']              = params.step
if (params.tools)               summary['Tools']             = "${params.tools}" 
if (params.skipQC)              summary['QC tools skip']     = skipQC.join(', ')

if ('haplotypecaller' in tools)              summary['GVCF']              = params.noGVCF ? 'No' : 'Yes'
if ('strelka' in tools && 'manta' in tools ) summary['Strelka BP']        = params.noStrelkaBP ? 'No' : 'Yes'
if (params.sequencing_center)                summary['Sequenced by']      = params.sequencing_center
if (params.pon && 'mutect2' in tools)        summary['Panel of normals']  = params.pon
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']        = "${params.outDir}/Reports"
summary['Launch dir']        = workflow.launchDir
summary['Working dir']       = workflow.workDir
summary['Script dir']        = workflow.projectDir
summary['User']              = workflow.userName

if (params.acLoci)              summary['acLoci']            = params.acLoci
if (params.acLociGC)            summary['acLociGC']          = params.acLociGC
if (params.chrDir)              summary['chrDir']            = params.chrDir
if (params.chrLength)           summary['chrLength']         = params.chrLength
if (params.dbsnp)               summary['dbsnp']             = params.dbsnp
if (params.fasta)               summary['fasta']             = params.fasta
if (params.germlineResource)    summary['germlineResource']  = params.germlineResource
if (params.intervals)           summary['intervals']         = params.intervals
if (params.knownIndels)         summary['knownIndels']       = params.knownIndels.join(', ')
if (params.snpeffDb)            summary['snpeffDb']          = params.snpeffDb
if (params.vepCacheVersion)     summary['vepCacheVersion']   = params.vepCacheVersion

if (workflow.profile == 'awsbatch') {
    summary['AWS Region']        = params.awsregion
    summary['AWS Queue']         = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description)  summary['Config Description']  = params.config_profile_description
if (params.config_profile_contact)      summary['Config Contact']      = params.config_profile_contact
if (params.config_profile_url)          summary['Config URL']          = params.config_profile_url
if (params.email) {
    summary['E-mail Address']        = params.email
    summary['MultiQC maxsize']       = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k, v -> "${k.padRight(18)}: $v" }.join("\n")
if (params.monochrome_logs) log.info "----------------------------------------------------"
else log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

process GetVersionAll {
  publishDir "${params.outDir}/Reports/MultiQC", mode: params.publishDirMode

  input:
    file(versions) from Channel.fromPath("${params.outDir}/Reports/ToolsVersion/*").collect().ifEmpty(null)

  output:
    file ("tool_versions_mqc.yaml") into versionsForMultiQC

  when: !params.noReports

  script:
  """
  bcftools version > v_bcftools.txt 2>&1 || true
  bwa &> v_bwa.txt 2>&1 || true
  configManta.py --version > v_manta.txt 2>&1 || true
  configureStrelkaGermlineWorkflow.py --version > v_strelka.txt 2>&1 || true
  echo "${workflow.manifest.version}" &> v_sarek.txt 2>&1 || true
  echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
  echo "SNPEFF version"\$(snpEff -h 2>&1) > v_snpeff.txt
  fastqc -v > v_fastqc.txt 2>&1 || true
  freebayes --version > v_freebayes.txt 2>&1 || true
  gatk ApplyBQSR --help 2>&1 | grep Version: > v_gatk.txt 2>&1 || true
  multiqc --version &> v_multiqc.txt 2>&1 || true
  qualimap --version &> v_qualimap.txt 2>&1 || true
  samtools --version &> v_samtools.txt 2>&1 || true
  vcftools --version &> v_vcftools.txt 2>&1 || true
  echo "VEP 95.2" &> v_vep.txt

  scrape_tool_versions.py &> tool_versions_mqc.yaml
  """
}

if (params.verbose && !params.noReports) versionsForMultiQC = versionsForMultiQC.view {
  "MultiQC tools version:\n\
  File  : [${it.fileName}]"
}

reportsForMultiQC = Channel.empty()
  .mix(
    Channel.fromPath("${params.outDir}/Reports/bamQC/*", type: 'dir'),
    Channel.fromPath("${params.outDir}/Reports/BCFToolsStats/*"),
    Channel.fromPath("${params.outDir}/Reports/FastQC/*/*"),
    Channel.fromPath("${params.outDir}/Reports/MarkDuplicates/*"),
    Channel.fromPath("${params.outDir}/Reports/SamToolsStats/*"),
    Channel.fromPath("${params.outDir}/Annotation/*/snpEff/*.csv"),
    Channel.fromPath("${params.outDir}/Reports/VCFTools/*"),
  ).collect()

process RunMultiQC {
  publishDir "${params.outDir}/Reports/MultiQC", mode: params.publishDirMode

  input:
    file (multiqcConfig) from createMultiQCconfig()
    file (reports) from reportsForMultiQC
    file (versions) from versionsForMultiQC

  output:
    set file("*multiqc_report.html"), file("*multiqc_data") into multiQCReport

  when: !params.noReports

  script:
  """
  multiqc -f -v .
  """
}

if (params.verbose) multiQCReport = multiQCReport.view {
  "MultiQC report:\n\
  File  : [${it[0].fileName}]\n\
  Dir   : [${it[1].fileName}]"
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkUppmaxProject() {
  // check if UPPMAX project number is specified
  return !(workflow.profile == 'slurm' && !params.project)
}

def createMultiQCconfig() {
  def file = workDir.resolve('multiqc_config.yaml')
  file.text  = """
  custom_logo: ${baseDir}/docs/images/Sarek_no_Border.png
  custom_logo_url: http://opensource.scilifelab.se/projects/sarek
  custom_logo_title: 'Sarek'
  report_header_info:
  - Contact Name: ${params.callName}
  - Contact E-mail: ${params.contactMail}
  - Genome: ${params.genome}
  top_modules:
  - 'fastqc'
  - 'picard'
  - 'samtools'
  - 'qualimap'
  - 'bcftools'
  - 'vcftools'
  - 'snpeff'
  """.stripIndent()

  return file
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def helpMessage() {
  // Display help message
  this.sarekMessage()
  log.info "    Usage:"
  log.info "       nextflow run runMultiQC.nf"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outDir
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Containers"
  if (params.repository != "") log.info "  Repository   : " + params.repository
  if (params.containerPath != "") log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def sarekMessage() {
  // Display Sarek message
  log.info "Sarek - Workflow For Somatic And Germline Variations ~ ${workflow.manifest.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def startMessage() {
  // Display start message
  SarekUtils.sarek_ascii()
  this.sarekMessage()
  this.minimalInformationMessage()
}
def endMessage() {
  // Display complete message
  this.nextflowMessage()
  this.sarekMessage()
  this.minimalInformationMessage()
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
def nfcoreHeader() {
    // Log colors ANSI codes
    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_dim    = params.monochrome_logs ? '' : "\033[2m";
    c_black  = params.monochrome_logs ? '' : "\033[0;30m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue   = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan   = params.monochrome_logs ? '' : "\033[0;36m";
    c_white  = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
        ${c_white}____${c_reset}
      ${c_white}.´ _  `.${c_reset}
     ${c_white}/  ${c_green}|\\${c_reset}`-_ \\${c_reset}     ${c_blue} __        __   ___     ${c_reset}
    ${c_white}|   ${c_green}| \\${c_reset}  `-|${c_reset}    ${c_blue}|__`  /\\  |__) |__  |__/${c_reset}
     ${c_white}\\ ${c_green}|   \\${c_reset}  /${c_reset}     ${c_blue}.__| /¯¯\\ |  \\ |___ |  \\${c_reset}
      ${c_white}`${c_green}|${c_reset}____${c_green}\\${c_reset}´${c_reset}

    ${c_purple}  nf-core/sarek v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/sarek] Successful: ${workflow.runName}"
    if (!workflow.success) subject = "[nf-core/sarek] FAILED: ${workflow.runName}"

    summary['Duration']     = workflow.duration
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
// Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    //def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]

    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" , mqcFile: mqc_report]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/sarek] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, params.email ].execute() << email_txt
            log.info "[nf-core/sarek] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "$baseDir/Reports/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File(output_d, "Sarek_${step}_pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "Sarek_${step}_pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es)${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt}${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt}${c_reset}"
    }

    if (workflow.success) log.info "${c_purple}[nf-core/sarek]${c_green} Pipeline completed successfully${c_reset}"
    else {
        checkHostname()
        log.info "${c_purple}[nf-core/sarek]${c_red} Pipeline completed with errors${c_reset}"
    }

  endMessage()
}

workflow.onError {
  // Display error message
  this.nextflowMessage()
  this.sarekMessage()
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
