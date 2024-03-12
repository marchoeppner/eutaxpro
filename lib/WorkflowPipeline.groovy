//
// This file holds several functions specific to the workflow/esga.nf in the nf-core/esga pipeline
//

class WorkflowPipeline {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (!params.run_name) {
            log.info 'Must provide a run_name (--run_name)'
            System.exit(1)
        }
        if (!params.input && !params.build_references) {
            log.info 'This pipeline requires a sample sheet as input (--input)'
            System.exit(1)
        }
        if (!params.reference_base && !params.build_references) {
            log.info 'No local taxonomy reference specified - downloading on-the-fly instead...'
            log.info 'Consider installing the reference(s) as specified in our documentation!'
        }
        if (!params.build_references) {
            if (params.primer_set && !params.references.primers.keySet().contains(params.primer_set)) {
                log.info "The primer set ${params.primer_set} is not currently configured."
                System.exit(1)
            }
            if (!params.primer_set && !params.gene) {
                log.info 'You have to specify which gene you are targeting if you do not use a built-in primer set (--gene)'
                System.exit(1)
            }
            if (!params.primer_set && !params.primers_txt && !params.primers_fa) {
                log.info 'No primer set (--primer_set) or custom primer configuration (--primers_txt) provided. Exiting...'
                System.exit(1)
            }
            if (!params.primer_set && !params.primers_txt && !params.primers_fa) {
                log.info 'No primer information provided, exiting...'
                System.exit(1)
            }
        }
    }

}
