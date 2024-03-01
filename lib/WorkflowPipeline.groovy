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
        if (!params.references.primers.keySet().contains(params.primer_set)) {
            log.info "The primer set ${params.primer_set} is not currently configured."
            System.exit(1)
        }
    }

}
