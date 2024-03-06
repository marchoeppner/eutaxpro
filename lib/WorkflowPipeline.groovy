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
        if (!params.reference_base && !params.build_references) {
            log.info "No local taxonomy reference specified - downloading on-the-fly instead..."
            log.info "Consider installing the reference(s) as specified in our documentation!"
        }
        if(!params.primer_set) {
            log.info "No primer set specified (--primer_set) - cannot work without one ..."
            System.exit(1)
        }
        
    }

}
