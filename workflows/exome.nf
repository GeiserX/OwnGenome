/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The --tools names an exome run leaves out
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    With --model_type WES the reads cover the capture regions, about 2% of
    the genome. The steps that read depth along the whole genome (CNVpytor,
    TelomereHunter, viral DNA in the unmapped reads), or weigh a locus's
    depth against a background spread over it (Cyrius, Parascopy,
    Locityper, MEGAnE), or read variants spread over it (ROH, archaic
    segments, MoChA) have no input to read, and their results would look
    like answers. Cyrius, Locityper and MEGAnE are whole-genome tools by
    their own documentation; Parascopy's background windows lie across
    GRCh38. They are not run, and the run log says "not run: exome"
    (main.nf). survivor_merge goes with them when CNVpytor was one of its
    two callers. docs/wes-guide.md lists them with the reason.

    main.nf takes them out of its tool list, and the workflows ask
    selected_tools(), so both read the same rule.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The --tools names a whole genome runs and an exome does not
def exome_skipped_tools() {
    return ['cnvpytor', 'roh', 'telomere_hunter', 'archaic', 'mocha', 'viral',
            'cyrius', 'parascopy', 'locityper', 'mei']
}

// The names of TOOLS (a --tools list) this run does not run because it is an exome
def exome_dropped(List tools) {
    if (params.model_type != 'WES') {
        return []
    }
    def dropped = tools.findAll { t -> t in exome_skipped_tools() }
    def callers = ['manta', 'delly', 'cnvpytor'].findAll { c -> tools.contains(c) && !dropped.contains(c) }
    if (tools.contains('survivor_merge') && dropped.contains('cnvpytor') && callers.size() < 2) {
        dropped << 'survivor_merge'
    }
    return dropped
}

// --tools as this run runs it: the names given, less the ones an exome leaves out
def selected_tools() {
    def tools = params.tools ? params.tools.toString().split(',').collect { t -> t.trim() }.findAll { t -> t } : []
    return tools - exome_dropped(tools)
}
