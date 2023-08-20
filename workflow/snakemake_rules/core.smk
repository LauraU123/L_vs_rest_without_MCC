'''
This part of the workflow expects input files
            sequences = "data/sequences.fasta"
            metadata = "data/metadata.tsv"
'''

L_or_rest = ["onlyNS1-M", "onlyL"]
L_OR_REST = ["onlyNS1-M", "onlyL"]
rule wrangle_metadata:
    input:
        metadata="data/{a_or_b}/metadata.tsv",
    output:
        metadata="data/{a_or_b}/metadata_by_accession.tsv"
    params:
        strain_id=lambda w: config.get("strain_id_field", "strain"),
    shell:
        """
        python3 scripts/wrangle_metadata.py --metadata {input.metadata} \
                    --strain-id {params.strain_id} \
                    --output {output.metadata}
        """

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = "data/{a_or_b}/sequences.fasta"
    output:
        sequence_index = build_dir + "/{a_or_b}/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule newreference:
    message:
        """
        Making new reference
        """
    input:
        oldreference = "config/{a_or_b}reference.gbk"
    output:
        newreferencegbk = build_dir + "/{a_or_b}/{gene}_reference.gbk",
        newreferencefasta = build_dir + "/{a_or_b}/{gene}_reference.fasta",
    params:
        gene = lambda w: w.gene,
    shell:
        """
        python scripts/newreference.py \
            --reference {input.oldreference} \
            --output-genbank {output.newreferencegbk} \
            --output-fasta {output.newreferencefasta} \
            --gene {params.gene}
        """

rule filter:
    message:
        """
        filtering sequences
        """
    input:
        sequences = "data/{a_or_b}/sequences.fasta",
        reference = "config/{a_or_b}reference.gbk",
        metadata = "data/{a_or_b}/metadata_by_accession.tsv",
        sequence_index = rules.index_sequences.output,
        exclude = config['exclude']
    output:
    	sequences = build_dir + "/{a_or_b}/filtered.fasta"
    params:
    	group_by = config["filter"]["group_by"],
    	min_coverage = lambda w: f'genome_coverage>{config["filter"]["min_coverage"].get("genome", 10000)}',
    	subsample_max_sequences = lambda w: config["filter"]["subsample_max_sequences"].get("genome", 1000),
        strains = "config/dropped_strains.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --query '{params.min_coverage}' \
            --exclude {params.strains}
        """

rule nextalign:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = build_dir + "/{a_or_b}/genome_reference.fasta"
    output:
        alignment = build_dir + "/{a_or_b}/sequences.aligned.fasta"
    threads: 4
    shell:
        """
        nextalign run -j {threads}\
            --reference {input.reference} \
            --output-fasta {output.alignment} \
            {input.sequences}
        """

rule cut:
    input:
        oldalignment = rules.nextalign.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        slicedalignment = build_dir + "/{a_or_b}/{gene}_slicedalignment.fasta"
    params:
        gene = lambda w: w.gene
    shell:
        """
        python scripts/cut.py \
            --oldalignment {input.oldalignment} \
            --slicedalignment {output.slicedalignment} \
            --reference {input.reference} \
            --gene {params.gene}
        """

rule realign:
    input:
        slicedalignment = rules.cut.output.slicedalignment,
        reference = build_dir + "/{a_or_b}/{gene}_reference.fasta"
    output:
        realigned = build_dir + "/{a_or_b}/{gene}gene_aligned.fasta"
    threads: 4
    shell:
        """
        augur align --nthreads {threads} \
            --sequences {input.slicedalignment} \
            --reference-sequence {input.reference} \
            --output {output.realigned}
        """


rule hybrid_align:
    input:
        original = rules.nextalign.output.alignment,
        G_alignment = build_dir + "/{a_or_b}/Ggene_aligned.fasta",
        reference = "config/areference.gbk"
    output:
        hybrid_alignment = build_dir + "/{a_or_b}/hybrid_alignment.fasta"
    params:
        gene = "genome"
    shell:
        """
        python scripts/align_for_tree.py \
            --realign {input.G_alignment} \
            --original {input.original} \
            --reference {input.reference} \
            --output {output.hybrid_alignment} \
            --gene {params.gene}
        """

#def get_alignment(w):
#    if w.build_name == "genome":
#        return rules.hybrid_align.output.hybrid_alignment
#    else:
#        return build_dir + f"/a/genome_aligned.fasta"


#rule genes:
#    input:
#        sequences = build_dir + "/a/hybrid_alignment.fasta"
#    output:
#        l = build_dir + "/a/onlySN1-M%GENE.fasta"
#    shell:
#        """
#        python "scripts/onlySN1-Msequences.py" \
#        --alignment {input.sequences} \
#        --output {output.l}
#        """

rule split:
    input:
        sequences = build_dir + "/{a_or_b}/hybrid_alignment.fasta",
        reference = "config/{a_or_b}reference.gbk",
    output:
        l = build_dir + "/{a_or_b}/onlyLGENE.fasta",
        rest = build_dir + "/{a_or_b}/onlyNS1-MGENE.fasta"
    shell:
        """
        python "scripts/split_aligned.py" \
        --oldalignment {input.sequences} \
        --loutput {output.l} \
        --restoutput {output.rest} \
        --reference {input.reference}
        """

rule trimL:
    """masking sequences"""
    input:
        alignment = build_dir + "/{a_or_b}/onlyLGENE.fasta"
    output:
        alignment = build_dir + "/{a_or_b}/onlyLGENE_trimmed.fasta"
    shell:
        """
        augur mask \
        --mask-from-end 150 \
        --sequences {input.alignment} \
        --output {output.alignment}
        """

rule trimRest:
    """masking sequences"""
    input:
        alignment = build_dir + "/{a_or_b}/onlyNS1-MGENE.fasta"
    output:
        alignment = build_dir + "/{a_or_b}/onlyNS1-MGENE_trimmed.fasta"
    shell:
        """
        augur mask \
        --mask-from-beginning 120 \
        --sequences {input.alignment} \
        --output {output.alignment}
        """

        
rule tree:
    message: "Building tree"
    input:
        alignment = build_dir + "/{a_or_b}/{L_or_rest}GENE_trimmed.fasta"
    output:
        tree = build_dir + "/{a_or_b}/{L_or_rest}tree_raw.nwk"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
        """

rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = rules.tree.input.alignment,
        metadata = rules.filter.input.metadata
    output:
        tree = build_dir + "/{a_or_b}/{L_or_rest}tree_refined.nwk",
    params:
        coalescent = config["refine"]["coalescent"],
        clock_filter_iqd = config["refine"]["clock_filter_iqd"],
        date_inference = config["refine"]["date_inference"]
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --timetree \
            --clock-filter-iqd {params.clock_filter_iqd}
        """



rule branch_lengths:
    message:
        """
        Refining tree
        """
    input:
        tree = build_dir + "/{a_or_b}/{L_or_rest}tree_refined.nwk",
    output:
        node_data = build_dir + "/{a_or_b}/{L_or_rest}branch_lengths.json"
    shell:
        """
        python3 scripts/make-branch-lengths.py \
            --tree {input.tree} \
            --output-node-data {output.node_data} \
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree =  build_dir + "/{a_or_b}/{L_or_rest}tree_refined.nwk",
        alignment = rules.tree.input.alignment
    output:
        node_data = build_dir + "/{a_or_b}/{L_or_rest}nt_muts.json"
    params:
    	inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = build_dir + "/{a_or_b}/{L_or_rest}tree_refined.nwk",
        node_data = rules.ancestral.output.node_data,
        reference = build_dir + "/{a_or_b}/genome_reference.gbk",
    output:
        node_data = build_dir + "/{a_or_b}/{L_or_rest}aa_muts.json"
    params:
    	alignment_file_mask = build_dir + "/{a_or_b}/aligned_{L_or_rest}%GENE.fasta"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
            --alignment-output {params.alignment_file_mask}
        """

rule traits:
    input:
        tree = build_dir + "/{a_or_b}/{L_or_rest}tree_refined.nwk",
        metadata = rules.filter.input.metadata
    output:
        node_data = build_dir + "/{a_or_b}/{L_or_rest}traits.json"
    log:
        "logs/{a_or_b}/{L_or_rest}traits_genome_rsv.txt"
    params:
    	columns = config["traits"]["columns"]
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = build_dir + "/{a_or_b}/{L_or_rest}tree_refined.nwk",
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = "config/clades_G_a.tsv"
    output:
        node_data = build_dir + "/{a_or_b}/{L_or_rest}clades_G.json"
    log:
        "logs/{a_or_b}/{L_or_rest}_clades_genome.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """
