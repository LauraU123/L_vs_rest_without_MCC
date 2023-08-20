configfile: "config/configfile.yaml"

build_dir = 'results_folder'
auspice_dir = 'auspice'

rule all:
    input:
        expand("auspice/rsv_{a_or_b}_onlyL.json",
        a_or_b = config.get("subtypes",['b'])),
        expand("auspice/rsv_{a_or_b}_onlyNS1-M.json",
        a_or_b = config["subtypes"])


include: "workflow/snakemake_rules/core.smk"

include: "workflow/snakemake_rules/export.smk"

include: "workflow/snakemake_rules/download.smk"


rule clean:
    params:
        targets = ["auspice", "results"]
    shell:
        """
        rm -rf {params.targets}
        """

rule clobber:
    params:
        targets = ["data", "auspice", "results"]
    shell:
        """
        rm -rf {params.targets}
        """
