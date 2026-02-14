# List the final desired files here
rule all:
    input:
        "results/tables/treatment_dea.csv",
        "results/plots/final_annotations.pdf"

rule get_data:
    output: "data/processed/scc_raw.rds"
    script: "scripts/01_get_data.R"

rule preprocess:
    input:  "data/processed/scc_raw.rds"
    output:
        rds="data/processed/scc_preprocessed.rds",
        plot="results/plots/qc_plots.pdf",
        vln_png="results/plots/qc_vln.png",
        scatter_png="results/plots/qc_scatter.png"
    script: "scripts/02_qc_and_preprocess.R"

rule clustering:
    input:  "data/processed/scc_preprocessed.rds"
    output:
        rds="data/processed/scc_clustered.rds",
        plots="results/plots/umap_plots.pdf",
        elbow_png="results/plots/elbow_plot.png",
        umap_png="results/plots/umap_treatment.png"
    script: "scripts/03_clustering_umap.R"

rule annotate_dea:
    input:  "data/processed/scc_clustered.rds"
    output:
        rds="data/processed/scc_final.rds",
        plots="results/plots/final_annotations.pdf",
        immgen_png="results/plots/umap_immgen.png",
        bcr_png="results/plots/feature_bcr.png",
        igkc_png="results/plots/feature_igkc.png",
        markers="results/tables/cluster_markers.csv",
        dea="results/tables/treatment_dea.csv"
    script: "scripts/04_annotation_and_dea.R"
