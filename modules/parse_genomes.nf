
// Import modules fetchFTP for genome FASTA files
include {
    fetchFTP;
} from './modules' addParams(
    ftp_threads: params.ftp_threads,
    ftp_suffix: "_genomic.fna.gz",
    local_suffix: ".fasta.gz",
)

// Parse the manifest and sanitize the fields
process sanitize_manifest {
    container "${params.container__pandas}"
    label 'io_limited'

    input:
        path "raw.manifest.csv"
    
    output:
        path "manifest.csv", emit: manifest
    
"""
#!/usr/bin/env python3

import pandas as pd
import re

df = pd.read_csv("raw.manifest.csv")

print("Subsetting to three columns")
df = df.reindex(
    columns = [
        "GenBank FTP",
        "#Organism Name",
        "uri"
    ]
)

# Remove rows where the "GenBank FTP" doesn't start with "ftp://"
input_count = df.shape[0]
df = df.loc[
    (df["GenBank FTP"].fillna(
        ""
    ).apply(
        lambda n: str(n).startswith("ftp://")
    )) | (
        df["uri"].fillna("").apply(len) > 0
    )
]
print("%d / %d rows have valid FTP or file paths" % (input_count, df.shape[0]))

# Force organism names to be alphanumeric
df = df.apply(
    lambda c: c.apply(lambda n: re.sub('[^0-9a-zA-Z .]+', '_', n)) if c.name == "#Organism Name" else c
)

df.to_csv("manifest.csv", index=None, sep=",")
"""
}

// Parse the local manifest and sanitize the fields
process sanitize_manifest_local {
    container "${params.container__pandas}"
    label 'io_limited'

    input:
        path "raw.manifest.csv"
    
    output:
        path "manifest_local.csv", emit: manifest
    
"""
#!/usr/bin/env python3

import pandas as pd
import re

df = pd.read_csv("raw.manifest.csv")

print("Subsetting to two columns")
df = df.reindex(
    columns = [
        "#Organism Name",
        "uri"
    ]
)

# Force organism names to be alphanumeric
df = df.apply(
    lambda c: c.apply(lambda n: re.sub('[^0-9a-zA-Z .]+', '_', str(n))) if c.name == "#Organism Name" else c
)

df.to_csv("manifest_local.csv", index=None, sep=",")
"""
}

workflow parse_genomes {

    // Parse the manifest to get a name and FTP prefix for each genome
    sanitize_manifest(
        file("${params.genomes}")
    )

    sanitize_manifest.out.map {
        r -> r.splitCsv(
            header: true
        )
    }.flatten(
    ).branch {
        remote: it["uri"] == ""
        local: true
    }.set {
        split_genome_ch
    }

    // Download the genomes from the NCBI FTP server
    fetchFTP(
        split_genome_ch.remote.map {
            r -> [
                r["GenBank FTP"].split("/")[-1].replaceAll(/"/, ""), // Unique genome ID
                r["#Organism Name"],                                 // Readable name
                r["GenBank FTP"]                                     // FTP folder containing the genome
            ]
        }
    )

    // Add both sets of files to the combined channel
    fetchFTP
        .out
        .mix(
            split_genome_ch.local.map {
                r -> [
                    r["uri"].split("/")[-1],
                    r["#Organism Name"],
                    file(r["uri"])
                ]
            }
        )
        .set{genomes_ch}

    // Parse the manifest to get a name and FTP prefix for each genome
    sanitize_manifest_local(
        file("${params.genomes_local}")
    )

    sanitize_manifest_local.out.map {
        r -> r.splitCsv(
            header: true
        )
    }.flatten(
    ).map {
        r -> [
                r["uri"].split("/")[-1],
                r["#Organism Name"],
                file(r["uri"])
            ]
    }
    .set{genomes_local_ch}

    genomes_ch
        .mix(genomes_local_ch)
        .set{joined_fasta_ch}

    emit:
    joined_fasta_ch
}