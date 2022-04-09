
import pandas as pd
from tqdm.notebook import tqdm
import yaml
import re
import os
from pathlib import Path
import pandas as pd
from genomic_benchmarks.utils.datasets import (
    _download_url,
    _fastagz2dict,
    _get_reference_name,
    _rev,
)
from genomic_benchmarks.utils.paths import REF_CACHE_PATH


def download_dataset(
    dest_file='rbp_encode_eclip.csv',
    cache_path=REF_CACHE_PATH,
    force_download=False
):
    """
    Transform an interval-list genomic dataset into a full-seq genomic dataset.
            Parameters:
                    cache_path (str or Path): Folder to store the downloaded references.
                    force_download (bool): If True, force downloading of references.
    """

    interval_list_dataset = Path('./csv') 
    with open(interval_list_dataset / 'metadata.yaml', "r") as fr:
        metadata = yaml.safe_load(fr)

    refs = _download_references(metadata, cache_path=cache_path, force=force_download)
    fastas = _load_fastas_into_memory(refs, cache_path=cache_path)

    mode = 'w'
    for c in tqdm(metadata["classes"]):
        dt_filename = interval_list_dataset / (c + ".csv.gz")
        dt = pd.read_csv(str(dt_filename), compression="gzip")

        ref_name = _get_reference_name(metadata["classes"][c]["url"])
        dt["seq"] = _fill_seq_column(fastas[ref_name], dt)

        dt.to_csv(dest_file, mode=mode, index=False)
        mode = 'a'


EXTRA_PREPROCESSING = {
    # known extra preprocessing steps
    "default": [None, None, lambda x: x],
    "ENSEMBL_HUMAN_GENOME": [24, "MT", lambda x: "chr" + x],  # use only chromosomes, not contigs, and add chr prefix
    "ENSEMBL_MOUSE_GENOME": [21, "MT", lambda x: "chr" + x],  # use only chromosomes, not contigs, and add chr prefix
    "ENSEMBL_HUMAN_TRANSCRIPTOME": [
        190_000,
        None,
        lambda x: re.sub("ENST([0-9]*)[.][0-9]*", "ENST\\1", x),
    ],  # remove the version number from the ensembl id
}


def _load_fastas_into_memory(refs, cache_path):
    # load all references into memory
    fastas = {}
    for ref in refs:
        ref_path = Path(cache_path) / _get_reference_name(ref[0])
        ref_type = ref[1]
        ref_extra_preprocessing = ref[2] if ref[2] is not None else "default"
        if ref_extra_preprocessing not in EXTRA_PREPROCESSING:
            raise ValueError(f"Unknown extra preprocessing: {ref_extra_preprocessing}")

        if ref_type == "fa.gz":
            fasta = _fastagz2dict(
                ref_path,
                fasta_total=EXTRA_PREPROCESSING[ref_extra_preprocessing][0],
                stop_id=EXTRA_PREPROCESSING[ref_extra_preprocessing][1],
                region_name_transform=EXTRA_PREPROCESSING[ref_extra_preprocessing][2],
            )
            fastas[_get_reference_name(ref[0])] = fasta
        else:
            raise ValueError(f"Unknown reference type {ref_type}")
    return fastas


def _download_references(metadata, cache_path, force=False):
    # download all references from the metadata into cache_path folder
    cache_path = Path(cache_path)
    if not cache_path.exists():
        cache_path.mkdir(parents=True)

    refs = {(c["url"], c["type"], c.get("extra_processing")) for c in metadata["classes"].values()}

    for ref in refs:
        ref_path = cache_path / _get_reference_name(ref[0])
        if not ref_path.exists() or force:
            _download_url(ref[0], ref_path)
        else:
            print(f"Reference {ref_path} already exists. Skipping.")

    return refs


def _fill_seq_column(fasta, df):
    # fill seq column in DataFrame tab
    if not all([r in fasta for r in df["chr"]]):
        missing_regions = list({r for r in df["chr"] if r not in fasta})
        if len(missing_regions) > 5:
            missing_regions = missing_regions[:6]
        raise ValueError("Some chromosomes not found in the reference, e.g. " + " ".join([str(r) for r in missing_regions]))
    output = pd.Series(
        [
            _rev(fasta[region][start:end], strand)
            for region, start, end, strand in zip(df["chr"], df["start"], df["end"], df["strand"])
        ]
    )
    return output

if __name__ == "__main__":

    try:
        os.remove('rbp_encode_eclip.csv')
    except OSError:
        pass

    download_dataset(
        dest_file='rbp_encode_eclip.csv',
        cache_path=REF_CACHE_PATH,
        force_download=False
    )