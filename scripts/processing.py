from xml.etree import ElementTree as et
from dataclasses import dataclass 
from pathlib import Path
import shutil
import gzip

import sys
sys.path.append('/home/hsarkar/code/slideseq-tools/src/')

from slideseq.metadata import Manifest
from typing import Dict, List, Tuple
import logging
log = logging.getLogger(__name__)
import csv
from slideseq.metadata import Manifest, split_sample_lanes, validate_run_df
from slideseq.config import get_config
import os

import slideseq.alignment_quality as alignment_quality
import slideseq.bead_matching as bead_matching
from slideseq.config import Config, get_config
from slideseq.library import Base, Library
from slideseq.metadata import Manifest
from slideseq.pipeline.write_matrix import write_sparse_matrix
from slideseq.plot.plot_library_metrics import make_library_plots
from slideseq.retag_bam import write_retagged_bam
from slideseq.util import give_group_access, rsync_to_google, run_command
from slideseq.util.logger import create_logger


def calc_alignment_metrics(
    config: Config,
    input_base: Base,
    library: Library,
    tmp_dir: Path,
    matched_barcodes: Path = None,
    cell_tag: str = "XC",
):
    # Bam tag histogram (cells)
    cmd = config.dropseq_cmd(
        "BamTagHistogram", input_base.bam, input_base.reads_per_cell(cell_tag), tmp_dir
    )
    cmd.extend(
        [
            f"TAG={cell_tag}",
            "FILTER_PCR_DUPLICATES=false",
            f"READ_MQ={library.base_quality}",
        ]
    )

    run_command(cmd, "BamTagHistogram (cells)", library)

    # Bam tag histogram (UMIs)
    cmd = config.dropseq_cmd(
        "BamTagHistogram", input_base.bam, input_base.reads_per_umi, tmp_dir
    )
    cmd.extend(
        [
            "TAG=XM",
            "FILTER_PCR_DUPLICATES=false",
            f"READ_MQ={library.base_quality}",
        ]
    )

    run_command(cmd, "BamTagHistogram (UMIs)", library)

    # Collect RnaSeq metrics
    cmd = config.picard_cmd("CollectRnaSeqMetrics", tmp_dir)
    cmd.extend(
        [
            "--INPUT",
            input_base.bam,
            "--REF_FLAT",
            library.reference.ref_flat,
            "--OUTPUT",
            input_base.frac_intronic_exonic,
            "--STRAND_SPECIFICITY",
            "NONE",
            "--RIBOSOMAL_INTERVALS",
            library.reference.ribosomal_intervals,
        ]
    )

    run_command(cmd, "CollectRnaSeqMetrics", library)

    # Base distribution at read position for raw cellular barcode
    cmd = config.dropseq_cmd(
        "BaseDistributionAtReadPosition",
        input_base.bam,
        input_base.xc_distribution,
        tmp_dir,
    )
    cmd.extend(["TAG=XC"])
    run_command(cmd, "BaseDistributionAtReadPosition (Cellular)", library)

    # Base distribution at read position for molecular barcode
    cmd = config.dropseq_cmd(
        "BaseDistributionAtReadPosition",
        input_base.bam,
        input_base.xm_distribution,
        tmp_dir,
    )
    cmd.extend(["TAG=XM"])
    run_command(cmd, "BaseDistributionAtReadPosition (Molecular)", library)

    # Gather read quality metrics
    cmd = config.dropseq_cmd(
        "GatherReadQualityMetrics",
        input_base.bam,
        input_base.read_quality_metrics,
        tmp_dir,
    )
    run_command(cmd, "GatherReadQualityMetrics", library)

    if matched_barcodes is not None:
        cmd = config.dropseq_cmd(
            "SingleCellRnaSeqMetricsCollector",
            input_base.bam,
            input_base.frac_intronic_exonic_per_cell,
            tmp_dir,
            compression=6,
        )
        cmd.extend(
            [
                f"CELL_BARCODE_TAG={cell_tag}",
                f"READ_MQ={library.base_quality}",
                f"CELL_BC_FILE={matched_barcodes}",
                f"RIBOSOMAL_INTERVALS={library.reference.ribosomal_intervals}",
                f"ANNOTATIONS_FILE={library.reference.annotations}",
            ]
        )
        run_command(cmd, "SingleCellRnaSeqMetricsCollector", library)

    # Select cells by num transcripts
    cmd = config.dropseq_cmd(
        "SelectCellsByNumTranscripts",
        input_base.bam,
        input_base.selected_cells,
        tmp_dir,
        compression=6,
    )
    cmd.extend(
        [
            f"CELL_BARCODE_TAG={cell_tag}",
            f"MIN_TRANSCRIPTS_PER_CELL={library.min_transcripts_per_cell}",
            f"READ_MQ={library.base_quality}",
        ]
    )
    if library.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif library.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])
    run_command(cmd, "SelectCellsByNumTranscripts", library)

    # Generate digital expression files for all Illumina barcodes
    cmd = config.dropseq_cmd(
        "DigitalExpression",
        input_base.bam,
        input_base.digital_expression,
        tmp_dir,
        compression=6,
    )
    cmd.extend(
        [
            f"CELL_BARCODE_TAG={cell_tag}",
            f"CELL_BC_FILE={input_base.selected_cells}",
            f"SUMMARY={input_base.digital_expression_summary}",
            f"READ_MQ={library.base_quality}",
            "OUTPUT_HEADER=false",
            f"UEI={library}",
        ]
    )
    if library.locus_function_list == "intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=null", "LOCUS_FUNCTION_LIST=INTRONIC"])
    elif library.locus_function_list == "exonic+intronic":
        cmd.extend(["LOCUS_FUNCTION_LIST=INTRONIC"])
    run_command(cmd, "DigitalExpression", library)


def get_flowcell(run_info: et.ElementTree) -> str:
    """
    Get the flowcell name from RunInfo.xml. This should be more reliable than
    trying to parse the run directory, which might have been renamed.

    :param run_info: ElementTree representing RunInfo.xml
    :return: The flowcell name for this run
    """
    flowcell = run_info.find("./Run/Flowcell").text
    return flowcell


def get_read_structure(run_info: et.ElementTree) -> str:
    """
    Get read structure from RunInfo.xml. Assumes one index sequence only,
    will warn if two are present and ignore the second.

    :param run_info: ElementTree representing RunInfo.xml
    :return: Formatting string representing the read structure
    """
    read_elems = run_info.findall("./Run/Reads/Read[@NumCycles][@Number]")
    read_elems.sort(key=lambda el: int(el.get("Number")))

    if len(read_elems) == 4:
        # two index reads. We will just ignore the second index
        log.warning(
            "This sequencing run has two index reads, we are ignoring the second one"
        )
        return "{}T{}B{}S{}T".format(*(el.get("NumCycles") for el in read_elems))
    elif len(read_elems) != 3:
        raise ValueError(f"Expected three reads, got {len(read_elems)}")

    return "{}T{}B{}T".format(*(el.get("NumCycles") for el in read_elems))

def get_lanes(run_info: et.ElementTree) -> range:
    """
    Return the lanes for this run, as a range
    :param run_info: ElementTree representing RunInfo.xml
    :return: range object for the lanes of the run
    """
    lane_count = int(run_info.find("./Run/FlowcellLayout[@LaneCount]").get("LaneCount"))
    return range(1, lane_count + 1)


@dataclass
class RunInfo:
    """A dataclass to represent RunInfo.xml for a sequencing run (a.k.a. flowcell)"""

    run_dir: Path
    flowcell: str
    lanes: range
    read_structure: str

    @property
    def demux_log(self):
        return f"demultiplex.{self.flowcell}.L00$TASK_ID.log"

    @property
    def basecall_dir(self):
        return self.run_dir / "Data" / "Intensities" / "BaseCalls"

    def alignment_log(self, lane: int):
        return f"alignment.{self.flowcell}.L{lane:03d}.$TASK_ID.log"


def get_run_info(run_dir: Path) -> RunInfo:
    # open the RunInfo.xml file and parse it with element tree
    with (run_dir / "RunInfo.xml").open() as f:
        run_info = et.parse(f)

    return RunInfo(
        run_dir,
        get_flowcell(run_info),
        get_lanes(run_info),
        get_read_structure(run_info),
    )


run_info = get_run_info(
    Path('/home/hsarkar/code/slide-seq-data/P25255/220329_A00689_0513_BHYK23DRXY')
)

from openpyxl import load_workbook
import pandas as pd
import slideseq.util.constants as constants
ws = load_workbook('/home/hsarkar/notebooks/slide_seq_analysis/2023/example_metadata.xlsx').active
from itertools import islice
data = ws.values
cols = next(data)
data = list(data)
#idx = [r[0] for r in data]
#data = (islice(r, 1, None) for r in data)
df = pd.DataFrame(data, columns=cols)
run_df = df[constants.METADATA_COLS]
run_df.columns = [c.lower() for c in constants.METADATA_COLS]
import numpy as np
run_df = run_df.fillna(value=0)
run_df = run_df.astype(constants.METADATA_TYPES)


run_name = '1'
output_dir = Path('/home/hsarkar/code/slide-seq-data/P25255/processed_data')/\
    Path(run_name)
flowcell_dirs = sorted(Path(fd) for fd in set(run_df.bclpath))
manifest_file = output_dir / "manifest.yaml"
metadata_file = output_dir / "metadata.csv"

run_info_list = [ run_info ]
config = get_config()

manifest = Manifest(
            run_name=run_name,
            flowcell_dirs=flowcell_dirs,
            workflow_dir=output_dir,
            library_dir=Path('/home/hsarkar/code/slide-seq-data/P25255/processed_data/libraries'),
            metadata_file=metadata_file,
            metadata=split_sample_lanes(run_df, run_info_list),
            email_addresses=sorted(
                set(e.strip() for v in run_df.email for e in v.split(","))
            ),
        )

md = split_sample_lanes(run_df, run_info_list)
n_libraries = len(list(manifest.libraries))

from slideseq.util import give_group_access, rsync_to_google, run_command, start_popen

import os
force = False

for library_index in range(n_libraries):

    library = manifest.get_library(library_index)
    if not(library.name == 'P25255_1004'):
        print('skipping ', library.name)
        continue
    barcode_matching_file = (
        library.barcode_matching_dir / f"{library}_barcode_matching.txt.gz"
    )
    barcode_coordinate_file = (
        library.barcode_matching_dir / f"{library}_barcode_xy.txt.gz"
    )
    matched_barcodes_file = (
        library.barcode_matching_dir / f"{library}_matched_barcodes.txt.gz"
    )

    # merge lanes
    cmd = config.picard_cmd("MergeSamFiles", manifest.tmp_dir)
    cmd.extend(
        [
            "--CREATE_INDEX",
            "true",
            "--CREATE_MD5_FILE",
            "false",
            "--OUTPUT",
            library.merged.bam,
            "--SORT_ORDER",
            "coordinate",
            "--ASSUME_SORTED",
            "true",
        ]
    )

    for bam_file in library.processed_bams:
        cmd.extend(["--INPUT", bam_file])

    run_command(cmd, "MergeBamFiles", library)
    # generate various metrics files, including digital expression matrix
    print('Create matrices')
    calc_alignment_metrics(config, library.merged, library, manifest.tmp_dir)

    library.barcode_matching_dir.mkdir(exist_ok=True, parents=True)
    shutil.copy(library.bead_barcodes, library.barcode_matching_dir)
    shutil.copy(library.bead_locations, library.barcode_matching_dir)

    barcode_list, barcode_mapping, bead_xy, _ = bead_matching.match_barcodes(
        library.merged.selected_cells, library.bead_barcodes, library.bead_locations
    )

    bead_matching.write_barcode_mapping(
        barcode_mapping, bead_xy, barcode_matching_file
    )
    bead_matching.write_barcode_xy(barcode_list, bead_xy, barcode_coordinate_file)

    with gzip.open(matched_barcodes_file, "wt") as out:
        for bead_bc in sorted(set(barcode_mapping.values())):
            print(bead_bc, file=out)

    # subset to the matched beads and add combined barcode as XB tag
    write_retagged_bam(library.merged.bam, library.matched.bam, barcode_mapping)

    # do it all again, but should be faster on the smaller file
    calc_alignment_metrics(
        config,
        library.matched,
        library,
        manifest.tmp_dir,
        matched_barcodes_file,
        "XB",
    )

    write_sparse_matrix(library)

    make_library_plots(library, bead_xy)
    log.info(f"Processing for {library} complete")