from xml.etree import ElementTree as et
from dataclasses import dataclass 
from pathlib import Path

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
from slideseq.alignment_quality import write_alignment_stats


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
# TagBamWithReadSequenceExtended (Cellular)
drop_seq_cmd_1 = {}
# TagBamWithReadSequenceExtended (Molecular)
drop_seq_cmd_2 = {}

# FilterBam
drop_seq_cmd_3 = {}
# TrimStartingSequence
drop_seq_cmd_4 = {}
# PolyATrimmer
drop_seq_cmd_5 = {}

# Sam2Fastq
picard_cmd_1 = {}
# star
star_cmd = {}
# SortSam
picard_cmd_2 = {}
# MergedBam
picard_cmd_3 = {}
# TagReadWithGeneFunction
drop_seq_cmd_6 = {}



for library_index in range(n_libraries):
    for run_info in run_info_list:
        for lane in run_info.lanes:
            print('[{},{},{}]'.format(
                library_index,run_info.flowcell, lane))
            sample = manifest.get_sample(
                library_index, run_info.flowcell, lane)

             # Don't run for 1001-L001
            if not(sample.name == 'P25255_1001' and sample.lane == 2):
                print('skipping ', sample.name, sample.lane)
                continue
            
            if sample is None:
                print('sample not present')
            else:
                for barcode_ubam in sample.barcode_ubams:
                    if sample.raw_ubam.exists():
                        pass
                    elif not barcode_ubam.exists():
                        raise ValueError(f"{barcode_ubam} does not exist")
            (sample.lane_dir / "alignment").mkdir(exist_ok=True)
            bead_structure = sample.get_bead_structure()
            xc_range = ":".join(f"{i}-{j}" for c, i, j in bead_structure if c == "C")
            xm_range = ":".join(f"{i}-{j}" for c, i, j in bead_structure if c == "M")

            procs = []
            
            cmd = config.picard_cmd("SortSam", manifest.tmp_dir, mem="24g")
            cmd.extend(
                [
                    "-I",
                    sample.aligned_bam,
                    "-O",
                    "/dev/stdout",
                    "--SORT_ORDER",
                    "queryname",
                ]
            )
            
            procs.append(start_popen(cmd, "SortSam", sample, lane))
            cmd = [str(arg) for arg in cmd]
            # cmd = ' '.join(cmd)
            picard_cmd_2[(library_index,run_info.flowcell,lane)] = cmd
            
            # Merge unmapped bam and aligned bam
            cmd = config.picard_cmd("MergeBamAlignment", manifest.tmp_dir, mem="24g")
            cmd.extend(
                [
                    "-R",
                    sample.reference.fasta,
                    "--UNMAPPED",
                    sample.polya_filtered_ubam,
                    "--ALIGNED",
                    "/dev/stdin",
                    "-O",
                    "/dev/stdout",
                    "--COMPRESSION_LEVEL",
                    "0",
                    "--INCLUDE_SECONDARY_ALIGNMENTS",
                    "false",
                    "--CLIP_ADAPTERS",
                    "false",
                ]
            )
            procs.append(start_popen(cmd, "MergeBamAlignment", sample, lane, procs[-1]))
            

             # Tag read with interval
            cmd = config.dropseq_cmd(
                "TagReadWithInterval", "/dev/stdin", "/dev/stdout", manifest.tmp_dir
            )
            cmd.extend([f"INTERVALS={sample.reference.intervals}", "TAG=XG"])

            procs.append(start_popen(cmd, "TagReadWithInterval", sample, lane, procs[-1]))

            # Tag read with gene function            
            cmd = config.dropseq_cmd(
                "TagReadWithGeneFunction",
                "/dev/stdin",
                sample.processed_bam,
                manifest.tmp_dir,
                compression=5,
            )
            cmd.extend(
                [f"ANNOTATIONS_FILE={sample.reference.annotations}", "CREATE_INDEX=false"]
            )

            procs.append(start_popen(cmd, "TagReadWithGeneFunction", sample, lane, procs[-1]))             
            
            # close intermediate streams
            for p in procs[:-1]:
                p.stdout.close()

            # wait for final process to finish
            procs[-1].communicate()
            log.debug("Finished with post-alignment processing")