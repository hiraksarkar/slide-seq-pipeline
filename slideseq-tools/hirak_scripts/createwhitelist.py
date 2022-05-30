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
from subprocess import run
from collections import defaultdict
import tqdm
import re

# util functions
## Extract position
@dataclass
class BeadBarcodeLocations:
    bi1: int = 0
    be1: int = 0
    si: int = 0
    se: int = 0
    bi2: int = 0
    be2: int = 0
    ui: int = 0
    ue: int = 0


def extract_beadmatrix(bead_structure):
    bcLocations = BeadBarcodeLocations()
    matches = re.findall(r'(\d+)([A-Z]{1})', bead_structure)
    see_C = False
    see_X = False
    cum_index = 0
    for m in matches:
        if m[1] == 'C':
            if not see_C:
                bcLocations.be1 = int(m[0])
                see_C = True
            else:
                bcLocations.bi2 = cum_index
                bcLocations.be2 = cum_index + int(m[0])
        elif m[1] == 'X':
            if not see_X:
                bcLocations.si = cum_index
                bcLocations.se = cum_index + int(m[0])
                see_X = True
        elif m[1] == 'M':
            bcLocations.ui = cum_index
            bcLocations.ue = cum_index + int(m[0])
        cum_index = cum_index + int(m[0])
    return (bcLocations)


def extract_frequent_spacer(lines, bcLocation):
    spacer_seq = defaultdict(int)
    si,se = bcLocation.si, bcLocation.se
    for line in tqdm.tqdm(lines):
        spacer_seq[line[si:se]] += 1
    freq_spacer = sorted(
        spacer_seq.items(), 
        key=lambda item: item[1], 
        reverse=True
    )[0][0]
    return(spacer_seq, freq_spacer)


def get_lines(file_name):
    fp = open(file_name)
    barcodes_list = [line.strip() for line in tqdm.tqdm(fp.readlines())]
    fp.close()
    return(barcodes_list)


def get_mutations(s):
    for d in ['A','T','C','G']:
        yield s+d



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

bcLocation = extract_beadmatrix('8C18X6C1X8M1X')
bi1,be1 = bcLocation.bi1, bcLocation.be1
bi2,be2 = bcLocation.bi2, bcLocation.be2
si,se = bcLocation.si, bcLocation.se
ui = bcLocation.ui

def merge_white_list(path_list):
    white_list_per_sample = {}
    whitelist_freq = defaultdict(int)
    for file_name in path_list:
        barcode_list = get_lines(file_name)
        spacer_seq, freq_spacer = extract_frequent_spacer(barcode_list, bcLocation)
        pure_lines = [i for i,b in enumerate(barcode_list) if ((b[si:se] == freq_spacer) and not('N' in b))]
        for i in tqdm.tqdm(pure_lines):
            whitelist_freq[ barcode_list[i][bi1:be1] + barcode_list[i][bi2:be2] ] += 1

    whitelist_sort = dict(
        sorted(
            whitelist_freq.items(), 
            key=lambda item: item[1], 
            reverse=True
    ))
    print(len(whitelist_sort))
    whitelist = dict(list(whitelist_sort.items())[:70000])
    whitelist_list = [k for k in whitelist]

    return(whitelist_list)

# Create dummy locations 
with open('/home/hsarkar/code/SlideSeqRecover/outdata/BeadLocations.txt') as fp:
    x = np.array(list(map(float,fp.readline().strip().split(','))))
    y = np.array(list(map(float,fp.readline().strip().split(','))))

for library_index in range(n_libraries):
    library = manifest.get_library(library_index)
    
    whitelist = merge_white_list(library.raw_barcode)
    library.bead_barcodes.parent.mkdir(exist_ok=True)
    with open(library.bead_barcodes, 'w') as fp:
        for wl in whitelist:
            fp.write("{}\n".format(",".join(list(wl))))    

    library.bead_locations.parent.mkdir(exist_ok=True)
    with open(library.bead_locations, 'w') as fp:
        fp.write("{}\n".format(",".join(list(map(str,x[:len(whitelist)])))))
        fp.write("{}\n".format(",".join(list(map(str,y[:len(whitelist)])))))