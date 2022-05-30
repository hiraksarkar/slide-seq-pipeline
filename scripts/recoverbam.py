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
import edlib

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

import random
import editdistance
def get_hamming_ball_seq(seq):
    return(
        set([seq[:i]+x+seq[i+1:] for x in ['A','T','G','C'] for i in range(len(seq))])
    )

def replace_N(seq):
    return([seq.replace('N',x) for x in ['A','T','G','C']])

def flatten(t):
    return [item for sublist in t for item in sublist]

import random
import editdistance
def get_hamming_ball_seq(seq):
    return(
        set([seq[:i]+x+seq[i+1:] for x in ['A','T','G','C'] for i in range(len(seq))])
    )

def replace_N(seq):
    return([seq.replace('N',x) for x in ['A','T','G','C']])

def flatten(t):
    return [item for sublist in t for item in sublist]

def get_recovery2(
    sample,
    whitelist_size = 70000,
    ref_end = 33,
    edit_distance_cutoff=2,
    total_sample = 10000,
    strict_matching=False
):

    file_name = sample.raw_barcode

    print('Reading ', file_name)
    print('Making whitelist')

    # Make whitelist
    
    whitelist_freq = defaultdict(int)
    barcode_list = get_lines(file_name)
    bcLocation = extract_beadmatrix(sample.data.bead_structure)
    bi1,be1 = bcLocation.bi1, bcLocation.be1
    bi2,be2 = bcLocation.bi2, bcLocation.be2
    si,se = bcLocation.si, bcLocation.se
    ui = bcLocation.ui
    k = be2 - bi2
    
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
    whitelist = dict(list(whitelist_sort.items())[ : whitelist_size ])
    whitelist_id = 0
    for kmer in tqdm.tqdm(whitelist):
        whitelist[ kmer ] = whitelist_id
        whitelist_id += 1
    whitelist_list = [kmer for kmer in whitelist]

    barcode_part_1_set = defaultdict(set)
    for kmer in tqdm.tqdm(whitelist): 
        barcode_part_1_set[kmer[:be1]].add(kmer[be1:])

    # store the associated second part of barcode
    barcode_part_1_list = {}
    for key in tqdm.tqdm(barcode_part_1_set): 
        barcode_part_1_list[key] = list(barcode_part_1_set[key])

    exact_match_lines = [i for i in range(len(barcode_list)) if (
        barcode_list[i][bi1:be1] + barcode_list[i][bi2:be2] in whitelist)
    ]
    to_recover = list(set(range(len(barcode_list))) - set(exact_match_lines))

    alignment_results = [
        edlib.align(
                freq_spacer,
                barcode_list[i],
                task ='locations',
                mode = "HW",
                k = 4
            ) for i in tqdm.tqdm(to_recover)
    ]
    alignment_results_dict = {}
    for i,a in tqdm.tqdm(enumerate(alignment_results)):
        alignment_results_dict[ to_recover[i] ] = a
    #print(len(to_recover))
    #to_recover_sample = random.sample(to_recover, total_sample)
    print(len(set(to_recover)))
    recovered = {}
    key_set = []
    recovered_map = {}
    exact_match = 0
    part1_em = 0
    part1_ham = 0
    part1_shift = 0
    part2_em = 0
    part2_alignment = 0
    match_found = 0
    bad_start = 0
    rematch = {}
    status = {}
    align_status = {}
    ref_lens = []
    num_kmers = []

    print('start recovering...')
    for i in tqdm.tqdm(to_recover):
        a = alignment_results_dict[i]
        if (a['editDistance'] == -1):
            continue
        start_pos = a['locations'][0][0]
        end_pos = a['locations'][0][1]
        part1 = barcode_list[i][:start_pos]
        orig_part1 = part1
        ref = barcode_list[i][end_pos:ref_end]
        #print(len(ref))
        # chuck out if start position is too early
        if (start_pos < 7 or start_pos > 8):
            bad_start += 1
            continue
        if(len(ref) < k):
            continue
        ref_lens.append(len(ref))
        if len(part1) == 7:
            part1 += 'N'
        # first do hamming-sensitive search
        exit_part1 = True
        exit_part2 = False
        sad_part_2 = False
        if not(part1 in barcode_part_1_set):
            edit_distance_cutoff_part2 = edit_distance_cutoff
            exit_part1 = False
            if 'N' in part1:
                part1_tmp = flatten([get_hamming_ball_seq(s) for s in replace_N(part1)])
            else:
                part1_tmp = get_hamming_ball_seq(part1)
            for kmer in part1_tmp:
                if kmer in barcode_part_1_set:
                    part1 = kmer
                    exit_part1 = True
                    # start second part
                    # search part2 
                    # kmer search
                    num_kmers += [len(ref) - k + 1]
                    for j in range(len(ref) - k + 1):
                        if(ref[j:j+k] in barcode_part_1_set[part1]):
                            recovered[ i ] = whitelist[ part1 + ref[j:j+k] ]
                            # recovered_map[i] = whitelist[ part1 + ref[j:j+k] ]
                            # rematch[barcode_list[i][:start_pos]] = part1
                            # key_set.append(i)
                            exit_part2 = True
                            part1_ham += 1
                            part2_em += 1
                            match_found += 1
                            # status[i] = ['ham_and_kmer']
                            #print(len(key_set), match_found)
                            break
                    # align
                    if exit_part2 == False:
                        if strict_matching:
                            ed = editdistance.eval(part1, orig_part1)
                            if ed <= edit_distance_cutoff:
                                edit_distance_cutoff_part2 = edit_distance_cutoff - ed
                            else:
                                sad_part_2 = True
                                continue
                        tmp_a_results = [
                            edlib.align(
                                p,
                                ref,
                                task = 'locations',
                                mode = "HW",
                                k = edit_distance_cutoff_part2
                            ) for p in barcode_part_1_list[part1]
                        ]
                        for j, sub_a in enumerate(tmp_a_results):
                            if 0 <= sub_a['editDistance'] <= edit_distance_cutoff_part2:
                                recovered[ i ] = \
                                    whitelist[ part1 + barcode_part_1_list[part1][j] ]
                                # recovered_map[i] = whitelist[ part1 + barcode_part_1_list[part1][j] ]
                                match_found += 1
                                # key_set.append(i) 
                                part1_ham += 1
                                part2_alignment += 1
                                
                                # status[i] = ['ham_and_align']
        
                                # align_status[i] = (barcode_part_1_list[part1][j],
                                #                    ref,
                                #                    sub_a
                                #                   )
                                exit_part2 = True
                                #print(len(key_set), match_found)
                                break
                if (exit_part1 & exit_part2):
                    break
        else:
            part1_em += 1
            # search part2 
            # kmer search
            num_kmers += [len(ref) - k + 1]
            for j in range(len(ref) - k + 1):
                if(ref[j:j+k] in barcode_part_1_set[part1]):
                    recovered[ i ] = whitelist[ part1 + ref[j:j+k] ]
                    # recovered_map[i] = whitelist[ part1 + ref[j:j+k] ]
                    exit_part2 = True
                    part2_em += 1
                    exact_match += 1
                    match_found += 1
                    # key_set.append(i)
                    # status[i] = ['match_and_kmer']
                    #print(len(key_set), match_found)
                    break
            # align
            if exit_part2 == False:
                edit_distance_cutoff_part2 = edit_distance_cutoff - 1
                tmp_a_results = [
                    edlib.align(
                        p,
                        ref,
                        task = 'locations',
                        mode = "HW",
                        k = edit_distance_cutoff_part2
                    ) for p in barcode_part_1_list[part1]
                ]
                
                for j, sub_a in enumerate(tmp_a_results):
                    if 0 <= sub_a['editDistance'] <= edit_distance_cutoff_part2:
                        recovered[ i ] = \
                            whitelist[ part1 + barcode_part_1_list[part1][j] ]
                        # recovered_map[i] = whitelist[ part1 + barcode_part_1_list[part1][j] ]
                        # align_status[i] = (barcode_part_1_list[part1][j],
                        #                            ref,
                        #                            sub_a
                        #                           )
                        part2_alignment += 1
                        match_found += 1
                        exit_part2 = True
                        # key_set.append(i) 
                        # status[i] = ['match_and_align']
                        #print(len(key_set), match_found)
                        break
                        
    print('total recovered {:.2f}%'.format(len(recovered)/ len(to_recover) *100))
                        
    return(
        {
            'recovered': recovered,
            'whitelist': whitelist_list
        }
    )

import pysam
import copy

for library_index in range(n_libraries):
    for run_info in run_info_list:
        for lane in run_info.lanes:

            sample = manifest.get_sample(
                library_index, run_info.flowcell, lane)

            if not(sample.name == 'P25255_1004'):
                print('skipping ', sample.name, sample.lane)
                continue

            print(sample.name, sample.lane)
            ret = get_recovery2(
                sample
            )

            recovered = ret['recovered']
            whitelist = ret['whitelist']

            print(list(recovered.keys())[:10])

            bam_file = sample.raw_ubam
            new_bam_file = sample.recovered_ubam

            print('Writing the new bam')
            line_index = 0
            save = pysam.set_verbosity(0)
            with pysam.AlignmentFile(bam_file, mode="rb", threads=8, check_sq=False) as r_bam:
                with pysam.AlignmentFile(
                    new_bam_file, mode="wb", template= r_bam, threads=8
                ) as n_bam:
                    for a in r_bam:
                        # even lines are barcode
                        r_index = line_index // 2
                        if ( (line_index % 2 != 0) or not(r_index in recovered)):
                            n_bam.write(a)
                        else:
                            # write modified version
                            cor = whitelist[ recovered[ r_index ] ]
                            q = copy.copy(a)
                            q.query_sequence = cor[bi1:be1] + a.query_sequence[si:se] + cor[be1:] + a.query_sequence[be2:]
                            q.query_qualities = a.query_qualities
                            assert( len(q.query_sequence) == len(a.query_sequence) )
                            n_bam.write(q)
                        line_index += 1
                        