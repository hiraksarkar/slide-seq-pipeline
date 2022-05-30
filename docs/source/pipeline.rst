Pipeline
========
Pre-requisite softwares
-------------------

To run the pipeline successfully one needs to install/download several third-party tools.

**Tools**

- `Picard tools <https://broadinstitute.github.io/picard/>`_
- `Dropseq tools <https://github.com/broadinstitute/Drop-seq/tree/v2.5.1>`_
- `STAR <https://github.com/alexdobin/STAR>`_


Before going into the installation procedure. We should go over the actual input format. 
The input folder contains the raw sequencing reads. 

The operation starts from the ``BASECALLS_DIR``, that contains the actual intensities and the 
``bcl`` files. For example the basic structure of such file will be 

::

        ➜slide-seq-pipeline git:(main) ls /home/hsarkar/code/slide-seq-data/P25255/220329_A00689_0513_BHYK23DRXY/
        Config            InterOp                 Recipe           RunInfo.xml        SequenceComplete.txt
        CopyComplete.txt  Logs                    RTA3.cfg         runParameters.xml  Thumbnail_Images
        Data              P25255_SampleSheet.txt  RTAComplete.txt  RunParameters.xml


The ``bcl`` files are present in ``Data/Intensities/BaseCalls/``.  The ``RunInfo.xml`` used to be the file
that has the information that is to be parsed before writing down the information.

Use the existing `package <https://github.com/MacoskoLab/slideseq-tools>`_ so that we can import existing structures 
from the ``slideseq-tools``. 

Install existing ``slideseq-tools`` package 
::

        >git clone https://github.com/MacoskoLab/slideseq-tools.git

Create a spread-sheet that contains some of the following fields that we need to feel. An example spread-sheet
is provided in `spreadsheet` folder of the repo. 

::

        Index(['library', 'date', 'flowcell', 'run_name', 'bclpath', 'lane',
        'sample_barcode', 'bead_structure', 'reference', 'run_barcodematching',
        'locus_function_list', 'start_sequence', 'base_quality',
        'min_transcripts_per_cell', 'email', 'puckcaller_path', 'bead_type',
        'gen_read1_plot', 'gen_downsampling'],
        dtype='object')


**Desctiption**


library: Name of the library (Each puck will have a specic name)
puckcaller_path: A location containing the files from the image analysis pipeline, ``BeasBarcodes.txt`` and ``BeadLocations.txt``
run_name: Depends on how many different runs we are having
BCLPath: Actual location to the BCL files
sample_barcod: It's a 8bp barcode that can be obtained from  ``RunInfo.xml`` file with the column name ``index`` (I am not sure about this yet)
reference: The actual reference file I used ``GRCh38.fasta``
bead_structure: Determines the actual length of the cell barcode and the corresponding UMI

**Cosntruction**


The construction of the dataframe with the help of the ``RunInfo.xml`` is depicted in the
next embedded notebook.