Pipeline
========
Pre-requisite softwares
-------------------

To run the pipeline successfully one needs to install/download several third-party tools.
Before going into the installation procedure. We should go over the actual input format. 
The input folder contains the raw sequencing reads. 

The operation starts from the `BASECALLS_DIR`, that contains the actual intensities and the 
`BCL` files. For example the basic structure of such file will be 

::
    ➜  slide-seq-pipeline git:(main) ls /home/hsarkar/code/slide-seq-data/P25255/220329_A00689_0513_BHYK23DRXY/
        Config            InterOp                 Recipe           RunInfo.xml        SequenceComplete.txt
        CopyComplete.txt  Logs                    RTA3.cfg         runParameters.xml  Thumbnail_Images
        Data              P25255_SampleSheet.txt  RTAComplete.txt  RunParameters.xml

