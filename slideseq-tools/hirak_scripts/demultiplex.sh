TMP_DIR="/home/hsarkar/code/slide-seq-data/P25255/processed_data/tmp"
PICARD_JAR="/home/hsarkar/code/slideseq-tools/soft/picard/build/libs/picard.jar"
SGE_TASK_ID="2"
BASECALLS_DIR="/home/hsarkar/code/slide-seq-data/P25255/220329_A00689_0513_BHYK23DRXY/Data/Intensities/BaseCalls"
READ_STRUCTURE="42T8B41T"
OUTPUT_DIR="/home/hsarkar/code/slide-seq-data/P25255/processed_data/1"
FLOWCELL="HYK23DRXY"

java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms8g -Xmx124g \
  -jar ${PICARD_JAR} CheckIlluminaDirectory \
  --TMP_DIR ${TMP_DIR} \
  --LANES ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE}

# extract barcodes
java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms64g -Xmx124g \
  -jar ${PICARD_JAR} ExtractIlluminaBarcodes \
  --TMP_DIR ${TMP_DIR} \
  --LANE ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE} \
  --OUTPUT_DIR ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/barcodes \
  --BARCODE_FILE ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/barcode_params.txt \
  --METRICS_FILE ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/${FLOWCELL}.barcode_metrics.txt \
  --COMPRESS_OUTPUTS true \
  --NUM_PROCESSORS 8

# create uBAM files
java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms64g -Xmx124g \
  -jar ${PICARD_JAR} IlluminaBasecallsToSam \
  --TMP_DIR ${TMP_DIR} \
  --LANE ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE} \
  --RUN_BARCODE ${FLOWCELL} \
  --BARCODES_DIR ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/barcodes \
  --LIBRARY_PARAMS ${OUTPUT_DIR}/${FLOWCELL}/L00${SGE_TASK_ID}/library_params.txt \
  --INCLUDE_NON_PF_READS false \
  --APPLY_EAMSS_FILTER false \
  --ADAPTERS_TO_CHECK null \
  --IGNORE_UNEXPECTED_BARCODES true \
  --SEQUENCING_CENTER BI

# give group access to all files and directories
chmod -R --silent g+rwX ${OUTPUT_DIR}