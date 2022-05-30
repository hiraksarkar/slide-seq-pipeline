TMP_DIR="/home/hsarkar/code/slide-seq-data/P25255/processed_data/tmp"
PICARD_JAR="/home/hsarkar/code/slideseq-tools/soft/picard/build/libs/picard.jar"
SGE_TASK_ID="1"
BASECALLS_DIR="/home/hsarkar/code/slide-seq-data/P25255/220329_A00689_0513_BHYK23DRXY/Data/Intensities/BaseCalls"
READ_STRUCTURE="42T8B41T"
OUTPUT_DIR="/home/hsarkar/code/slide-seq-data/P25255/processed_data/1"
FLOWCELL="HYK23DRXY"

echo "java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms8g -Xmx124g \
  -jar ${PICARD_JAR} CheckIlluminaDirectory \
  --TMP_DIR ${TMP_DIR} \
  --LANES ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE}"

java -Djava.io.tmpdir=${TMP_DIR} -XX:+UseParallelGC \
  -XX:GCTimeLimit=20 -XX:GCHeapFreeLimit=10 -Xms8g -Xmx124g \
  -jar ${PICARD_JAR} CheckIlluminaDirectory \
  --TMP_DIR ${TMP_DIR} \
  --LANES ${SGE_TASK_ID} \
  --BASECALLS_DIR ${BASECALLS_DIR} \
  --READ_STRUCTURE ${READ_STRUCTURE}

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