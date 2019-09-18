bam=$1
output=$2

# make directory to save bam-readcount output, sinvict requires that it is only file in folder
mkdir /home/analysis_in

# run bam-readcount, suppress warnings to 3 (otherwise will be one per read)
/home/bam-readcount/bin/bam-readcount -w 3 $bam > /home/analysis_in/readcount.txt

# run sinvict
/home/sinvict/sinvict -t /home/analysis_in -o /home

# TODO process sinvict output
cp /home/calls_level1.sinvict $output
