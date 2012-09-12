#
# Usage:
#
# generate_data_files.sh [keyword for time extraction]
#
#
for i in 001 002 004 008 016 032 064 100 128 256 512; do
    echo "${i}	$(./extract_min_time.sh out_${i}.txt $1)" >> $1.dat;
done
