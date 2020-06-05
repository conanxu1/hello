#! /bin/sh
fdupes -r -f . | grep -v ^$ | tee duplicate.txt cat duplicate.txt | \
while read file; do rm -v "$file"; done
