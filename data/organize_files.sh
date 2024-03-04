#!/bin/bash

# Loop over all the gzipped files in the current directory
for file in *.gz; do
    # Extract the full prefix (e.g., GSM6738782_WD_4539 from GSM6738782_WD_4539_barcodes.tsv.gz)
    prefix=$(echo "$file" | grep -oP '.*(?=_)')
    
    # Skip if no prefix is found or it's a directory
    if [ -z "$prefix" ] || [ -d "$prefix" ]; then
        continue
    fi
    
    # Create a directory for the full prefix if it doesn't exist
    if [ ! -d "$prefix" ]; then
        mkdir "$prefix"
    fi
    
    # Move the file to the new directory with the new name
    mv "$file" "${prefix}/"
done

echo "Files have been reorganized and moved to respective directories."
