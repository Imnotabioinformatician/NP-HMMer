#!/bin/bash

# Check if a directory has been passed as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/your/fasta/files"
    exit 1
fi

# The directory with FASTA files as the first argument
DIR=$1

# Output file for concatenated sequences
OUTPUT="concatenated_sequences.fasta"

# Clear or create the output file
> "$OUTPUT"

# Function to process each file
process_file() {
    local FILE=$1
    local NAME=$(basename "$FILE" .fasta)
    
    while IFS= read -r line || [ -n "$line" ]; do
        if [[ $line == '>'* ]]; then
            # When the line starts with '>', it's a header line
            # Append filename to the header
            echo ">$NAME | ${line:1}" >> "$OUTPUT"
        else
            # Otherwise, it's a sequence line
            echo "$line" >> "$OUTPUT"
        fi
    done < "$FILE"
}

# Loop through each FASTA file in the specified directory
for FILE in "$DIR"/*.fasta; do
    process_file "$FILE"
    
    # Optionally, add a newline after each file's sequences for better readability
    echo "" >> "$OUTPUT"
done

echo "Concatenation complete. Output available in $OUTPUT"
