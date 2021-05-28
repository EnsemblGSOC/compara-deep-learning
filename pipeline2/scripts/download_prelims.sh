#!/bin/sh


# set up the required directories
mkdir -p downloads/gtfs/
mkdir -p downloads/pep/
mkdir -p downloads/cds/
mkdir -p logs


echo "Downloading gtfs"
wget -i download_links/gtf_link.txt --show-progress -o logs/gtf_download_log.txt -P downloads/gtfs/
echo "Download complete"

echo "Downloading protein files"
wget -i download_links/protein_seq.txt --show-progress -o logs/pep_download_log.txt -P downloads/pep/
echo "Download complete"

echo "Downloading cds"
wget -i download_links/seq_links.txt --show-progress -o logs/seq_download_log.txt -P downloads/cds/
echo "Download complete"

echo "ALL downloads completed!"