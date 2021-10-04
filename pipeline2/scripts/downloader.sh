#!/bin/sh


# Script takes 3 arguments to download files using wget.
# Run script with command line as:
# bash downloader.sh download_links.txt out_dir_path out_log_path


# The text file containing download links as first script argument
links_text_file=$1
# The output directory
out_dir=$2
# File path to output log of the download
log_path=$3
# Make the output directory
mkdir -p $out_dir


echo "Download started"
wget -i $links_text_file -N --show-progress -o $log_path -P $out_dir
echo "Download complete"
echo "ALL downloads completed!"


# Running the same line again should check whether the downloaded files downloaded successfully
# If there's a mismatch or missing file, it will download again
echo "Checking Downloads"
wget -i $links_text_file -N --show-progress -q -P $out_dir
echo "Downloads checked"

