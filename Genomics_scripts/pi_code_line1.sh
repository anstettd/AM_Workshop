#!/bin/bash

# create a new file called "output.txt"
touch output.txt

# loop through all .txt files in the directory
for file in *.txt; do
  # use the tail command to remove the first line from each file
  # and write the resulting output to "output.txt"
  tail -n +2 "$file" >> output.txt
done

# prepend the first line from the first .txt file to "output.txt"
head -n 1 *.txt >> output.txt