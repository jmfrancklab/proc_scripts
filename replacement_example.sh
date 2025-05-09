#!/usr/bin/bash
file_list=$(find . -path './doc*' -prune -o -iname '*.py' -exec grep {} -lie '\<HH_weighted_integral\>' \;)
echo "$file_list"
sed -i 's/\<HH_weighted_integral\>/Heaviside_time_domain/g' $file_list
