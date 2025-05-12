#!/usr/bin/bash
file_list=$(find . -path './doc*' -prune -o -iname '*.py' -exec grep {} -lie '\<Heaviside_time_domain\>' \;)
echo "$file_list"
sed -i 's/\<Heaviside_time_domain\>/heaviside_time_domain/g' $file_list
