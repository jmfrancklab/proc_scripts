#!/usr/bin/bash
file_list=$(find . -path './doc*' -prune -o -iname '*.py' -exec grep {} -lie '\<integral_w_errors\>' \;)
echo "$file_list"
sed -i 's/\<integral_w_errors\>/frequency_domain_integral/g' $file_list
