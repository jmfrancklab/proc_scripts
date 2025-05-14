#!/usr/bin/bash
file_list=$(find . -path './doc*' -prune -o -iname '*.py' -exec grep {} -lie '\<coherence_mask\>' \;)
echo "$file_list"
sed -i 's/\<coherence_mask\>/coherence_unmask_fn/g' $file_list
