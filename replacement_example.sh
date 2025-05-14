#!/usr/bin/bash
file_list=$(find . -path './doc*' -prune -o -iname '*.py' -exec grep {} -lie '\<Delta_p_mask_fn\>' \;)
echo "$file_list"
sed -i 's/\<Delta_p_mask_fn\>/coherence_mask/g' $file_list
