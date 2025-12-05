#!/bin/bash

# Output file
OUTPUT="/Users/apple/Research/FlashDeconv/docs/merged_papers.md"

# Clear output file
> "$OUTPUT"

# Add header
echo "# Merged Reference Papers" >> "$OUTPUT"
echo "" >> "$OUTPUT"
echo "Generated on: $(date)" >> "$OUTPUT"
echo "" >> "$OUTPUT"
echo "---" >> "$OUTPUT"
echo "" >> "$OUTPUT"

# Paper 1
echo "# Paper 1: s41467-023-37168-7" >> "$OUTPUT"
echo "" >> "$OUTPUT"
cat "/Users/apple/Research/FlashDeconv/docs/s41467-023-37168-7 (1).md" >> "$OUTPUT"
echo "" >> "$OUTPUT"
echo "---" >> "$OUTPUT"
echo "" >> "$OUTPUT"

# Paper 2
echo "# Paper 2: s41467-023-43600-9" >> "$OUTPUT"
echo "" >> "$OUTPUT"
cat "/Users/apple/Research/FlashDeconv/docs/s41467-023-43600-9 (1).md" >> "$OUTPUT"
echo "" >> "$OUTPUT"
echo "---" >> "$OUTPUT"
echo "" >> "$OUTPUT"

# Paper 3
echo "# Paper 3: s41592-022-01480-9" >> "$OUTPUT"
echo "" >> "$OUTPUT"
cat "/Users/apple/Research/FlashDeconv/docs/s41592-022-01480-9.md" >> "$OUTPUT"
echo "" >> "$OUTPUT"
echo "---" >> "$OUTPUT"
echo "" >> "$OUTPUT"

# Paper 4
echo "# Paper 4: SPOTlight" >> "$OUTPUT"
echo "" >> "$OUTPUT"
cat "/Users/apple/Research/FlashDeconv/docs/spotless.md" >> "$OUTPUT"
echo "" >> "$OUTPUT"

echo "Merged papers saved to: $OUTPUT"
wc -l "$OUTPUT"
