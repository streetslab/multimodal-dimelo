import gzip
from collections import defaultdict

# Group CDS by transcript
transcripts = defaultdict(list)

with gzip.open('/clusterfs/nilah/oberon/datasets/annotations/hs1.ncbiRefSeq.gtf.gz', 'rt') as gtf:
    for line in gtf:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9 or fields[2] != 'CDS':
            continue
        
        chrom = fields[0]
        start = int(fields[3]) - 1  # Convert to 0-based
        end = int(fields[4])
        strand = fields[6]
        
        # Extract transcript_id
        transcript_id = None
        for attr in fields[8].split(';'):
            if 'transcript_id' in attr:
                transcript_id = attr.split('"')[1]
                break
        
        # Filter for curated transcripts only (NM_ and NR_ prefixes)
        # Skip predicted models (XM_, XR_)
        if transcript_id and (transcript_id.startswith('NM_') or transcript_id.startswith('NR_')):
            transcripts[transcript_id].append((chrom, start, end, strand))

# Find stop codon positions
with open('/clusterfs/nilah/oberon/datasets/annotations/hs1_stop_codons_high_confidence.bed', 'w') as out:
    for transcript_id, cds_list in transcripts.items():
        if not cds_list:
            continue
        
        cds_list.sort(key=lambda x: x[1])
        chrom, _, _, strand = cds_list[0]
        
        if strand == '+':
            stop_start = cds_list[-1][2]
            stop_end = stop_start + 3
        else:
            stop_end = cds_list[0][1]
            stop_start = stop_end - 3
        
        out.write(f"{chrom}\t{stop_start}\t{stop_end}\t{transcript_id}\t.\t{strand}\n")

print(f"Done! Found {len(transcripts)} high-confidence transcripts with stop codons")