#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import pyjaspar
from Bio import motifs
import numpy as np
from Bio.Seq import Seq

# List of well-known, human TF IDs from JASPAR
KNOWN_TF_IDS = [
    'MA0139.1',  # CTCF
    'MA0148.4',  # FOXA1
    'MA0037.3',  # GATA3
    'MA0079.3',  # SP1
    'MA0147.3',  # MYC
    'MA0104.4',  # MYCN
    'MA0098.3',  # ETS1
    'MA0062.2',  # GABPA
    'MA0114.4',  # HNF4A
    'MA0106.3',  # TP53
    'MA0105.4',  # NFKB1
    'MA0080.5',  # SPI1
    'MA0265.1',  # JUN
    'MA0059.3',  # MYC::MAX
    'MA0099.3',  # FOS::JUN
    'MA0137.3',  # REST
    'MA0107.1',  # RELA
    'MA0043.2',  # HLF
    'MA0517.1',  # STAT1::STAT2
    'MA0142.1',  # POUX5
]

# Define a custom JSON encoder for BioPython objects
class BioEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Seq):
            return str(obj)
        return super().default(obj)

def download_tfs():
    """Download specific TFs from JASPAR and save them to a file."""
    jaspar_db = pyjaspar.jaspardb()
    
    # Dictionary to store TF motifs
    tf_motifs = {}
    
    print(f"Downloading {len(KNOWN_TF_IDS)} TFs from JASPAR...")
    
    for tf_id in KNOWN_TF_IDS:
        try:
            motif = jaspar_db.fetch_motif_by_id(tf_id)
            if motif:
                # Verify that the motif has a valid PSSM
                try:
                    pssm = motif.pssm
                    # Test calculating on a simple sequence
                    test_seq = Seq("ACTGACTGACTG")
                    if len(test_seq) >= len(motif):
                        score = pssm.calculate(test_seq[:len(motif)])
                        if np.isnan(score) or np.isinf(score):
                            print(f"  Skipping {tf_id} ({motif.name}): Invalid PSSM calculation result")
                            continue
                    
                    print(f"  Downloaded {tf_id} ({motif.name})")
                    
                    # Create a serializable representation of the position frequency matrix
                    pfm = {}
                    for base in motif.counts:
                        pfm[base] = [float(x) for x in motif.counts[base]]
                    
                    # Store the motif details
                    tf_motifs[tf_id] = {
                        'name': motif.name,
                        'pfm': pfm,  # Position frequency matrix
                        'consensus': str(motif.consensus),
                        'length': len(motif),
                        'id': motif.matrix_id
                    }
                except Exception as e:
                    print(f"  Error with PSSM for {tf_id}: {e}")
            else:
                print(f"  Could not find motif {tf_id}")
        except Exception as e:
            print(f"  Error fetching {tf_id}: {e}")
    
    # Save the results to a JSON file
    with open('known_tfs.json', 'w') as f:
        json.dump(tf_motifs, f, indent=2, cls=BioEncoder)
    
    print(f"Downloaded {len(tf_motifs)} out of {len(KNOWN_TF_IDS)} TFs")
    return tf_motifs

if __name__ == "__main__":
    download_tfs() 