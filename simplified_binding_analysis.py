#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import json
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from Bio import motifs
import matplotlib.pyplot as plt
import seaborn as sns

class SimplifiedBindingAnalyzer:
    def __init__(self, genes_dir='genes', results_dir='results'):
        self.genes_dir = genes_dir
        self.results_dir = results_dir
        self.gene_sequences = {}
        self.binding_results = {}
        
        # Create results directory if it doesn't exist
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
    
    def load_gene_sequences(self):
        """Load all gene sequences from FASTA files in the genes directory."""
        gene_files = glob.glob(os.path.join(self.genes_dir, '*.fa'))
        
        for gene_file in gene_files:
            gene_name = os.path.basename(gene_file).split('.')[0]
            
            try:
                for record in SeqIO.parse(gene_file, "fasta"):
                    self.gene_sequences[gene_name] = {
                        'id': record.id,
                        'description': record.description,
                        'sequence': str(record.seq)
                    }
                    break  # Take only the first record
            except Exception as e:
                print(f"Error loading gene {gene_name}: {e}")
        
        print(f"Loaded {len(self.gene_sequences)} gene sequences")
        return self.gene_sequences
    
    def load_motifs(self, json_file='known_tfs.json'):
        """Load motifs from a JSON file."""
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        motifs_dict = {}
        
        for tf_id, tf_data in data.items():
            try:
                # Create a position weight matrix from the PFM
                pwm = {}
                for base, counts in tf_data['pfm'].items():
                    pwm[base] = counts
                
                # Create a motif object
                tf_motif = motifs.Motif(alphabet='ACGT', counts=pwm)
                tf_motif.matrix_id = tf_id
                tf_motif.name = tf_data['name']
                
                motifs_dict[tf_id] = tf_motif
            except Exception as e:
                print(f"Error creating motif for {tf_id}: {e}")
        
        print(f"Loaded {len(motifs_dict)} motifs")
        return motifs_dict
    
    def find_binding_sites(self, motif, sequence, threshold=0.75):
        """Find binding sites for a given motif in a DNA sequence."""
        try:
            motif_len = len(motif)
            seq_len = len(sequence)
            
            # Create position-specific scoring matrix
            pssm = motif.pssm
            
            sites = []
            scores = []
            
            # Get PSSM min and max - handle NaN values
            pssm_min = pssm.min if not np.isnan(pssm.min) else 0
            pssm_max = pssm.max if not np.isnan(pssm.max) else 1
            
            # Avoid division by zero
            if pssm_max == pssm_min:
                pssm_max = pssm_min + 1
            
            # Scan the sequence for potential binding sites
            for i in range(seq_len - motif_len + 1):
                subseq = sequence[i:i+motif_len]
                
                # Skip if the sequence contains invalid nucleotides
                if not set(subseq.upper()).issubset({'A', 'C', 'G', 'T'}):
                    continue
                
                try:
                    # Score for forward strand
                    score = pssm.calculate(Seq(subseq))
                    
                    # Calculate relative score (0-1 scale) - handle NaN values
                    rel_score = (score - pssm_min) / (pssm_max - pssm_min)
                    
                    # If score passes threshold, record the site
                    if not np.isnan(rel_score) and rel_score >= threshold:
                        # Position relative to TSS (assuming -2000 to +500)
                        position = i - 2000
                        
                        sites.append({
                            'position': position,
                            'sequence': subseq,
                            'score': float(score),
                            'relative_score': float(rel_score),
                            'strand': '+'
                        })
                        scores.append(float(rel_score))
                    
                    # Check reverse complement
                    rev_seq = str(Seq(subseq).reverse_complement())
                    rev_pssm = pssm.reverse_complement()
                    rev_score = rev_pssm.calculate(Seq(rev_seq))
                    
                    # Calculate reverse score
                    rev_rel_score = (rev_score - pssm_min) / (pssm_max - pssm_min)
                    
                    if not np.isnan(rev_rel_score) and rev_rel_score >= threshold:
                        sites.append({
                            'position': position,
                            'sequence': subseq,
                            'score': float(rev_score),
                            'relative_score': float(rev_rel_score),
                            'strand': '-'
                        })
                        scores.append(float(rev_rel_score))
                except Exception as e:
                    # Skip this site if there's an error
                    continue
            
            return {
                'sites': sites,
                'total_sites': len(sites),
                'avg_score': float(np.mean(scores)) if scores else 0,
                'max_score': float(np.max(scores)) if scores else 0
            }
        except Exception as e:
            print(f"Error in find_binding_sites: {e}")
            return {
                'sites': [],
                'total_sites': 0,
                'avg_score': 0,
                'max_score': 0
            }
    
    def analyze_binding_sites(self, threshold=0.75):
        """Analyze binding sites for all genes using all motifs."""
        # Load gene sequences if not already loaded
        if not self.gene_sequences:
            self.load_gene_sequences()
        
        # Load motifs
        motifs_dict = self.load_motifs()
        
        # Process each motif
        for tf_id, motif in motifs_dict.items():
            print(f"Processing TF: {tf_id} ({motif.name})")
            tf_results = {}
            
            # Analyze each gene sequence
            for gene_name, gene_data in self.gene_sequences.items():
                binding_data = self.find_binding_sites(motif, gene_data['sequence'], threshold)
                tf_results[gene_name] = binding_data
                
                print(f"  Gene {gene_name}: Found {binding_data['total_sites']} binding sites")
            
            # Store results
            self.binding_results[tf_id] = {
                'tf_name': motif.name,
                'matrix_id': tf_id,
                'binding_data': tf_results
            }
        
        return self.binding_results
    
    def save_results(self, filename=None):
        """Save binding analysis results to a JSON file."""
        if filename is None:
            filename = os.path.join(self.results_dir, 'simplified_binding_results.json')
        
        # Convert results to JSON-serializable format
        json_results = {}
        for tf_id, tf_data in self.binding_results.items():
            json_results[tf_id] = {
                'tf_name': tf_data['tf_name'],
                'matrix_id': tf_data['matrix_id'],
                'binding_data': {}
            }
            
            for gene_name, binding_data in tf_data['binding_data'].items():
                # Convert sites to serializable format
                sites_json = []
                for site in binding_data['sites']:
                    sites_json.append({
                        'position': int(site['position']),
                        'sequence': str(site['sequence']),
                        'score': float(site['score']),
                        'relative_score': float(site['relative_score']),
                        'strand': str(site['strand'])
                    })
                
                json_results[tf_id]['binding_data'][gene_name] = {
                    'total_sites': int(binding_data['total_sites']),
                    'avg_score': float(binding_data['avg_score']),
                    'max_score': float(binding_data['max_score']),
                    'sites': sites_json
                }
        
        try:
            with open(filename, 'w') as f:
                json.dump(json_results, f, indent=2)
        except TypeError as e:
            print(f"Error serializing to JSON: {e}")
            # Try a more robust approach with default converter
            with open(filename, 'w') as f:
                json.dump(json_results, f, indent=2, default=str)
        
        print(f"Saved binding results to {filename}")
        return filename
    
    def generate_summary_dataframe(self):
        """Generate a summary dataframe of binding sites per gene per TF."""
        if not self.binding_results:
            print("No binding results available.")
            return None
        
        summary_data = []
        
        for tf_id, tf_data in self.binding_results.items():
            tf_name = tf_data['tf_name']
            
            for gene_name, binding_data in tf_data['binding_data'].items():
                summary_data.append({
                    'TF_ID': tf_id,
                    'TF_Name': tf_name,
                    'Gene': gene_name,
                    'Binding_Sites': binding_data['total_sites'],
                    'Avg_Score': binding_data['avg_score'],
                    'Max_Score': binding_data['max_score']
                })
        
        df = pd.DataFrame(summary_data)
        
        # Save to CSV
        csv_file = os.path.join(self.results_dir, 'binding_sites_summary.csv')
        df.to_csv(csv_file, index=False)
        print(f"Saved summary to {csv_file}")
        
        return df
    
    def generate_heatmap(self, output_file=None):
        """Generate a heatmap of binding sites per gene per TF."""
        df = self.generate_summary_dataframe()
        
        if df is None or df.empty:
            print("No data available for heatmap.")
            return None
        
        # Pivot the dataframe to create a matrix of genes vs TFs
        heatmap_data = df.pivot_table(
            values='Binding_Sites', 
            index='Gene', 
            columns='TF_Name', 
            fill_value=0
        )
        
        # Create the heatmap
        plt.figure(figsize=(12, 10))
        sns.heatmap(heatmap_data, cmap='viridis', annot=True, fmt='.0f', linewidths=0.5)
        plt.title('Transcription Factor Binding Sites per Gene')
        plt.tight_layout()
        
        if output_file is None:
            output_file = os.path.join(self.results_dir, 'binding_heatmap.png')
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Generated heatmap at {output_file}")
        return output_file
    
    def generate_binding_site_distribution(self, output_file=None):
        """Generate a histogram of binding site positions relative to TSS."""
        if not self.binding_results:
            print("No binding results available.")
            return None
        
        all_positions = []
        
        for tf_id, tf_data in self.binding_results.items():
            for gene_name, binding_data in tf_data['binding_data'].items():
                for site in binding_data['sites']:
                    all_positions.append(site['position'])
        
        plt.figure(figsize=(12, 6))
        
        if all_positions:
            plt.hist(all_positions, bins=50, alpha=0.7, color='navy')
            plt.title('Distribution of Binding Sites Relative to Transcription Start Site')
            plt.xlabel('Position Relative to TSS (bp)')
            plt.ylabel('Number of Binding Sites')
            plt.axvline(x=0, color='red', linestyle='--', label='TSS')
            plt.legend()
            plt.grid(True, alpha=0.3)
        else:
            plt.text(0.5, 0.5, "No binding sites found", 
                    horizontalalignment='center', verticalalignment='center',
                    transform=plt.gca().transAxes, fontsize=14)
            plt.title('No Binding Sites Found')
        
        if output_file is None:
            output_file = os.path.join(self.results_dir, 'binding_site_distribution.png')
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Generated binding site distribution plot at {output_file}")
        return output_file

if __name__ == "__main__":
    analyzer = SimplifiedBindingAnalyzer()
    analyzer.load_gene_sequences()
    analyzer.analyze_binding_sites(threshold=0.75)
    analyzer.save_results()
    analyzer.generate_heatmap()
    analyzer.generate_binding_site_distribution() 