#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import pyjaspar
from Bio import motifs
import numpy as np
from Bio.Seq import Seq
import argparse

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
    'MA0265.1',  # ABF1
    'MA0059.3',  # MYC::MAX
    'MA0099.3',  # FOS::JUN
    'MA0137.3',  # REST
    'MA0107.1',  # RELA
    'MA0043.2',  # HLF
    'MA0517.1',  # STAT1::STAT2
    'MA0142.1',  # Pou5f1::Sox2
    # Добавленные ТФ:
    'MA0138.2',  # REST
    'MA0144.2',  # STAT3
    'MA0145.3',  # Tcf7l2
    'MA0154.4',  # EBF1
    'MA0163.1',  # PLAG1
    'MA0258.2',  # ESR1
    'MA0259.1',  # HIF1A::ARNT
    'MA0442.2',  # SOX10
    'MA0491.2',  # JUND
    'MA0736.1',  # NR2F2
    'MA0778.1',  # NRF1
    'MA0778.1',  # NRF1
    'MA0867.2',  # E2F4
    'MA1124.1',  # ZNF24
    'MA1416.1',  # FOXP3
    'MA1513.1',  # BCL6
    'MA1535.1',  # FOXL1
    'MA1634.1',  # LHX1
    'MA1635.1',  # RUNX3
    'MA1684.1',  # NKX2-2
]

# Define a custom JSON encoder for BioPython objects
class BioEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Seq):
            return str(obj)
        return super().default(obj)

def download_all_human_tfs():
    """Download all human TFs from JASPAR database."""
    jaspar_db = pyjaspar.jaspardb()
    
    # Получаем все матрицы с фильтром для человека (species_id=9606)
    all_human_tfs = jaspar_db.fetch_motifs(species=[9606])
    
    # Формируем список ID для всех найденных ТФ
    human_tf_ids = [motif.matrix_id for motif in all_human_tfs]
    print(f"Найдено {len(human_tf_ids)} человеческих транскрипционных факторов в JASPAR")
    
    return human_tf_ids

def download_tfs(use_extended=False, max_tfs=None):
    """Download specific TFs from JASPAR and save them to a file.
    
    Args:
        use_extended: Если True, скачивает все человеческие ТФ из JASPAR
        max_tfs: Максимальное количество ТФ для загрузки (None - без ограничения)
    """
    jaspar_db = pyjaspar.jaspardb()
    
    # Dictionary to store TF motifs
    tf_motifs = {}
    
    # Определяем список ТФ для загрузки
    tf_ids_to_download = KNOWN_TF_IDS
    
    if use_extended:
        try:
            tf_ids_to_download = download_all_human_tfs()
            if max_tfs and max_tfs > 0:
                tf_ids_to_download = tf_ids_to_download[:max_tfs]
        except Exception as e:
            print(f"Ошибка при получении списка всех человеческих ТФ: {e}")
            print("Используем стандартный список ТФ")
    
    print(f"Загрузка {len(tf_ids_to_download)} транскрипционных факторов из JASPAR...")
    
    for tf_id in tf_ids_to_download:
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
                            print(f"  Пропускаем {tf_id} ({motif.name}): неверный результат расчета PSSM")
                            continue
                    
                    print(f"  Загружен {tf_id} ({motif.name})")
                    
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
                    print(f"  Ошибка с PSSM для {tf_id}: {e}")
            else:
                print(f"  Не удалось найти мотив {tf_id}")
        except Exception as e:
            print(f"  Ошибка при загрузке {tf_id}: {e}")
    
    # Save the results to a JSON file
    with open('known_tfs.json', 'w') as f:
        json.dump(tf_motifs, f, indent=2, cls=BioEncoder)
    
    print(f"Загружено {len(tf_motifs)} из {len(tf_ids_to_download)} транскрипционных факторов")
    return tf_motifs

def parse_args():
    parser = argparse.ArgumentParser(description='Загрузка транскрипционных факторов из базы данных JASPAR')
    parser.add_argument('--extended', action='store_true', help='Загрузить все человеческие ТФ из JASPAR')
    parser.add_argument('--max-tfs', type=int, default=None, help='Максимальное количество ТФ для загрузки')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    download_tfs(use_extended=args.extended, max_tfs=args.max_tfs) 