import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, motifs
from Bio.Seq import Seq
from collections import Counter
import logomaker

# Функция для чтения FASTA файла
def read_fasta(file_path):
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            return str(record.seq)
    return None

# Функция для извлечения последовательностей вокруг сайта связывания
def extract_binding_site_sequences(sequence, positions, site_length, context_length=5):
    extracted_seqs = []
    for pos, pattern in positions:
        pattern_len = len(pattern.split(' ')[0])  # Учитываем только сам паттерн, без "(rev)"
        start = max(0, pos - context_length)
        end = min(len(sequence), pos + pattern_len + context_length)
        context_seq = sequence[start:end]
        extracted_seqs.append((pos, pattern, context_seq))
    return extracted_seqs

# Функция для построения PWM (Position Weight Matrix) и логотипа мотива
def create_motif_logo(sequences, tf_name, output_dir='tf_logos'):
    if not sequences:
        print(f"Нет последовательностей для создания логотипа {tf_name}")
        return
    
    # Убедимся, что все последовательности одинаковой длины
    seq_length = len(sequences[0])
    seqs = [seq for seq in sequences if len(seq) == seq_length]
    
    if not seqs:
        print(f"Нет последовательностей одинаковой длины для {tf_name}")
        return
    
    # Создаем каунт-матрицу
    counts_mat = logomaker.alignment_to_matrix(seqs)
    
    # Создаем логотип
    os.makedirs(output_dir, exist_ok=True)
    plt.figure(figsize=(10, 3))
    logo = logomaker.Logo(counts_mat, color_scheme='classic')
    
    # Настраиваем внешний вид
    logo.style_spines(visible=False)
    logo.style_xticks(rotation=0, fmt='%d', anchor=0)
    
    plt.title(f'Motif binding site for {tf_name}')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{tf_name}_logo.png')
    plt.close()
    
    return counts_mat

# Функция для создания тепловой карты частот нуклеотидов в мотиве
def create_heatmap(matrix, tf_name, output_dir='tf_heatmaps'):
    if matrix is None:
        return
    
    os.makedirs(output_dir, exist_ok=True)
    plt.figure(figsize=(10, 4))
    sns.heatmap(matrix, cmap='viridis', annot=True, fmt='.2f')
    plt.title(f'Nucleotide frequencies in {tf_name} motif')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{tf_name}_heatmap.png')
    plt.close()

# Функция для анализа и визуализации
def analyze_tf_binding_sites(sequence_file, positions_file):
    # Чтение последовательности
    sequence = read_fasta(sequence_file)
    if not sequence:
        print(f"Could not read file {sequence_file}")
        return
    
    # Чтение позиций сайтов связывания
    df_positions = pd.read_csv(positions_file)
    
    # Группировка по транскрипционным факторам
    tf_groups = df_positions.groupby('Transcription Factor')
    
    # Для каждого транскрипционного фактора
    for tf_name, group in tf_groups:
        print(f"Analyzing motifs for {tf_name}...")
        
        # Получаем позиции и паттерны
        positions = [(row['Position'], row['Pattern']) for _, row in group.iterrows()]
        
        # Определяем среднюю длину сайта связывания для этого TF
        pattern_lengths = [len(pattern.split(' ')[0]) for _, pattern in positions]
        avg_length = int(np.mean(pattern_lengths)) if pattern_lengths else 0
        
        # Извлекаем последовательности вокруг сайтов связывания
        binding_sites = extract_binding_site_sequences(sequence, positions, avg_length)
        
        # Готовим последовательности для анализа мотива
        context_sequences = [site for _, _, site in binding_sites]
        
        if len(context_sequences) > 5:  # Минимум 5 последовательностей для надежного анализа
            # Обрезаем до одинаковой длины (минимальной)
            min_length = min(len(seq) for seq in context_sequences)
            aligned_seqs = [seq[:min_length] for seq in context_sequences]
            
            # Создаем мотив
            matrix = create_motif_logo(aligned_seqs, tf_name)
            
            # Создаем тепловую карту
            create_heatmap(matrix, tf_name)
            
            # Сохраняем примеры сайтов связывания
            save_example_sites(tf_name, binding_sites[:10])

# Функция для сохранения примеров сайтов связывания
def save_example_sites(tf_name, sites, output_dir='tf_examples'):
    os.makedirs(output_dir, exist_ok=True)
    
    with open(f'{output_dir}/{tf_name}_examples.txt', 'w', encoding='utf-8') as f:
        f.write(f"Binding site examples for {tf_name}:\n\n")
        for pos, pattern, context in sites:
            # Определяем индекс начала паттерна в контексте
            pattern_clean = pattern.split(' ')[0]
            pattern_start = context.find(pattern_clean)
            
            # Если паттерн найден в контексте
            if pattern_start != -1:
                # Подготавливаем строку с выделенным паттерном
                before = context[:pattern_start]
                after = context[pattern_start + len(pattern_clean):]
                
                f.write(f"Position {pos}: {before}[{pattern_clean}]{after}\n")
            else:
                f.write(f"Position {pos}: {context} (pattern: {pattern})\n")
        
# Основная функция
def main():
    # Пути к файлам
    sequence_file = 'Homo_sapiens_ATF3_sequence (1).fa'
    positions_file = 'tf_binding_sites_positions.csv'
    
    # Анализ и визуализация
    analyze_tf_binding_sites(sequence_file, positions_file)
    
    print("\nAnalysis completed. Results saved in directories:")
    print("- tf_logos/ - motif logos")
    print("- tf_heatmaps/ - nucleotide frequency heatmaps")
    print("- tf_examples/ - binding site examples")

if __name__ == "__main__":
    main() 