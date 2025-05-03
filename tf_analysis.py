import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

# Функция для чтения FASTA файла
def read_fasta(file_path):
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            return str(record.seq)
    return None

# Определение основных мотивов связывания транскрипционных факторов
tf_motifs = {
    'TATA-box': ['TATAAA', 'TATATA'],
    'CAAT-box': ['CCAAT', 'ATTGG'],
    'GC-box': ['GGGCGG', 'CCGCCC'],
    'E-box': ['CANNTG'],
    'CRE': ['TGACGTCA'],
    'AP-1': ['TGAGTCA'],
    'NF-kB': ['GGGACTTTCC'],
    'SP1': ['GGGCGG'],
    'GATA': ['WGATAR'],
    'HRE': ['RCGTG'],
    'STAT': ['TTCNNNNGAA'],
    'NFAT': ['GGAAA'],
    'ETS': ['GGAA', 'TTCC'],
    'IRF': ['GAAANN'],
    'CREB': ['TGACG']
}

# Функция для поиска мотивов в последовательности
def find_motifs(sequence, motifs_dict):
    results = {}
    
    for tf_name, patterns in motifs_dict.items():
        results[tf_name] = []
        
        for pattern in patterns:
            if 'N' in pattern:
                # Обработка паттернов с N (любой нуклеотид)
                pattern_regex = pattern.replace('NN', '..')
                pattern_regex = pattern_regex.replace('N', '.')
                
                # Обработка IUPAC кодов
                pattern_regex = pattern_regex.replace('R', '[AG]')
                pattern_regex = pattern_regex.replace('Y', '[CT]')
                pattern_regex = pattern_regex.replace('W', '[AT]')
                
                matches = re.finditer(pattern_regex, sequence)
                for match in matches:
                    results[tf_name].append((match.start(), match.group()))
            else:
                # Прямой поиск точных совпадений
                pos = 0
                while True:
                    pos = sequence.find(pattern, pos)
                    if pos == -1:
                        break
                    results[tf_name].append((pos, pattern))
                    pos += 1
                
                # Поиск в обратной комплементарной последовательности
                rev_comp = str(Seq(pattern).reverse_complement())
                pos = 0
                while True:
                    pos = sequence.find(rev_comp, pos)
                    if pos == -1:
                        break
                    results[tf_name].append((pos, rev_comp + ' (rev)'))
                    pos += 1
    
    return results

# Функция для визуализации результатов
def visualize_results(motifs_results, sequence_length):
    # Считаем количество сайтов для каждого TF
    tf_counts = {tf: len(sites) for tf, sites in motifs_results.items()}
    
    # Создаем DataFrame для графиков
    df_counts = pd.DataFrame(list(tf_counts.items()), columns=['Transcription Factor', 'Binding Sites Count'])
    df_counts = df_counts.sort_values('Binding Sites Count', ascending=False)
    
    # Создаем DataFrame с позициями сайтов
    positions_data = []
    for tf, sites in motifs_results.items():
        for pos, pattern in sites:
            positions_data.append({
                'Transcription Factor': tf,
                'Position': pos,
                'Pattern': pattern
            })
    
    df_positions = pd.DataFrame(positions_data)
    
    # Отображаем график количества сайтов
    plt.figure(figsize=(12, 6))
    sns.barplot(x='Transcription Factor', y='Binding Sites Count', data=df_counts)
    plt.xticks(rotation=45, ha='right')
    plt.title('Количество сайтов связывания для каждого транскрипционного фактора')
    plt.tight_layout()
    plt.savefig('tf_binding_sites_count.png')
    
    # Отображаем распределение сайтов вдоль последовательности
    plt.figure(figsize=(14, 8))
    sns.scatterplot(x='Position', y='Transcription Factor', data=df_positions)
    plt.title('Расположение сайтов связывания транскрипционных факторов')
    plt.xlim(0, sequence_length)
    plt.tight_layout()
    plt.savefig('tf_binding_sites_distribution.png')
    
    # Создаем тепловую карту плотности сайтов связывания
    plt.figure(figsize=(14, 8))
    tf_density = {}
    for tf in motifs_results:
        positions = [pos for pos, _ in motifs_results[tf]]
        density = Counter([pos // 100 for pos in positions])  # Разбиваем на участки по 100 нуклеотидов
        tf_density[tf] = [density.get(i, 0) for i in range((sequence_length // 100) + 1)]
    
    df_density = pd.DataFrame(tf_density).transpose()
    sns.heatmap(df_density, cmap='viridis')
    plt.title('Плотность сайтов связывания транскрипционных факторов')
    plt.xlabel('Позиция (x100 нуклеотидов)')
    plt.tight_layout()
    plt.savefig('tf_binding_sites_heatmap.png')
    
    return df_counts, df_positions

# Основная функция
def main():
    # Путь к файлу FASTA
    fasta_path = 'Homo_sapiens_ATF3_sequence (1).fa'
    
    # Чтение последовательности
    sequence = read_fasta(fasta_path)
    if not sequence:
        print(f"Не удалось прочитать файл {fasta_path}")
        return
    
    print(f"Загружена последовательность длиной {len(sequence)} нуклеотидов")
    
    # Поиск мотивов транскрипционных факторов
    motifs_results = find_motifs(sequence, tf_motifs)
    
    # Выводим результаты поиска для каждого фактора
    for tf, sites in motifs_results.items():
        if sites:
            print(f"{tf}: найдено {len(sites)} сайтов связывания")
            # Выводим первые 5 сайтов для примера
            for i, (pos, pattern) in enumerate(sites[:5]):
                print(f"  - Позиция {pos}: {pattern}")
            if len(sites) > 5:
                print(f"  ... и еще {len(sites) - 5} сайтов")
        else:
            print(f"{tf}: сайты связывания не найдены")
    
    # Визуализация результатов
    df_counts, df_positions = visualize_results(motifs_results, len(sequence))
    
    # Сохраняем результаты в CSV файлы
    df_counts.to_csv('tf_binding_sites_counts.csv', index=False)
    df_positions.to_csv('tf_binding_sites_positions.csv', index=False)
    
    print("\nРезультаты сохранены в файлы:")
    print("- tf_binding_sites_counts.csv")
    print("- tf_binding_sites_positions.csv")
    print("- tf_binding_sites_count.png")
    print("- tf_binding_sites_distribution.png")
    print("- tf_binding_sites_heatmap.png")

if __name__ == "__main__":
    main() 