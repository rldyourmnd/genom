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
    'CREB': ['TGACG'],
    # Добавляем дополнительные транскрипционные факторы, характерные для ATF
    'ATF/CREB': ['TGACGTCA', 'TGACGTCA'],
    'AP-2': ['GCCNNNGGC'],
    'OCT': ['ATGCAAAT'],
    'YY1': ['CCATNTT'],
    'NRF-1': ['GCGCATGCGC']
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
def visualize_results(motifs_results, sequence_length, output_prefix="atf2500"):
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
    plt.title('Количество сайтов связывания для каждого транскрипционного фактора в ATF2500')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_tf_binding_sites_count.png')
    
    # Отображаем распределение сайтов вдоль последовательности
    plt.figure(figsize=(14, 8))
    
    # Ограничиваемся топ-10 TF для лучшей визуализации
    top_tfs = df_counts.head(10)['Transcription Factor'].tolist()
    df_top_positions = df_positions[df_positions['Transcription Factor'].isin(top_tfs)]
    
    sns.scatterplot(x='Position', y='Transcription Factor', data=df_top_positions)
    plt.title('Расположение сайтов связывания транскрипционных факторов в ATF2500')
    plt.xlim(0, sequence_length)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_tf_binding_sites_distribution.png')
    
    # Создаем тепловую карту плотности сайтов связывания
    plt.figure(figsize=(14, 8))
    tf_density = {}
    bin_size = 50  # Размер бина (окна) для анализа плотности
    
    for tf in top_tfs:
        if tf in motifs_results:
            positions = [pos for pos, _ in motifs_results[tf]]
            density = Counter([pos // bin_size for pos in positions])
            bins = sequence_length // bin_size + 1
            tf_density[tf] = [density.get(i, 0) for i in range(bins)]
    
    df_density = pd.DataFrame(tf_density).transpose()
    sns.heatmap(df_density, cmap='viridis')
    plt.title('Плотность сайтов связывания транскрипционных факторов в ATF2500')
    plt.xlabel(f'Позиция (x{bin_size} нуклеотидов)')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_tf_binding_sites_heatmap.png')
    
    return df_counts, df_positions

# Основная функция
def main():
    # Путь к файлу FASTA
    fasta_path = 'ATF2500.fa'
    
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
    df_counts.to_csv('atf2500_tf_binding_sites_counts.csv', index=False)
    df_positions.to_csv('atf2500_tf_binding_sites_positions.csv', index=False)
    
    print("\nРезультаты анализа ATF2500 сохранены в файлы:")
    print("- atf2500_tf_binding_sites_counts.csv")
    print("- atf2500_tf_binding_sites_positions.csv")
    print("- atf2500_tf_binding_sites_count.png")
    print("- atf2500_tf_binding_sites_distribution.png")
    print("- atf2500_tf_binding_sites_heatmap.png")
    
    # Создание краткого отчета в Markdown
    create_markdown_report(df_counts, len(sequence))

def create_markdown_report(df_counts, sequence_length):
    """Создает отчет в формате Markdown с результатами анализа"""
    
    with open('atf2500_analysis_report.md', 'w', encoding='utf-8') as f:
        f.write("# Анализ транскрипционных факторов в последовательности ATF2500\n\n")
        
        f.write("## Общая информация\n\n")
        f.write(f"- Длина последовательности: {sequence_length} нуклеотидов\n")
        f.write(f"- Проанализировано: {len(df_counts)} типов транскрипционных факторов\n")
        f.write(f"- Всего найдено сайтов связывания: {df_counts['Binding Sites Count'].sum()}\n\n")
        
        f.write("## Топ-10 транскрипционных факторов по количеству сайтов связывания\n\n")
        f.write("| Транскрипционный фактор | Количество сайтов |\n")
        f.write("|--------------------------|-------------------|\n")
        
        for _, row in df_counts.head(10).iterrows():
            f.write(f"| {row['Transcription Factor']} | {row['Binding Sites Count']} |\n")
        
        f.write("\n## Визуализации\n\n")
        f.write("1. **atf2500_tf_binding_sites_count.png** - Гистограмма количества сайтов связывания\n")
        f.write("2. **atf2500_tf_binding_sites_distribution.png** - Распределение сайтов вдоль последовательности\n")
        f.write("3. **atf2500_tf_binding_sites_heatmap.png** - Тепловая карта плотности сайтов связывания\n\n")
        
        f.write("## Выводы\n\n")
        
        # Наиболее представленный фактор
        top_tf = df_counts.iloc[0]['Transcription Factor']
        top_count = df_counts.iloc[0]['Binding Sites Count']
        
        f.write(f"1. Наиболее представленным транскрипционным фактором является **{top_tf}** ({top_count} сайтов связывания)\n")
        f.write("2. Для данной последовательности характерно наличие следующих регуляторных элементов:\n")
        
        # Проверка наличия базовых элементов промотора
        core_elements = ["TATA-box", "CAAT-box", "GC-box"]
        for element in core_elements:
            if element in df_counts['Transcription Factor'].values:
                count = df_counts[df_counts['Transcription Factor'] == element]['Binding Sites Count'].values[0]
                f.write(f"   - **{element}**: {count} сайтов\n")
        
        # Проверка наличия факторов, связанных с ATF
        atf_elements = ["CRE", "CREB", "AP-1", "ATF/CREB"]
        for element in atf_elements:
            if element in df_counts['Transcription Factor'].values:
                count = df_counts[df_counts['Transcription Factor'] == element]['Binding Sites Count'].values[0]
                if count > 0:
                    f.write(f"   - **{element}**: {count} сайтов\n")
        
        f.write("\n3. Данный анализ предоставляет информацию о потенциальных сайтах связывания транскрипционных факторов.")
        f.write(" Для подтверждения функциональной значимости выявленных сайтов рекомендуется проведение экспериментальных исследований.\n")
    
    print("Создан отчет: atf2500_analysis_report.md")

if __name__ == "__main__":
    main() 