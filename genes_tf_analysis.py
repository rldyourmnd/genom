import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from io import StringIO
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams

# Устанавливаем шрифт с поддержкой кириллицы
rcParams['font.family'] = 'DejaVu Sans'

# Определяем мотивы транскрипционных факторов
motifs = {
    # Базовые элементы промотора
    'TATA-box': ['TATAAA', 'TATAAT', 'TATAAG', 'TATAAC'],
    'CAAT-box': ['CCAAT', 'ATTGG'],
    'GC-box': ['GGGCGG', 'CCGCCC'],
    
    # Транскрипционные факторы
    'E-box': ['CANNTG', 'CACGTG', 'CAGCTG', 'CATGTG', 'CAACTG', 'CAGGTG', 'CACATG', 'CAAGTG', 'CATCTG', 'CAATTG', 'CACCTG', 'CAAATG'],
    'CRE': ['TGACGTCA'],
    'AP-1': ['TGAGTCA', 'TGACTCA'],
    'NF-kB': ['GGGACTTTCC', 'GGGRNNYYCC'],
    'SP1': ['GGGCGG', 'CCGCCC'],
    'GATA': ['WGATAR', 'GATA', 'GATAA', 'GATAG'],
    'HRE': ['RCGTG'],
    'STAT': ['TTCNNNGAA'],
    'NFAT': ['GGAAA', 'TTTCC'],
    'ETS': ['GGAA', 'TTCC'],
    'IRF': ['GAAANN', 'AANNNGAAA', 'NAANNNGAAA'],
    'CREB': ['TGACG', 'CGTCA'],
    'ATF/CREB': ['TGACGTCA'],
    'OCT': ['ATGCAAAT', 'ATTTGCAT'],
    'YY1': ['CCATNTT', 'AANGAT', 'AANGATGG'],
    'NRF-1': ['GCGCATGCGC'],
    'AP-2': ['GCCNNNGGC', 'GCCNNGGC'],
    'P53': ['RRRCWWGYYY', 'RRRCATGYYY'],
    'MYC': ['CACGTG'],
    'HSF': ['AGAAN', 'NTTCT'],
    'SRF': ['CCATATGG']
}

# Функция для чтения FASTA-файла
def read_fasta(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            # Проверяем формат FASTA
            if not content.startswith('>'):
                return None
            
            # Берем последовательность (все кроме заголовка)
            lines = content.split('\n')
            header = lines[0]
            sequence = ''.join(lines[1:]).upper()
            return {'header': header, 'sequence': sequence}
    except Exception as e:
        print(f"Ошибка при чтении файла {file_path}: {e}")
        return None

# Функция для создания обратно-комплементарной последовательности
def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 
                  'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 
                  'D': 'H', 'H': 'D', 'V': 'B'}
    return ''.join(complement.get(base, base) for base in reversed(sequence.upper()))

# Функция для поиска мотивов в последовательности
def find_motifs(sequence, motifs_dict):
    results = {}
    positions = {}
    
    # Преобразуем последовательность в верхний регистр
    sequence = sequence.upper()
    
    for tf, patterns in motifs_dict.items():
        results[tf] = 0
        positions[tf] = []
        
        for pattern in patterns:
            # Преобразуем IUPAC-код в регулярные выражения
            regex_pattern = pattern
            regex_pattern = regex_pattern.replace('R', '[AG]')
            regex_pattern = regex_pattern.replace('Y', '[CT]')
            regex_pattern = regex_pattern.replace('W', '[AT]')
            regex_pattern = regex_pattern.replace('S', '[GC]')
            regex_pattern = regex_pattern.replace('M', '[AC]')
            regex_pattern = regex_pattern.replace('K', '[GT]')
            regex_pattern = regex_pattern.replace('B', '[CGT]')
            regex_pattern = regex_pattern.replace('D', '[AGT]')
            regex_pattern = regex_pattern.replace('H', '[ACT]')
            regex_pattern = regex_pattern.replace('V', '[ACG]')
            regex_pattern = regex_pattern.replace('N', '[ACGT]')
            
            # Поиск прямого соответствия
            try:
                for match in re.finditer(regex_pattern, sequence):
                    results[tf] += 1
                    positions[tf].append({
                        'Position': match.start() + 1,
                        'Pattern': match.group(),
                        'Direction': 'forward'
                    })
            except re.error as e:
                print(f"Ошибка регулярного выражения для паттерна {pattern}: {e}")
                continue
            
            # Поиск обратного соответствия
            rev_comp_pattern = reverse_complement(pattern)
            rev_regex_pattern = rev_comp_pattern
            rev_regex_pattern = rev_regex_pattern.replace('R', '[AG]')
            rev_regex_pattern = rev_regex_pattern.replace('Y', '[CT]')
            rev_regex_pattern = rev_regex_pattern.replace('W', '[AT]')
            rev_regex_pattern = rev_regex_pattern.replace('S', '[GC]')
            rev_regex_pattern = rev_regex_pattern.replace('M', '[AC]')
            rev_regex_pattern = rev_regex_pattern.replace('K', '[GT]')
            rev_regex_pattern = rev_regex_pattern.replace('B', '[CGT]')
            rev_regex_pattern = rev_regex_pattern.replace('D', '[AGT]')
            rev_regex_pattern = rev_regex_pattern.replace('H', '[ACT]')
            rev_regex_pattern = rev_regex_pattern.replace('V', '[ACG]')
            rev_regex_pattern = rev_regex_pattern.replace('N', '[ACGT]')
            
            try:
                for match in re.finditer(rev_regex_pattern, sequence):
                    results[tf] += 1
                    positions[tf].append({
                        'Position': match.start() + 1,
                        'Pattern': match.group(),
                        'Direction': 'reverse'
                    })
            except re.error as e:
                print(f"Ошибка регулярного выражения для обратного паттерна {rev_comp_pattern}: {e}")
                continue
    
    return results, positions

# Функция для визуализации результатов
def visualize_results(gene_name, motifs_results, sequence_length, output_dir='results'):
    # Создаем директорию для результатов, если её нет
    os.makedirs(output_dir, exist_ok=True)
    gene_dir = os.path.join(output_dir, gene_name)
    os.makedirs(gene_dir, exist_ok=True)
    
    # Преобразуем результаты в DataFrame
    df = pd.DataFrame({
        'Transcription Factor': list(motifs_results.keys()),
        'Binding Sites Count': list(motifs_results.values())
    })
    
    # Сортируем по количеству сайтов
    df = df.sort_values(by='Binding Sites Count', ascending=False)
    
    # Добавляем нормализованные значения на 1000 bp
    df['Sites per 1000 bp'] = df['Binding Sites Count'] / sequence_length * 1000
    
    # Сохраняем результаты в CSV
    df.to_csv(os.path.join(gene_dir, f"{gene_name}_tf_binding_sites_counts.csv"), index=False)
    
    # Строим гистограмму
    plt.figure(figsize=(12, 8))
    sns.barplot(x='Transcription Factor', y='Binding Sites Count', data=df.head(15))
    plt.title(f'Количество сайтов связывания транскрипционных факторов в {gene_name}')
    plt.xlabel('Транскрипционный фактор')
    plt.ylabel('Количество сайтов связывания')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(gene_dir, f"{gene_name}_tf_binding_sites_count.png"), dpi=300)
    plt.close()
    
    # Строим гистограмму для нормализованных значений
    plt.figure(figsize=(12, 8))
    sns.barplot(x='Transcription Factor', y='Sites per 1000 bp', data=df.head(15))
    plt.title(f'Плотность сайтов связывания транскрипционных факторов в {gene_name}')
    plt.xlabel('Транскрипционный фактор')
    plt.ylabel('Сайтов на 1000 нуклеотидов')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(gene_dir, f"{gene_name}_tf_binding_sites_density.png"), dpi=300)
    plt.close()
    
    return df

# Функция для сравнения генов по транскрипционным факторам
def compare_genes(all_results, output_dir='results'):
    # Создаем директорию для сравнительных результатов
    comparison_dir = os.path.join(output_dir, 'comparison')
    os.makedirs(comparison_dir, exist_ok=True)
    
    # Собираем данные из всех генов
    all_data = []
    for gene_name, data in all_results.items():
        gene_df = pd.DataFrame({
            'Transcription Factor': list(data['motifs_results'].keys()),
            'Binding Sites Count': list(data['motifs_results'].values()),
            'Sites per 1000 bp': [count / data['sequence_length'] * 1000 for count in data['motifs_results'].values()],
            'Gene': gene_name
        })
        all_data.append(gene_df)
    
    # Объединяем данные
    combined_df = pd.concat(all_data)
    
    # Сохраняем полные данные
    combined_df.to_csv(os.path.join(comparison_dir, "all_genes_tf_comparison.csv"), index=False)
    
    # Создаем сводную таблицу для плотности сайтов
    pivot_df = combined_df.pivot(index='Transcription Factor', columns='Gene', values='Sites per 1000 bp')
    pivot_df = pivot_df.fillna(0)
    
    # Сохраняем сводную таблицу
    pivot_df.to_csv(os.path.join(comparison_dir, "tf_density_by_gene.csv"))
    
    # Строим тепловую карту для топ-15 транскрипционных факторов
    top_tfs = combined_df.groupby('Transcription Factor')['Binding Sites Count'].sum().nlargest(15).index
    pivot_subset = pivot_df.loc[top_tfs]
    
    plt.figure(figsize=(15, 10))
    sns.heatmap(pivot_subset, annot=True, cmap='YlGnBu', fmt='.2f')
    plt.title('Плотность сайтов связывания транскрипционных факторов по генам')
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, "tf_density_heatmap.png"), dpi=300)
    plt.close()
    
    # Строим кластерную тепловую карту для всех транскрипционных факторов
    plt.figure(figsize=(18, 14))
    sns.clustermap(pivot_df, cmap='YlGnBu', figsize=(18, 14), 
                  dendrogram_ratio=(.1, .2), standard_scale=1)
    plt.savefig(os.path.join(comparison_dir, "tf_density_clustermap.png"), dpi=300)
    plt.close()
    
    # Сравнение топ-5 транскрипционных факторов для каждого гена
    top5_by_gene = {}
    for gene_name, data in all_results.items():
        gene_df = pd.DataFrame({
            'Transcription Factor': list(data['motifs_results'].keys()),
            'Binding Sites Count': list(data['motifs_results'].values()),
            'Sites per 1000 bp': [count / data['sequence_length'] * 1000 for count in data['motifs_results'].values()]
        })
        top5_by_gene[gene_name] = gene_df.sort_values(by='Binding Sites Count', ascending=False).head(5)
    
    # Создаем отчет в markdown
    with open(os.path.join(comparison_dir, "genes_comparison_report.md"), 'w', encoding='utf-8') as f:
        f.write("# Сравнительный анализ транскрипционных факторов в промоторах генов\n\n")
        
        f.write("## Общая информация\n\n")
        f.write("| Ген | Длина последовательности | Всего сайтов связывания | Количество различных факторов |\n")
        f.write("|-----|--------------------------|--------------------------|--------------------------------|\n")
        
        for gene_name, data in all_results.items():
            num_factors = sum(1 for count in data['motifs_results'].values() if count > 0)
            total_sites = sum(data['motifs_results'].values())
            f.write(f"| {gene_name} | {data['sequence_length']} | {total_sites} | {num_factors} |\n")
        
        f.write("\n## Топ-5 транскрипционных факторов для каждого гена\n\n")
        
        for gene_name, top5_df in top5_by_gene.items():
            f.write(f"### {gene_name}\n\n")
            f.write("| Транскрипционный фактор | Количество сайтов | Сайтов на 1000 нуклеотидов |\n")
            f.write("|--------------------------|-------------------|------------------------------|\n")
            
            for _, row in top5_df.iterrows():
                tf = row['Transcription Factor']
                count = row['Binding Sites Count']
                density = row['Sites per 1000 bp']
                f.write(f"| {tf} | {count} | {density:.2f} |\n")
            
            f.write("\n")
        
        # Анализ общих транскрипционных факторов
        all_tfs = set()
        for data in all_results.values():
            present_tfs = {tf for tf, count in data['motifs_results'].items() if count > 0}
            all_tfs.update(present_tfs)
        
        common_tfs = set.intersection(*[{tf for tf, count in data['motifs_results'].items() if count > 0} 
                                      for data in all_results.values()])
        
        f.write("\n## Общие транскрипционные факторы\n\n")
        f.write(f"Все гены содержат сайты связывания для следующих {len(common_tfs)} транскрипционных факторов:\n\n")
        f.write(", ".join(sorted(common_tfs)))
        
        # Уникальные транскрипционные факторы для каждого гена
        f.write("\n\n## Уникальные транскрипционные факторы\n\n")
        
        for gene_name, data in all_results.items():
            present_tfs = {tf for tf, count in data['motifs_results'].items() if count > 0}
            other_genes_tfs = set.union(*[{tf for tf, count in other_data['motifs_results'].items() if count > 0} 
                                        for other_gene, other_data in all_results.items() 
                                        if other_gene != gene_name])
            
            unique_tfs = present_tfs - other_genes_tfs
            
            f.write(f"### {gene_name}\n\n")
            if unique_tfs:
                f.write(f"Уникальные транскрипционные факторы ({len(unique_tfs)}):\n\n")
                f.write(", ".join(sorted(unique_tfs)))
            else:
                f.write("Нет уникальных транскрипционных факторов")
            f.write("\n\n")
        
        # Ключевые выводы
        f.write("\n## Ключевые выводы\n\n")
        
        # Находим самые распространенные факторы
        tf_counts = {tf: sum(1 for data in all_results.values() 
                            if data['motifs_results'].get(tf, 0) > 0) 
                    for tf in all_tfs}
        
        most_common_tfs = sorted(tf_counts.items(), key=lambda x: x[1], reverse=True)[:5]
        
        f.write("1. Наиболее представленные транскрипционные факторы во всех генах:\n")
        for tf, count in most_common_tfs:
            f.write(f"   - {tf}: присутствует в {count} из {len(all_results)} генов\n")
        
        # Анализируем распределение плотности сайтов связывания
        avg_densities = pivot_df.mean(axis=1).sort_values(ascending=False)
        top_density_tfs = avg_densities.head(5)
        
        f.write("\n2. Транскрипционные факторы с наибольшей средней плотностью сайтов связывания:\n")
        for tf, density in top_density_tfs.items():
            f.write(f"   - {tf}: в среднем {density:.2f} сайтов на 1000 нуклеотидов\n")
        
        # Гены с наибольшим разнообразием факторов
        gene_tf_diversity = {gene: sum(1 for count in data['motifs_results'].values() if count > 0) 
                           for gene, data in all_results.items()}
        
        most_diverse_genes = sorted(gene_tf_diversity.items(), key=lambda x: x[1], reverse=True)[:3]
        
        f.write("\n3. Гены с наибольшим разнообразием транскрипционных факторов:\n")
        for gene, count in most_diverse_genes:
            f.write(f"   - {gene}: {count} различных транскрипционных факторов\n")
        
        f.write("\n4. Анализ кластеризации показывает, что гены группируются по паттернам транскрипционных факторов, что может указывать на функциональные связи между ними.\n")
        
        f.write("\n5. Для более детального понимания функциональных различий между этими генами рекомендуется экспериментальная валидация обнаруженных сайтов связывания.\n")

    # Создаем визуализацию сравнения топ-10 транскрипционных факторов по плотности
    top10_tfs = combined_df.groupby('Transcription Factor')['Sites per 1000 bp'].mean().nlargest(10).index
    top10_data = combined_df[combined_df['Transcription Factor'].isin(top10_tfs)]
    
    plt.figure(figsize=(15, 10))
    sns.barplot(x='Transcription Factor', y='Sites per 1000 bp', hue='Gene', data=top10_data)
    plt.title('Сравнение плотности сайтов связывания топ-10 транскрипционных факторов')
    plt.xlabel('Транскрипционный фактор')
    plt.ylabel('Сайтов на 1000 нуклеотидов')
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Ген', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, "top10_tf_density_comparison.png"), dpi=300)
    plt.close()
    
    return combined_df

# Основная функция для анализа всех генов
def analyze_all_genes(genes_dir='genes', output_dir='results'):
    # Создаем директорию для результатов
    os.makedirs(output_dir, exist_ok=True)
    
    # Получаем список FASTA-файлов
    fasta_files = [f for f in os.listdir(genes_dir) if f.endswith('.fa') or f.endswith('.fasta')]
    
    all_results = {}
    
    # Анализируем каждый ген
    for fasta_file in fasta_files:
        file_path = os.path.join(genes_dir, fasta_file)
        gene_name = os.path.splitext(fasta_file)[0]
        
        print(f"Анализ гена {gene_name}...")
        
        # Читаем FASTA-файл
        fasta_data = read_fasta(file_path)
        if not fasta_data:
            print(f"Ошибка при чтении файла {file_path}, пропускаем...")
            continue
        
        sequence = fasta_data['sequence']
        
        # Ищем мотивы транскрипционных факторов
        motifs_results, positions = find_motifs(sequence, motifs)
        
        # Визуализируем результаты
        results_df = visualize_results(gene_name, motifs_results, len(sequence), output_dir)
        
        # Сохраняем позиции сайтов в CSV
        positions_df = pd.DataFrame([
            {'Transcription Factor': tf, 'Position': pos['Position'], 
             'Pattern': pos['Pattern'], 'Direction': pos['Direction']}
            for tf, pos_list in positions.items() for pos in pos_list
        ])
        
        if not positions_df.empty:
            positions_df.to_csv(os.path.join(output_dir, gene_name, f"{gene_name}_tf_binding_sites_positions.csv"), index=False)
        
        # Сохраняем результаты
        all_results[gene_name] = {
            'motifs_results': motifs_results,
            'positions': positions,
            'sequence_length': len(sequence)
        }
        
        # Создаем отчет в markdown
        create_gene_report(gene_name, motifs_results, positions, len(sequence), output_dir)
    
    # Сравниваем гены
    if all_results:
        comparison_df = compare_genes(all_results, output_dir)
        print("Анализ завершен. Результаты сохранены в директории", output_dir)
    else:
        print("Не удалось проанализировать ни один ген")

# Функция для создания отчета по отдельному гену
def create_gene_report(gene_name, motifs_results, positions, sequence_length, output_dir='results'):
    gene_dir = os.path.join(output_dir, gene_name)
    os.makedirs(gene_dir, exist_ok=True)
    
    with open(os.path.join(gene_dir, f"{gene_name}_analysis_report.md"), 'w', encoding='utf-8') as f:
        f.write(f"# Анализ сайтов связывания транскрипционных факторов в промоторе гена {gene_name}\n\n")
        
        f.write("## Общая информация\n\n")
        f.write(f"* **Длина последовательности**: {sequence_length} нуклеотидов\n")
        
        # Общее количество сайтов
        total_sites = sum(motifs_results.values())
        f.write(f"* **Всего сайтов связывания**: {total_sites}\n")
        
        # Количество различных транскрипционных факторов
        present_tfs = sum(1 for count in motifs_results.values() if count > 0)
        f.write(f"* **Количество различных транскрипционных факторов**: {present_tfs}\n\n")
        
        f.write("## Распределение транскрипционных факторов\n\n")
        f.write("| Транскрипционный фактор | Количество сайтов | Сайтов на 1000 нуклеотидов |\n")
        f.write("|--------------------------|-------------------|------------------------------|\n")
        
        # Сортируем по количеству сайтов
        sorted_motifs = sorted(motifs_results.items(), key=lambda x: x[1], reverse=True)
        
        for tf, count in sorted_motifs:
            if count > 0:
                density = count / sequence_length * 1000
                f.write(f"| {tf} | {count} | {density:.2f} |\n")
        
        f.write("\n")
        
        # Добавляем примеры сайтов связывания для топ-5 факторов
        f.write("## Примеры сайтов связывания\n\n")
        
        top5_tfs = [tf for tf, _ in sorted_motifs[:5] if motifs_results[tf] > 0]
        
        for tf in top5_tfs:
            f.write(f"### {tf}\n\n")
            tf_positions = positions[tf]
            
            if tf_positions:
                f.write("| Позиция | Паттерн | Направление |\n")
                f.write("|---------|---------|-------------|\n")
                
                # Показываем до 10 примеров
                for pos_data in tf_positions[:10]:
                    f.write(f"| {pos_data['Position']} | {pos_data['Pattern']} | {pos_data['Direction']} |\n")
                
                if len(tf_positions) > 10:
                    f.write(f"\n... и еще {len(tf_positions) - 10} сайтов\n")
            else:
                f.write("Нет данных о позициях\n")
            
            f.write("\n")
        
        f.write("## Визуализации\n\n")
        f.write(f"![Количество сайтов связывания]({gene_name}_tf_binding_sites_count.png)\n\n")
        f.write(f"![Плотность сайтов связывания]({gene_name}_tf_binding_sites_density.png)\n\n")
        
        f.write("## Выводы\n\n")
        
        if total_sites > 0:
            # Топ-3 транскрипционных фактора
            top3_tfs = [tf for tf, count in sorted_motifs[:3] if count > 0]
            
            f.write(f"1. В промоторе гена {gene_name} обнаружено {total_sites} потенциальных сайтов связывания транскрипционных факторов.\n")
            
            if top3_tfs:
                f.write(f"2. Наиболее представленные факторы: {', '.join(top3_tfs)}.\n")
            
            # Проверяем наличие базовых элементов промотора
            core_elements = []
            if motifs_results.get('TATA-box', 0) > 0:
                core_elements.append("TATA-box")
            if motifs_results.get('CAAT-box', 0) > 0:
                core_elements.append("CAAT-box")
            if motifs_results.get('GC-box', 0) > 0:
                core_elements.append("GC-box")
            
            if core_elements:
                f.write(f"3. Обнаружены базовые элементы промотора: {', '.join(core_elements)}.\n")
            else:
                f.write("3. Базовые элементы промотора (TATA-box, CAAT-box, GC-box) не обнаружены.\n")
            
            # Особенности промотора
            if motifs_results.get('TATA-box', 0) > 0:
                f.write("4. Наличие TATA-box указывает на классический тип промотора.\n")
            elif motifs_results.get('GC-box', 0) > 0 or motifs_results.get('SP1', 0) > 0:
                f.write("4. Обогащение GC-box/SP1 сайтами характерно для CpG-островков и конститутивно экспрессирующихся генов.\n")
            
            # Регуляторные особенности
            regulatory_features = []
            if motifs_results.get('P53', 0) > 0:
                regulatory_features.append("p53-зависимая регуляция")
            if motifs_results.get('NF-kB', 0) > 0:
                regulatory_features.append("воспалительный ответ")
            if motifs_results.get('CREB', 0) > 0 or motifs_results.get('ATF/CREB', 0) > 0:
                regulatory_features.append("cAMP-зависимая регуляция")
            if motifs_results.get('AP-1', 0) > 0:
                regulatory_features.append("реакция на клеточный стресс")
            
            if regulatory_features:
                f.write(f"5. Паттерн сайтов связывания указывает на возможную {', '.join(regulatory_features)}.\n")
            
            f.write("6. Для более детального понимания функциональной значимости обнаруженных сайтов рекомендуется экспериментальная валидация.\n")
        else:
            f.write(f"В промоторе гена {gene_name} не обнаружены сайты связывания из анализируемого набора транскрипционных факторов.\n")

if __name__ == "__main__":
    analyze_all_genes() 