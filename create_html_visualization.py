#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import base64
from io import BytesIO
import numpy as np
from datetime import datetime
import shutil
from tf_gene_database import get_tf_info, get_gene_info, get_tf_description_html, get_gene_description_html

class HTMLVisualizer:
    def __init__(self, results_dir='results', output_dir='html_report'):
        self.results_dir = results_dir
        self.output_dir = output_dir
        
        # Создаем директорию, если она не существует
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    def load_binding_results(self, filename=None):
        """Загрузка результатов анализа сайтов связывания из JSON файла."""
        if filename is None:
            # Сначала пробуем новое имя файла, затем старое
            if os.path.exists(os.path.join(self.results_dir, 'simplified_binding_results.json')):
                filename = os.path.join(self.results_dir, 'simplified_binding_results.json')
            else:
                filename = os.path.join(self.results_dir, 'binding_analysis_results.json')
        
        if not os.path.exists(filename):
            print(f"Ошибка: Файл результатов {filename} не найден.")
            return None
        
        with open(filename, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    def plot_to_base64(self, plt):
        """Конвертация графика matplotlib в base64 для встраивания в HTML."""
        buffer = BytesIO()
        plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        buffer.seek(0)
        img_str = base64.b64encode(buffer.read()).decode('utf-8')
        plt.close()
        return img_str
    
    def generate_gene_profile(self, gene_name, binding_results):
        """Генерация графика профиля связывания для конкретного гена."""
        # Собираем данные о связывании для этого гена по всем ТФ
        tf_data = []
        
        for tf_id, tf_info in binding_results.items():
            tf_name = tf_info['tf_name']
            
            if gene_name in tf_info['binding_data']:
                binding_data = tf_info['binding_data'][gene_name]
                
                if binding_data['total_sites'] > 0:
                    tf_data.append({
                        'TF_ID': tf_id,
                        'TF_Name': tf_name,
                        'Total_Sites': binding_data['total_sites'],
                        'Avg_Score': binding_data['avg_score'],
                        'Max_Score': binding_data['max_score']
                    })
        
        # Сортировка по количеству сайтов связывания
        tf_data = sorted(tf_data, key=lambda x: x['Total_Sites'], reverse=True)
        
        # Создание графика
        plt.figure(figsize=(10, 6))
        
        if tf_data:
            # Берем топ-20 ТФ
            top_tfs = tf_data[:20]
            
            # Столбчатая диаграмма сайтов связывания
            x = [f"{d['TF_Name']} ({d['TF_ID']})" for d in top_tfs]
            y = [d['Total_Sites'] for d in top_tfs]
            
            plt.bar(x, y, color='skyblue')
            plt.xticks(rotation=90)
            plt.ylabel('Количество сайтов связывания')
            plt.title(f'Основные транскрипционные факторы, связывающиеся с {gene_name}')
            plt.tight_layout()
        else:
            plt.text(0.5, 0.5, "Сайтов связывания не найдено", 
                    horizontalalignment='center', verticalalignment='center')
            plt.title(f'Сайты связывания ТФ для {gene_name} не найдены')
        
        img_str = self.plot_to_base64(plt)
        return img_str, tf_data
    
    def generate_tf_profile(self, tf_id, tf_name, binding_results):
        """Генерация графика профиля связывания для конкретного ТФ."""
        # Проверяем, есть ли этот ТФ в результатах
        if tf_id not in binding_results:
            return None, []
        
        # Собираем данные о связывании этого ТФ по всем генам
        gene_data = []
        
        for gene_name, binding_data in binding_results[tf_id]['binding_data'].items():
            if binding_data['total_sites'] > 0:
                gene_data.append({
                    'Gene': gene_name,
                    'Total_Sites': binding_data['total_sites'],
                    'Avg_Score': binding_data['avg_score'],
                    'Max_Score': binding_data['max_score']
                })
        
        # Сортировка по количеству сайтов связывания
        gene_data = sorted(gene_data, key=lambda x: x['Total_Sites'], reverse=True)
        
        # Создание графика
        plt.figure(figsize=(10, 6))
        
        if gene_data:
            # Берем топ-20 генов
            top_genes = gene_data[:20]
            
            # Столбчатая диаграмма сайтов связывания
            x = [d['Gene'] for d in top_genes]
            y = [d['Total_Sites'] for d in top_genes]
            
            plt.bar(x, y, color='lightgreen')
            plt.xticks(rotation=90)
            plt.ylabel('Количество сайтов связывания')
            plt.title(f'Основные гены, с которыми связывается {tf_name} ({tf_id})')
            plt.tight_layout()
        else:
            plt.text(0.5, 0.5, "Сайтов связывания не найдено", 
                    horizontalalignment='center', verticalalignment='center')
            plt.title(f'Сайты связывания для {tf_name} ({tf_id}) не найдены')
        
        img_str = self.plot_to_base64(plt)
        return img_str, gene_data
    
    def check_interactive_heatmap(self):
        """Проверка наличия интерактивной тепловой карты или ее создание при необходимости."""
        interactive_heatmap_path = os.path.join(self.results_dir, 'interactive_heatmap.html')
        
        if not os.path.exists(interactive_heatmap_path):
            try:
                from interactive_heatmap import InteractiveHeatmap
                heatmap_creator = InteractiveHeatmap(self.results_dir)
                heatmap_creator.create_detailed_heatmap(interactive_heatmap_path)
            except Exception as e:
                print(f"Предупреждение: Не удалось создать интерактивную тепловую карту: {e}")
                return False
        
        # Копирование в выходную директорию
        target_path = os.path.join(self.output_dir, 'interactive_heatmap.html')
        shutil.copy(interactive_heatmap_path, target_path)
        return True
    
    def generate_html_report(self):
        """Генерация полного HTML-отчета по анализу сайтов связывания."""
        binding_results = self.load_binding_results()
        
        if binding_results is None:
            print("Ошибка: Не удалось загрузить результаты анализа.")
            return None
        
        # Пробуем создать/использовать интерактивную тепловую карту
        has_interactive_heatmap = self.check_interactive_heatmap()
        
        # Подготовка аналитических данных
        # Подсчет общего количества сайтов связывания для каждого ТФ
        tf_binding_counts = {}
        # Список генов для каждого ТФ с хотя бы 1 сайтом связывания
        tf_target_genes = {}
        # Общее количество сайтов связывания по всем ТФ для каждого гена
        gene_binding_counts = {}
        # Позиционное распределение сайтов связывания
        position_distribution = {}
        # Самые сильные сайты связывания для каждого ТФ
        strongest_binding_sites = {}
        
        if binding_results:
            # Получаем уникальные имена генов из данных первого ТФ
            first_tf = list(binding_results.values())[0]
            gene_names = list(first_tf['binding_data'].keys())
            
            # Инициализируем счетчики для генов
            for gene_name in gene_names:
                gene_binding_counts[gene_name] = 0
                position_distribution[gene_name] = {}
            
            # Заполняем данные для аналитики
            for tf_id, tf_data in binding_results.items():
                total_sites = 0
                tf_target_genes[tf_id] = []
                strongest_binding_sites[tf_id] = {
                    'gene': '',
                    'score': 0,
                    'sequence': '',
                    'position': 0
                }
                
                for gene_name, binding_sites in tf_data['binding_data'].items():
                    if isinstance(binding_sites, dict) and 'sites' in binding_sites:
                        # Новый формат с детальной информацией о сайтах
                        sites_count = len(binding_sites['sites']) if 'sites' in binding_sites else 0
                        
                        # Находим самый сильный сайт связывания
                        if sites_count > 0:
                            for site in binding_sites['sites']:
                                if 'score' in site and site['score'] > strongest_binding_sites[tf_id]['score']:
                                    strongest_binding_sites[tf_id] = {
                                        'gene': gene_name,
                                        'score': site['score'],
                                        'sequence': site.get('sequence', ''),
                                        'position': site.get('position', 0)
                                    }
                                
                                # Собираем информацию о позициях сайтов связывания
                                if 'position' in site:
                                    pos = site['position']
                                    # Группируем позиции по 100 bp
                                    pos_bin = pos // 100 * 100
                                    if pos_bin not in position_distribution[gene_name]:
                                        position_distribution[gene_name][pos_bin] = 0
                                    position_distribution[gene_name][pos_bin] += 1
                    else:
                        # Старый формат или формат с общим числом сайтов
                        sites_count = binding_sites.get('total_sites', 0)
                        if sites_count == 0 and isinstance(binding_sites, list):
                            sites_count = len(binding_sites)
                    
                    total_sites += sites_count
                    gene_binding_counts[gene_name] += sites_count
                    
                    if sites_count > 0:
                        tf_target_genes[tf_id].append(gene_name)
                
                tf_binding_counts[tf_id] = total_sites
        
        # Сортировка ТФ по количеству сайтов связывания (от большего к меньшему)
        top_tfs = sorted(tf_binding_counts.items(), key=lambda x: x[1], reverse=True)
        # Сортировка генов по количеству сайтов связывания (от большего к меньшему)
        top_genes = sorted(gene_binding_counts.items(), key=lambda x: x[1], reverse=True)
        # ТФ с наибольшим количеством целевых генов
        tf_by_gene_count = sorted(tf_target_genes.items(), key=lambda x: len(x[1]), reverse=True)
        
        # Подсчет распределения ТФ по количеству сайтов связывания
        high_activity_tfs = [tf_id for tf_id, count in tf_binding_counts.items() if count > 50]
        medium_activity_tfs = [tf_id for tf_id, count in tf_binding_counts.items() if 10 <= count <= 50]
        low_activity_tfs = [tf_id for tf_id, count in tf_binding_counts.items() if 1 <= count < 10]
        inactive_tfs = [tf_id for tf_id, count in tf_binding_counts.items() if count == 0]
        
        # Анализ позиционного распределения сайтов связывания
        combined_position_dist = {}
        for gene_positions in position_distribution.values():
            for pos_bin, count in gene_positions.items():
                if pos_bin not in combined_position_dist:
                    combined_position_dist[pos_bin] = 0
                combined_position_dist[pos_bin] += count
        
        # Анализ семейств ТФ
        tf_families = {}
        for tf_id in tf_binding_counts.keys():
            tf_info = get_tf_info(tf_id)
            if tf_info:
                family = tf_info.family
                if family not in tf_families:
                    tf_families[family] = {
                        'count': 0,
                        'total_sites': 0,
                        'tfs': []
                    }
                tf_families[family]['count'] += 1
                tf_families[family]['total_sites'] += tf_binding_counts[tf_id]
                tf_families[family]['tfs'].append(tf_id)
        
        # Сортировка семейств ТФ по количеству сайтов связывания
        top_families = sorted(tf_families.items(), key=lambda x: x[1]['total_sites'], reverse=True)
        
        # Создание index.html
        html = f"""<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Отчет по анализу сайтов связывания транскрипционных факторов</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            color: #333;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .summary {{
            background-color: #f9f9f9;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .card {{
            background-color: white;
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 15px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 20px;
        }}
        th, td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        .nav {{
            background-color: #2c3e50;
            padding: 10px;
            margin-bottom: 20px;
        }}
        .nav a {{
            color: white;
            text-decoration: none;
            margin-right: 15px;
        }}
        img {{
            max-width: 100%;
            height: auto;
        }}
        .button {{
            display: inline-block;
            background-color: #3498db;
            color: white;
            padding: 10px 15px;
            text-decoration: none;
            border-radius: 4px;
            font-weight: bold;
            margin-top: 10px;
        }}
        .button:hover {{
            background-color: #2980b9;
        }}
        .analytics-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 20px;
        }}
        .progress-bar {{
            background-color: #e9ecef;
            border-radius: 0.25rem;
            height: 20px;
            margin-bottom: 5px;
        }}
        .progress-bar-fill {{
            height: 100%;
            border-radius: 0.25rem;
            background-color: #3498db;
        }}
        .family-card {{
            margin-bottom: 15px;
            padding: 10px;
            border-left: 4px solid #3498db;
            background-color: #f8f9fa;
        }}
        .highlight {{
            background-color: #ffffcc;
            font-weight: bold;
        }}
        .small-text {{
            font-size: 0.85em;
            color: #666;
        }}
    </style>
</head>
<body>
    <div class="nav">
        <a href="index.html">Обзор</a>
        <a href="genes.html">Гены</a>
        <a href="tfs.html">Транскрипционные факторы</a>
        {f'<a href="interactive_heatmap.html">Интерактивная тепловая карта</a>' if has_interactive_heatmap else ''}
    </div>
    
    <div class="container">
        <h1>Отчет по анализу сайтов связывания транскрипционных факторов</h1>
        
        <div class="summary">
            <h2>Сводка по анализу</h2>
            <p>Отчет сгенерирован: {datetime.now().strftime("%d.%m.%Y %H:%M:%S")}</p>
            <p>Количество проанализированных транскрипционных факторов: {len(binding_results)}</p>
            <p>Количество проанализированных генов: {len(gene_names)}</p>
            <p>Общее количество обнаруженных сайтов связывания: {sum(tf_binding_counts.values())}</p>
            
            <p>Этот анализ исследует сайты связывания транскрипционных факторов (ТФ) в промоторных областях генов
            от -2000 до +500 относительно сайта начала транскрипции (TSS).</p>
            
            {f'<a href="interactive_heatmap.html" class="button">Просмотреть интерактивную тепловую карту</a>' if has_interactive_heatmap else ''}
        </div>
        
        <div class="card">
            <h2>Общая аналитика</h2>
            <div class="analytics-grid">
                <div>
                    <h3>Топ-5 ТФ по количеству сайтов связывания</h3>
                    <table>
                        <tr>
                            <th>Транскрипционный фактор</th>
                            <th>Количество сайтов</th>
                        </tr>
"""

        # Добавляем топ-5 ТФ по количеству сайтов связывания
        max_sites = top_tfs[0][1] if top_tfs else 0
        for i, (tf_id, count) in enumerate(top_tfs[:5]):
            if count > 0:
                tf_info = get_tf_info(tf_id)
                tf_name = tf_info.name if tf_info else tf_id
                percent = (count / max_sites) * 100 if max_sites > 0 else 0
                
                html += f"""
                        <tr>
                            <td>{tf_name}</td>
                            <td>
                                <div class="progress-bar">
                                    <div class="progress-bar-fill" style="width: {percent}%;"></div>
                                </div>
                                {count}
                            </td>
                        </tr>"""
                
        html += """
                    </table>
                </div>
                <div>
                    <h3>Топ-5 генов по количеству сайтов связывания</h3>
                    <table>
                        <tr>
                            <th>Ген</th>
                            <th>Количество сайтов</th>
                        </tr>
"""

        # Добавляем топ-5 генов по количеству сайтов связывания
        max_gene_sites = top_genes[0][1] if top_genes else 0
        for i, (gene_name, count) in enumerate(top_genes[:5]):
            if count > 0:
                percent = (count / max_gene_sites) * 100 if max_gene_sites > 0 else 0
                
                html += f"""
                        <tr>
                            <td>{gene_name}</td>
                            <td>
                                <div class="progress-bar">
                                    <div class="progress-bar-fill" style="width: {percent}%;"></div>
                                </div>
                                {count}
                            </td>
                        </tr>"""
                
        html += """
                    </table>
                </div>
                <div>
                    <h3>Топ-5 ТФ по количеству регулируемых генов</h3>
                    <table>
                        <tr>
                            <th>Транскрипционный фактор</th>
                            <th>Количество генов</th>
                        </tr>
"""

        # Добавляем топ-5 ТФ по количеству регулируемых генов
        max_gene_count = len(tf_by_gene_count[0][1]) if tf_by_gene_count else 0
        for i, (tf_id, genes) in enumerate(tf_by_gene_count[:5]):
            if len(genes) > 0:
                tf_info = get_tf_info(tf_id)
                tf_name = tf_info.name if tf_info else tf_id
                percent = (len(genes) / max_gene_count) * 100 if max_gene_count > 0 else 0
                
                html += f"""
                        <tr>
                            <td>{tf_name}</td>
                            <td>
                                <div class="progress-bar">
                                    <div class="progress-bar-fill" style="width: {percent}%;"></div>
                                </div>
                                {len(genes)}
                            </td>
                        </tr>"""
                
        html += """
                    </table>
                </div>
                <div>
                    <h3>Распределение сайтов связывания по ТФ</h3>
                    <table>
                        <tr>
                            <th>Категория</th>
                            <th>Количество ТФ</th>
                        </tr>
"""

        # Отображаем распределение ТФ по категориям активности
        total_tfs = len(tf_binding_counts)
        
        html += f"""
                        <tr>
                            <td>Высокая активность (>50 сайтов)</td>
                            <td>
                                <div class="progress-bar">
                                    <div class="progress-bar-fill" style="width: {len(high_activity_tfs)/total_tfs*100 if total_tfs else 0}%;"></div>
                                </div>
                                {len(high_activity_tfs)}
                            </td>
                        </tr>
                        <tr>
                            <td>Средняя активность (10-50 сайтов)</td>
                            <td>
                                <div class="progress-bar">
                                    <div class="progress-bar-fill" style="width: {len(medium_activity_tfs)/total_tfs*100 if total_tfs else 0}%;"></div>
                                </div>
                                {len(medium_activity_tfs)}
                            </td>
                        </tr>
                        <tr>
                            <td>Низкая активность (1-9 сайтов)</td>
                            <td>
                                <div class="progress-bar">
                                    <div class="progress-bar-fill" style="width: {len(low_activity_tfs)/total_tfs*100 if total_tfs else 0}%;"></div>
                                </div>
                                {len(low_activity_tfs)}
                            </td>
                        </tr>
                        <tr>
                            <td>Неактивные (0 сайтов)</td>
                            <td>
                                <div class="progress-bar">
                                    <div class="progress-bar-fill" style="width: {len(inactive_tfs)/total_tfs*100 if total_tfs else 0}%;"></div>
                                </div>
                                {len(inactive_tfs)}
                            </td>
                        </tr>
                    </table>
                </div>
            </div>
        </div>
        
        <div class="card">
            <h2>Самые сильные сайты связывания</h2>
            <table>
                <tr>
                    <th>Транскрипционный фактор</th>
                    <th>Ген</th>
                    <th>Последовательность</th>
                    <th>Позиция</th>
                    <th>Оценка</th>
                </tr>
"""

        # Добавляем информацию о самых сильных сайтах связывания
        top_binding_sites = sorted([(tf_id, data) for tf_id, data in strongest_binding_sites.items() if data['score'] > 0], 
                                 key=lambda x: x[1]['score'], reverse=True)[:10]
        
        for tf_id, site_data in top_binding_sites:
            tf_info = get_tf_info(tf_id)
            tf_name = tf_info.name if tf_info else tf_id
            
            html += f"""
                <tr>
                    <td>{tf_name}</td>
                    <td>{site_data['gene']}</td>
                    <td><span class="highlight">{site_data['sequence']}</span></td>
                    <td>{site_data['position']}</td>
                    <td>{site_data['score']:.3f}</td>
                </tr>"""
        
        html += """
            </table>
        </div>
        
        <div class="card">
            <h2>Распределение ТФ по семействам</h2>
"""

        # Добавляем информацию о семействах ТФ
        for family, family_data in top_families[:10]:
            tf_count = family_data['count']
            total_sites = family_data['total_sites']
            
            html += f"""
            <div class="family-card">
                <h3>{family} ({tf_count} ТФ)</h3>
                <p>Общее количество сайтов связывания: {total_sites}</p>
                <p>Основные ТФ: """
            
            # Добавляем основные ТФ этого семейства
            family_tfs = []
            for tf_id in family_data['tfs']:
                tf_info = get_tf_info(tf_id)
                tf_name = tf_info.name if tf_info else tf_id
                family_tfs.append(f'{tf_name} ({tf_binding_counts[tf_id]} сайтов)')
            
            html += ", ".join(family_tfs[:5])
            if len(family_tfs) > 5:
                html += f" и {len(family_tfs) - 5} других"
            
            html += """</p>
            </div>"""
        
        html += """
        </div>
        
        <div class="card">
            <h2>Общее распределение сайтов связывания</h2>
            <img src="binding_site_distribution.png" alt="Распределение сайтов связывания">
        </div>
        
        <div class="card">
            <h2>Тепловая карта сайтов связывания ТФ по генам</h2>
            <img src="binding_heatmap.png" alt="Тепловая карта ТФ">
            {f'<br><a href="interactive_heatmap.html" class="button">Просмотреть интерактивную версию</a>' if has_interactive_heatmap else ''}
        </div>
    </div>
</body>
</html>"""
        
        with open(os.path.join(self.output_dir, 'index.html'), 'w', encoding='utf-8') as f:
            f.write(html)
        
        # Копирование графиков в выходную директорию
        shutil.copy(os.path.join(self.results_dir, 'binding_site_distribution.png'), 
                   os.path.join(self.output_dir, 'binding_site_distribution.png'))
        shutil.copy(os.path.join(self.results_dir, 'binding_heatmap.png'), 
                   os.path.join(self.output_dir, 'binding_heatmap.png'))
        
        # Создание genes.html
        genes_html = f"""<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Анализ генов | Отчет по сайтам связывания ТФ</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            color: #333;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .summary {{
            background-color: #f9f9f9;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .card {{
            background-color: white;
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 15px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 20px;
        }}
        th, td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        .nav {{
            background-color: #2c3e50;
            padding: 10px;
            margin-bottom: 20px;
        }}
        .nav a {{
            color: white;
            text-decoration: none;
            margin-right: 15px;
        }}
        .gene-profile {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-bottom: 20px;
        }}
        .search-container {{
            margin-bottom: 20px;
        }}
        #geneSearch {{
            padding: 10px;
            width: 300px;
            border: 1px solid #ddd;
            border-radius: 4px;
        }}
        .progress-bar {{
            background-color: #e9ecef;
            border-radius: 0.25rem;
            height: 20px;
            margin-bottom: 5px;
        }}
        .progress-bar-fill {{
            height: 100%;
            border-radius: 0.25rem;
            background-color: #3498db;
        }}
        .highlight {{
            background-color: #ffffcc;
            font-weight: bold;
        }}
        .info-box {{
            background-color: #f8f9fa;
            border-left: 4px solid #3498db;
            padding: 10px;
            margin-bottom: 15px;
        }}
        .site-distributions {{
            display: flex;
            margin-top: 10px;
        }}
        .site-distributions > div {{
            flex: 1;
            padding: 10px;
        }}
    </style>
</head>
<body>
    <div class="nav">
        <a href="index.html">Обзор</a>
        <a href="genes.html">Гены</a>
        <a href="tfs.html">Транскрипционные факторы</a>
        {f'<a href="interactive_heatmap.html">Интерактивная тепловая карта</a>' if has_interactive_heatmap else ''}
    </div>
    
    <div class="container">
        <h1>Анализ генов</h1>
        
        <div class="search-container">
            <input type="text" id="geneSearch" placeholder="Введите название гена для поиска...">
        </div>
        
        <div class="summary">
            <h2>Сводка по генам</h2>
            <p>Количество проанализированных генов: {len(gene_names)}</p>
            <p>Общее количество обнаруженных сайтов связывания: {sum(gene_binding_counts.values())}</p>
            <p>Средний промотор гена содержит {sum(gene_binding_counts.values()) / len(gene_names):.1f} сайтов связывания</p>
        </div>
        
        <div class="card">
            <h2>Топ гены по количеству сайтов связывания</h2>
            <table>
                <tr>
                    <th>Позиция</th>
                    <th>Ген</th>
                    <th>Количество сайтов</th>
                    <th>Количество ТФ</th>
                    <th>Основные ТФ</th>
                </tr>
"""

        # Добавляем топ-20 генов по количеству сайтов связывания
        for i, (gene_name, count) in enumerate(top_genes[:20]):
            # Подсчет количества ТФ, имеющих хотя бы 1 сайт связывания
            tf_count = 0
            main_tfs = []
            
            # Собираем данные о ТФ для этого гена
            gene_tf_sites = {}
            for tf_id, tf_data in binding_results.items():
                binding_sites = tf_data['binding_data'].get(gene_name, {})
                sites_count = 0
                
                if isinstance(binding_sites, dict):
                    if 'sites' in binding_sites:
                        sites_count = len(binding_sites['sites'])
                    elif 'total_sites' in binding_sites:
                        sites_count = binding_sites['total_sites']
                
                if sites_count > 0:
                    tf_count += 1
                    gene_tf_sites[tf_id] = sites_count
            
            # Сортируем ТФ по количеству сайтов связывания
            sorted_tfs = sorted(gene_tf_sites.items(), key=lambda x: x[1], reverse=True)
            for tf_id, sites in sorted_tfs[:3]:
                tf_info = get_tf_info(tf_id)
                tf_name = tf_info.name if tf_info else tf_id
                main_tfs.append(f"{tf_name} ({sites})")
            
            genes_html += f"""
                <tr>
                    <td>{i+1}</td>
                    <td><a href="#gene-{gene_name}">{gene_name}</a></td>
                    <td>{count}</td>
                    <td>{tf_count}</td>
                    <td>{', '.join(main_tfs) if main_tfs else 'Нет данных'}</td>
                </tr>"""
        
        genes_html += """
            </table>
        </div>
        
        <h2>Профили генов</h2>
"""
        
        # Добавляем профили для каждого гена
        for gene_name in gene_names:
            gene_info = get_gene_info(gene_name)
            gene_desc = gene_info.description if gene_info else "Нет информации"
            gene_func = gene_info.function if gene_info else "Нет информации"
            
            # Собираем информацию о ТФ, связывающихся с этим геном
            gene_tfs = []
            for tf_id, tf_data in binding_results.items():
                binding_sites = tf_data['binding_data'].get(gene_name, {})
                sites_count = 0
                
                if isinstance(binding_sites, dict):
                    if 'sites' in binding_sites:
                        sites_count = len(binding_sites['sites'])
                    elif 'total_sites' in binding_sites:
                        sites_count = binding_sites['total_sites']
                
                if sites_count > 0:
                    gene_tfs.append((tf_id, sites_count))
            
            # Сортируем ТФ по количеству сайтов связывания
            gene_tfs.sort(key=lambda x: x[1], reverse=True)
            
            genes_html += f"""
        <div id="gene-{gene_name}" class="card gene-card">
            <h3>{gene_name}</h3>
            <div class="gene-profile">
                <div>
                    <div class="info-box">
                        <h4>Описание</h4>
                        <p>{gene_desc}</p>
                        <h4>Функция</h4>
                        <p>{gene_func}</p>
                    </div>
                    <h4>Статистика сайтов связывания</h4>
                    <p>Общее количество сайтов: <strong>{gene_binding_counts.get(gene_name, 0)}</strong></p>
                    <p>Количество транскрипционных факторов: <strong>{len(gene_tfs)}</strong></p>
                </div>
                <div>
                    <h4>Транскрипционные факторы</h4>
                    <table>
                        <tr>
                            <th>ТФ</th>
                            <th>Количество сайтов</th>
                        </tr>
"""
            
            # Добавляем топ-10 ТФ для этого гена
            max_sites = gene_tfs[0][1] if gene_tfs else 0
            for tf_id, sites_count in gene_tfs[:10]:
                tf_info = get_tf_info(tf_id)
                tf_name = tf_info.name if tf_info else tf_id
                
                percent = (sites_count / max_sites) * 100 if max_sites > 0 else 0
                
                genes_html += f"""
                        <tr>
                            <td>{tf_name}</td>
                            <td>
                                <div class="progress-bar">
                                    <div class="progress-bar-fill" style="width: {percent}%;"></div>
                                </div>
                                {sites_count}
                            </td>
                        </tr>"""
            
            genes_html += """
                    </table>
                </div>
            </div>
        </div>
"""
        
        genes_html += """
        <script>
            // Скрипт для поиска генов
            document.getElementById('geneSearch').addEventListener('keyup', function() {
                const searchTerm = this.value.toLowerCase();
                const geneCards = document.querySelectorAll('.gene-card');
                
                geneCards.forEach(card => {
                    const geneName = card.querySelector('h3').textContent.toLowerCase();
                    if (geneName.includes(searchTerm)) {
                        card.style.display = '';
                    } else {
                        card.style.display = 'none';
                    }
                });
            });
        </script>
    </div>
</body>
</html>"""
        
        with open(os.path.join(self.output_dir, 'genes.html'), 'w', encoding='utf-8') as f:
            f.write(genes_html)
        
        # Создание tfs.html
        tfs_html = f"""<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Анализ транскрипционных факторов | Отчет по сайтам связывания ТФ</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            color: #333;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .summary {{
            background-color: #f9f9f9;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .card {{
            background-color: white;
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 15px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 20px;
        }}
        th, td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        .nav {{
            background-color: #2c3e50;
            padding: 10px;
            margin-bottom: 20px;
        }}
        .nav a {{
            color: white;
            text-decoration: none;
            margin-right: 15px;
        }}
        .search-container {{
            margin-bottom: 20px;
        }}
        #tfSearch {{
            padding: 10px;
            width: 300px;
            border: 1px solid #ddd;
            border-radius: 4px;
        }}
        .progress-bar {{
            background-color: #e9ecef;
            border-radius: 0.25rem;
            height: 20px;
            margin-bottom: 5px;
        }}
        .progress-bar-fill {{
            height: 100%;
            border-radius: 0.25rem;
            background-color: #3498db;
        }}
        .tf-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 20px;
        }}
        .highlight {{
            background-color: #ffffcc;
            font-weight: bold;
        }}
        .info-box {{
            background-color: #f8f9fa;
            border-left: 4px solid #3498db;
            padding: 10px;
            margin-bottom: 15px;
        }}
        .filter-group {{
            margin-bottom: 15px;
        }}
        .filter-btn {{
            padding: 5px 10px;
            margin-right: 5px;
            border: none;
            border-radius: 3px;
            cursor: pointer;
            background-color: #f0f0f0;
        }}
        .filter-btn.active {{
            background-color: #3498db;
            color: white;
        }}
        .status-indicator {{
            display: inline-block;
            width: 12px;
            height: 12px;
            border-radius: 50%;
            margin-right: 5px;
        }}
        .status-high {{
            background-color: #27ae60;
        }}
        .status-medium {{
            background-color: #f39c12;
        }}
        .status-low {{
            background-color: #3498db;
        }}
        .status-inactive {{
            background-color: #e74c3c;
        }}
        .tf-card {{
            border-left: 4px solid #3498db;
        }}
        .gene-list {{
            height: 200px;
            overflow-y: auto;
            border: 1px solid #ddd;
            padding: 10px;
            margin-top: 10px;
        }}
    </style>
</head>
<body>
    <div class="nav">
        <a href="index.html">Обзор</a>
        <a href="genes.html">Гены</a>
        <a href="tfs.html">Транскрипционные факторы</a>
        {f'<a href="interactive_heatmap.html">Интерактивная тепловая карта</a>' if has_interactive_heatmap else ''}
    </div>
    
    <div class="container">
        <h1>Анализ транскрипционных факторов</h1>
        
        <div class="search-container">
            <input type="text" id="tfSearch" placeholder="Введите название ТФ для поиска...">
        </div>
        
        <div class="filter-group">
            <span>Фильтр по активности: </span>
            <button class="filter-btn active" data-filter="all">Все</button>
            <button class="filter-btn" data-filter="high"><span class="status-indicator status-high"></span> Высокая (>50)</button>
            <button class="filter-btn" data-filter="medium"><span class="status-indicator status-medium"></span> Средняя (10-50)</button>
            <button class="filter-btn" data-filter="low"><span class="status-indicator status-low"></span> Низкая (1-9)</button>
            <button class="filter-btn" data-filter="inactive"><span class="status-indicator status-inactive"></span> Неактивные</button>
        </div>
        
        <div class="summary">
            <h2>Сводка по транскрипционным факторам</h2>
            <p>Количество проанализированных ТФ: {len(binding_results)}</p>
            <p>Общее количество обнаруженных сайтов связывания: {sum(tf_binding_counts.values())}</p>
            <p>Распределение ТФ по уровню активности:</p>
            <ul>
                <li><span class="status-indicator status-high"></span> Высокая активность (>50 сайтов): {len(high_activity_tfs)}</li>
                <li><span class="status-indicator status-medium"></span> Средняя активность (10-50 сайтов): {len(medium_activity_tfs)}</li>
                <li><span class="status-indicator status-low"></span> Низкая активность (1-9 сайтов): {len(low_activity_tfs)}</li>
                <li><span class="status-indicator status-inactive"></span> Неактивные (0 сайтов): {len(inactive_tfs)}</li>
            </ul>
        </div>
        
        <div class="card">
            <h2>Топ транскрипционные факторы</h2>
            <table>
                <tr>
                    <th>Позиция</th>
                    <th>ТФ</th>
                    <th>Семейство</th>
                    <th>Количество сайтов</th>
                    <th>Количество генов</th>
                    <th>Статус</th>
                </tr>
"""

        # Добавляем топ-20 ТФ по количеству сайтов связывания
        for i, (tf_id, count) in enumerate(top_tfs[:20]):
            tf_info = get_tf_info(tf_id)
            tf_name = tf_info.name if tf_info else tf_id
            tf_family = tf_info.family if tf_info else "Неизвестно"
            
            # Определяем статус активности
            status_class = ""
            status_text = ""
            if tf_id in high_activity_tfs:
                status_class = "status-high"
                status_text = "Высокая"
            elif tf_id in medium_activity_tfs:
                status_class = "status-medium"
                status_text = "Средняя"
            elif tf_id in low_activity_tfs:
                status_class = "status-low"
                status_text = "Низкая"
            else:
                status_class = "status-inactive"
                status_text = "Неактивный"
            
            tfs_html += f"""
                <tr class="tf-row" data-status="{status_text.lower()}">
                    <td>{i+1}</td>
                    <td><a href="#tf-{tf_id}">{tf_name}</a></td>
                    <td>{tf_family}</td>
                    <td>{count}</td>
                    <td>{len(tf_target_genes.get(tf_id, []))}</td>
                    <td><span class="status-indicator {status_class}"></span> {status_text}</td>
                </tr>"""
        
        tfs_html += """
            </table>
        </div>
        
        <h2>Профили транскрипционных факторов</h2>
        <div class="tf-grid">
"""
        
        # Добавляем профили для каждого ТФ
        for tf_id, count in tf_binding_counts.items():
            tf_info = get_tf_info(tf_id)
            tf_name = tf_info.name if tf_info else tf_id
            tf_desc = tf_info.description if tf_info else "Нет информации"
            tf_func = tf_info.function if tf_info else "Нет информации"
            tf_family = tf_info.family if tf_info else "Неизвестно"
            
            # Определяем статус активности
            status_class = ""
            status_text = ""
            data_status = ""
            if tf_id in high_activity_tfs:
                status_class = "status-high"
                status_text = "Высокая"
                data_status = "high"
            elif tf_id in medium_activity_tfs:
                status_class = "status-medium"
                status_text = "Средняя"
                data_status = "medium"
            elif tf_id in low_activity_tfs:
                status_class = "status-low"
                status_text = "Низкая"
                data_status = "low"
            else:
                status_class = "status-inactive"
                status_text = "Неактивный"
                data_status = "inactive"
            
            # Собираем список целевых генов
            target_genes = tf_target_genes.get(tf_id, [])
            
            tfs_html += f"""
            <div id="tf-{tf_id}" class="card tf-card" data-status="{data_status}">
                <h3>{tf_name} <span class="tf-id">({tf_id})</span></h3>
                <p><span class="status-indicator {status_class}"></span> <strong>Статус активности:</strong> {status_text}</p>
                <div class="info-box">
                    <p><strong>Семейство:</strong> {tf_family}</p>
                    <p><strong>Описание:</strong> {tf_desc}</p>
                    <p><strong>Функция:</strong> {tf_func}</p>
                </div>
                <p><strong>Количество сайтов связывания:</strong> {count}</p>
                <p><strong>Количество регулируемых генов:</strong> {len(target_genes)}</p>
"""
            
            # Добавляем информацию о наиболее сильном сайте связывания
            strongest_site = strongest_binding_sites.get(tf_id, {})
            if strongest_site.get('score', 0) > 0:
                tfs_html += f"""
                <div class="info-box">
                    <h4>Наиболее сильный сайт связывания</h4>
                    <p><strong>Ген:</strong> {strongest_site['gene']}</p>
                    <p><strong>Последовательность:</strong> <span class="highlight">{strongest_site['sequence']}</span></p>
                    <p><strong>Позиция:</strong> {strongest_site['position']}</p>
                    <p><strong>Оценка:</strong> {strongest_site['score']:.3f}</p>
                </div>
"""
            
            # Добавляем список целевых генов
            if target_genes:
                tfs_html += f"""
                <h4>Целевые гены ({len(target_genes)})</h4>
                <div class="gene-list">
"""
                
                # Сортируем гены по количеству сайтов связывания
                gene_counts = []
                for gene_name in target_genes:
                    binding_sites = binding_results[tf_id]['binding_data'].get(gene_name, {})
                    sites_count = 0
                    
                    if isinstance(binding_sites, dict):
                        if 'sites' in binding_sites:
                            sites_count = len(binding_sites['sites'])
                        elif 'total_sites' in binding_sites:
                            sites_count = binding_sites['total_sites']
                    
                    gene_counts.append((gene_name, sites_count))
                
                gene_counts.sort(key=lambda x: x[1], reverse=True)
                
                for gene_name, sites_count in gene_counts:
                    tfs_html += f"""
                    <p><strong>{gene_name}</strong> - {sites_count} сайтов</p>
"""
                
                tfs_html += """
                </div>
"""
            
            tfs_html += """
            </div>
"""
        
        tfs_html += """
        </div>
        
        <script>
            // Скрипт для поиска ТФ
            document.getElementById('tfSearch').addEventListener('keyup', function() {
                const searchTerm = this.value.toLowerCase();
                const tfCards = document.querySelectorAll('.tf-card');
                
                tfCards.forEach(card => {
                    const tfName = card.querySelector('h3').textContent.toLowerCase();
                    if (tfName.includes(searchTerm)) {
                        card.style.display = '';
                    } else {
                        card.style.display = 'none';
                    }
                });
            });
            
            // Скрипт для фильтрации по активности
            document.querySelectorAll('.filter-btn').forEach(button => {
                button.addEventListener('click', function() {
                    // Убираем класс active со всех кнопок
                    document.querySelectorAll('.filter-btn').forEach(btn => {
                        btn.classList.remove('active');
                    });
                    
                    // Добавляем класс active нажатой кнопке
                    this.classList.add('active');
                    
                    // Фильтруем карточки ТФ
                    const filter = this.getAttribute('data-filter');
                    const tfCards = document.querySelectorAll('.tf-card');
                    
                    tfCards.forEach(card => {
                        if (filter === 'all' || card.getAttribute('data-status') === filter) {
                            card.style.display = '';
                        } else {
                            card.style.display = 'none';
                        }
                    });
                });
            });
        </script>
    </div>
</body>
</html>"""
        
        with open(os.path.join(self.output_dir, 'tfs.html'), 'w', encoding='utf-8') as f:
            f.write(tfs_html)
            
        # Копируем интерактивную тепловую карту, если она доступна
        if has_interactive_heatmap:
            shutil.copy(os.path.join(self.results_dir, 'interactive_heatmap.html'), 
                       os.path.join(self.output_dir, 'interactive_heatmap.html'))
            
        return os.path.join(self.output_dir, 'index.html')

if __name__ == "__main__":
    visualizer = HTMLVisualizer()
    visualizer.generate_html_report() 