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
            
            <p>Этот анализ исследует сайты связывания транскрипционных факторов (ТФ) в промоторных областях генов
            от -2000 до +500 относительно сайта начала транскрипции (TSS).</p>
            
            {f'<a href="interactive_heatmap.html" class="button">Просмотреть интерактивную тепловую карту</a>' if has_interactive_heatmap else ''}
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
    <title>Профили генов - Анализ сайтов связывания ТФ</title>
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
        .gene-info {{
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 15px;
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
        <h1>Профили генов</h1>
"""
        
        # Получаем уникальные имена генов из данных первого ТФ
        if binding_results:
            first_tf = list(binding_results.values())[0]
            gene_names = list(first_tf['binding_data'].keys())
            
            for gene_name in gene_names:
                img_str, gene_tf_data = self.generate_gene_profile(gene_name, binding_results)
                
                # Получаем информацию о гене из нашей базы знаний
                gene_info = get_gene_info(gene_name)
                
                genes_html += f"""
        <div class="card">
            <h2>{gene_name}</h2>
"""
                
                if gene_info:
                    genes_html += f"""
            <div class="gene-info">
                <p><strong>Полное название:</strong> {gene_info.description}</p>
                <p><strong>Функция:</strong> {gene_info.function}</p>
                <p><strong>Сигнальный путь:</strong> {gene_info.pathway}</p>
            </div>
"""
                
                genes_html += f"""
            <img src="data:image/png;base64,{img_str}" alt="Профиль {gene_name}">
            
            <h3>Топ-10 ТФ, связывающихся с {gene_name}</h3>
            <table>
                <tr>
                    <th>ID ТФ</th>
                    <th>Название ТФ</th>
                    <th>Сайты связывания</th>
                    <th>Средний скор</th>
                    <th>Макс. скор</th>
                </tr>
"""
                
                for tf in gene_tf_data[:10]:  # Показываем топ-10
                    genes_html += f"""
                <tr>
                    <td>{tf['TF_ID']}</td>
                    <td>{tf['TF_Name']}</td>
                    <td>{tf['Total_Sites']}</td>
                    <td>{tf['Avg_Score']:.3f}</td>
                    <td>{tf['Max_Score']:.3f}</td>
                </tr>
"""
                
                genes_html += """
            </table>
        </div>
"""
        
        genes_html += """
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
    <title>Профили транскрипционных факторов - Анализ сайтов связывания ТФ</title>
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
        .tf-info {{
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 15px;
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
        <h1>Профили транскрипционных факторов</h1>
"""
        
        # Генерируем профили для ТФ с сайтами связывания
        for tf_id, tf_info in binding_results.items():
            tf_name = tf_info['tf_name']
            
            img_str, tf_gene_data = self.generate_tf_profile(tf_id, tf_name, binding_results)
            
            # Получаем информацию о ТФ из нашей базы знаний
            tf_detailed_info = get_tf_info(tf_id)
            
            if img_str:
                tfs_html += f"""
        <div class="card">
            <h2>{tf_name} ({tf_id})</h2>
            
            <div class="tf-info">
"""
                
                if tf_detailed_info:
                    tfs_html += f"""
                <p><strong>Полное название:</strong> {tf_detailed_info.description}</p>
                <p><strong>Семейство:</strong> {tf_detailed_info.family}</p>
                <p><strong>Функция:</strong> {tf_detailed_info.function}</p>
"""
                else:
                    tfs_html += f"""
                <p><strong>ID матрицы:</strong> {tf_id}</p>
                <p><strong>Название:</strong> {tf_name}</p>
"""
                
                tfs_html += f"""
            </div>
            
            <img src="data:image/png;base64,{img_str}" alt="Профиль {tf_name}">
            
            <h3>Топ-10 генов, связываемых {tf_name}</h3>
            <table>
                <tr>
                    <th>Ген</th>
                    <th>Сайты связывания</th>
                    <th>Средний скор</th>
                    <th>Макс. скор</th>
                    <th>Описание гена</th>
                </tr>
"""
                
                for gene in tf_gene_data[:10]:  # Показываем топ-10
                    gene_name = gene['Gene']
                    gene_detailed_info = get_gene_info(gene_name)
                    gene_desc = gene_detailed_info.description if gene_detailed_info else ""
                    
                    tfs_html += f"""
                <tr>
                    <td>{gene_name}</td>
                    <td>{gene['Total_Sites']}</td>
                    <td>{gene['Avg_Score']:.3f}</td>
                    <td>{gene['Max_Score']:.3f}</td>
                    <td>{gene_desc}</td>
                </tr>
"""
                
                tfs_html += """
            </table>
        </div>
"""
        
        tfs_html += """
    </div>
</body>
</html>"""
        
        with open(os.path.join(self.output_dir, 'tfs.html'), 'w', encoding='utf-8') as f:
            f.write(tfs_html)
        
        print(f"Сгенерирован HTML-отчет в директории {self.output_dir}")
        return self.output_dir

if __name__ == "__main__":
    visualizer = HTMLVisualizer()
    visualizer.generate_html_report() 