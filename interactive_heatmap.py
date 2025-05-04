#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import os
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from tf_gene_database import get_tf_info, get_gene_info, get_tf_description_html, get_gene_description_html, TF_DATABASE, GENE_DATABASE

class InteractiveHeatmap:
    def __init__(self, results_dir='results'):
        self.results_dir = results_dir
        
    def load_binding_results(self, filename=None):
        """Загрузка результатов анализа сайтов связывания из JSON файла."""
        if filename is None:
            # Сначала попробуем новое имя файла, затем старое
            if os.path.exists(os.path.join(self.results_dir, 'simplified_binding_results.json')):
                filename = os.path.join(self.results_dir, 'simplified_binding_results.json')
            else:
                filename = os.path.join(self.results_dir, 'binding_analysis_results.json')
        
        if not os.path.exists(filename):
            print(f"Ошибка: Файл результатов {filename} не найден.")
            return None
        
        with open(filename, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    def create_summary_dataframe(self, binding_results):
        """Создание сводной таблицы из результатов анализа."""
        summary_data = []
        
        for tf_id, tf_data in binding_results.items():
            tf_name = tf_data['tf_name']
            
            for gene_name, binding_data in tf_data['binding_data'].items():
                summary_data.append({
                    'TF_ID': tf_id,
                    'TF_Name': tf_name,
                    'Gene': gene_name,
                    'Binding_Sites': binding_data['total_sites'],
                    'Avg_Score': binding_data.get('avg_score', 0),
                    'Max_Score': binding_data.get('max_score', 0)
                })
        
        return pd.DataFrame(summary_data)
    
    def create_interactive_heatmap(self, output_file=None):
        """Генерация интерактивной тепловой карты с использованием Plotly."""
        binding_results = self.load_binding_results()
        
        if binding_results is None:
            print("Результаты анализа не найдены.")
            return None
        
        # Создание сводной таблицы
        df = self.create_summary_dataframe(binding_results)
        
        if df.empty:
            print("Нет данных для тепловой карты.")
            return None
        
        # Создание матрицы для тепловой карты: гены vs ТФ
        heatmap_data = df.pivot_table(
            values='Binding_Sites', 
            index='Gene', 
            columns='TF_Name', 
            fill_value=0
        )
        
        # Создание интерактивной тепловой карты
        fig = go.Figure(data=go.Heatmap(
            z=heatmap_data.values,
            x=heatmap_data.columns,
            y=heatmap_data.index,
            colorscale='Viridis',
            text=[[f"Ген: {gene}<br>ТФ: {tf}<br>Сайты связывания: {value}" 
                  for tf, value in zip(heatmap_data.columns, row)]
                 for gene, row in zip(heatmap_data.index, heatmap_data.values)],
            hoverinfo='text',
            colorbar=dict(title='Сайты связывания')
        ))
        
        fig.update_layout(
            title='Сайты связывания транскрипционных факторов по генам',
            xaxis=dict(title='Транскрипционный фактор'),
            yaxis=dict(title='Ген'),
            hoverlabel=dict(
                bgcolor="white",
                font_size=12,
                font_family="Arial"
            )
        )
        
        # Сохранение в виде интерактивного HTML
        if output_file is None:
            output_file = os.path.join(self.results_dir, 'interactive_heatmap.html')
        
        fig.write_html(
            output_file,
            full_html=True,
            include_plotlyjs='cdn',
        )
        
        print(f"Интерактивная тепловая карта сохранена в {output_file}")
        return output_file
    
    def create_detailed_heatmap(self, output_file=None):
        """Генерация детальной интерактивной тепловой карты с подсказками и кликабельными ячейками."""
        binding_results = self.load_binding_results()
        
        if binding_results is None:
            print("Результаты анализа не найдены.")
            return None
        
        # Создание сводной таблицы
        df = self.create_summary_dataframe(binding_results)
        
        if df.empty:
            print("Нет данных для тепловой карты.")
            return None
        
        # Создание матрицы для тепловой карты
        heatmap_data = df.pivot_table(
            values='Binding_Sites', 
            index='Gene', 
            columns='TF_Name', 
            fill_value=0
        )
        
        # Получение максимальных значений для масштабирования цвета
        max_scores = df.pivot_table(
            values='Max_Score', 
            index='Gene', 
            columns='TF_Name', 
            fill_value=0
        )
        
        # Получение средних значений для дополнительной информации
        avg_scores = df.pivot_table(
            values='Avg_Score', 
            index='Gene', 
            columns='TF_Name', 
            fill_value=0
        )
        
        # Создание текста для подсказок с подробной информацией
        hovertext = []
        for i, gene in enumerate(heatmap_data.index):
            row_text = []
            for j, tf in enumerate(heatmap_data.columns):
                tf_id = df[(df['Gene'] == gene) & (df['TF_Name'] == tf)]['TF_ID'].values[0] if not df[(df['Gene'] == gene) & (df['TF_Name'] == tf)].empty else "Неизвестно"
                sites = heatmap_data.iloc[i, j]
                avg = avg_scores.iloc[i, j]
                max_score = max_scores.iloc[i, j]
                
                # Получаем информацию из базы данных
                tf_info = get_tf_info(tf_id)
                gene_info = get_gene_info(gene)
                
                text = f"<b>Ген:</b> {gene}<br>"
                if gene_info:
                    text += f"<b>Описание гена:</b> {gene_info.description}<br>"
                    text += f"<b>Путь:</b> {gene_info.pathway}<br>"
                
                text += f"<b>ТФ:</b> {tf} ({tf_id})<br>"
                if tf_info:
                    text += f"<b>Семейство ТФ:</b> {tf_info.family}<br>"
                
                text += f"<b>Сайты связывания:</b> {sites}<br>"
                text += f"<b>Средний скор:</b> {avg:.3f}<br>"
                text += f"<b>Макс. скор:</b> {max_score:.3f}"
                
                # Добавление позиций сайтов связывания, если они доступны и есть сайты
                if sites > 0:
                    try:
                        site_positions = [site['position'] for site in binding_results[tf_id]['binding_data'][gene]['sites']]
                        if len(site_positions) > 10:
                            pos_str = ", ".join(str(pos) for pos in site_positions[:10]) + f"... (и ещё {len(site_positions)-10})"
                        else:
                            pos_str = ", ".join(str(pos) for pos in site_positions)
                        text += f"<br><b>Позиции:</b> {pos_str}"
                    except (KeyError, TypeError):
                        pass
                
                row_text.append(text)
            hovertext.append(row_text)
        
        # Создание графика с подграфиками для дополнительных элементов управления
        fig = make_subplots(rows=1, cols=1)
        
        # Добавление тепловой карты
        fig.add_trace(go.Heatmap(
            z=heatmap_data.values,
            x=heatmap_data.columns,
            y=heatmap_data.index,
            colorscale='Viridis',
            text=hovertext,
            hoverinfo='text',
            colorbar=dict(title='Сайты связывания'),
            hoverongaps=False
        ))
        
        # Обновление макета
        fig.update_layout(
            title={
                'text': 'Интерактивная тепловая карта сайтов связывания транскрипционных факторов',
                'y':0.95,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            xaxis=dict(
                title='Транскрипционный фактор',
                tickangle=45,
            ),
            yaxis=dict(title='Ген'),
            hoverlabel=dict(
                bgcolor="white",
                font_size=12,
                font_family="Arial"
            ),
            height=800,
            margin=dict(l=100, r=20, t=100, b=100),
            autosize=True,
        )
        
        # Добавление обработчика кликов для отображения подробной информации с использованием JavaScript
        js_code = """
        <script>
        var heatmapData = %s;
        var tfDatabase = %s;
        var geneDatabase = %s;
        
        function showDetails(gene, tf, tfId) {
            var detailsDiv = document.getElementById('details-panel');
            if (!detailsDiv) {
                detailsDiv = document.createElement('div');
                detailsDiv.id = 'details-panel';
                detailsDiv.style.position = 'fixed';
                detailsDiv.style.right = '20px';
                detailsDiv.style.top = '100px';
                detailsDiv.style.width = '500px';
                detailsDiv.style.maxHeight = '80%%';
                detailsDiv.style.overflow = 'auto';
                detailsDiv.style.backgroundColor = 'white';
                detailsDiv.style.border = '1px solid #ccc';
                detailsDiv.style.borderRadius = '5px';
                detailsDiv.style.padding = '15px';
                detailsDiv.style.boxShadow = '0 2px 10px rgba(0,0,0,0.2)';
                detailsDiv.style.zIndex = '1000';
                document.body.appendChild(detailsDiv);
            }
            
            // Поиск данных о сайтах связывания
            try {
                var bindingData = heatmapData[tfId]['binding_data'][gene];
                var sites = bindingData['sites'];
                
                var html = '<div style="text-align:right"><button onclick="closeDetails()" style="border:none;background:none;cursor:pointer;font-size:16px;">✖</button></div>';
                
                // Информация о гене
                html += '<h3 style="color:#2c3e50;margin-top:10px;">' + gene + '</h3>';
                if (geneDatabase[gene]) {
                    var geneInfo = geneDatabase[gene];
                    html += '<div class="gene-info" style="background-color:#f8f9fa;padding:10px;border-radius:5px;margin-bottom:15px;">';
                    html += '<p><strong>Описание:</strong> ' + geneInfo.description + '</p>';
                    html += '<p><strong>Функция:</strong> ' + geneInfo.function + '</p>';
                    html += '<p><strong>Сигнальный путь:</strong> ' + geneInfo.pathway + '</p>';
                    html += '</div>';
                }
                
                // Информация о ТФ
                html += '<h3 style="color:#2c3e50;margin-top:10px;">' + tf + ' (' + tfId + ')</h3>';
                if (tfDatabase[tfId]) {
                    var tfInfo = tfDatabase[tfId];
                    html += '<div class="tf-info" style="background-color:#f8f9fa;padding:10px;border-radius:5px;margin-bottom:15px;">';
                    html += '<p><strong>Полное название:</strong> ' + tfInfo.description + '</p>';
                    html += '<p><strong>Семейство:</strong> ' + tfInfo.family + '</p>';
                    html += '<p><strong>Функция:</strong> ' + tfInfo.function + '</p>';
                    html += '</div>';
                }
                
                // Информация о сайтах связывания
                html += '<div class="binding-info" style="margin-top:15px;">';
                html += '<h3 style="color:#2c3e50;">Информация о связывании</h3>';
                html += '<p><strong>Всего сайтов связывания:</strong> ' + bindingData['total_sites'] + '</p>';
                html += '<p><strong>Средний скор:</strong> ' + bindingData['avg_score'].toFixed(3) + '</p>';
                html += '<p><strong>Максимальный скор:</strong> ' + bindingData['max_score'].toFixed(3) + '</p>';
                html += '</div>';
                
                if (sites && sites.length > 0) {
                    html += '<h4 style="color:#2c3e50;margin-top:15px;">Сайты связывания</h4>';
                    html += '<table style="width:100%%; border-collapse:collapse;">';
                    html += '<tr style="background-color:#f2f2f2;"><th style="text-align:left;padding:5px;border:1px solid #ddd;">Позиция</th>';
                    html += '<th style="text-align:left;padding:5px;border:1px solid #ddd;">Последовательность</th>';
                    html += '<th style="text-align:left;padding:5px;border:1px solid #ddd;">Скор</th>';
                    html += '<th style="text-align:left;padding:5px;border:1px solid #ddd;">Цепь</th></tr>';
                    
                    for (var i = 0; i < sites.length; i++) {
                        var site = sites[i];
                        html += '<tr style="' + (i %% 2 === 0 ? 'background-color:#f9f9f9;' : '') + '">';
                        html += '<td style="padding:5px;border:1px solid #ddd;">' + site['position'] + '</td>';
                        html += '<td style="padding:5px;border:1px solid #ddd;font-family:monospace;">' + site['sequence'] + '</td>';
                        html += '<td style="padding:5px;border:1px solid #ddd;">' + site['relative_score'].toFixed(3) + '</td>';
                        html += '<td style="padding:5px;border:1px solid #ddd;">' + (site['strand'] === '+' ? 'Прямая' : 'Обратная') + '</td>';
                        html += '</tr>';
                    }
                    html += '</table>';
                }
                
                detailsDiv.innerHTML = html;
                detailsDiv.style.display = 'block';
            } catch (e) {
                console.error("Ошибка при отображении деталей", e);
            }
        }
        
        function closeDetails() {
            var detailsDiv = document.getElementById('details-panel');
            if (detailsDiv) {
                detailsDiv.style.display = 'none';
            }
        }
        
        document.addEventListener('click', function(e) {
            var target = e.target;
            // Если клик на ячейке тепловой карты
            if (target.tagName === 'rect' && target.parentNode.classList.contains('heatmaplayer')) {
                // Получение координат - это хак, но работает
                try {
                    var gd = document.querySelector('.js-plotly-plot');
                    var plotData = gd._fullData[0];
                    
                    // Найти индексы строки и столбца по точке клика
                    var xaxis = gd._fullLayout.xaxis;
                    var yaxis = gd._fullLayout.yaxis;
                    
                    var x = xaxis.p2c(e.clientX - gd.getBoundingClientRect().left);
                    var y = yaxis.p2c(e.clientY - gd.getBoundingClientRect().top);
                    
                    var col = Math.floor(x + 0.5);
                    var row = Math.floor(y + 0.5);
                    
                    if (row >= 0 && row < plotData.y.length && col >= 0 && col < plotData.x.length) {
                        var gene = plotData.y[row];
                        var tf = plotData.x[col];
                        
                        // Найти ID ТФ из данных
                        for (var tfId in heatmapData) {
                            if (heatmapData[tfId]['tf_name'] === tf) {
                                showDetails(gene, tf, tfId);
                                break;
                            }
                        }
                    }
                } catch (e) {
                    console.error("Ошибка при обработке клика", e);
                }
            }
        });
        </script>
        """
        
        # Подготовка JSON данных для использования в JS
        tf_data_js = {}
        for tf_id, tf_obj in TF_DATABASE.items():
            tf_data_js[tf_id] = {
                'name': tf_obj.name,
                'description': tf_obj.description,
                'family': tf_obj.family,
                'function': tf_obj.function
            }
        
        gene_data_js = {}
        for gene_name, gene_obj in GENE_DATABASE.items():
            gene_data_js[gene_name] = {
                'name': gene_obj.name,
                'description': gene_obj.description,
                'function': gene_obj.function,
                'pathway': gene_obj.pathway
            }
        
        # Формирование JS с параметрами
        js_code = js_code % (json.dumps(binding_results), json.dumps(tf_data_js), json.dumps(gene_data_js))
        
        if output_file is None:
            output_file = os.path.join(self.results_dir, 'interactive_heatmap.html')
        
        # Запись HTML с внедренным JavaScript
        with open(output_file, 'w', encoding='utf-8') as f:
            raw_html = pio.to_html(fig, full_html=True, include_plotlyjs='cdn')
            # Вставка пользовательского JavaScript перед </body>
            final_html = raw_html.replace('</body>', js_code + '</body>')
            f.write(final_html)
        
        print(f"Подробная интерактивная тепловая карта сохранена в {output_file}")
        return output_file

if __name__ == "__main__":
    heatmap = InteractiveHeatmap()
    heatmap.create_detailed_heatmap() 