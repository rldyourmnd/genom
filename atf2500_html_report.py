import os
import pandas as pd
import base64
from io import BytesIO
import matplotlib.pyplot as plt
from PIL import Image

# Функция для конвертации изображений в base64 для встраивания в HTML
def img_to_base64(img_path):
    try:
        with open(img_path, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
            return encoded_string
    except Exception as e:
        print(f"Error converting image {img_path} to base64: {e}")
        return ""

# Функция для создания HTML отчета
def create_html_report():
    # Загрузка данных
    tf_counts = pd.read_csv('atf2500_tf_binding_sites_counts.csv')
    tf_positions = pd.read_csv('atf2500_tf_binding_sites_positions.csv')
    
    # Сортировка факторов по количеству сайтов
    tf_counts_sorted = tf_counts.sort_values('Binding Sites Count', ascending=False)
    
    # Создание HTML документа
    html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>ATF2500 Promoter Analysis</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                line-height: 1.6;
                margin: 0;
                padding: 20px;
                color: #333;
            }
            h1, h2, h3 {
                color: #2c3e50;
            }
            .container {
                max-width: 1200px;
                margin: 0 auto;
            }
            .summary {
                background-color: #f5f5f5;
                padding: 15px;
                border-radius: 5px;
                margin-bottom: 20px;
            }
            table {
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
            }
            th, td {
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }
            th {
                background-color: #f2f2f2;
            }
            tr:nth-child(even) {
                background-color: #f9f9f9;
            }
            .visualization {
                margin-bottom: 30px;
            }
            .visualization img {
                max-width: 100%;
                height: auto;
                display: block;
                margin: 0 auto;
                border: 1px solid #ddd;
                border-radius: 5px;
            }
            .tf-section {
                margin-bottom: 20px;
                border: 1px solid #ddd;
                padding: 15px;
                border-radius: 5px;
            }
            .binding-sites {
                font-family: monospace;
                max-height: 200px;
                overflow-y: auto;
                background-color: #f8f9fa;
                padding: 10px;
                border-radius: 5px;
            }
            .site-position {
                font-weight: bold;
                color: #2980b9;
            }
            .site-pattern {
                color: #27ae60;
            }
            .comparison {
                display: flex;
                justify-content: space-between;
                margin-top: 20px;
            }
            .comparison-chart {
                width: 48%;
            }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>ATF2500 Promoter Transcription Factor Analysis</h1>
            
            <div class="summary">
                <h2>Analysis Summary</h2>
                <p>This report presents the analysis of potential transcription factor binding sites in the ATF2500 promoter sequence.</p>
                <ul>
                    <li><strong>Sequence Length:</strong> 2501 nucleotides</li>
                    <li><strong>Analyzed TFs:</strong> 20 different transcription factors</li>
                    <li><strong>Total Binding Sites:</strong> 146 potential binding sites</li>
                </ul>
            </div>
            
            <h2>Transcription Factor Binding Sites Overview</h2>
            
            <table>
                <tr>
                    <th>Transcription Factor</th>
                    <th>Binding Sites Count</th>
                    <th>Description</th>
                </tr>
    """
    
    # Словарь с описаниями транскрипционных факторов
    tf_descriptions = {
        'ETS': 'ETS family factors play key roles in cell differentiation, cell cycle control, and cell migration',
        'IRF': 'Interferon regulatory factors involved in immune response and cell cycle regulation',
        'E-box': 'Binding sites for bHLH transcription factors, important for tissue-specific gene expression',
        'NFAT': 'Nuclear factor of activated T-cells, involved in immune response',
        'CAAT-box': 'Promoter element important for RNA polymerase binding and transcription initiation',
        'TATA-box': 'Core promoter element required for accurate transcription initiation',
        'GC-box': 'Promoter element often found in housekeeping genes, binds SP1',
        'SP1': 'Specificity Protein 1, constitutively expressed transcription factor',
        'CREB': 'cAMP response element-binding protein, regulates transcription in response to cAMP levels',
        'STAT': 'Signal Transducer and Activator of Transcription, activated by cytokines and growth factors',
        'AP-1': 'Activator Protein 1, regulates gene expression in response to various stimuli',
        'CRE': 'cAMP Response Element, regulates expression in response to cAMP levels',
        'NF-kB': 'Nuclear Factor kappa B, involved in immune and inflammatory responses',
        'GATA': 'GATA binding proteins, important for development and differentiation',
        'HRE': 'Hormone Response Element, mediates transcriptional responses to hormones',
        'ATF/CREB': 'ATF/CREB (Activating Transcription Factor/cAMP Response Element-Binding protein) family members, involved in cellular stress response',
        'AP-2': 'Activator Protein 2, involved in development and differentiation',
        'OCT': 'Octamer-binding transcription factor, important for embryonic development',
        'YY1': 'Yin Yang 1, multifunctional transcription factor with roles in development and cancer',
        'NRF-1': 'Nuclear Respiratory Factor 1, regulates mitochondrial biogenesis and function'
    }
    
    # Добавление строк в таблицу
    for _, row in tf_counts_sorted.iterrows():
        tf_name = row['Transcription Factor']
        count = row['Binding Sites Count']
        description = tf_descriptions.get(tf_name, "")
        
        html += f"""
        <tr>
            <td>{tf_name}</td>
            <td>{count}</td>
            <td>{description}</td>
        </tr>
        """
    
    html += """
            </table>
            
            <h2>Binding Sites Visualization</h2>
            
            <div class="visualization">
                <h3>Transcription Factor Binding Sites Count</h3>
                <img src="data:image/png;base64,{0}" alt="TF Binding Sites Count">
            </div>
            
            <div class="visualization">
                <h3>Distribution Along Sequence</h3>
                <img src="data:image/png;base64,{1}" alt="TF Distribution">
            </div>
            
            <div class="visualization">
                <h3>Binding Sites Density Heatmap</h3>
                <img src="data:image/png;base64,{2}" alt="TF Density Heatmap">
            </div>
            
            <h2>Top Transcription Factors Details</h2>
    """.format(
        img_to_base64('atf2500_tf_binding_sites_count.png'),
        img_to_base64('atf2500_tf_binding_sites_distribution.png'),
        img_to_base64('atf2500_tf_binding_sites_heatmap.png')
    )
    
    # Детальная информация о топ-5 факторах
    top_tfs = tf_counts_sorted.head(5)['Transcription Factor'].tolist()
    
    for tf_name in top_tfs:
        # Получаем все сайты связывания для данного TF
        tf_sites = tf_positions[tf_positions['Transcription Factor'] == tf_name]
        
        html += f"""
            <div class="tf-section">
                <h3>{tf_name}</h3>
                <p>{tf_descriptions.get(tf_name, "")}</p>
                <p><strong>Number of binding sites:</strong> {len(tf_sites)}</p>
                
                <h4>Binding Sites Examples:</h4>
                <div class="binding-sites">
        """
        
        # Добавляем примеры сайтов связывания (максимум 10)
        for _, site in tf_sites.head(10).iterrows():
            html += f"""
                    <div>
                        <span class="site-position">Position {site['Position']}:</span> 
                        <span class="site-pattern">{site['Pattern']}</span>
                    </div>
            """
        
        if len(tf_sites) > 10:
            html += f"<div>... and {len(tf_sites) - 10} more sites</div>"
        
        html += """
                </div>
            </div>
        """
    
    # Сравнение с ATF3
    html += """
        <h2>Comparison with ATF3 Promoter Analysis</h2>
        
        <p>Below is a comparison of the transcription factor binding sites found in the ATF2500 promoter versus the previously analyzed ATF3 promoter:</p>
        
        <div class="comparison">
            <div class="comparison-chart">
                <h3>ATF2500 Promoter (2,501 bp)</h3>
                <ul>
                    <li>Total binding sites: 146</li>
                    <li>Top factor: ETS (64 sites)</li>
                    <li>High GC-box content (24 sites)</li>
                    <li>Low TATA-box presence (1 site)</li>
                    <li>No CRE, AP-1, NF-kB, GATA, HRE, STAT sites found</li>
                </ul>
            </div>
            
            <div class="comparison-chart">
                <h3>ATF3 Promoter (57,944 bp)</h3>
                <ul>
                    <li>Total binding sites: > 2,600</li>
                    <li>Top factor: ETS (1,268 sites)</li>
                    <li>High E-box content (291 sites)</li>
                    <li>Significant TATA-box presence (146 sites)</li>
                    <li>Contains CRE (2 sites) and AP-1 (12 sites)</li>
                </ul>
            </div>
        </div>
        
        <h2>Conclusions</h2>
        
        <ol>
            <li>The ETS family transcription factors are the most abundant in both promoters, suggesting their importance in ATF regulation.</li>
            <li>The ATF2500 promoter is characterized by a high GC content, with numerous GC-box and SP1 binding sites, indicating it may be a CpG island promoter.</li>
            <li>The ATF2500 promoter has only one TATA-box, suggesting it may rely more on other core promoter elements for transcription initiation.</li>
            <li>The ATF2500 sequence lacks CRE and AP-1 sites that are present in the ATF3 promoter, which may indicate different stress-response mechanisms.</li>
            <li>The binding site pattern suggests the ATF2500 promoter may be regulated differently than the ATF3 promoter, despite both being related to the ATF family.</li>
        </ol>
        
        <p><em>Note: This analysis provides information about potential transcription factor binding sites. Experimental validation is recommended to confirm the functional significance of these sites.</em></p>
    """
    
    # Завершение HTML документа
    html += """
        </div>
    </body>
    </html>
    """
    
    # Запись HTML в файл
    with open('atf2500_analysis_report.html', 'w', encoding='utf-8') as f:
        f.write(html)
    
    print("HTML report generated: atf2500_analysis_report.html")

if __name__ == "__main__":
    create_html_report() 