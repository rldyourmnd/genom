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

# Функция для загрузки примеров сайтов связывания
def load_examples(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        return content
    except Exception as e:
        print(f"Error loading examples from {file_path}: {e}")
        return ""

# Функция создания HTML отчета
def create_html_report():
    # Загрузка данных
    tf_counts = pd.read_csv('tf_binding_sites_counts.csv')
    
    # Создание HTML документа
    html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Transcription Factor Analysis for ATF3 Promoter</title>
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
            .tf-section {
                margin-bottom: 40px;
                border: 1px solid #ddd;
                padding: 15px;
                border-radius: 5px;
            }
            .logo-heatmap {
                display: flex;
                flex-wrap: wrap;
                gap: 20px;
                margin-top: 15px;
            }
            .image-container {
                max-width: 45%;
            }
            .image-container img {
                max-width: 100%;
                height: auto;
            }
            .examples {
                background-color: #f8f9fa;
                padding: 15px;
                border-radius: 5px;
                font-family: monospace;
                white-space: pre;
                overflow-x: auto;
            }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Transcription Factor Analysis for ATF3 Promoter</h1>
            
            <div class="summary">
                <h2>Analysis Summary</h2>
                <p>This report presents the analysis of the ATF3 promoter sequence for potential transcription factor binding sites.
                The analysis identified binding motifs for various transcription factors and their distribution along the sequence.</p>
                <p>Length of analyzed sequence: 57,944 nucleotides</p>
            </div>
            
            <h2>Overview of Transcription Factor Binding Sites</h2>
            
            <table>
                <tr>
                    <th>Transcription Factor</th>
                    <th>Binding Sites Count</th>
                    <th>Description</th>
                </tr>
    """
    
    # Dictionary with descriptions of transcription factors
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
        'HRE': 'Hormone Response Element, mediates transcriptional responses to hormones'
    }
    
    # Sort TF counts in descending order
    tf_counts_sorted = tf_counts.sort_values('Binding Sites Count', ascending=False)
    
    # Add rows to the table
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
            
            <h2>Overall Distribution Visualizations</h2>
            
            <div class="logo-heatmap">
                <div class="image-container">
                    <h3>Binding Sites Count</h3>
                    <img src="data:image/png;base64,{0}" alt="TF Binding Sites Count">
                </div>
                <div class="image-container">
                    <h3>Distribution Along Sequence</h3>
                    <img src="data:image/png;base64,{1}" alt="TF Distribution">
                </div>
                <div class="image-container">
                    <h3>Density Heatmap</h3>
                    <img src="data:image/png;base64,{2}" alt="TF Density Heatmap">
                </div>
            </div>
            
            <h2>Detailed Analysis by Transcription Factor</h2>
    """.format(
        img_to_base64('tf_binding_sites_count.png'),
        img_to_base64('tf_binding_sites_distribution.png'),
        img_to_base64('tf_binding_sites_heatmap.png')
    )
    
    # Add sections for each transcription factor
    for tf_name in tf_counts_sorted['Transcription Factor']:
        # Skip factors with no binding sites
        if tf_counts_sorted[tf_counts_sorted['Transcription Factor'] == tf_name]['Binding Sites Count'].values[0] == 0:
            continue
            
        logo_path = f'tf_logos/{tf_name}_logo.png'
        heatmap_path = f'tf_heatmaps/{tf_name}_heatmap.png'
        examples_path = f'tf_examples/{tf_name}_examples.txt'
        
        # Check if files exist
        has_logo = os.path.exists(logo_path)
        has_heatmap = os.path.exists(heatmap_path)
        examples = load_examples(examples_path) if os.path.exists(examples_path) else ""
        
        html += f"""
            <div class="tf-section">
                <h3>{tf_name}</h3>
                <p>{tf_descriptions.get(tf_name, "")}</p>
                
                <div class="logo-heatmap">
        """
        
        if has_logo:
            html += f"""
                    <div class="image-container">
                        <h4>Motif Logo</h4>
                        <img src="data:image/png;base64,{img_to_base64(logo_path)}" alt="{tf_name} Logo">
                    </div>
            """
            
        if has_heatmap:
            html += f"""
                    <div class="image-container">
                        <h4>Nucleotide Frequencies</h4>
                        <img src="data:image/png;base64,{img_to_base64(heatmap_path)}" alt="{tf_name} Heatmap">
                    </div>
            """
            
        html += """
                </div>
        """
        
        if examples:
            html += f"""
                <h4>Example Binding Sites</h4>
                <div class="examples">{examples}</div>
            """
            
        html += """
            </div>
        """
    
    # Close the HTML document
    html += """
        </div>
    </body>
    </html>
    """
    
    # Write the HTML to a file
    with open('tf_analysis_report.html', 'w', encoding='utf-8') as f:
        f.write(html)
    
    print("HTML report generated: tf_analysis_report.html")

if __name__ == "__main__":
    create_html_report() 