import os
import base64
import re
import markdown
from jinja2 import Template
import glob

# Функция для преобразования изображения в base64
def img_to_base64(image_path):
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except Exception as e:
        print(f"Ошибка при конвертации изображения {image_path}: {e}")
        return ""

# Функция для чтения markdown файла и преобразования его в HTML
def markdown_to_html(markdown_path):
    try:
        with open(markdown_path, 'r', encoding='utf-8') as md_file:
            md_content = md_file.read()
            # Обработка путей к изображениям
            md_content = re.sub(r'!\[(.*?)\]\((.*?)\)', r'![\1](../\2)', md_content)
            return markdown.markdown(md_content, extensions=['tables', 'fenced_code'])
    except Exception as e:
        print(f"Ошибка при чтении markdown файла {markdown_path}: {e}")
        return "<p>Ошибка при загрузке отчета</p>"

# Функция для создания HTML-отчета
def create_html_report(results_dir='results'):
    # Получаем список всех генов
    gene_dirs = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d)) and d != 'comparison']
    
    # Загружаем отчеты для каждого гена
    gene_reports = {}
    for gene in gene_dirs:
        report_path = os.path.join(results_dir, gene, f"{gene}_analysis_report.md")
        if os.path.exists(report_path):
            gene_reports[gene] = markdown_to_html(report_path)
    
    # Загружаем сравнительный отчет
    comparison_path = os.path.join(results_dir, 'comparison', 'genes_comparison_report.md')
    comparison_report = markdown_to_html(comparison_path) if os.path.exists(comparison_path) else ""
    
    # Загружаем изображения для сравнительного анализа
    comparison_images = {}
    for img_path in glob.glob(os.path.join(results_dir, 'comparison', '*.png')):
        img_name = os.path.basename(img_path)
        comparison_images[img_name] = img_to_base64(img_path)
    
    # HTML шаблон
    template_str = """
    <!DOCTYPE html>
    <html lang="ru">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Анализ транскрипционных факторов в промоторах генов</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                line-height: 1.6;
                margin: 0;
                padding: 0;
                color: #333;
                background-color: #f8f9fa;
            }
            .container {
                max-width: 1200px;
                margin: 0 auto;
                padding: 20px;
            }
            h1, h2, h3, h4 {
                color: #2c3e50;
            }
            h1 {
                text-align: center;
                margin-bottom: 30px;
                border-bottom: 2px solid #3498db;
                padding-bottom: 10px;
            }
            h2 {
                margin-top: 30px;
                border-bottom: 1px solid #e0e0e0;
                padding-bottom: 5px;
            }
            h3 {
                color: #2980b9;
            }
            table {
                border-collapse: collapse;
                width: 100%;
                margin: 20px 0;
            }
            th, td {
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }
            th {
                background-color: #f2f2f2;
                font-weight: bold;
            }
            tr:nth-child(even) {
                background-color: #f9f9f9;
            }
            img {
                max-width: 100%;
                height: auto;
                margin: 10px 0;
                border: 1px solid #ddd;
                border-radius: 4px;
                padding: 5px;
            }
            .tabs {
                overflow: hidden;
                border: 1px solid #ccc;
                background-color: #f1f1f1;
                border-radius: 5px 5px 0 0;
                margin-top: 20px;
            }
            .tab-button {
                background-color: inherit;
                float: left;
                border: none;
                outline: none;
                cursor: pointer;
                padding: 14px 16px;
                transition: 0.3s;
                font-size: 17px;
            }
            .tab-button:hover {
                background-color: #ddd;
            }
            .tab-button.active {
                background-color: #3498db;
                color: white;
            }
            .tab-content {
                display: none;
                padding: 20px;
                border: 1px solid #ccc;
                border-top: none;
                border-radius: 0 0 5px 5px;
                background-color: white;
            }
            .tab-content.active {
                display: block;
            }
            .comparison-tab {
                margin-top: 30px;
            }
            .navbar {
                overflow: hidden;
                background-color: #2c3e50;
                position: sticky;
                top: 0;
                z-index: 100;
            }
            .navbar a {
                float: left;
                display: block;
                color: white;
                text-align: center;
                padding: 14px 20px;
                text-decoration: none;
            }
            .navbar a:hover {
                background-color: #3498db;
            }
            .content {
                margin-top: 20px;
            }
            .comparison-images {
                display: flex;
                flex-wrap: wrap;
                justify-content: center;
                gap: 20px;
                margin-top: 20px;
            }
            .comparison-image {
                max-width: 45%;
                box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                margin-bottom: 20px;
            }
            .gene-selector {
                margin: 20px 0;
                padding: 10px;
                background-color: #f1f1f1;
                border-radius: 5px;
            }
            select {
                padding: 8px;
                border-radius: 4px;
                border: 1px solid #ddd;
                width: 200px;
            }
            button {
                padding: 8px 12px;
                background-color: #3498db;
                color: white;
                border: none;
                border-radius: 4px;
                cursor: pointer;
            }
            button:hover {
                background-color: #2980b9;
            }
            footer {
                text-align: center;
                margin-top: 50px;
                padding: 20px;
                background-color: #2c3e50;
                color: white;
            }
        </style>
    </head>
    <body>
        <div class="navbar">
            <a href="#home">Главная</a>
            <a href="#genes">Гены</a>
            <a href="#comparison">Сравнительный анализ</a>
        </div>
        
        <div class="container">
            <section id="home">
                <h1>Анализ транскрипционных факторов в промоторах генов</h1>
                <div class="content">
                    <p>В данном отчете представлены результаты анализа сайтов связывания транскрипционных факторов в промоторных последовательностях различных генов. Анализ включает в себя:</p>
                    <ul>
                        <li>Идентификацию потенциальных сайтов связывания различных транскрипционных факторов</li>
                        <li>Количественный анализ обнаруженных сайтов</li>
                        <li>Визуализацию распределения сайтов в промоторах</li>
                        <li>Сравнительный анализ паттернов связывания между генами</li>
                    </ul>
                    <p>Для просмотра результатов анализа отдельных генов перейдите во вкладку "Гены". Для просмотра сравнительного анализа всех генов перейдите во вкладку "Сравнительный анализ".</p>
                </div>
            </section>
            
            <section id="genes">
                <h2>Анализ отдельных генов</h2>
                <div class="gene-selector">
                    <label for="gene-select">Выберите ген: </label>
                    <select id="gene-select">
                        {% for gene in gene_reports %}
                        <option value="{{ gene }}">{{ gene }}</option>
                        {% endfor %}
                    </select>
                    <button onclick="showGeneReport()">Показать отчет</button>
                </div>
                
                <div id="gene-reports">
                    {% for gene, report in gene_reports.items() %}
                    <div id="{{ gene }}-report" class="gene-report" style="display: none;">
                        {{ report|safe }}
                    </div>
                    {% endfor %}
                </div>
            </section>
            
            <section id="comparison">
                <h2>Сравнительный анализ генов</h2>
                <div class="comparison-content">
                    {{ comparison_report|safe }}
                    
                    <h3>Визуализации сравнительного анализа</h3>
                    <div class="comparison-images">
                        {% for img_name, img_data in comparison_images.items() %}
                        <div class="comparison-image">
                            <h4>{{ img_name.replace('_', ' ').replace('.png', '') }}</h4>
                            <img src="data:image/png;base64,{{ img_data }}" alt="{{ img_name }}">
                        </div>
                        {% endfor %}
                    </div>
                </div>
            </section>
        </div>
        
        <footer>
            <p>© 2025 Анализ транскрипционных факторов в промоторах генов</p>
        </footer>
        
        <script>
            function showGeneReport() {
                const geneSelect = document.getElementById('gene-select');
                const selectedGene = geneSelect.value;
                
                // Скрываем все отчеты
                const geneReports = document.getElementsByClassName('gene-report');
                for (let i = 0; i < geneReports.length; i++) {
                    geneReports[i].style.display = 'none';
                }
                
                // Показываем выбранный отчет
                document.getElementById(selectedGene + '-report').style.display = 'block';
                
                // Заменяем относительные пути к изображениям
                fixImagePaths(selectedGene);
            }
            
            function fixImagePaths(gene) {
                const report = document.getElementById(gene + '-report');
                const images = report.getElementsByTagName('img');
                for (let i = 0; i < images.length; i++) {
                    const src = images[i].getAttribute('src');
                    if (src && src.includes('../')) {
                        images[i].src = src.replace('../', '');
                    }
                }
            }
            
            // Показываем первый ген при загрузке страницы
            window.onload = function() {
                const geneSelect = document.getElementById('gene-select');
                if (geneSelect && geneSelect.options.length > 0) {
                    const firstGene = geneSelect.options[0].value;
                    document.getElementById(firstGene + '-report').style.display = 'block';
                    fixImagePaths(firstGene);
                }
            };
            
            // Плавная прокрутка к секциям
            document.querySelectorAll('.navbar a').forEach(anchor => {
                anchor.addEventListener('click', function(e) {
                    e.preventDefault();
                    const targetId = this.getAttribute('href');
                    document.querySelector(targetId).scrollIntoView({
                        behavior: 'smooth'
                    });
                });
            });
        </script>
    </body>
    </html>
    """
    
    # Создаем HTML-отчет с использованием шаблона Jinja2
    template = Template(template_str)
    html_content = template.render(
        gene_reports=gene_reports,
        comparison_report=comparison_report,
        comparison_images=comparison_images
    )
    
    # Сохраняем HTML-отчет
    output_path = 'genes_tf_analysis_report.html'
    with open(output_path, 'w', encoding='utf-8') as html_file:
        html_file.write(html_content)
    
    print(f"HTML-отчет успешно создан: {output_path}")
    return output_path

if __name__ == "__main__":
    create_html_report() 