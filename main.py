#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import argparse
from simplified_binding_analysis import SimplifiedBindingAnalyzer
from create_html_visualization import HTMLVisualizer
from interactive_heatmap import InteractiveHeatmap

def parse_args():
    parser = argparse.ArgumentParser(description='Анализ сайтов связывания транскрипционных факторов')
    parser.add_argument('--genes-dir', type=str, default='genes', 
                        help='Директория с FASTA файлами генов')
    parser.add_argument('--results-dir', type=str, default='results', 
                        help='Директория для сохранения результатов')
    parser.add_argument('--html-dir', type=str, default='html_report', 
                        help='Директория для HTML визуализации')
    parser.add_argument('--threshold', type=float, default=0.75, 
                        help='Порог для определения сайтов связывания (0.0-1.0)')
    parser.add_argument('--download-motifs', action='store_true', 
                        help='Загрузить мотивы JASPAR перед анализом')
    parser.add_argument('--skip-analysis', action='store_true', 
                        help='Пропустить анализ генов, только сгенерировать визуализацию')
    parser.add_argument('--interactive-only', action='store_true',
                        help='Только создать интерактивную визуализацию')
    return parser.parse_args()

def main():
    # Разбор аргументов командной строки
    args = parse_args()
    
    # Создание директории результатов, если она не существует
    if not os.path.exists(args.results_dir):
        os.makedirs(args.results_dir)
    
    # Создание HTML директории, если она не существует
    if not os.path.exists(args.html_dir):
        os.makedirs(args.html_dir)
    
    # Проверка, нужно ли загрузить мотивы
    if args.download_motifs:
        print("=== Шаг 1: Загрузка мотивов JASPAR ===")
        try:
            from download_known_tfs import download_tfs
            download_tfs()
        except Exception as e:
            print(f"Ошибка при загрузке мотивов: {e}")
            print("Используем существующие мотивы...")
    
    if args.interactive_only:
        print("=== Генерация только интерактивной визуализации ===")
        interactive_heatmap = InteractiveHeatmap(results_dir=args.results_dir)
        interactive_heatmap.create_detailed_heatmap()
        print("Интерактивная тепловая карта создана. Открываем HTML отчет...")
        try:
            import webbrowser
            webbrowser.open(os.path.join(args.html_dir, 'index.html'))
        except Exception:
            pass
        return
    
    if not args.skip_analysis:
        # Шаг 2: Анализ последовательностей генов на сайты связывания ТФ
        print("=== Шаг 2: Анализ последовательностей генов на сайты связывания ТФ ===")
        analyzer = SimplifiedBindingAnalyzer(genes_dir=args.genes_dir, results_dir=args.results_dir)
        analyzer.load_gene_sequences()
        analyzer.analyze_binding_sites(threshold=args.threshold)
        
        # Шаг 3: Сохранение результатов и генерация визуализаций
        print("=== Шаг 3: Сохранение результатов и генерация визуализаций ===")
        results_file = analyzer.save_results(os.path.join(args.results_dir, 'binding_analysis_results.json'))
        analyzer.generate_heatmap()
        analyzer.generate_binding_site_distribution()
        
        # Создание интерактивной тепловой карты
        print("=== Шаг 4: Создание интерактивной тепловой карты ===")
        interactive_heatmap = InteractiveHeatmap(results_dir=args.results_dir)
        interactive_heatmap.create_detailed_heatmap()
    else:
        print("=== Пропуск анализа, используем существующие результаты ===")
        # Все равно создаем интерактивную тепловую карту, если она не существует
        interactive_heatmap = InteractiveHeatmap(results_dir=args.results_dir)
        interactive_heatmap.create_detailed_heatmap()
    
    # Шаг 5: Генерация HTML визуализации
    print("=== Шаг 5: Генерация HTML визуализации ===")
    try:
        visualizer = HTMLVisualizer(results_dir=args.results_dir, output_dir=args.html_dir)
        html_dir = visualizer.generate_html_report()
        
        print("=== Анализ завершен ===")
        print(f"Результаты сохранены в: {args.results_dir}")
        print(f"HTML визуализация доступна в: {args.html_dir}")
        print(f"Откройте {os.path.join(args.html_dir, 'index.html')} в веб-браузере для просмотра отчета")
        
        # Пробуем открыть отчет в веб-браузере
        try:
            import webbrowser
            webbrowser.open(os.path.join(args.html_dir, 'index.html'))
        except Exception:
            pass
    except Exception as e:
        print(f"Ошибка при создании HTML отчета: {e}")
        print("=== Анализ завершен с ошибками ===")
        print(f"Результаты сохранены в: {args.results_dir}")

if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    print(f"Общее время выполнения: {elapsed_time:.2f} секунд") 