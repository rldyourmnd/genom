import os
import pandas as pd
import base64
from io import BytesIO
from PIL import Image
import markdown
import re

def img_to_base64(img_path):
    """Convert an image to base64 for embedding in HTML"""
    try:
        with open(img_path, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
            return encoded_string
    except Exception as e:
        print(f"Error converting image {img_path} to base64: {e}")
        return ""

def md_to_html(md_file):
    """Convert markdown to HTML with image embedding"""
    with open(md_file, 'r') as f:
        md_content = f.read()
    
    # Convert markdown to HTML
    html_content = markdown.markdown(md_content, extensions=['tables'])
    
    # Replace image links with base64 embedded images
    img_tags = re.findall(r'!\[.*?\]\((.*?)\)', md_content)
    for img_path in img_tags:
        base64_img = img_to_base64(img_path)
        if base64_img:
            img_tag = f'<img src="data:image/png;base64,{base64_img}" alt="{os.path.basename(img_path)}" style="max-width:100%;">'
            html_content = html_content.replace(f'<img src="{img_path}" alt', img_tag.replace('alt', ' alt'))
    
    return html_content

def create_html_report():
    """Create HTML report with enhanced styling and image embedding"""
    # Get HTML content from markdown conversion
    md_html = md_to_html('atf_comparison_report.md')
    
    # Create HTML with styling
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>ATF Sequences Comparison Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                line-height: 1.6;
                margin: 0;
                padding: 20px;
                color: #333;
                max-width: 1200px;
                margin: 0 auto;
            }}
            h1, h2, h3 {{
                color: #2c3e50;
            }}
            h1 {{
                border-bottom: 2px solid #3498db;
                padding-bottom: 10px;
            }}
            h2 {{
                border-bottom: 1px solid #3498db;
                padding-bottom: 5px;
                margin-top: 30px;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
            }}
            th, td {{
                border: 1px solid #ddd;
                padding: 12px;
                text-align: left;
            }}
            th {{
                background-color: #f2f2f2;
                font-weight: bold;
            }}
            tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}
            tr:hover {{
                background-color: #f1f1f1;
            }}
            img {{
                max-width: 100%;
                height: auto;
                display: block;
                margin: 20px auto;
                border: 1px solid #ddd;
                border-radius: 5px;
                box-shadow: 0 3px 10px rgba(0,0,0,0.1);
            }}
            .conclusion {{
                background-color: #f8f9fa;
                padding: 15px;
                border-left: 5px solid #3498db;
                margin: 20px 0;
            }}
            .conclusion ol {{
                margin: 10px 0;
                padding-left: 20px;
            }}
            .conclusion li {{
                margin-bottom: 8px;
            }}
            .side-by-side {{
                display: flex;
                justify-content: space-between;
                flex-wrap: wrap;
            }}
            .side-by-side > div {{
                width: 48%;
            }}
            @media (max-width: 768px) {{
                .side-by-side > div {{
                    width: 100%;
                    margin-bottom: 20px;
                }}
            }}
            .details-section {{
                margin: 30px 0;
            }}
            .alignment-results {{
                background-color: #ebf5fb;
                padding: 15px;
                border-radius: 5px;
            }}
            .alignment-results ul {{
                list-style-type: none;
                padding-left: 0;
            }}
            .alignment-results li {{
                margin-bottom: 8px;
                padding-left: 20px;
                position: relative;
            }}
            .alignment-results li:before {{
                content: "•";
                position: absolute;
                left: 0;
                color: #3498db;
                font-weight: bold;
            }}
            .tf-images {{
                display: flex;
                flex-wrap: wrap;
                justify-content: center;
                gap: 20px;
                margin-top: 20px;
            }}
            .tf-image {{
                flex: 0 0 calc(50% - 20px);
                max-width: calc(50% - 20px);
                margin-bottom: 20px;
            }}
            .tf-image img {{
                width: 100%;
            }}
            .tf-image h4 {{
                text-align: center;
                margin: 10px 0;
            }}
            .header {{
                display: flex;
                align-items: center;
                justify-content: space-between;
                margin-bottom: 20px;
            }}
            .header-left {{
                flex: 1;
            }}
            .header-right {{
                text-align: right;
                color: #7f8c8d;
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <div class="header-left">
                <h1>ATF Sequences Comparative Analysis</h1>
            </div>
            <div class="header-right">
                <p>Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d')}</p>
            </div>
        </div>
        
        {md_html}
        
        <div class="details-section">
            <h2>Transcription Factor Distribution Patterns</h2>
            <div class="tf-images">
    """
    
    # Add TF position comparison images
    tf_imgs = [f for f in os.listdir('.') if f.startswith('tf_position_comparison_') and f.endswith('.png')]
    for img in tf_imgs:
        tf_name = img.replace('tf_position_comparison_', '').replace('.png', '')
        html += f"""
                <div class="tf-image">
                    <h4>{tf_name} Distribution</h4>
                    <img src="data:image/png;base64,{img_to_base64(img)}" alt="{tf_name} Distribution">
                </div>
        """
    
    html += """
            </div>
        </div>
        
        <div class="conclusion">
            <h2>Key Insights</h2>
            <p>Based on the comprehensive analysis of both ATF2500 and ATF3 sequences, we can observe several important patterns:</p>
            <ol>
                <li>Both promoters show strong enrichment for ETS family binding sites, suggesting a conserved regulatory mechanism.</li>
                <li>The ATF2500 sequence (2,501 bp) contains proportionally more GC-box and SP1 binding sites than the longer ATF3 sequence (57,944 bp).</li>
                <li>The sequence alignment shows 75.45% identity in the aligned regions, indicating significant conservation despite the size difference.</li>
                <li>The distribution of binding sites along both sequences shows clustering patterns that may be biologically significant.</li>
                <li>While ATF2500 has 5 unique transcription factor binding sites not found in ATF3, all transcription factors found in ATF3 are also present in ATF2500.</li>
            </ol>
            <p>These findings provide valuable insights for further experimental investigation of the regulatory mechanisms of these promoters.</p>
        </div>
    </body>
    </html>
    """
    
    # Write HTML to file
    with open('atf_comparison_report.html', 'w', encoding='utf-8') as f:
        f.write(html)
    
    print("HTML comparison report generated: atf_comparison_report.html")

if __name__ == "__main__":
    create_html_report() 