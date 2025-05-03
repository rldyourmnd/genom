import sys

def remove_duplicates(input_file, output_file=None):
    seen_lines = set()
    output_lines = []
    
    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f:
            line_stripped = line.strip()
            if line_stripped not in seen_lines:
                seen_lines.add(line_stripped)
                output_lines.append(line)
    
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.writelines(output_lines)
    else:
        with open(input_file, 'w', encoding='utf-8') as f:
            f.writelines(output_lines)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python remove_duplicates.py input_file [output_file]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    remove_duplicates(input_file, output_file)